/*
* CIS 441 Project 1C.
* 
* My version of Visual Studio won't let me use a C++ compiler version prior to 2014, so I'll do my best to avoid using newer syntax that doesn't work with your compiler.
* 
* Thomas Mitchell
*/

#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>


double ceil__441(double f)
{
    return ceil(f - 0.00001);
}

double floor__441(double f)
{
    return floor(f + 0.00001);
}


vtkImageData*
NewImage(int height, int width)
{
    vtkImageData* img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData* img, const char* filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter* writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

class Triangle
{
public:
    double         X[3];
    double         Y[3];
    unsigned char color[3];

    //these are used to simplify returning the lower/upper rows for each column queried
    double bottomSlope, bottomStartPt[2];
    double topSlope, topStartPt[2];

    void sortPointsByX() //this will sort X[] values (ascending) for this triangle and shuffle the Y[] as well so it matches
                         //using vectors and built-in sorting methods was much too slow for the volume of triangles we're processing, so I had to make this
    {
        double points[3][2];
        double tmp[2];
        
        //load our X's and Y's into our working space
        for (int i = 0; i < 3; i++)
        {
            points[i][0] = X[i];
            points[i][1] = Y[i];
        }

        //sorting method
        for (int c = 0; c < 2; c++)
        {
            for (int i = 0; i < 2; i++)
            {
                if (points[i][0] > points[i + 1][0])
                {
                    tmp[0] = points[i][0];
                    tmp[1] = points[i][1];
                    points[i][0] = points[i + 1][0];
                    points[i][1] = points[i + 1][1];
                    points[i + 1][0] = tmp[0];
                    points[i + 1][1] = tmp[1];
                }
            }
        }

        //set the class member variables appropriately
        for (int i = 0; i < 3; i++)
        {
            X[i] = points[i][0];
            Y[i] = points[i][1];
        }
    }

    void prepare() //run this before using getRowBounds()
    {
        //determine the points with matching x values, as well as the point with a unique x


        if (X[0] != X[1] && X[0] != X[2]) //left-pointing triangle
        {
            bottomStartPt[0] = X[0];
            bottomStartPt[1] = Y[0];

            bottomSlope = (std::min(Y[1], Y[2]) - Y[0]) / (X[1] - X[0]);

            topStartPt[0] = X[0];
            topStartPt[1] = Y[0];

            topSlope = (std::max(Y[1], Y[2]) - Y[0]) / (X[1] - X[0]);
        }
        else //(X[2] != X[0] && X[2] != X[1]) //right-pointing triangle
        {
            bottomStartPt[0] = X[0];
            bottomStartPt[1] = std::min(Y[0], Y[1]);

            bottomSlope = (Y[2] - std::min(Y[0], Y[1])) / (X[2] - X[0]);

            topStartPt[0] = X[0];
            topStartPt[1] = std::max(Y[0], Y[1]);

            topSlope = (Y[2] - std::max(Y[0], Y[1])) / (X[2] - X[0]);
        }
    }

    void getRowBounds(int column, int *rowLimits) //get [lower bound, upper bound] for given column. Make sure to run the above funtion first
    {
        rowLimits[0] = ceil__441(bottomStartPt[1] + bottomSlope * (column - bottomStartPt[0]));
        rowLimits[1] = floor__441(topStartPt[1] + topSlope * (column - topStartPt[0])) + 1;

    }
};

class Screen
{
public:
    unsigned char* buffer;
    int width, height;
};

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader* rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1c_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData* pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints* pts = pd->GetPoints();
    vtkCellArray* cells = pd->GetPolys();
    vtkFloatArray* colors = (vtkFloatArray*)pd->GetPointData()->GetArray("color_nodal");
    float* color_ptr = colors->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType* ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal(); cells->GetNextCell(npts, ptIds); idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
        tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
        tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
        tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
        tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
        tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
        tris[idx].color[0] = (unsigned char)color_ptr[4 * ptIds[0] + 0];
        tris[idx].color[1] = (unsigned char)color_ptr[4 * ptIds[0] + 1];
        tris[idx].color[2] = (unsigned char)color_ptr[4 * ptIds[0] + 2];
    }
    cerr << "Done reading" << endl;

    return tris;
}

int main()
{
    vtkImageData* image = NewImage(1344, 1786);
    unsigned char* buffer =
        (unsigned char*)image->GetScalarPointer(0, 0, 0);
    int npixels = 1344 * 1786;
    for (int i = 0; i < npixels * 3; i++)
        buffer[i] = 0;

    std::vector<Triangle> readTriangles = GetTriangles(); //triangles read from file
    std::vector<Triangle> triangles; //triangles we're actually gonna use

    //loop though each triangle read from file and add two generated triangles to triangles
    for (int t = 0; t < readTriangles.size(); t++)
    {
        readTriangles[t].sortPointsByX(); //have the points in this triangle sorted by their X values

        double slopeBetweenExtremeXs = (readTriangles[t].Y[2] - readTriangles[t].Y[0]) / (readTriangles[t].X[2] - readTriangles[t].X[0]); //the slope between the points of highest and lowest x values

        //the new point has the x of our midrange existing point.
        //y of new pt = y of leftmost point + slope between the leftmost and rightmost points * num of steps between leftmost's x and new pt's x
        double newPoint[] = { readTriangles[t].X[1], readTriangles[t].Y[0] + slopeBetweenExtremeXs * (readTriangles[t].X[1] - readTriangles[t].X[0])};
        
        //putting together and pushing the two triangles before iterating through the loop again
        Triangle leftTriangle;
        leftTriangle.X[0] = readTriangles[t].X[0];
        leftTriangle.Y[0] = readTriangles[t].Y[0];
        leftTriangle.X[1] = newPoint[0];
        leftTriangle.Y[1] = newPoint[1];
        leftTriangle.X[2] = readTriangles[t].X[1];
        leftTriangle.Y[2] = readTriangles[t].Y[1];
        leftTriangle.color[0] = readTriangles[t].color[0];
        leftTriangle.color[1] = readTriangles[t].color[1];
        leftTriangle.color[2] = readTriangles[t].color[2];
        triangles.push_back(leftTriangle);

        Triangle rightTriangle;
        rightTriangle.X[0] = newPoint[0];
        rightTriangle.Y[0] = newPoint[1];
        rightTriangle.X[1] = readTriangles[t].X[1];
        rightTriangle.Y[1] = readTriangles[t].Y[1];
        rightTriangle.X[2] = readTriangles[t].X[2];
        rightTriangle.Y[2] = readTriangles[t].Y[2];
        rightTriangle.color[0] = readTriangles[t].color[0];
        rightTriangle.color[1] = readTriangles[t].color[1];
        rightTriangle.color[2] = readTriangles[t].color[2];
        triangles.push_back(rightTriangle);
    }

    Screen screen;
    screen.buffer = buffer;
    screen.width = 1786;
    screen.height = 1344;

    int* rowLimits = (int*)malloc(sizeof(int) * 2);

    //looping through our generated triangles
    for (int t = 0; t < triangles.size(); t++)
    {
        triangles[t].prepare(); //compute the slopes/other info necessary to get the row bounds for each scanline
        
        int initialScanPosition = ceil__441(triangles[t].X[0]); //starting value for the column scan

        int maxScanPosition = floor__441(triangles[t].X[2]); //ending value for column scan

        //each scanline (columns)
        for (int scanPosition = initialScanPosition; scanPosition < maxScanPosition + 1; scanPosition++)
        {
            //oob check
            if (scanPosition >= screen.width
                || scanPosition < 0)
                continue;

            triangles[t].getRowBounds(scanPosition, rowLimits); //get the min and max rows we need to add to the screen buffer

            //fill in the pixels by row (for each column scan)
            for (int currentRow = rowLimits[0];
                currentRow < rowLimits[1];
                currentRow++)
            {
                //oob check
                if (currentRow >= screen.height
                    || currentRow < 0)
                    continue;

                //set this pixel's color
                buffer[(int)(currentRow * screen.width * 3 + scanPosition * 3)] = triangles[t].color[0];
                buffer[(int)(currentRow * screen.width * 3 + scanPosition * 3 + 1)] = triangles[t].color[1];
                buffer[(int)(currentRow * screen.width * 3 + scanPosition * 3 + 2)] = triangles[t].color[2];

            }
        }
    }

    free(rowLimits);

    WriteImage(image, "goDucks");

	return 0;
}