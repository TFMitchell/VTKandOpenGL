/*
* CIS 441 Project 1D
* 
* I'll try my best to avoid using syntax Hank's compiler doesn't support.
* 
* The output file makes the same image as the example, but is a few KB bigger. This shouldn't be a problem.
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

class Triangle
{
public:
    double         X[3];
    double         Y[3];
    double         Z[3];
    double colors[3][3];

    //store these as member variables to reduce the number of calculations needed
    double bottomSlope, bottomRGBGradient[3], bottomZGradient, bottomStartPt[6]; //[6] is { X, Y, Z, R, G, B }
    double topSlope, topRGBGradient[3], topZGradient, topStartPt[6];

    void sortPointsByX() //this will sort X[] values (ascending) for this triangle and shuffle the Y[] as well so it matches
                         //using vectors and built-in sorting methods was much too slow for the volume of triangles we're processing, so I had to make this
    {
        double pointsAndColor[3][6]; //3 points of x, y, z, r, g, b
        double tmp[6];

        //load our X's and Y's into our working space
        for (int i = 0; i < 3; i++)
        {
            pointsAndColor[i][0] = X[i];
            pointsAndColor[i][1] = Y[i];
            pointsAndColor[i][2] = Z[i];
            pointsAndColor[i][3] = colors[i][0];
            pointsAndColor[i][4] = colors[i][1];
            pointsAndColor[i][5] = colors[i][2];
        }

        //sorting method
        for (int c = 0; c < 2; c++)
        {
            for (int i = 0; i < 2; i++)
            {
                if (pointsAndColor[i][0] > pointsAndColor[i + 1][0])
                {
                    tmp[0] = pointsAndColor[i][0];
                    tmp[1] = pointsAndColor[i][1];
                    tmp[2] = pointsAndColor[i][2];
                    tmp[3] = pointsAndColor[i][3];
                    tmp[4] = pointsAndColor[i][4];
                    tmp[5] = pointsAndColor[i][5];

                    pointsAndColor[i][0] = pointsAndColor[i + 1][0];
                    pointsAndColor[i][1] = pointsAndColor[i + 1][1];
                    pointsAndColor[i][2] = pointsAndColor[i + 1][2];
                    pointsAndColor[i][3] = pointsAndColor[i + 1][3];
                    pointsAndColor[i][4] = pointsAndColor[i + 1][4];
                    pointsAndColor[i][5] = pointsAndColor[i + 1][5];
                    pointsAndColor[i + 1][0] = tmp[0];
                    pointsAndColor[i + 1][1] = tmp[1];
                    pointsAndColor[i + 1][2] = tmp[2];
                    pointsAndColor[i + 1][3] = tmp[3];
                    pointsAndColor[i + 1][4] = tmp[4];
                    pointsAndColor[i + 1][5] = tmp[5];
                }
            }
        }

        //set the class member variables appropriately
        for (int i = 0; i < 3; i++)
        {
            X[i] = pointsAndColor[i][0];
            Y[i] = pointsAndColor[i][1];
            Z[i] = pointsAndColor[i][2];
            colors[i][0] = pointsAndColor[i][3];
            colors[i][1] = pointsAndColor[i][4];
            colors[i][2] = pointsAndColor[i][5];
        }
    }

    void prepare() //run this before using getRowBounds()
    {
        //determine the points with matching x values, as well as the point with a unique x

        if (X[0] != X[1] && X[0] != X[2]) //left-pointing triangle
        {
            bottomStartPt[0] = X[0];
            bottomStartPt[1] = Y[0];
            bottomStartPt[2] = Z[0];
            bottomStartPt[3] = colors[0][0];
            bottomStartPt[4] = colors[0][1];
            bottomStartPt[5] = colors[0][2];

            //the point on the right with a lower y value is needed
            int idxLowYPt = 1;
            if (Y[2] < Y[1])
                idxLowYPt = 2;

            bottomSlope = (Y[idxLowYPt] - Y[0]) / (X[idxLowYPt] - X[0]);
            bottomZGradient = (Z[idxLowYPt] - Z[0]) / (X[idxLowYPt] - X[0]);
            
            bottomRGBGradient[0] = (colors[idxLowYPt][0] - colors[0][0]) / (X[idxLowYPt] - X[0]);
            bottomRGBGradient[1] = (colors[idxLowYPt][1] - colors[0][1]) / (X[idxLowYPt] - X[0]);
            bottomRGBGradient[2] = (colors[idxLowYPt][2] - colors[0][2]) / (X[idxLowYPt] - X[0]);

            topStartPt[0] = X[0];
            topStartPt[1] = Y[0];
            topStartPt[2] = Z[0];
            topStartPt[3] = colors[0][0];
            topStartPt[4] = colors[0][1];
            topStartPt[5] = colors[0][2];

            //the point of the highest y
            int idxHighPt = 1;
            if (idxLowYPt == 1)
                idxHighPt = 2;

            topSlope = (Y[idxHighPt] - Y[0]) / (X[idxHighPt] - X[0]);
            topZGradient = (Z[idxHighPt] - Z[0]) / (X[idxHighPt] - X[0]);

            topRGBGradient[0] = (colors[idxHighPt][0] - colors[0][0]) / (X[idxHighPt] - X[0]);
            topRGBGradient[1] = (colors[idxHighPt][1] - colors[0][1]) / (X[idxHighPt] - X[0]);
            topRGBGradient[2] = (colors[idxHighPt][2] - colors[0][2]) / (X[idxHighPt] - X[0]);
        }

        else //(X[2] != X[0] && X[2] != X[1]) //right-pointing triangle
        {
            //the point on the left with a lower y value is needed
            int idxLowYPt = 0;
            if (Y[1] < Y[0])
                idxLowYPt = 1;

            bottomStartPt[0] = X[idxLowYPt];
            bottomStartPt[1] = Y[idxLowYPt];
            bottomStartPt[2] = Z[idxLowYPt];
            bottomStartPt[3] = colors[idxLowYPt][0];
            bottomStartPt[4] = colors[idxLowYPt][1];
            bottomStartPt[5] = colors[idxLowYPt][2];

            bottomSlope = (Y[2] - Y[idxLowYPt]) / (X[2] - X[idxLowYPt]);
            bottomZGradient = (Z[2] - Z[idxLowYPt]) / (X[2] - X[idxLowYPt]);

            bottomRGBGradient[0] = (colors[2][0] - colors[idxLowYPt][0]) / (X[2] - X[idxLowYPt]);
            bottomRGBGradient[1] = (colors[2][1] - colors[idxLowYPt][1]) / (X[2] - X[idxLowYPt]);
            bottomRGBGradient[2] = (colors[2][2] - colors[idxLowYPt][2]) / (X[2] - X[idxLowYPt]);


            int idxHighPt = 0;
            if (idxLowYPt == 0)
                idxHighPt = 1;

            topStartPt[0] = X[idxHighPt];
            topStartPt[1] = Y[idxHighPt];
            topStartPt[2] = Z[idxHighPt];
            topStartPt[3] = colors[idxHighPt][0];
            topStartPt[4] = colors[idxHighPt][1];
            topStartPt[5] = colors[idxHighPt][2];

            topSlope = (Y[2] - Y[idxHighPt]) / (X[2] - X[idxHighPt]);
            topZGradient = (Z[2] - Z[idxHighPt]) / (X[2] - X[idxHighPt]);

            topRGBGradient[0] = (colors[2][0] - colors[idxHighPt][0]) / (X[2] - X[idxHighPt]);
            topRGBGradient[1] = (colors[2][1] - colors[idxHighPt][1]) / (X[2] - X[idxHighPt]);
            topRGBGradient[2] = (colors[2][2] - colors[idxHighPt][2]) / (X[2] - X[idxHighPt]);
        }
    }


    void getRowBounds(int column, double* rowLimits, double** colorsAtLimits, double* zAtLimits) //for a given column to scan, get back the rows between which we need to fill in, the colors at those endpoints,
                                                                                                  //and the z values at the endpoints. Make sure to run the above funtion first
    {
        rowLimits[0] = bottomStartPt[1] + bottomSlope * (column - bottomStartPt[0]);
        rowLimits[1] = topStartPt[1] + topSlope * (column - topStartPt[0]);

        colorsAtLimits[0][0] = bottomStartPt[3] + bottomRGBGradient[0] * ((double)column - bottomStartPt[0]);
        colorsAtLimits[0][1] = bottomStartPt[4] + bottomRGBGradient[1] * ((double)column - bottomStartPt[0]);
        colorsAtLimits[0][2] = bottomStartPt[5] + bottomRGBGradient[2] * ((double)column - bottomStartPt[0]);

        colorsAtLimits[1][0] = topStartPt[3] + topRGBGradient[0] * ((double)column - topStartPt[0]);
        colorsAtLimits[1][1] = topStartPt[4] + topRGBGradient[1] * ((double)column - topStartPt[0]);
        colorsAtLimits[1][2] = topStartPt[5] + topRGBGradient[2] * ((double)column - topStartPt[0]);

        zAtLimits[0] = bottomStartPt[2] + bottomZGradient * ((double)column - bottomStartPt[0]);
        zAtLimits[1] = topStartPt[2] + topZGradient * ((double)column - topStartPt[0]);
    }
};

class Screen
{
public:
    unsigned char* buffer;
    int width, height;
};


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

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader* rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1d_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData* pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints* pts = pd->GetPoints();
    vtkCellArray* cells = pd->GetPolys();
    vtkFloatArray* var = (vtkFloatArray*)pd->GetPointData()->GetArray("hardyglobal");
    float* color_ptr = var->GetPointer(0);
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
        tris[idx].Z[0] = pts->GetPoint(ptIds[0])[2];
        tris[idx].Z[1] = pts->GetPoint(ptIds[1])[2];
        tris[idx].Z[2] = pts->GetPoint(ptIds[2])[2];
        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 },
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 }
        };
        for (int j = 0; j < 3; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0; r < 7; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val - mins[r]) / (maxs[r] - mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0] + proportion * (RGB[r + 1][0] - RGB[r][0])) / 255.0;
            tris[idx].colors[j][1] = (RGB[r][1] + proportion * (RGB[r + 1][1] - RGB[r][1])) / 255.0;
            tris[idx].colors[j][2] = (RGB[r][2] + proportion * (RGB[r + 1][2] - RGB[r][2])) / 255.0;
        }
    }
    return tris;
}


int main()
{
    vtkImageData* image = NewImage(1000, 1000);

    unsigned char* buffer = (unsigned char*)image->GetScalarPointer(0, 0, 0);
    for (int i = 0; i < 1000 * 1000 * 3; i++)
        buffer[i] = 0;
    
    double *zBuffer = (double*) malloc (sizeof(double) * 1000 * 1000);
    for (int i = 0; i < 1000 * 1000; i++)
        zBuffer[i] = -1.f;
       
    std::vector<Triangle> readTriangles = GetTriangles(); //triangles read from file
    std::vector<Triangle> triangles; //triangles we're actually gonna use

    //loop though each triangle read from file and add two generated triangles to triangles
    for (int t = 0; t < readTriangles.size(); t++)
    {
        readTriangles[t].sortPointsByX(); //have the points in this triangle sorted by their X values

        double slopeBetweenExtremeXs = (readTriangles[t].Y[2] - readTriangles[t].Y[0]) / (readTriangles[t].X[2] - readTriangles[t].X[0]); //the slope between the points of highest and lowest x values
        double zGradientBetweenExtremeXs = (readTriangles[t].Z[2] - readTriangles[t].Z[0]) / (readTriangles[t].X[2] - readTriangles[t].X[0]);
        double colorGradientBetweenExtremeXs[3] =
        {
            (readTriangles[t].colors[2][0] - readTriangles[t].colors[0][0]) / (readTriangles[t].X[2] - readTriangles[t].X[0]),
            (readTriangles[t].colors[2][1] - readTriangles[t].colors[0][1]) / (readTriangles[t].X[2] - readTriangles[t].X[0]),
            (readTriangles[t].colors[2][2] - readTriangles[t].colors[0][2]) / (readTriangles[t].X[2] - readTriangles[t].X[0])
        };

        //the new point has the x of our midrange existing point.
        //y of new pt = y of leftmost point + slope between the leftmost and rightmost points * num of steps between leftmost's x and new pt's x
        double newPoint[] =
        {
            readTriangles[t].X[1],
            readTriangles[t].Y[0] + slopeBetweenExtremeXs * (readTriangles[t].X[1] - readTriangles[t].X[0]),
            readTriangles[t].Z[0] + zGradientBetweenExtremeXs * (readTriangles[t].X[1] - readTriangles[t].X[0]),
            readTriangles[t].colors[0][0] + colorGradientBetweenExtremeXs[0] * (readTriangles[t].X[1] - readTriangles[t].X[0]),
            readTriangles[t].colors[0][1] + colorGradientBetweenExtremeXs[1] * (readTriangles[t].X[1] - readTriangles[t].X[0]),
            readTriangles[t].colors[0][2] + colorGradientBetweenExtremeXs[2] * (readTriangles[t].X[1] - readTriangles[t].X[0])
         };


        //putting together and pushing the two triangles before iterating through the loop again
        Triangle leftTriangle;
        leftTriangle.X[0] = readTriangles[t].X[0];
        leftTriangle.Y[0] = readTriangles[t].Y[0];
        leftTriangle.Z[0] = readTriangles[t].Z[0];
        leftTriangle.colors[0][0] = readTriangles[t].colors[0][0];
        leftTriangle.colors[0][1] = readTriangles[t].colors[0][1];
        leftTriangle.colors[0][2] = readTriangles[t].colors[0][2];

        leftTriangle.X[1] = readTriangles[t].X[1];
        leftTriangle.Y[1] = readTriangles[t].Y[1];
        leftTriangle.Z[1] = readTriangles[t].Z[1];
        leftTriangle.colors[1][0] = readTriangles[t].colors[1][0];
        leftTriangle.colors[1][1] = readTriangles[t].colors[1][1];
        leftTriangle.colors[1][2] = readTriangles[t].colors[1][2];

        leftTriangle.X[2] = newPoint[0];
        leftTriangle.Y[2] = newPoint[1];
        leftTriangle.Z[2] = newPoint[2];
        leftTriangle.colors[2][0] = newPoint[3];
        leftTriangle.colors[2][1] = newPoint[4];
        leftTriangle.colors[2][2] = newPoint[5];

        triangles.push_back(leftTriangle);

        Triangle rightTriangle;
        rightTriangle.X[0] = readTriangles[t].X[1];
        rightTriangle.Y[0] = readTriangles[t].Y[1];
        rightTriangle.Z[0] = readTriangles[t].Z[1];
        rightTriangle.colors[0][0] = readTriangles[t].colors[1][0];
        rightTriangle.colors[0][1] = readTriangles[t].colors[1][1];
        rightTriangle.colors[0][2] = readTriangles[t].colors[1][2];

        rightTriangle.X[1] = newPoint[0];
        rightTriangle.Y[1] = newPoint[1];
        rightTriangle.Z[1] = newPoint[2];
        rightTriangle.colors[1][0] = newPoint[3];
        rightTriangle.colors[1][1] = newPoint[4];
        rightTriangle.colors[1][2] = newPoint[5];

        rightTriangle.X[2] = readTriangles[t].X[2];
        rightTriangle.Y[2] = readTriangles[t].Y[2];
        rightTriangle.Z[2] = readTriangles[t].Z[2];
        rightTriangle.colors[2][0] = readTriangles[t].colors[2][0];
        rightTriangle.colors[2][1] = readTriangles[t].colors[2][1];
        rightTriangle.colors[2][2] = readTriangles[t].colors[2][2];

        triangles.push_back(rightTriangle);

    }

    Screen screen;
    screen.buffer = buffer;
    screen.width = 1000;
    screen.height = 1000;

    double* rowLimits = (double*)malloc(sizeof(double) * 2);
    double* zAtLimits = (double*)malloc(sizeof(double) * 2);
    double** colorsAtLimits = (double**)malloc(sizeof(double*) * 2); //2x3 array
    for (int i = 0; i < 2; i++)
        colorsAtLimits[i] = (double*)malloc(sizeof(double) * 3);

    //looping through our generated triangles
    for (int t = 0; t < triangles.size(); t++)
    {
        triangles[t].prepare(); //compute the slopes/other info necessary to get the row bounds for each scanline

        //each scanline (columns)
        for (int scanPosition = ceil__441(triangles[t].X[0]); scanPosition <= floor__441(triangles[t].X[2]); scanPosition++)
        {
            //oob check
            if (scanPosition >= screen.width
                || scanPosition < 0)
                continue;

            triangles[t].getRowBounds(scanPosition, rowLimits, colorsAtLimits, zAtLimits); //get the min and max rows we need to add to the screen buffer

            //fill in the pixels by row (for each column scan)
            for (int currentRow = ceil__441(rowLimits[0]);
                currentRow <= floor__441(rowLimits[1]);
                currentRow++)
            {
                //oob check
                if (currentRow >= screen.height
                    || currentRow < 0)
                    continue;

                //stop if this point is hidden
                double zValue = zAtLimits[0] + (zAtLimits[1] - zAtLimits[0]) * ((double)currentRow - rowLimits[0]) / ((double)rowLimits[1] - rowLimits[0]);
                if (zValue < zBuffer[currentRow * screen.width + scanPosition])
                    continue;

                zBuffer[currentRow * screen.width + scanPosition] = zValue; //set z buffer

                //set this pixel's color
                buffer[currentRow * screen.width * 3 + scanPosition * 3 + 0] = ceil__441((colorsAtLimits[0][0] + (colorsAtLimits[1][0] - colorsAtLimits[0][0]) * ((double)currentRow - rowLimits[0]) / (rowLimits[1] - rowLimits[0])) * 255.f);
                buffer[currentRow * screen.width * 3 + scanPosition * 3 + 1] = ceil__441((colorsAtLimits[0][1] + (colorsAtLimits[1][1] - colorsAtLimits[0][1]) * ((double)currentRow - rowLimits[0]) / (rowLimits[1] - rowLimits[0])) * 255.f);
                buffer[currentRow * screen.width * 3 + scanPosition * 3 + 2] = ceil__441((colorsAtLimits[0][2] + (colorsAtLimits[1][2] - colorsAtLimits[0][2]) * ((double)currentRow - rowLimits[0]) / (rowLimits[1] - rowLimits[0])) * 255.f);           
            }
        }
    }

    free(rowLimits);
    free(zAtLimits);
    for (int i = 0; i < 2; i++)
        free(colorsAtLimits[i]);
    free(colorsAtLimits);

    WriteImage(image, "1D");

    return 0;
}