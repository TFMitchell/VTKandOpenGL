/*
* CIS 441 Project 1F
* 
* I'll try my best to avoid using syntax Hank's compiler doesn't support.
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
#include <vtkDoubleArray.h>
#define M_PI 3.14159265358979323846 //this was necessary for my program to compile on a Windows machine. You might need to comment it out
#define NORMALS

struct LightingParameters
{
    LightingParameters(void)
    {
        lightDir[0] = -0.6;
        lightDir[1] = 0;
        lightDir[2] = -0.8;
        Ka = 0.3;
        Kd = 0.7;
        Ks = 2.8;
        alpha = 50.5;
    };


    double lightDir[3]; // The direction of the light source
    double Ka;          // The coefficient for ambient lighting
    double Kd;          // The coefficient for diffuse lighting
    double Ks;          // The coefficient for specular lighting
    double alpha;       // The exponent term for specular lighting


};

class Matrix
{
public:
    double          A[4][4];  // A[i][j] means row i, column j

    void            TransformPoint(const double* ptIn, double* ptOut);
    static Matrix   ComposeMatrices(const Matrix&, const Matrix&);
    void            Print(ostream& o);
};

static void normalize(double* myVector) //normalizes a vector in the form of a double with three elements
{
    double magnitude = sqrt(pow(myVector[0], 2) + pow(myVector[1], 2) + pow(myVector[2], 2));

    for (int i = 0; i < 3; i++)
        myVector[i] = myVector[i] / magnitude;
}


double CalculateDiffuseShading(double* normal, LightingParameters lp)
{

    return (lp.lightDir[0] * normal[0] + lp.lightDir[1] * normal[1] + lp.lightDir[2] * normal[2]) * lp.Kd;
}



class Camera
{
public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    Matrix ViewTransform(void) //get the viewTransform matrix
    { 
        Matrix rv;

        rv.A[0][0] = 1.f / tan(angle / 2);
        rv.A[0][1] = 0.f;
        rv.A[0][2] = 0.f;
        rv.A[0][3] = 0.f;

        rv.A[1][0] = 0.f;
        rv.A[1][1] = 1.f / tan(angle / 2);
        rv.A[1][2] = 0.f;
        rv.A[1][3] = 0.f;

        rv.A[2][0] = 0.f;
        rv.A[2][1] = 0.f;
        rv.A[2][2] = (far + near) / (far - near);
        rv.A[2][3] = -1.f;

        rv.A[3][0] = 0.f;
        rv.A[3][1] = 0.f;
        rv.A[3][2] = 2.f * far * near / (far - near);
        rv.A[3][3] = 0.f;

        return rv;
    };

    

    Matrix CameraTransform(void) //get the cameraTransform matrix
    { 
        //we use (position - focus) a lot, so declaring it here makes the code easier to understand
        double posMinusFocus[3];
        for (int i = 0; i < 3; i++)
            posMinusFocus[i] = position[i] - focus[i];

        //multiply up by position - focus
        double U[] =
        {
            up[1] * posMinusFocus[2] - up[2] * posMinusFocus[1],
            up[2] * posMinusFocus[0] - up[0] * posMinusFocus[2],
            up[0] * posMinusFocus[1] - up[1] * posMinusFocus[0]
        };
        normalize(U);

        //multiply position - focus by U
        double V[] =
        {
            posMinusFocus[1] * U[2] - posMinusFocus[2] * U[1],
            posMinusFocus[2] * U[0] - posMinusFocus[0] * U[2],
            posMinusFocus[0] * U[1] - posMinusFocus[1] * U[0]
        };
        normalize(V);

        //just position - focus normalized
        double W[] =
        {
            posMinusFocus[0],
            posMinusFocus[1],
            posMinusFocus[2]
        };
        normalize(W);

        double t[] =
        {
            0.f - position[0],
            0.f - position[1],
            0.f - position[2]
        };

        double dotProducts[] =
        {
            U[0] * t[0] + U[1] * t[1] + U[2] * t[2],
            V[0] * t[0] + V[1] * t[1] + V[2] * t[2],
            W[0] * t[0] + W[1] * t[1] + W[2] * t[2]
        };

        Matrix rv;

        rv.A[0][0] = U[0];
        rv.A[1][0] = U[1]; 
        rv.A[2][0] = U[2];
        rv.A[3][0] = dotProducts[0];

        rv.A[0][1] = V[0];
        rv.A[1][1] = V[1];
        rv.A[2][1] = V[2];
        rv.A[3][1] = dotProducts[1];

        rv.A[0][2] = W[0];
        rv.A[1][2] = W[1];
        rv.A[2][2] = W[2];
        rv.A[3][2] = dotProducts[2];

        rv.A[0][3] = 0.f;
        rv.A[1][3] = 0.f;
        rv.A[2][3] = 0.f;
        rv.A[3][3] = 1.f;

        return rv;
    };
};


double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames - 2 * ramp;
    double height = 1. / (nNonRamp + 4 * ramp / M_PI);
    if (curFrame < ramp)
    {
        double factor = 2 * height * ramp / M_PI;
        double eval = cos(M_PI / 2 * ((double)curFrame) / ramp);
        return (1. - eval) * factor;
    }
    else if (curFrame > nFrames - ramp)
    {
        int amount_left = nFrames - curFrame;
        double factor = 2 * height * ramp / M_PI;
        double eval = cos(M_PI / 2 * ((double)amount_left / ramp));
        return 1. - (1 - eval) * factor;
    }
    double amount_in_quad = ((double)curFrame - ramp);
    double quad_part = amount_in_quad * height;
    double curve_part = height * (2 * ramp) / M_PI;
    return quad_part + curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes / 10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI / 6;
    c.position[0] = 40 * sin(2 * M_PI * t);
    c.position[1] = 40 * cos(2 * M_PI * t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;

    return c;
}

LightingParameters
GetLighting(Camera c)
{
    LightingParameters lp;
    lp.lightDir[0] = c.position[0] - c.focus[0];
    lp.lightDir[1] = c.position[1] - c.focus[1];
    lp.lightDir[2] = c.position[2] - c.focus[2];
    double mag = sqrt(lp.lightDir[0] * lp.lightDir[0]
        + lp.lightDir[1] * lp.lightDir[1]
        + lp.lightDir[2] * lp.lightDir[2]);
    if (mag > 0)
    {
        lp.lightDir[0] /= mag;
        lp.lightDir[1] /= mag;
        lp.lightDir[2] /= mag;
    }

    return lp;
}


void
Matrix::Print(ostream& o)
{
    for (int i = 0; i < 4; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix& M1, const Matrix& M2)
{
    Matrix rv;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0; k < 4; k++)
                rv.A[i][j] += M1.A[i][k] * M2.A[k][j];
        }

    return rv;
}

class Triangle
{
public:
    double         X[3];
    double         Y[3];
    double         Z[3];
    double colors[3][3];
    double normals[3][3];

    //store these as member variables to reduce the number of calculations needed
    double bottomSlope, bottomRGBGradient[3], bottomDiffuseGradient, bottomSpecularGradient, bottomZGradient, bottomStartPt[8]; //[8] is { X, Y, Z, R, G, B, diffuse, specular }
    double topSlope, topRGBGradient[3], topDiffuseGradient, topSpecularGradient, topZGradient, topStartPt[8];

    void transform(Matrix transformation)
    {
        for (int i = 0; i < 3; i++)
        {
            double newX, newY, newZ;
            newX = transformation.A[0][0] * X[i] + transformation.A[1][0] * Y[i] + transformation.A[2][0] * Z[i] + transformation.A[3][0];
            newY = transformation.A[0][1] * X[i] + transformation.A[1][1] * Y[i] + transformation.A[2][1] * Z[i] + transformation.A[3][1];
            newZ = transformation.A[0][2] * X[i] + transformation.A[1][2] * Y[i] + transformation.A[2][2] * Z[i] + transformation.A[3][2];
            double W = transformation.A[0][3] * X[i] + transformation.A[1][3] * Y[i] + transformation.A[2][3] * Z[i] + transformation.A[3][3];

            X[i] = newX / W;
            Y[i] = newY / W;
            Z[i] = newZ / W;
        }     
    }

    void sortPointsByX() //this will sort X[] values (ascending) for this triangle and shuffle the Y[] as well so it matches
                         //using vectors and built-in sorting methods was much too slow for the volume of triangles we're processing, so I had to make this
    {
        double pointsAndColor[3][9]; //3 points of x, y, z, r, g, b, normals
        double tmp[9];

        //load our X's and Y's into our working space
        for (int i = 0; i < 3; i++)
        {
            pointsAndColor[i][0] = X[i];
            pointsAndColor[i][1] = Y[i];
            pointsAndColor[i][2] = Z[i];
            pointsAndColor[i][3] = colors[i][0];
            pointsAndColor[i][4] = colors[i][1];
            pointsAndColor[i][5] = colors[i][2];
            pointsAndColor[i][6] = normals[i][0];
            pointsAndColor[i][7] = normals[i][1];
            pointsAndColor[i][8] = normals[i][2];
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
                    tmp[6] = pointsAndColor[i][6];
                    tmp[7] = pointsAndColor[i][7];
                    tmp[8] = pointsAndColor[i][8];

                    pointsAndColor[i][0] = pointsAndColor[i + 1][0];
                    pointsAndColor[i][1] = pointsAndColor[i + 1][1];
                    pointsAndColor[i][2] = pointsAndColor[i + 1][2];
                    pointsAndColor[i][3] = pointsAndColor[i + 1][3];
                    pointsAndColor[i][4] = pointsAndColor[i + 1][4];
                    pointsAndColor[i][5] = pointsAndColor[i + 1][5];
                    pointsAndColor[i][6] = pointsAndColor[i + 1][6];
                    pointsAndColor[i][7] = pointsAndColor[i + 1][7];
                    pointsAndColor[i][8] = pointsAndColor[i + 1][8];
                    pointsAndColor[i + 1][0] = tmp[0];
                    pointsAndColor[i + 1][1] = tmp[1];
                    pointsAndColor[i + 1][2] = tmp[2];
                    pointsAndColor[i + 1][3] = tmp[3];
                    pointsAndColor[i + 1][4] = tmp[4];
                    pointsAndColor[i + 1][5] = tmp[5];
                    pointsAndColor[i + 1][6] = tmp[6];
                    pointsAndColor[i + 1][7] = tmp[7];
                    pointsAndColor[i + 1][8] = tmp[8];
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
            normals[i][0] = pointsAndColor[i][6];
            normals[i][1] = pointsAndColor[i][7];
            normals[i][2] = pointsAndColor[i][8];
        }
    }

    void prepare(LightingParameters lp) //run this before using getRowBounds()
    {
        //determine the points with matching x values, as well as the point with a unique x

        if (X[0] != X[1] && X[0] != X[2]) //left-pointing triangle
        {
            //the point on the right with a lower y value is needed
            int idxLowYPt = 1;
            if (Y[2] < Y[1])
                idxLowYPt = 2;

            bottomStartPt[0] = X[0];
            bottomStartPt[1] = Y[0];
            bottomStartPt[2] = Z[0];
            bottomStartPt[3] = colors[0][0];
            bottomStartPt[4] = colors[0][1];
            bottomStartPt[5] = colors[0][2];
            bottomStartPt[6] = CalculateDiffuseShading(normals[0], lp);
            bottomStartPt[7] = 0.f;


            bottomSlope = (Y[idxLowYPt] - Y[0]) / (X[idxLowYPt] - X[0]);
            bottomZGradient = (Z[idxLowYPt] - Z[0]) / (X[idxLowYPt] - X[0]);

            bottomRGBGradient[0] = (colors[idxLowYPt][0] - colors[0][0]) / (X[idxLowYPt] - X[0]);
            bottomRGBGradient[1] = (colors[idxLowYPt][1] - colors[0][1]) / (X[idxLowYPt] - X[0]);
            bottomRGBGradient[2] = (colors[idxLowYPt][2] - colors[0][2]) / (X[idxLowYPt] - X[0]);

            bottomDiffuseGradient = (CalculateDiffuseShading(normals[idxLowYPt], lp) - bottomStartPt[6]) / (X[idxLowYPt] - X[0]);

            //the point of the highest y
            int idxHighPt = 1;
            if (idxLowYPt == 1)
                idxHighPt = 2;

            topStartPt[0] = X[0];
            topStartPt[1] = Y[0];
            topStartPt[2] = Z[0];
            topStartPt[3] = colors[0][0];
            topStartPt[4] = colors[0][1];
            topStartPt[5] = colors[0][2];
            topStartPt[6] = CalculateDiffuseShading(normals[0], lp);
            topStartPt[7] = 0.f;


            topSlope = (Y[idxHighPt] - Y[0]) / (X[idxHighPt] - X[0]);
            topZGradient = (Z[idxHighPt] - Z[0]) / (X[idxHighPt] - X[0]);

            topRGBGradient[0] = (colors[idxHighPt][0] - colors[0][0]) / (X[idxHighPt] - X[0]);
            topRGBGradient[1] = (colors[idxHighPt][1] - colors[0][1]) / (X[idxHighPt] - X[0]);
            topRGBGradient[2] = (colors[idxHighPt][2] - colors[0][2]) / (X[idxHighPt] - X[0]);

            topDiffuseGradient = (CalculateDiffuseShading(normals[idxHighPt], lp) - topStartPt[6]) / (X[idxHighPt] - X[0]);
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
            bottomStartPt[6] = CalculateDiffuseShading(normals[idxLowYPt], lp);
            bottomStartPt[7] = 0.f;

            bottomSlope = (Y[2] - Y[idxLowYPt]) / (X[2] - X[idxLowYPt]);
            bottomZGradient = (Z[2] - Z[idxLowYPt]) / (X[2] - X[idxLowYPt]);

            bottomRGBGradient[0] = (colors[2][0] - colors[idxLowYPt][0]) / (X[2] - X[idxLowYPt]);
            bottomRGBGradient[1] = (colors[2][1] - colors[idxLowYPt][1]) / (X[2] - X[idxLowYPt]);
            bottomRGBGradient[2] = (colors[2][2] - colors[idxLowYPt][2]) / (X[2] - X[idxLowYPt]);

            bottomDiffuseGradient = (CalculateDiffuseShading(normals[2], lp) - bottomStartPt[6]) / (X[2] - X[idxLowYPt]);

            int idxHighPt = 0;
            if (idxLowYPt == 0)
                idxHighPt = 1;

            topStartPt[0] = X[idxHighPt];
            topStartPt[1] = Y[idxHighPt];
            topStartPt[2] = Z[idxHighPt];
            topStartPt[3] = colors[idxHighPt][0];
            topStartPt[4] = colors[idxHighPt][1];
            topStartPt[5] = colors[idxHighPt][2];
            topStartPt[6] = CalculateDiffuseShading(normals[idxHighPt], lp);
            topStartPt[7] = 0.f;

            topSlope = (Y[2] - Y[idxHighPt]) / (X[2] - X[idxHighPt]);
            topZGradient = (Z[2] - Z[idxHighPt]) / (X[2] - X[idxHighPt]);

            topRGBGradient[0] = (colors[2][0] - colors[idxHighPt][0]) / (X[2] - X[idxHighPt]);
            topRGBGradient[1] = (colors[2][1] - colors[idxHighPt][1]) / (X[2] - X[idxHighPt]);
            topRGBGradient[2] = (colors[2][2] - colors[idxHighPt][2]) / (X[2] - X[idxHighPt]);

            topDiffuseGradient = (CalculateDiffuseShading(normals[2], lp) - topStartPt[6]) / (X[2] - X[idxHighPt]);
        }
    }

    void getRowBounds(double column, double* rowLimits, double** colorsAtLimits, double* diffuseShaderAtLimits, double* zAtLimits) //for a given column to scan, get back the rows between which we need to fill in, the colors at those endpoints,
                                                                                                  //and the z values at the endpoints. Make sure to run the above funtion first
    {
        rowLimits[0] = bottomStartPt[1] + bottomSlope * (column - bottomStartPt[0]);
        rowLimits[1] = topStartPt[1] + topSlope * (column - topStartPt[0]);

        colorsAtLimits[0][0] = bottomStartPt[3] + bottomRGBGradient[0] * (column - bottomStartPt[0]);
        colorsAtLimits[0][1] = bottomStartPt[4] + bottomRGBGradient[1] * (column - bottomStartPt[0]);
        colorsAtLimits[0][2] = bottomStartPt[5] + bottomRGBGradient[2] * (column - bottomStartPt[0]);

        colorsAtLimits[1][0] = topStartPt[3] + topRGBGradient[0] * (column - topStartPt[0]);
        colorsAtLimits[1][1] = topStartPt[4] + topRGBGradient[1] * (column - topStartPt[0]);
        colorsAtLimits[1][2] = topStartPt[5] + topRGBGradient[2] * (column - topStartPt[0]);

        diffuseShaderAtLimits[0] = bottomStartPt[6] + bottomDiffuseGradient * (column - bottomStartPt[0]);

        diffuseShaderAtLimits[1] = topStartPt[6] + topDiffuseGradient * (column - topStartPt[0]);

        zAtLimits[0] = bottomStartPt[2] + bottomZGradient * (column - bottomStartPt[0]);
        zAtLimits[1] = topStartPt[2] + topZGradient * (column - topStartPt[0]);
    }
};


std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader* rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1f_geometry.vtk");
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
    vtkDoubleArray* var = (vtkDoubleArray*)pd->GetPointData()->GetArray("hardyglobal");
    double* color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray* n = (vtkFloatArray*)pd->GetPointData()->GetNormals();
    float* normals = n->GetPointer(0);
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
        double* pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3 * ptIds[0] + 0];
        tris[idx].normals[0][1] = normals[3 * ptIds[0] + 1];
        tris[idx].normals[0][2] = normals[3 * ptIds[0] + 2];
#endif
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS
        tris[idx].normals[1][0] = normals[3 * ptIds[1] + 0];
        tris[idx].normals[1][1] = normals[3 * ptIds[1] + 1];
        tris[idx].normals[1][2] = normals[3 * ptIds[1] + 2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3 * ptIds[2] + 0];
        tris[idx].normals[2][1] = normals[3 * ptIds[2] + 1];
        tris[idx].normals[2][2] = normals[3 * ptIds[2] + 2];
#endif

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


double ceil__441(double f)
{
    return ceil(f - 0.00001);
}

double floor__441(double f)
{
    return floor(f + 0.00001);
}


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


int main()
{
    //set up the image and screen. We'll fill with the 0's/-1's later
    vtkImageData* image = NewImage(1000, 1000);
    unsigned char* buffer = (unsigned char*)image->GetScalarPointer(0, 0, 0); 
    double *zBuffer = (double*) malloc (sizeof(double) * 1000 * 1000);

    Screen screen;
    screen.buffer = buffer;
    screen.width = 1000;
    screen.height = 1000;

    //allocating some arrays we'll need later
    double* rowLimits = (double*)malloc(sizeof(double) * 2);
    double* zAtLimits = (double*)malloc(sizeof(double) * 2);
    double** colorsAtLimits = (double**)malloc(sizeof(double*) * 2); //2x3 array
    for (int i = 0; i < 2; i++)
        colorsAtLimits[i] = (double*)malloc(sizeof(double) * 3);
    double* diffuseShaderAtLimits = (double*)malloc(sizeof(double) * 2); 



    //for each of the camera's four perspectives
    for (int f = 0; f < 4; f++)
    {
        std::vector<Triangle> readTriangles = GetTriangles(); //triangles read from file. Re-reading adds extra io, so this could be improved
        std::vector<Triangle> triangles; //triangles we're actually gonna use. 

        Camera cam = GetCamera(f * 250, 1000);

        LightingParameters lp = GetLighting(cam);

        Matrix camTrans = cam.CameraTransform();
        Matrix viewTrans = cam.ViewTransform();
        Matrix viewSpace = Matrix::ComposeMatrices(camTrans, viewTrans); //multiply the three matrices

        //loop though each triangle read from file and add two generated triangles to triangles
        for (int t = 0; t < readTriangles.size(); t++)
        {
            readTriangles[t].transform(viewSpace); //adjust the points for this triangle according to the total matrix

            readTriangles[t].sortPointsByX(); //have the points in this triangle sorted by their X values

            double slopeBetweenExtremeXs = (readTriangles[t].Y[2] - readTriangles[t].Y[0]) / (readTriangles[t].X[2] - readTriangles[t].X[0]); //the slope between the points of highest and lowest x values
            double zGradientBetweenExtremeXs = (readTriangles[t].Z[2] - readTriangles[t].Z[0]) / (readTriangles[t].X[2] - readTriangles[t].X[0]);
   
            double colorGradientBetweenExtremeXs[] =
            {
                (readTriangles[t].colors[2][0] - readTriangles[t].colors[0][0]) / (readTriangles[t].X[2] - readTriangles[t].X[0]),
                (readTriangles[t].colors[2][1] - readTriangles[t].colors[0][1]) / (readTriangles[t].X[2] - readTriangles[t].X[0]),
                (readTriangles[t].colors[2][2] - readTriangles[t].colors[0][2]) / (readTriangles[t].X[2] - readTriangles[t].X[0])
            };

            double normalGradientBetweenExtremeXs[] =
            {
                (readTriangles[t].normals[2][0] - readTriangles[t].normals[0][0]) / (readTriangles[t].X[2] - readTriangles[t].X[0]),
                (readTriangles[t].normals[2][1] - readTriangles[t].normals[0][1]) / (readTriangles[t].X[2] - readTriangles[t].X[0]),
                (readTriangles[t].normals[2][2] - readTriangles[t].normals[0][2]) / (readTriangles[t].X[2] - readTriangles[t].X[0])
            };

            //the new point has the x of our midrange existing point.
            //y of new pt = y of leftmost point + slope between the leftmost and rightmost points * num of steps between leftmost's x and new pt's x
            double newPoint[] = //form x, y, z, r, g, b, normals for each three colors
            {
                readTriangles[t].X[1],
                readTriangles[t].Y[0] + slopeBetweenExtremeXs * (readTriangles[t].X[1] - readTriangles[t].X[0]),
                readTriangles[t].Z[0] + zGradientBetweenExtremeXs * (readTriangles[t].X[1] - readTriangles[t].X[0]),
                readTriangles[t].colors[0][0] + colorGradientBetweenExtremeXs[0] * (readTriangles[t].X[1] - readTriangles[t].X[0]),
                readTriangles[t].colors[0][1] + colorGradientBetweenExtremeXs[1] * (readTriangles[t].X[1] - readTriangles[t].X[0]),
                readTriangles[t].colors[0][2] + colorGradientBetweenExtremeXs[2] * (readTriangles[t].X[1] - readTriangles[t].X[0]),
                readTriangles[t].normals[0][0] + normalGradientBetweenExtremeXs[0] * (readTriangles[t].X[1] - readTriangles[t].X[0]),
                readTriangles[t].normals[0][1] + normalGradientBetweenExtremeXs[1] * (readTriangles[t].X[1] - readTriangles[t].X[0]),
                readTriangles[t].normals[0][2] + normalGradientBetweenExtremeXs[2] * (readTriangles[t].X[1] - readTriangles[t].X[0])
            };


            //putting together and pushing the two triangles before iterating through the loop again
            Triangle leftTriangle;
            leftTriangle.X[0] = readTriangles[t].X[0];
            leftTriangle.Y[0] = readTriangles[t].Y[0];
            leftTriangle.Z[0] = readTriangles[t].Z[0];
            leftTriangle.colors[0][0] = readTriangles[t].colors[0][0];
            leftTriangle.colors[0][1] = readTriangles[t].colors[0][1];
            leftTriangle.colors[0][2] = readTriangles[t].colors[0][2];
            leftTriangle.normals[0][0] = readTriangles[t].normals[0][0];
            leftTriangle.normals[0][1] = readTriangles[t].normals[0][1];
            leftTriangle.normals[0][2] = readTriangles[t].normals[0][2];

            leftTriangle.X[1] = readTriangles[t].X[1];
            leftTriangle.Y[1] = readTriangles[t].Y[1];
            leftTriangle.Z[1] = readTriangles[t].Z[1];
            leftTriangle.colors[1][0] = readTriangles[t].colors[1][0];
            leftTriangle.colors[1][1] = readTriangles[t].colors[1][1];
            leftTriangle.colors[1][2] = readTriangles[t].colors[1][2];
            leftTriangle.normals[1][0] = readTriangles[t].normals[1][0];
            leftTriangle.normals[1][1] = readTriangles[t].normals[1][1];
            leftTriangle.normals[1][2] = readTriangles[t].normals[1][2];

            leftTriangle.X[2] = newPoint[0];
            leftTriangle.Y[2] = newPoint[1];
            leftTriangle.Z[2] = newPoint[2];
            leftTriangle.colors[2][0] = newPoint[3];
            leftTriangle.colors[2][1] = newPoint[4];
            leftTriangle.colors[2][2] = newPoint[5];
            leftTriangle.normals[2][0] = newPoint[6];
            leftTriangle.normals[2][1] = newPoint[7];
            leftTriangle.normals[2][2] = newPoint[8];

            triangles.push_back(leftTriangle);

            Triangle rightTriangle;
            rightTriangle.X[0] = readTriangles[t].X[1];
            rightTriangle.Y[0] = readTriangles[t].Y[1];
            rightTriangle.Z[0] = readTriangles[t].Z[1];
            rightTriangle.colors[0][0] = readTriangles[t].colors[1][0];
            rightTriangle.colors[0][1] = readTriangles[t].colors[1][1];
            rightTriangle.colors[0][2] = readTriangles[t].colors[1][2];
            rightTriangle.normals[0][0] = readTriangles[t].normals[1][0];
            rightTriangle.normals[0][1] = readTriangles[t].normals[1][1];
            rightTriangle.normals[0][2] = readTriangles[t].normals[1][2];

            rightTriangle.X[1] = newPoint[0];
            rightTriangle.Y[1] = newPoint[1];
            rightTriangle.Z[1] = newPoint[2];
            rightTriangle.colors[1][0] = newPoint[3];
            rightTriangle.colors[1][1] = newPoint[4];
            rightTriangle.colors[1][2] = newPoint[5];
            rightTriangle.normals[1][0] = newPoint[6];
            rightTriangle.normals[1][1] = newPoint[7];
            rightTriangle.normals[1][2] = newPoint[8];

            rightTriangle.X[2] = readTriangles[t].X[2];
            rightTriangle.Y[2] = readTriangles[t].Y[2];
            rightTriangle.Z[2] = readTriangles[t].Z[2];
            rightTriangle.colors[2][0] = readTriangles[t].colors[2][0];
            rightTriangle.colors[2][1] = readTriangles[t].colors[2][1];
            rightTriangle.colors[2][2] = readTriangles[t].colors[2][2];
            rightTriangle.normals[2][0] = readTriangles[t].normals[2][0];
            rightTriangle.normals[2][1] = readTriangles[t].normals[2][1];
            rightTriangle.normals[2][2] = readTriangles[t].normals[2][2];

            triangles.push_back(rightTriangle);

        }

        //"zeroing" our buffers
        for (int i = 0; i < 1000 * 1000 * 3; i++)
            buffer[i] = 0;

        for (int i = 0; i < 1000 * 1000; i++)
            zBuffer[i] = -1.f;

        //looping through our generated triangles
        for (int t = 0; t < triangles.size(); t++)
        {
            triangles[t].prepare(lp); //compute the slopes/other info necessary to get the row bounds for each scanline

            for (int scanPosition = ceil__441(500.f + 500.f * (triangles[t].X[0]));
                scanPosition <= floor__441(500.f + 500.f * triangles[t].X[2]);
                scanPosition++)
            {
                //oob check
                if (scanPosition >= screen.width
                    || scanPosition < 0)
                    continue;

                triangles[t].getRowBounds(((double)scanPosition - 500.f) / 500.f, rowLimits, colorsAtLimits, diffuseShaderAtLimits, zAtLimits); //get the min and max rows we need to add to the screen buffer

                //fill in the pixels by row (for each column scan)
                for (int currentRow = ceil__441(500.f + 500.f * rowLimits[0]);
                    currentRow <= floor__441(500.f + 500.f * rowLimits[1]);
                    currentRow++)
                {
                    //oob check
                    if (currentRow >= screen.height
                        || currentRow < 0)
                        continue;

                    //stop if this point is hidden
                    double zValue = zAtLimits[0] + (zAtLimits[1] - zAtLimits[0]) * (((double)currentRow - 500.f) / 500.f - rowLimits[0]) / (rowLimits[1] - rowLimits[0]);
                    if (zValue < zBuffer[currentRow * screen.width + scanPosition])
                        continue;

                    zBuffer[currentRow * screen.width + scanPosition] = zValue; //set z buffer         

                    double shade = std::max((double) 0.f, diffuseShaderAtLimits[0] + (diffuseShaderAtLimits[1] - diffuseShaderAtLimits[0]) * (((double)currentRow - 500.f) / 500.f - rowLimits[0]) / (rowLimits[1] - rowLimits[0]));

                    //set this pixel's color
                    buffer[currentRow * screen.width * 3 + scanPosition * 3 + 0] = ceil__441(std::min(shade * (colorsAtLimits[0][0] + (colorsAtLimits[1][0] - colorsAtLimits[0][0]) * (((double)currentRow - 500.f) / 500.f - rowLimits[0]) / (rowLimits[1] - rowLimits[0])), (double)1.f) * 255.f);
                    buffer[currentRow * screen.width * 3 + scanPosition * 3 + 1] = ceil__441(std::min(shade * (colorsAtLimits[0][1] + (colorsAtLimits[1][1] - colorsAtLimits[0][1]) * (((double)currentRow - 500.f) / 500.f - rowLimits[0]) / (rowLimits[1] - rowLimits[0])), (double)1.f) * 255.f);
                    buffer[currentRow * screen.width * 3 + scanPosition * 3 + 2] = ceil__441(std::min(shade * (colorsAtLimits[0][2] + (colorsAtLimits[1][2] - colorsAtLimits[0][2]) * (((double)currentRow - 500.f) / 500.f - rowLimits[0]) / (rowLimits[1] - rowLimits[0])), (double)1.f) * 255.f);
                }
            }
        }

        //set up the filename and write
        char outputName[10];
        std::sprintf(outputName, "frame%03d", f * 250);
        WriteImage(image, outputName);
    }

    free(rowLimits);
    free(zAtLimits);
    for (int i = 0; i < 2; i++)
        free(colorsAtLimits[i]);
    free(colorsAtLimits);  
    free(diffuseShaderAtLimits);
    free(zBuffer);
    
    return 0;
}