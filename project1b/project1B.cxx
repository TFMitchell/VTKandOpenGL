/*
* Project 1B - CIS 441
* 
* I'll probably end up adding extra methods to future projects. But I only made changes to main for this one.
* 
* Thomas Mitchell
*/

#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;

double ceil__441(double f)
{
    return ceil(f-0.00001);
}

double floor__441(double f)
{
    return floor(f+0.00001);
}


vtkImageData *
NewImage(int height, int width)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
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
};

class Screen
{
  public:
      unsigned char   *buffer;
      int width, height;
};


std::vector<Triangle>
GetTriangles(void)
{
   std::vector<Triangle> rv(100);

   unsigned char colors[6][3] = { {255,128,0}, {255, 0, 127}, {0,204,204}, 
                                  {76,153,0}, {255, 204, 204}, {204, 204, 0}};
   for (int i = 0 ; i < 100 ; i++)
   {
       int idxI = i%10;
       int posI = idxI*100;
       int idxJ = i/10;
       int posJ = idxJ*100;
       int firstPt = (i%3);
       rv[i].X[firstPt] = posI;
       rv[i].Y[firstPt] = posJ;
       rv[i].X[(firstPt+1)%3] = posI;
       rv[i].Y[(firstPt+1)%3] = posJ+99;
       rv[i].X[(firstPt+2)%3] = posI+10*(idxI+1);
       rv[i].Y[(firstPt+2)%3] = posJ+20*(idxJ+1)-50;
       if (i == 5)
           rv[i].Y[firstPt] = -10;
       if (i == 49)
          rv[i].X[(firstPt+2)%3] = posI+20*(idxI+1);
       rv[i].color[0] = colors[i%6][0];
       rv[i].color[1] = colors[i%6][1];
       rv[i].color[2] = colors[i%6][2];
   }

   return rv;
}

int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = 1000*1000;
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;

   std::vector<Triangle> triangles = GetTriangles();
   
   Screen screen;
   screen.buffer = buffer;
   screen.width = 1000;
   screen.height = 1000;

   for each (Triangle triangle in triangles)
   {
       //this is the starting value for the column scan
       int initialScanPosition = std::min(ceil__441(triangle.X[0]),
           std::min(ceil__441(triangle.X[1]), ceil__441(triangle.X[2])));

       //ending value for column scan
       int maxScanPosition = std::max(floor__441(triangle.X[0]),
           std::max(floor__441(triangle.X[1]), floor__441(triangle.X[2])));


       double rightY = 0.f; //Y value of the rightmost point on the triangle
       double rightX = 0.f; //X value of rightmost point
       std::vector<double> leftYs; //vector of the two left Y values of the triangle
       double leftX = 0.f; //X value of either left point

       //for the three points, figure out if each is the rightmost or a left one
       for (int i = 0; i < 3; i++)
       {
           if (floor__441(triangle.X[i]) == maxScanPosition)
           {
               rightY = triangle.Y[i];
               rightX = triangle.X[i];
           }
           else
           {
               leftYs.push_back(triangle.Y[i]);

               leftX = triangle.X[i];
           }
       }

       //slopes of the bottom and top sides of the triangle
       double bottomSlope = (rightY - std::min(leftYs[0], leftYs[1])) / (rightX - leftX);
       double topSlope = (rightY - std::max(leftYs[0], leftYs[1])) / (rightX - leftX);

       //each scanline (columns)
       for (int scanPosition = initialScanPosition; scanPosition < maxScanPosition + 1; scanPosition++)
       {
           //oob check
           if (scanPosition >= 1000
               || scanPosition < 0)
               continue;

           //fill in the pixels by row (for each column scan)
           for (int currentRow = ceil__441( std::min(leftYs[0], leftYs[1]) + bottomSlope * (scanPosition - leftX)); 
               currentRow < floor__441( std::max(leftYs[0], leftYs[1]) + topSlope * (scanPosition - leftX) + 1); 
               currentRow++)
           {
               //oob check
               if (currentRow >= 1000
                   || currentRow < 0)
                   continue;

               //set this pixel's color
               buffer[(int)(currentRow * 1000 * 3 + scanPosition * 3)] = triangle.color[0];
               buffer[(int)(currentRow * 1000 * 3 + scanPosition * 3 + 1)] = triangle.color[1];
               buffer[(int)(currentRow * 1000 * 3 + scanPosition * 3 + 2)] = triangle.color[2];

           }
       }
   }

   WriteImage(image, "allTriangles");
}
