/*
* Project 1A
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


int main()
{
   std::cerr << "In main!" << endl;
   vtkImageData *image = NewImage(1350, 1024);

   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);

   unsigned char pixel[] = { 0, 0, 0 }; //the temp RGB holder
   int xthStrip = 0; // which strip we're dealing with for math purposes

   for (int r = 0; r < 1350; r++) 
   {
       xthStrip = r / 50; //get strip #

       //given in instructions
       pixel[0] = xthStrip / 9 * 255 / 2;
       pixel[1] = xthStrip / 3 % 3 * 255 / 2;
       pixel[2] = xthStrip % 3 * 255 / 2;

       //assign whole row the same
       for (int c = 0; c < 1024; c++)
       {
           buffer[r * 1024 * 3 + c * 3] = pixel[0];
           buffer[r * 1024 * 3 + c * 3 + 1] = pixel[1];
           buffer[r * 1024 * 3 + c * 3 + 2] = pixel[2];
       }
   }

   WriteImage(image, "proj1A");
}
