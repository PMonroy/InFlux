#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include "VectorXYZ.hpp"

void MakeVTKStructuredPointsFloat(vectorXYZ origin, vectorXYZ spacing, unsigned int dimx, unsigned int dimy, unsigned int dimz, vector<double> F,string filename)
{

  ofstream ofile(filename.c_str());
  unsigned int q;
  unsigned int dim=(dimx)*(dimy);
  
  ofile<<"# vtk DataFile Version 3.0"<<endl;
  ofile<<"HOLA CARACOLA "<<endl; 
  ofile<<"ASCII"<<endl;
  ofile<<"DATASET STRUCTURED_POINTS"<<endl;
  ofile<<"DIMENSIONS "<<dimx<<" "<<dimy<<" "<<1<<endl;
  ofile<<"ORIGIN "<<origin.x<<" "<<origin.y<<" "<<origin.z <<endl;
  ofile<<"SPACING "<<spacing.x<<" "<<spacing.y<<" "<<spacing.z <<endl;
  ofile<<"POINT_DATA "<<dim<<endl;
  ofile<<"SCALARS nose float 1"<<endl;
  ofile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<F.size(); q++) {
      ofile<<F[q]<<endl;
  }
  ofile.close();
}
void MakeVTKStructuredPointsInt(vectorXYZ origin, vectorXYZ spacing, unsigned int dimx, unsigned int dimy, unsigned int dimz, vector<int> F,string filename)
{

  ofstream ofile(filename.c_str());
  unsigned int q;
  unsigned int dim=(dimx)*(dimy);
  
  ofile<<"# vtk DataFile Version 3.0"<<endl;
  ofile<<"HOLA CARACOLA "<<endl; 
  ofile<<"ASCII"<<endl;
  ofile<<"DATASET STRUCTURED_POINTS"<<endl;
  ofile<<"DIMENSIONS "<<dimx<<" "<<dimy<<" "<<1<<endl;
  ofile<<"ORIGIN "<<origin.x<<" "<<origin.y<<" "<<origin.z <<endl;
  ofile<<"SPACING "<<spacing.x<<" "<<spacing.y<<" "<<spacing.z <<endl;
  ofile<<"POINT_DATA "<<dim<<endl;
  ofile<<"SCALARS noseint int 1"<<endl;
  ofile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<F.size(); q++) {
      ofile<<F[q]<<endl;
  }
  ofile.close();
}
