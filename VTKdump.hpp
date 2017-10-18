#ifndef VTKDUMP
#define VTKDUMP

#include "VectorXYZ.hpp"
#include "Constants.hpp"

void MakeVTKStructuredPointsFloat(vectorXYZ origin, vectorXYZ spacing, unsigned int dimx, unsigned int dimy, unsigned int dimz, vector<double> F,string filename);
void MakeVTKStructuredPointsInt(vectorXYZ origin, vectorXYZ spacing, unsigned int dimx, unsigned int dimy, unsigned int dimz, vector<int> F,string filename);

#endif
