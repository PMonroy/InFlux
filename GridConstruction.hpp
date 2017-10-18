#ifndef GRIDCONSTRUCTION
#define GRIDCONSTRUCTION

#include "VectorXYZ.hpp"
 
int MakeRegular2DGrid(vectorXYZ domainmin, vectorXYZ intergrid, vectorXYZ domainmax, int KcompNormal, vector <vectorXYZ> *point, int *dim0, int *dim1);
int MakeRandom2DGrid(vectorXYZ domainmin, vectorXYZ domainmax, int KcompNormal, int npoints, vector<vectorXYZ> *point);

int SatelliteGrid(double InterSpacing, vector<vectorXYZ> *point, vector<vectorXYZ> *satpoint);

#endif
