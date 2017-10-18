#ifndef VELOCITY
#define VELOCITY

#include "EqDate.hpp"
#include "VectorXYZ.hpp"



/* FUNCTIONS */
int SetupVflow(const string & vflowParamsFileName);
int CalcStdVz(eqdate rdate, double *stdVz);
int LoadVLatLonGrid(eqdate rdate);
int LoadVelocities(eqdate startdate, int ntau);


void FreeMemoryVelocities(int ntau);

int GetVelocity(double t,vectorXYZ point, vectorXYZ *vint);
int GradientVz(double t,vectorXYZ point, vectorXYZ *grad);
int TimeDerivativeV(double t,vectorXYZ point, vectorXYZ *dvdt);

#endif
