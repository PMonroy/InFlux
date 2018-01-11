#ifndef LENGINE
#define LENGINE

#include "mt19937arParallel.hpp"

int SetupLagrangianEngine(const string & LagEngParamsFileName, int random);
int GetVflowplusVsink(double t,vectorXYZ point, vectorXYZ *vint);
int RK4(double t0, double intstep, vectorXYZ *point, int (*velocity)(double t,vectorXYZ point, vectorXYZ *vint));
int Heun(double t0, double intstep, vectorXYZ *point, int (*velocity)(double ,vectorXYZ , vectorXYZ* ), StateRN *state);
#endif
