#ifndef LENGINE
#define LENGINE

int SetupLagrangianEngine(const string & LagEngParamsFileName);
int GetVflowplusVsink(double t,vectorXYZ point, vectorXYZ *vint);
int RK4(double t0, double intstep, vectorXYZ *point, int (*velocity)(double t,vectorXYZ point, vectorXYZ *vint));

#endif
