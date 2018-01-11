#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>

using namespace std;

#include "VFlow.hpp" // Functions to read velocities 
#include "Constants.hpp"
#include "ParseParameters.hpp"
#include "mt19937arParallel.hpp"

double vsink;
vectorXYZ sigma;
double Dh;
double Dz;

int SetupLagrangianEngine(const string & LagEngParamsFileName, int random) {
  vsink=getDoubleParam(LagEngParamsFileName, "Vsink");

  if(random==1){
  vectorXYZ diffusivity;
  diffusivity=getVectorXYZParam(LagEngParamsFileName, "Diffusivity");

  sigma.x = sqrt(2.0*diffusivity.x*secondsday);
  sigma.y = sqrt(2.0*diffusivity.y*secondsday);
  sigma.z = sqrt(2.0*diffusivity.z*secondsday);
  
#ifdef DEBUG
  cout << "Lagragian  parameters";
  cout << " Vsink = "<< vsink <<endl;
  cout << " Diffusivity = "<< diffusivity <<endl;
  cout << " Sigma = "<< sigma <<endl;
  cout << endl;
#endif
  }
    
#ifdef DEBUG
  cout << "Lagragian  parameters";
  cout << " Vsink = "<< vsink <<endl;
  cout << endl;
#endif  
  
  return 0;
}

int RK4(double t0, double intstep, vectorXYZ *point, int (*velocity)(double ,vectorXYZ , vectorXYZ* )) {
  vectorXYZ point2,point3,point4;
  vectorXYZ v1,v2,v3,v4;
  
  double tstep2,tstep6;
  double h;//scale factor for spherical coordinates
  double t;

  /* Time increments */
  tstep2 = intstep*0.5;
  tstep6 = intstep/6.0;
  
  /* Calculate V1: */
  if(velocity(t0,*point, &v1))
    return 1;
  h = rearth * cos(rads*(point->y));
  v1.x = degrees*(v1.x / h );
  v1.y = degrees*(v1.y / rearth);


  /* Calculate V2: */
  t = t0 + tstep2;
  point2 = *point + (tstep2 * v1);

  if(velocity(t,point2, &v2))
    return 1;

  h = rearth * cos(rads*(point2.y));
  v2.x = degrees*(v2.x / h); // rads velocity
  v2.y = degrees*(v2.y / rearth); // rads velocity


  /* Calculate V3: */
  point3 = *point + (tstep2 * v2);

  if(velocity(t,point3, &v3))
    return 1;

  h = rearth * cos(rads*(point3.y));
  v3.x = degrees*(v3.x / h);
  v3.y = degrees*(v3.y / rearth);

  
  /* Calculate V4: */
  t = t0 + intstep;
  point4 = *point + (intstep * v3);
  
  if(velocity(t,point4, &v4))
    return 1;

  h = rearth * cos(rads*(point4.y));
  v4.x = degrees*(v4.x / h);
  v4.y = degrees*(v4.y / rearth);

  /* Calculate Final point */  
  *point += (tstep6 * (v1 + v4 + 2.0*v2 + 2.0*v3));

  if(velocity(t0+intstep,*point, &v1))//Check if the new position is on the domain vflow
    return 1;
  
  
  return 0;
}

int GetVflowplusVsink(double t,vectorXYZ point, vectorXYZ *vint) {
  if(GetVelocity( t, point, vint))
    return 1;
  
  vint->z = vint->z - vsink;
  return 0;
}

int Heun(double t0, double intstep, vectorXYZ *point, int (*velocity)(double ,vectorXYZ , vectorXYZ* ),StateRN *state)
{
  vectorXYZ point1,point2;
  vectorXYZ v1,v2;
  vectorXYZ g1,g2;
  vectorXYZ u, uh, aux;
  vectorXYZ k,l;
  double h;

  u.x = sgaussrand(state); 
  u.y = sgaussrand(state); 
  u.z = sgaussrand(state); 
  
  /* Calculate V1: */
  point1 = *point;
  if(velocity(t0, point1, &v1))
    return 1;
  h = rearth * cos(rads*(point1.y));
  v1.x = degrees*(v1.x / h ); 
  v1.y = degrees*(v1.y / rearth);

  /* Calculate G1*/
  g1 = sigma;
  g1.x = degrees*(g1.x / h );
  g1.y = degrees*(g1.y / rearth);
 

  /* Calculate V2: */
  uh=sqrt(intstep)*u;
  aux = intstep * v1 + uh * g1;
  point2 = point1 + aux;

  if(velocity(t0+intstep,point2, &v2))
    return 1;
  h = rearth * cos(rads*(point2.y));
  v2.x = degrees*(v2.x / h ); // rads velocity
  v2.y = degrees*(v2.y / rearth); // rads velocity

  /* Calculate G2*/
  g2 = sigma;
  g2.x = degrees*(g2.x / h ); // rads velocity
  g2.y = degrees*(g2.y / rearth); // rads velocity

  /* Compute the new position */
  *point += 0.5*(aux + intstep*v2+ uh*g2);

  return 0;
}

