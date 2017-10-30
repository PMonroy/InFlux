#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

#include "VectorXYZ.hpp"
#include "Constants.hpp"
#include "VFlow.hpp"
#include "mt19937ar.hpp"

double RangeRand(double fmin, double fmax);

int MakeRegular2DGrid(vectorXYZ domainmin, vectorXYZ intergrid, vectorXYZ domainmax, int KcompNormal,
		      vector<vectorXYZ> *point, int *dim0, int *dim1){  
  double x,y,z;
  int i,j;

  int nx, ny, numpoints;

  // INTEGER SPATIAL DIMENSIONS
  ny = floor((domainmax.y-domainmin.y)/intergrid.y);
  nx = floor((domainmax.x-domainmin.x)/intergrid.x);
  
  numpoints=nx*ny;
  
  if(numpoints<1){
    cout <<"Error: Number of grid points is negative"<<endl;
    return 1;
  }

  (*point).reserve(numpoints);

  if(KcompNormal<0){
    z=domainmin.z; 
  }
  else if(KcompNormal>0){
    z=domainmax.z;
  }
  else{
    return 1;
  }  
    
  for(j=0; j<ny; j++){
    y=domainmin.y+intergrid.y*double(j);
    for(i=0; i<nx; i++){
      x=domainmin.x + intergrid.x*double(i);
      (*point).push_back(vectorXYZ(x,y,z));
    }
  }

  
  *dim0=nx;
  *dim1=ny;
  
  return 0;
}

int SatelliteGrid(double InterSpacing, vector<vectorXYZ> *point, vector<vectorXYZ> *satpoint){

  vectorXYZ Delta0,Delta1,Delta2,Delta3;
  
  Delta1.x=0.0;
  Delta1.y=degrees*(InterSpacing/rearth);
  Delta1.z=0.0;
  
  Delta3.x=0.0;
  Delta3.y=-1.0*degrees*(InterSpacing/rearth);
  Delta3.z=0.0;
  
  for(unsigned int q=0; q<(*point).size(); q++){
    
    Delta0.x=degrees*(InterSpacing/(rearth*cos(rads*(*point)[q].y)));
    Delta0.y=0.0;
    Delta0.z=0.0;
    
    Delta2.x=-1.0*Delta0.x;
    Delta2.y=0.0;
    Delta2.z=0.0;

    (*satpoint).push_back((*point)[q]+Delta0);
    (*satpoint).push_back((*point)[q]+Delta1);
    (*satpoint).push_back((*point)[q]+Delta2);
    (*satpoint).push_back((*point)[q]+Delta3);
    
  }

  
  return 0;
}


int MakeRandom2DGrid(vectorXYZ domainmin, vectorXYZ domainmax, int KcompNormal, int npoints, vector<vectorXYZ> *point){

  vectorXYZ position;
  
  (*point).reserve(npoints);

  if(KcompNormal<0){
    position.z=domainmin.z; 
  }else if(KcompNormal>0){
    position.z=domainmax.z;
  }else{
    return 1;
  }
  unsigned long s=1984;
  init_genrand(s);
  
  for(int q=0; q<npoints; q++){
    position.x=RangeRand(domainmin.x, domainmax.x);
    position.y=degrees*asin(RangeRand(sin(rads*domainmin.y), sin(rads*domainmax.y)));
    if(IsLand(position)==0) (*point).push_back(position);
  }
  return 0;
}
double RangeRand(double fmin, double fmax){
  double f=genrand_real1();
  return fmin+f*(fmax-fmin);
}
