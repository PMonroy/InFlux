#include <cstdio>
#include <cstdlib>
#include <iostream> // cout, endl...  
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <iomanip> // setfill, setw
#include <string>
#include <omp.h>

using namespace std;

#include "ParseParameters.hpp"
#include "GridConstruction.hpp"
#include "LagrangianEngine.hpp"
#include "Constants.hpp"
#include "VFlow.hpp"
#include "EqDate.hpp"
#include "VectorXYZ.hpp"
#include "VTKdump.hpp"
#include "mt19937arParallel.hpp"

string numprintf(int ndigits, int ndecimals, double number);
struct InFluxParameters {

  const vectorXYZ DomainBL;
  const vectorXYZ DomainTR;
  const vectorXYZ InterSpacing;
  const eqdate DepositionDate;
  const double TimeStep;
  const int Random;
  const double Vsink;
  
  // Define a constructor that will load stuff from a configuration file.
  InFluxParameters(const string & InFluxParamsFileName)
  :DomainBL(getVectorXYZParam(InFluxParamsFileName, "DomainBL"))
  ,DomainTR(getVectorXYZParam(InFluxParamsFileName, "DomainTR"))
  ,InterSpacing(getVectorXYZParam(InFluxParamsFileName, "InterSpacing"))
  ,DepositionDate(getEqDateParam(InFluxParamsFileName, "DepositionDate"))
  ,TimeStep(getDoubleParam(InFluxParamsFileName, "TimeStep"))
  ,Random(getIntParam(InFluxParamsFileName, "Random"))
  ,Vsink(getDoubleParam(InFluxParamsFileName, "Vsink"))
{}
};

int main(int argc, char **argv){

  /********************************************************************************
   * READ PARAMETERS FROM FILE
   ********************************************************************************/

  string SConfigurationFile; // File name that stores the parameters
  string me=argv[0];
  double radius;
  if(GetcmdlineCGParameters(argc, argv, &SConfigurationFile,&radius)) {//Get cmd line parameters 
    cout << me << ": Error getting parameters from command line." <<endl;
    return 1;
  }
  
#ifdef DEBUG
  cout << endl;
  cout << "InFlux Parameters from command line: "<<endl;
  cout <<" Configuraton File "<< SConfigurationFile<<endl;
#endif
  
  const InFluxParameters InFluxParams(SConfigurationFile);

#ifdef DEBUG
  cout<<endl;
  cout << "InFlux Parameters from file "<<SConfigurationFile<<" :"<<endl; 
  cout << " DomainBL "<< InFluxParams.DomainBL<<endl;
  cout << " DomainTR "<< InFluxParams.DomainTR<<endl;
  cout << " InterSpacing " <<InFluxParams.InterSpacing<<endl;
  cout << " Deposition Date "<< InFluxParams.DepositionDate.GetMday()<<"-";
  cout<<InFluxParams.DepositionDate.GetMonth()<<"-";
  cout<<InFluxParams.DepositionDate.GetYear()<<endl;
  cout << " TimeStep "<<InFluxParams.TimeStep<<endl;
  cout << " Random "<<InFluxParams.Random<<endl;
  cout << " Vsink "<<InFluxParams.Vsink<<endl;
  cout << " Radius "<< radius<<endl;
#endif

  string rawfilename;;
  size_t lastdot = SConfigurationFile.find_last_of(".");
  if(lastdot == string::npos){
    rawfilename=SConfigurationFile;
  } else {
    rawfilename=SConfigurationFile.substr(0,lastdot);
  }

  // Position FILE
  string posfilename=rawfilename+".pos";  
  ifstream ifile(posfilename.c_str());
  if(!ifile.is_open()){
    cout<<" Skipping unreadable file " << posfilename.c_str() <<" "<<endl; 
    return 1;
  }

  vector<vectorXYZ> tracer;
  vector<int> Flagfdepth;
  vector<int> FlagRK4;
  vector<double> FactorDensity;

  vectorXYZ tracerBuffer;
  int FlagfdepthBuffer;
  int FlagRK4Buffer;
  double FactorDensityBuffer;
  vectorXYZ UnitNormalBuffer;
  vectorXYZ VfinalBuffer;
  double muBuffer;
  
  
  while(ifile>>tracerBuffer
	>>FlagfdepthBuffer
	>>FlagRK4Buffer
	>>FactorDensityBuffer
	>>UnitNormalBuffer
	>>VfinalBuffer
	>>muBuffer){
    tracer.push_back(tracerBuffer); 
    Flagfdepth.push_back(FlagfdepthBuffer); 
    FlagRK4.push_back(FlagRK4Buffer); 
    FactorDensity.push_back(FactorDensityBuffer); 
  }
  ifile.close();

  // FILTER

  vector<double> FactorDensityCG(tracer.size(),0.0);
  vector<int> n(tracer.size(),0);

  for (unsigned int q=0; q<tracer.size(); q++) {
    
    if(FlagRK4[q]==1 || Flagfdepth[q]==0) continue;
    
    vectorXYZ delta;
    
    delta.x=radius/(rearth*cos(tracer[q].y*rads));
    delta.y=radius/rearth;
    delta.z=0.0;
    
    vectorXYZ max,min;
    
    max=tracer[q]+(degrees*delta);
    min=tracer[q]-(degrees*delta);

    vectorXYZ alpha;
    
    alpha.x=cos(tracer[q].y*rads)*rads;
    alpha.y=rads;
    alpha.z=0.0;
    
    for(unsigned int l=0; l<tracer.size(); l++) {
      if(FlagRK4[l]==0 &&
	 Flagfdepth[l]==1 &&
	 tracer[l].x<=max.x &&
	 tracer[l].y<=max.y &&
	 tracer[l].x>=min.x &&
	 tracer[l].y>=min.y){
	
	vectorXYZ beta;	
	beta=tracer[q]-tracer[l];
	beta*=alpha;
	double norm2beta=scalar(beta,beta);
	if(rearth*sqrt(norm2beta)<=radius){
	  FactorDensityCG[q]+=FactorDensity[l];
	  n[q]++;
	}
      }
    }
    if(n[q]>0) FactorDensityCG[q]/=double(n[q]);
  }

  /**********************************************
   * WRITE RESULTS
   **********************************************/
  string outputfilename=rawfilename+"CG"+numprintf(1,0,radius)+"m.pos";  
  ofstream posfile(outputfilename.c_str());
  for(unsigned int q=0; q<tracer.size(); q++){
    posfile<<tracer[q]<<" ";
    posfile<<Flagfdepth[q]<<" ";
    posfile<<FlagRK4[q]<<" ";
    posfile<<FactorDensityCG[q]<<" ";
    posfile<<n[q]<<endl;
  }
  posfile.close();

  // VTK FILE
  string VTKffilename=rawfilename+"CG"+numprintf(1,0,radius)+"m.vtk";  
  ofstream offile(VTKffilename.c_str());
  
  offile<<"# vtk DataFile Version 3.0"<<endl;
  offile<<"Final grid"<<endl; 
  offile<<"ASCII"<<endl;
  offile<<"DATASET POLYDATA"<<endl;
  offile<<"POINTS "<<tracer.size()<<" float"<<endl;
  for(unsigned int q=0; q<tracer.size(); q++){
    offile<<tracer[q].x<<" "<<tracer[q].y<<" "<< InFluxParams.DomainTR.z <<endl;
  }
  offile<<"VERTICES "<<tracer.size()<<" "<<tracer.size()*2<<endl;
  for(unsigned int q=0; q<tracer.size(); q++){
    offile<<"1 "<<q<<endl;
  }
  offile<<"POINT_DATA "<<tracer.size()<<endl;
  
  offile<<"SCALARS flagDepth int 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<Flagfdepth.size(); q++) {
    offile<<Flagfdepth[q]<<endl;
  }  

  offile<<"SCALARS flagRK4 int 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<FlagRK4.size(); q++) {
    offile<<FlagRK4[q]<<endl;
  }

  offile<<"SCALARS n int 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<n.size(); q++) {
    offile<<n[q]<<endl;
  }
 
  offile<<"SCALARS FactorDensity float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<FactorDensityCG.size(); q++) {
    offile<<FactorDensityCG[q]<<endl;
  }

  offile.close();

  return 0;
}

string numprintf(int ndigits, int ndecimals, double number) {
  ostringstream stream;// it needs to include <sstream>
  double factor=pow(10,ndecimals); // it needs to include <cmath>
  stream << fixed; //it needs to include <iostream>
  stream << setfill('0') << setw(ndigits);
  stream << setprecision(0) << internal<< (number*factor);
  return stream.str();
}
