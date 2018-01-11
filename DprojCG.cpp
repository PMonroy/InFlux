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

string numprintf(int ndigits, int ndecimals, double number);
struct DprojParameters {

  const vectorXYZ DomainBL;
  const vectorXYZ DomainTR;
  const vectorXYZ InterSpacing;
  const double DistSat;
  const eqdate DepositionDate;
  const double TimeStep;
  const double Period;
  const int Random;  
  const double Vsink;
  
  // Define a constructor that will load stuff from a configuration file.
  DprojParameters(const string & DprojParamsFileName)
  :DomainBL(getVectorXYZParam(DprojParamsFileName, "DomainBL"))
  ,DomainTR(getVectorXYZParam(DprojParamsFileName, "DomainTR"))
  ,InterSpacing(getVectorXYZParam(DprojParamsFileName, "InterSpacing"))
  ,DistSat(getDoubleParam(DprojParamsFileName, "DistSat"))
  ,DepositionDate(getEqDateParam(DprojParamsFileName, "DepositionDate"))
  ,TimeStep(getDoubleParam(DprojParamsFileName, "TimeStep"))
  ,Period(getDoubleParam(DprojParamsFileName, "Period"))
  ,Random(getIntParam(DprojParamsFileName, "Random"))
  ,Vsink(getDoubleParam(DprojParamsFileName, "Vsink"))
{}
};

int main(int argc, char **argv){

  /********************************************************************************
   * READ PARAMETERS FROM FILE
   ********************************************************************************/

  string SConfigurationFile; // File name that stores the parameters
  string me=argv[0];
  double radius;

  if(GetcmdlineCGParameters(argc, argv, &SConfigurationFile, &radius)) {//Get cmd line parameters 
    cout << me << ": Error getting parameters from command line." <<endl;
    return 1;
  }
  
#ifdef DEBUG
  cout << endl;
  cout << "Dproj Parameters from command line: "<<endl;
  cout <<" Configuraton File "<< SConfigurationFile<<endl;
#endif
  
  const DprojParameters DprojParams(SConfigurationFile);

#ifdef DEBUG
  cout<<endl;
  cout << "Dproj Parameters from file "<<SConfigurationFile<<" :"<<endl; 
  cout << " DomainBL "<< DprojParams.DomainBL<<endl;
  cout << " DomainTR "<< DprojParams.DomainTR<<endl;
  cout << " InterSpacing " <<DprojParams.InterSpacing<<endl;
  cout << " DistSat " <<DprojParams.DistSat<<endl;
  cout << " Deposition Date "<< DprojParams.DepositionDate.GetMday()<<"-";
  cout<<DprojParams.DepositionDate.GetMonth()<<"-";
  cout<<DprojParams.DepositionDate.GetYear()<<endl;
  cout << " TimeStep "<<DprojParams.TimeStep<<endl;
  cout << " Period "<<DprojParams.Period<<endl;
  cout << " Vsink "<<DprojParams.Vsink<<endl;
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
  vector<double> Projection;
  vector<double> Stretching;
  vector<double> FactorDensity;

  vectorXYZ tracerBuffer;
  int FlagfdepthBuffer;
  int FlagRK4Buffer;
  double ProjectionBuffer;
  double StretchingBuffer;
  double FactorDensityBuffer;
    
  while(ifile>>tracerBuffer
	>>FlagfdepthBuffer
	>>FlagRK4Buffer
	>>ProjectionBuffer
	>>StretchingBuffer
	>>FactorDensityBuffer){
    tracer.push_back(tracerBuffer); 
    Flagfdepth.push_back(FlagfdepthBuffer); 
    FlagRK4.push_back(FlagRK4Buffer);
    Projection.push_back(ProjectionBuffer); 
    Stretching.push_back(StretchingBuffer);  
    FactorDensity.push_back(FactorDensityBuffer); 
  }
  ifile.close();

    // FILTER

  vector<double> FactorDensityCG(tracer.size(),0.0);
  vector<double> FactorDensityHCG(tracer.size(),0.0);
  vector<double> ProjectionCG(tracer.size(),0.0);
  vector<double> StretchingCG(tracer.size(),0.0);
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
	  StretchingCG[q]+=Stretching[l];
	  ProjectionCG[q]+=Projection[l];
	  FactorDensityCG[q]+=FactorDensity[l];
	  FactorDensityHCG[q]+=(1.0/FactorDensity[l]);
	  n[q]++;
	}
      }
    }
    if(n[q]>0) StretchingCG[q]/=double(n[q]);
    if(n[q]>0) ProjectionCG[q]/=double(n[q]);
    if(n[q]>0) FactorDensityCG[q]/=double(n[q]);
    if(n[q]>0){
      FactorDensityHCG[q]=double(n[q])/FactorDensityHCG[q];
    }
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
    posfile<<StretchingCG[q]<<" ";
    posfile<<ProjectionCG[q]<<" ";
    posfile<<FactorDensityCG[q]<<" ";
    posfile<<FactorDensityHCG[q]<<" ";   
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
    offile<<tracer[q].x<<" "<<tracer[q].y<<" "<< DprojParams.DomainTR.z <<endl;
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
 
  offile<<"SCALARS Projection float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<ProjectionCG.size(); q++) {
    offile<<ProjectionCG[q]<<endl;
  }

  offile<<"SCALARS Stretching float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<StretchingCG.size(); q++) {
    offile<<StretchingCG[q]<<endl;
  }

  offile<<"SCALARS FactorDensity float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<FactorDensityCG.size(); q++) {
    offile<<FactorDensityCG[q]<<endl;
  }

  offile<<"SCALARS FactorDensityHarmonic float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<FactorDensityHCG.size(); q++) {
    offile<<FactorDensityHCG[q]<<endl;
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

