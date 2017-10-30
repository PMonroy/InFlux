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
struct CHistParameters {

  const vectorXYZ DomainBL;
  const vectorXYZ DomainTR;
  const vectorXYZ InterSpacing;
  const eqdate DepositionDate;
  const double TimeStep;
  const double Nparticles;
  const double Radius;
  const double Vsink;
  
  // Define a constructor that will load stuff from a configuration file.
  CHistParameters(const string & CHistParamsFileName)
  :DomainBL(getVectorXYZParam(CHistParamsFileName, "DomainBL"))
  ,DomainTR(getVectorXYZParam(CHistParamsFileName, "DomainTR"))
  ,InterSpacing(getVectorXYZParam(CHistParamsFileName, "InterSpacing"))
  ,DepositionDate(getEqDateParam(CHistParamsFileName, "DepositionDate"))
  ,TimeStep(getDoubleParam(CHistParamsFileName, "TimeStep"))
  ,Nparticles(getDoubleParam(CHistParamsFileName, "Nparticles"))
  ,Radius(getDoubleParam(CHistParamsFileName, "Radius"))
  ,Vsink(getDoubleParam(CHistParamsFileName, "Vsink"))
{}
};

int main(int argc, char **argv){

  /********************************************************************************
   * READ PARAMETERS FROM FILE
   ********************************************************************************/

  string SConfigurationFile; // File name that stores the parameters
  string me=argv[0];
  
  if(GetcmdlineParameters(argc, argv, &SConfigurationFile)) {//Get cmd line parameters 
    cout << me << ": Error getting parameters from command line." <<endl;
    return 1;
  }
  
#ifdef DEBUG
  cout << endl;
  cout << "CHist Parameters from command line: "<<endl;
  cout <<" Configuraton File "<< SConfigurationFile<<endl;
#endif
  
  const CHistParameters CHistParams(SConfigurationFile);

#ifdef DEBUG
  cout<<endl;
  cout << "CHist Parameters from file "<<SConfigurationFile<<" :"<<endl; 
  cout << " DomainBL "<< CHistParams.DomainBL<<endl;
  cout << " DomainTR "<< CHistParams.DomainTR<<endl;
  cout << " InterSpacing "<< CHistParams.InterSpacing<<endl;  
  cout << " Deposition Date "<< CHistParams.DepositionDate.GetMday()<<"-";
  cout<<CHistParams.DepositionDate.GetMonth()<<"-";
  cout<<CHistParams.DepositionDate.GetYear()<<endl;  
  cout << " Nparticles "<<CHistParams.Nparticles<<endl;
  cout << " Radius "<<CHistParams.Radius<<endl;
  cout << " TimeStep "<<CHistParams.TimeStep<<endl;
  cout << " Vsink "<<CHistParams.Vsink<<endl;
#endif

  // Total number of particles seeded

  double idensity;

  idensity=CHistParams.Nparticles/(pi*CHistParams.Radius*CHistParams.Radius);

  double DomainArea;

  DomainArea=((CHistParams.DomainTR.x-CHistParams.DomainBL.x)*rearth*rearth*rads)*(sin(CHistParams.DomainTR.y*rads)-sin(CHistParams.DomainBL.y*rads));
    
  int Npseeded;

  Npseeded=int(idensity*DomainArea);
  
#ifdef DEBUG
  cout<<endl;
  cout << "CHist derived parameters (km^2):"<<endl; 
  cout << " idensity "<< idensity*1e6 <<endl;
  cout << " DomainArea "<< DomainArea*1e-6 <<endl;  
  cout << " Npseeded "<< Npseeded <<endl;
#endif  

  /********************************************************************************
   * GRID CONSTRUCTION
   ********************************************************************************/

  SetupVflow(SConfigurationFile);// SETUP VELOCITY MODEL: read parameters from file

#ifdef DEBUG
  cout << "LOAD VELOCITY FIELD" << endl; //Verbose: loading vel. grid  
  cout << " *Loading velocity grid "; //Verbose: loading vel. grid
#endif
  if((LoadVLatLonGrid(CHistParams.DepositionDate))!=0){//Load velocity grid
    cout<<" [Fail] "<<endl;
    return 1;
  }
#ifdef DEBUG
  cout<<" [OK] " << endl;//Verbose: loading velocities from model 
#endif

  vector<vectorXYZ> itracer;
  int KcompNormal=2*(CHistParams.TimeStep>0)-1;
  int GridDimx, GridDimy;

  if(MakeRegular2DGrid(CHistParams.DomainBL, 
		       CHistParams.InterSpacing, 
		       CHistParams.DomainTR,
		       KcompNormal,
		       &itracer,
		       &GridDimx,
		       &GridDimy)){//Grid construction 
    cout << "Regular grid initial positions: [Fail]" << endl;
    return 1;
  }
  unsigned int NumGridTracers=itracer.size();
#ifdef DEBUG
  cout<<endl;
  cout << "Grid Initial Position: "<<endl; 
  cout << " Num. grid points "<< itracer.size() <<endl;
  cout << " Dim x "<< GridDimx << endl;
  cout << " Dim y "<< GridDimy << endl;
  cout << " depth " << itracer[0].z<< endl;
#endif

  if(MakeRandom2DGrid(CHistParams.DomainBL,
		      CHistParams.DomainTR,
		      KcompNormal,
		      Npseeded,
		      &itracer)){//Grid construction 
    cout << "Regular tracer initial positions: [Fail]" << endl;
    return 1;
  }

  unsigned int NumTracers=itracer.size();
  unsigned int NumSampleTracers=NumTracers-NumGridTracers;
    
#ifdef DEBUG
  cout<<endl;
  cout << "Initial number of tracers:"<<endl; 
  cout << " Number of tracers "<< NumTracers <<endl;
  cout << " Number of Sample tracers "<< NumSampleTracers <<endl;
#endif

  /********************************************************************************
   * SET UP TIME PARAMETERS
   ********************************************************************************/

  double StdVz;
  
  SetupVflow(SConfigurationFile);// SETUP VELOCITY MODEL: read parameters from file
  
  if(CalcStdVz(CHistParams.DepositionDate, &StdVz)){
    cout << " Calculation of Std Vz: [Fail]" << endl;
    return 1;
  }
  
#ifdef DEBUG
  cout << endl;
  cout << "SET UP TIME PARAMETERS: "<<endl; 
  cout << " Standard deviation of Vz at deposition time "<<StdVz<<endl;
#endif


  int tau=3*int((CHistParams.DomainTR.z-CHistParams.DomainBL.z)/(CHistParams.Vsink));
  
  eqdate inidate;
  double tstart;
  double tend;
  int ascnd = CHistParams.TimeStep > 0;

  if(ascnd){
    unsigned int initime = CHistParams.DepositionDate.GetDay();
    inidate.SetDay(initime-2);
    tstart = 2.0;
    tend = (double) tau;
  }
  else{
    unsigned int initime = CHistParams.DepositionDate.GetDay()-tau;
    inidate.SetDay(initime-2); 
    tend = 2.0;
    tstart = (double) tau;
  }

#ifdef DEBUG
  cout << " inidate "<< inidate.GetMday() <<"-";
  cout << inidate.GetMonth() <<"-";
  cout << inidate.GetYear() <<endl;
  cout << " tstart "<< tstart <<endl;
  cout << " tend " << tend <<endl;
  cout << " tau " << tau <<endl;
  cout << endl;
#endif

  /**********************************************
   * LOAD VELOCITY FIELD
   **********************************************/

#ifdef DEBUG
  cout<<" *Loading velocities from model " << endl;//Verbose: loading velocities from model 
#endif

  if(LoadVelocities(inidate,tau+4)!=0){// Load velocity field
    cout << "[Fail]"<< endl;
    return 1;
  }
  
  /**********************************************
   * LANGRANGIAN ENGINE
   **********************************************/

#ifdef DEBUG
  cout << "LANGRANGIAN ENGINE:" << endl;
#endif

  SetupLagrangianEngine(SConfigurationFile);
  
  vector<vectorXYZ> ftracer,tracerBuffer;
  ftracer=itracer;
  tracerBuffer=itracer;
 
  vector<double> ftime(NumTracers,0.0);
  vector<int> FlagRK4(NumTracers,0);

  double fdepth=(CHistParams.TimeStep<0)?CHistParams.DomainTR.z:CHistParams.DomainBL.z;
  vector<int> Flagfdepth(NumTracers,0);

#ifdef DEBUG
  cout << " fdepth " << fdepth<< endl;
#endif
  
  /****************************************************
   * TRACER LOOP
   ****************************************************/

  double t;
  
  omp_set_num_threads(20);
  #pragma omp parallel for default(shared) private(t)// Parallelizing the code for computing trajectories
  for (unsigned int q=0; q<NumTracers; q++) {
    for(t=tstart; ascnd==1?(t<tend):(t>=tend); t+=CHistParams.TimeStep) {
      
      //COMPUTE NEW POSITION
      tracerBuffer[q]=ftracer[q];// it stores previous position in tracerBuffer[q]
      if(RK4(t, CHistParams.TimeStep, &ftracer[q], GetVflowplusVsink)==1){
	ftime[q]=t;
	FlagRK4[q]=1;
	break;
      }
      ftime[q]=t+CHistParams.TimeStep;

      double f=(fdepth-tracerBuffer[q].z);
      double fmid=(fdepth-ftracer[q].z);
      double discriminant=f*fmid;   

      if(discriminant<=0.0) {
	double rtb,dt;
	rtb=0.0;
	dt=CHistParams.TimeStep;
	double tmid;
	for (int j=0;j<20;j++) {
	  
	  if(Flagfdepth[q]==1) break;
	  ftracer[q]=tracerBuffer[q];
	  tmid=rtb+(dt*=0.5);
	  if(RK4(t, tmid, &ftracer[q], GetVflowplusVsink)==1){
	    FlagRK4[q]=1;
	    break;
	  }
	  fmid=(fdepth-ftracer[q].z);
	  if (fmid <= 0.0) rtb=tmid;
	  if (dt < (CHistParams.TimeStep/16.0) || fmid == 0.0){
	    ftime[q]=t+tmid;
	    Flagfdepth[q]=1;
	  }
	}
      }    
      if(Flagfdepth[q]==1) break; 

    }
  }

  /* Free Velocities*/
  FreeMemoryVelocities(tau+4);

  /****************************************************
   * COMPUTATION NUMERICAL DENSITY
   ****************************************************/
  vector<int> density(NumGridTracers,0);
  vector<int> densitySquare(NumGridTracers,0);
  unsigned int l;
 
omp_set_num_threads(20);
#pragma omp parallel for default(shared) private(l)// Parallelizing the code for computing trajectories  
  for (unsigned int q=0; q<NumGridTracers; q++) {

    if(FlagRK4[q]==1 || Flagfdepth[q]==0) continue;
    
    vectorXYZ delta;

    delta.x=CHistParams.Radius/(rearth*cos(ftracer[q].y*rads));
    delta.y=CHistParams.Radius/rearth;
    delta.z=0.0;

    vectorXYZ max,min;

    max=ftracer[q]+(degrees*delta);
    min=ftracer[q]-(degrees*delta);
    
    for(l=NumGridTracers; l<NumTracers; l++) {
      if(FlagRK4[l]==0 &&
	 Flagfdepth[l]==1 &&
	 ftracer[l].x<=max.x &&
	 ftracer[l].y<=max.y &&
	 ftracer[l].x>=min.x &&
	 ftracer[l].y>=min.y){
	densitySquare[q]++;
	vectorXYZ alpha;

	alpha.x=cos(ftracer[q].y*rads)*rads;
	alpha.y=rads;
	alpha.z=0.0;
	
	vectorXYZ beta;

	beta=ftracer[q]-ftracer[l];
	beta*=alpha;
	
	if(rearth*sqrt(scalar(beta,beta))<=CHistParams.Radius) 	density[q]++;
      }
    }
  }


  /**********************************************
   * WRITE RESULTS
   **********************************************/

  string rawfilename;;
  size_t lastdot = SConfigurationFile.find_last_of(".");
  if(lastdot == string::npos){
    rawfilename=SConfigurationFile;
  } else {
    rawfilename=SConfigurationFile.substr(0,lastdot);
  }


  // Position FILE
  string posfilename=rawfilename+".pos";  
  ofstream posfile(posfilename.c_str());
  for(unsigned int q=0; q<NumGridTracers; q++){
    posfile<<itracer[q]<<" ";
    posfile<<ftracer[q]<<" ";
    posfile<<Flagfdepth[q]<<" ";
    posfile<<FlagRK4[q]<<" ";
    posfile<<density[q]/CHistParams.Nparticles<<" ";
    posfile<<density[q]<<endl;    
  }
  posfile.close();

  // FINAL GRID
  string VTKffilename=rawfilename+".vtk";  
  ofstream offile(VTKffilename.c_str());
  
  offile<<"# vtk DataFile Version 3.0"<<endl;
  offile<<"downward grid numerical density"<<endl; 
  offile<<"ASCII"<<endl;
  offile<<"DATASET POLYDATA"<<endl;
  offile<<"POINTS "<<NumGridTracers<<" float"<<endl;
  for(unsigned int q=0; q<NumGridTracers; q++){
    offile<<ftracer[q].x<<" "<<ftracer[q].y<<" "<< fdepth <<endl;
  }
  offile<<"VERTICES "<<NumGridTracers<<" "<<NumGridTracers*2<<endl;
  for(unsigned int q=0; q<NumGridTracers; q++){
    offile<<"1 "<<q<<endl;
  }
  offile<<"POINT_DATA "<<NumGridTracers<<endl;
  
  offile<<"SCALARS density int 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<density.size(); q++) {
      offile<<density[q]<<endl;
  }  
  offile<<"SCALARS densitySquare int 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<densitySquare.size(); q++) {
      offile<<densitySquare[q]<<endl;
  }  


  return 0;

}

string numprintf(int ndigits, int ndecimals, double number){
  ostringstream stream;// it needs to include <sstream>
  double factor=pow(10,ndecimals); // it needs to include <cmath>
  stream << fixed; //it needs to include <iostream>
  stream << setfill('0') << setw(ndigits);
  stream << setprecision(0) << internal<< (number*factor);
  return stream.str();
}
