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
#include "mt19937ar.hpp"

string numprintf(int ndigits, int ndecimals, double number);
struct CHistParameters {

  const vectorXYZ DomainBL;
  const vectorXYZ DomainTR;
  const vectorXYZ InterSpacing;
  const eqdate DepositionDate;
  const double TimeStep;
  const int Random;
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
  ,Random(getIntParam(CHistParamsFileName, "Random"))
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
  cout << CHistParams.DepositionDate.GetMonth()<<"-";
  cout << CHistParams.DepositionDate.GetYear()<<endl;  
  cout << " Nparticles "<<CHistParams.Nparticles<<endl;
  cout << " Radius "<<CHistParams.Radius<<endl;
  cout << " TimeStep "<<CHistParams.TimeStep<<endl;
  cout << " Random "<<CHistParams.Random<<endl;
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

  vector<vectorXYZ> Grid;
  int KcompNormal=2*(CHistParams.TimeStep>0)-1;
  int GridDimx, GridDimy;

  if(MakeRegular2DGrid(CHistParams.DomainBL, 
		       CHistParams.InterSpacing, 
		       CHistParams.DomainTR,
		       KcompNormal,
		       &Grid,
		       &GridDimx,
		       &GridDimy)){//Grid construction 
    cout << "Regular grid initial positions: [Fail]" << endl;
    return 1;
  }
  unsigned int NumGrid=Grid.size();
  
#ifdef DEBUG
  cout<<endl;
  cout << "Grid Initial Position: "<<endl; 
  cout << " Num. grid points "<< Grid.size() <<endl;
  cout << " Dim x "<< GridDimx << endl;
  cout << " Dim y "<< GridDimy << endl;
  cout << " depth " << Grid[0].z<< endl;
#endif
  
  vector<vectorXYZ> itracer;
  if(MakeRandom2DGrid(CHistParams.DomainBL,
		      CHistParams.DomainTR,
		      KcompNormal,
		      Npseeded,
		      &itracer)){//Grid construction 
    cout << "Regular tracer initial positions: [Fail]" << endl;
    return 1;
  }
  unsigned int NumTracers=itracer.size();
    
#ifdef DEBUG
  cout<<endl;
  cout << "Initial number of tracers:"<<endl; 
  cout << " Number of tracers "<< NumTracers <<endl;
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

  SetupLagrangianEngine(SConfigurationFile,CHistParams.Random);
  
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

#ifdef DEBUG
  cout << " Tracer loop " << endl;
#endif
  
  double t;

  int id;
  unsigned long seed;
  StateRN stateMT;
    
  omp_set_num_threads(20);
  
  #pragma omp parallel default(shared) private(t,id,seed,stateMT)// Parallelizing the code for computing trajectories
  {
    id=omp_get_thread_num();
    seed = 12345678 + (2*id+1);
    init_sgenrand(seed, &stateMT);
#pragma omp for
    for (unsigned int q=0; q<NumTracers; q++) {
      for(t=tstart; ascnd==1?(t<tend):(t>=tend); t+=CHistParams.TimeStep) {
      
	//COMPUTE NEW POSITION
	tracerBuffer[q]=ftracer[q];// it stores previous position in tracerBuffer[q]
	if(CHistParams.Random==1){
	  FlagRK4[q]=Heun(t, CHistParams.TimeStep, &ftracer[q], GetVflowplusVsink, &stateMT);
	} else {
	  FlagRK4[q]=RK4(t, CHistParams.TimeStep, &ftracer[q], GetVflowplusVsink);
	}
	
	if(FlagRK4[q]==1){
	  ftime[q]=t;
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
	    if(CHistParams.Random==1){
	      FlagRK4[q]=Heun(t, tmid, &ftracer[q], GetVflowplusVsink, &stateMT);
	    } else {
	      FlagRK4[q]=RK4(t, tmid, &ftracer[q], GetVflowplusVsink);
	    }
	    if(FlagRK4[q]==1) break;
	    
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
  }

  /****************************************************
   * COMPUTATION NUMERICAL DENSITY
   ****************************************************/
#ifdef DEBUG
  cout << " Numerical density " << endl;
#endif


  vector<int> numpoints(NumGrid,0);
  vector<double> density(NumGrid,0);
  int l;
  
omp_set_num_threads(20);
#pragma omp parallel for default(shared) private(l)// Parallelizing the code for computing trajectories
    
  for (int q=0; q<(int)NumTracers; q++) {
    
    if(FlagRK4[q]==1 || Flagfdepth[q]==0) continue;
    
    vectorXYZ delta;    
    delta.x=CHistParams.Radius/(rearth*cos(ftracer[q].y*rads));
    delta.y=CHistParams.Radius/rearth;
    delta.z=0.0;

    vectorXYZ min,max;

    max=ftracer[q]+(degrees*delta);
    min=ftracer[q]-(degrees*delta);
    
    int imin,imax,jmin,jmax;

    imin=((min.x-CHistParams.DomainBL.x)/CHistParams.InterSpacing.x)-5;
    if(imin<0) imin=0;
    if(imin>=GridDimx) continue;

    imax=((max.x-CHistParams.DomainBL.x)/CHistParams.InterSpacing.x)+5;
    if(imax>=GridDimx) imax=GridDimx-1;
    if(imax<0) continue;
    
    jmin=((min.y-CHistParams.DomainBL.y)/CHistParams.InterSpacing.y)-5;
    if(jmin<0) jmin=0;
    if(jmin>=GridDimy) continue;
    
    jmax=((max.y-CHistParams.DomainBL.y)/CHistParams.InterSpacing.y)+5;
    if(jmax>=GridDimy) jmax=GridDimy-1;
    if(jmax<0) continue;
    
    for(int j=jmin; j<jmax; j++){
      for(int i=imin; i<imax; i++){
	l=i+j*GridDimx;
	vectorXYZ alpha;
	
	alpha.x=cos(Grid[l].y*rads)*rads;
	alpha.y=rads;
	alpha.z=0.0;
	
	vectorXYZ beta;
	
	beta=ftracer[q]-Grid[l];
	beta*=alpha;
	double norm2beta=scalar(beta,beta);
	if(rearth*sqrt(norm2beta)<=CHistParams.Radius){
	  numpoints[l]++;
	}
      }
    }
  }
  
  // LAND RATIO

#ifdef DEBUG
  cout << " computation Land Ratio " << endl;
#endif

  int land;
  vector<double> LR(NumGrid,0);
  vectorXYZ testpoint;
  unsigned long s=123456789;
  init_genrand(s);
  for (unsigned int q=0; q<NumGrid; q++) {
    vectorXYZ delta;
    
    delta.x=CHistParams.Radius/(rearth*cos(Grid[q].y*rads));
    delta.y=CHistParams.Radius/rearth;
    delta.z=0.0;
    
    vectorXYZ max,min;
    
    max=Grid[q]+(degrees*delta);
    min=Grid[q]-(degrees*delta);
    
    vectorXYZ alpha;
	
    alpha.x=cos(Grid[q].y*rads)*rads;
    alpha.y=rads;
    alpha.z=0.0;
    land=0;
    l=0;
    while(l<CHistParams.Nparticles) {
      testpoint.x = RangeRand(min.x, max.x);
      testpoint.y = RangeRand(min.y, max.y);
      testpoint.z = fdepth;
    	
      vectorXYZ beta;
    
      beta=Grid[q]-testpoint;
      beta*=alpha;
      double norm2beta=scalar(beta,beta);
      if(rearth*sqrt(norm2beta)<=CHistParams.Radius){
	l++;
	land+=IsLand(testpoint);
      }
    }
    LR[q]=double(land)/CHistParams.Nparticles;
  }

  /* Free Velocities*/
  FreeMemoryVelocities(tau+4);

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
  for(unsigned int q=0; q<NumGrid; q++){
    density[q]=numpoints[q]/CHistParams.Nparticles;
    posfile<<Grid[q]<<" ";
    posfile<<LR[q]<<" ";
    posfile<<numpoints[q]<<" ";
    posfile<<density[q]<<endl;    
  }
  posfile.close();

  // REGULAR GRID Histogram
  string VTKffilename=rawfilename+".vtk";  
  ofstream offile(VTKffilename.c_str());
  
  offile<<"# vtk DataFile Version 3.0"<<endl;
  offile<<"Regular grid histogram"<<endl; 
  offile<<"ASCII"<<endl;
  offile<<"DATASET STRUCTURED_POINTS"<<endl;
  offile<<"DIMENSIONS "<<GridDimx<<" "<<GridDimy<<" 1"<<endl;
  offile<<"ORIGIN "<<CHistParams.DomainBL.x<<" "<<CHistParams.DomainBL.y<<" "<<CHistParams.DomainBL.z<<endl;
  offile<<"SPACINGS "<<CHistParams.InterSpacing.x<<" "<<CHistParams.InterSpacing.y<<" "<<CHistParams.InterSpacing.z <<endl;

  offile<<"POINT_DATA "<<NumGrid<<endl; 
  offile<<"SCALARS FactorDensity float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<density.size(); q++) {
      offile<<density[q]<<endl;
  }  
  offile<<"SCALARS NumericalDensity float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<numpoints.size(); q++) {
      offile<<numpoints[q]<<endl;
  }  
  offile<<"SCALARS LandRatio float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<LR.size(); q++) {
      offile<<LR[q]<<endl;
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
