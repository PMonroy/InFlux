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
//#include "VTKdump.hpp"

string numprintf(int ndigits, int ndecimals, double number);
struct InFluxNumParameters {

  const vectorXYZ DomainBL;
  const vectorXYZ DomainTR;
  const vectorXYZ InterSpacing;
  const eqdate DepositionDate;
  const double TimeStep;
  const double Vsink;
  const int Random;
  const int NumPoints;
  
  // Define a constructor that will load stuff from a configuration file.
  InFluxNumParameters(const string & InFluxNumParamsFileName)
  :DomainBL(getVectorXYZParam(InFluxNumParamsFileName, "DomainBL"))
  ,DomainTR(getVectorXYZParam(InFluxNumParamsFileName, "DomainTR"))
  ,InterSpacing(getVectorXYZParam(InFluxNumParamsFileName, "InterSpacing"))
  ,DepositionDate(getEqDateParam(InFluxNumParamsFileName, "DepositionDate"))
  ,TimeStep(getDoubleParam(InFluxNumParamsFileName, "TimeStep"))
  ,Vsink(getDoubleParam(InFluxNumParamsFileName, "Vsink"))
  ,Random(getIntParam(InFluxNumParamsFileName, "Random"))
  ,NumPoints(getIntParam(InFluxNumParamsFileName, "NumPoints"))
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
  cout << "InFluxNum Parameters from command line: "<<endl;
  cout <<" Configuraton File "<< SConfigurationFile<<endl;
#endif
  
  const InFluxNumParameters InFluxNumParams(SConfigurationFile);

#ifdef DEBUG
  cout<<endl;
  cout << "InFluxNum Parameters from file "<<SConfigurationFile<<" :"<<endl; 
  cout << " DomainBL "<< InFluxNumParams.DomainBL<<endl;
  cout << " DomainTR "<< InFluxNumParams.DomainTR<<endl;
  cout << " InterSpacing " <<InFluxNumParams.InterSpacing<<endl;
  cout << " Deposition Date "<< InFluxNumParams.DepositionDate.GetMday()<<"-";
  cout<<InFluxNumParams.DepositionDate.GetMonth()<<"-";
  cout<<InFluxNumParams.DepositionDate.GetYear()<<endl;
  cout << " TimeStep "<<InFluxNumParams.TimeStep<<endl;
  cout << " Vsink "<<InFluxNumParams.Vsink<<endl;
#endif

  /********************************************************************************
   * GRID CONSTRUCTION
   ********************************************************************************/

  vector<vectorXYZ> itracer;
  int KcompNormal=2*(InFluxNumParams.TimeStep>0)-1;
  
  /*if(MakeRandom2DGrid(InFluxNumParams.DomainBL, 
		       InFluxNumParams.DomainTR,
		       KcompNormal,
		       InFluxNumParams.NumPoints,
		       &itracer)){//Grid construction 
    cout << "Random grid initial positions: [Fail]" << endl;
    return 1;
    }*/
  int TracerDimx, TracerDimy;
  if(MakeRegular2DGrid(InFluxNumParams.DomainBL,
		       InFluxNumParams.InterSpacing*(1.0/20),
		       InFluxNumParams.DomainTR,
		       KcompNormal,
		       &itracer,
		       &TracerDimx,
		       &TracerDimy)){//Grid construction 
    cout << "Random grid initial positions: [Fail]" << endl;
    return 1;
  }
 

#ifdef DEBUG
  cout<<endl;
  cout << "Grid Initial Position: "<<endl; 
  cout << " Num. Tracers "<< itracer.size() <<endl;
#endif

  /********************************************************************************
   * SET UP TIME PARAMETERS
   ********************************************************************************/

  double StdVz;
  
  SetupVflow(SConfigurationFile);// SETUP VELOCITY MODEL: read parameters from file
  
  if(CalcStdVz(InFluxNumParams.DepositionDate, &StdVz)){
    cout << " Calculation of Std Vz: [Fail]" << endl;
    return 1;
  }
  
#ifdef DEBUG
  cout << endl;
  cout << "SET UP TIME PARAMETERS: "<<endl; 
  cout << " Standard deviation of Vz at deposition time "<<StdVz<<endl;
#endif

  //int tau=1+int((InFluxParams.DomainTR.z-InFluxParams.DomainBL.z)/(InFluxParams.Vsink-StdVz));
  //if(tau<0 ){
  //tau=360;
  //}

  int tau=2*int((InFluxNumParams.DomainTR.z-InFluxNumParams.DomainBL.z)/(InFluxNumParams.Vsink));
  
  eqdate inidate;
  double tstart;
  double tend;
  int ascnd = InFluxNumParams.TimeStep > 0;

  if(ascnd){
    unsigned int initime = InFluxNumParams.DepositionDate.GetDay();
    inidate.SetDay(initime-2);
    tstart = 2.0;
    tend = (double) tau;
  }
  else{
    unsigned int initime = InFluxNumParams.DepositionDate.GetDay()-tau;
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
  cout << "LOAD VELOCITY FIELD" << endl; //Verbose: loading vel. grid  
  cout << " *Loading velocity grid "; //Verbose: loading vel. grid
#endif

  if((LoadVLatLonGrid(InFluxNumParams.DepositionDate))!=0){//Load velocity grid
      cout<<"[Fail]"<<endl;
      return 1;
  }

#ifdef DEBUG
  cout<<"[OK]" << endl;//Success loading vel. grid
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

  SetupLagrangianEngine(SConfigurationFile, InFluxNumParams.Random);
  
  vector<vectorXYZ> buffertracer,ftracer;

  buffertracer=itracer;
  ftracer=itracer;
  unsigned int numtracers=itracer.size();
  double fdepth=(InFluxNumParams.TimeStep<0)?InFluxNumParams.DomainTR.z:InFluxNumParams.DomainBL.z;
  vector<int> Flagfdepth(numtracers,0);
  vector<double> ftime(numtracers,0.0);
  vector<int> FlagRK4(numtracers,0);

  /****************************************************
   * TRACER LOOP
   ****************************************************/

  double t;
  
omp_set_num_threads(20);
#pragma omp parallel for default(shared) private(t)// Parallelizing the code for computing trajectories
 
  for (unsigned int q=0; q<numtracers; q++) {
    for(t=tstart; ascnd==1?(t<tend):(t>=tend); t+=InFluxNumParams.TimeStep) {
      //COMPUTE NEW POSITION
      buffertracer[q]=ftracer[q];// it stores previous position in tracer1[q]
      if(RK4(t, InFluxNumParams.TimeStep, &ftracer[q], GetVflowplusVsink)==1){
	ftracer[q]=buffertracer[q];
	ftime[q]=t;
	FlagRK4[q]=1;
	break;
      }
      
      // CHECK IF IT ARRIVE DEPTH 
      double f=(fdepth-buffertracer[q].z);
      double fmid=(fdepth-ftracer[q].z);
      double discriminant=f*fmid;   
      double dt;

      dt=InFluxNumParams.TimeStep;      
      if(discriminant<=0.0) {	
	double rtb;
	rtb=0.0;
	double tmid;
	for (int j=0;j<20;j++) {	  
	  ftracer[q]=buffertracer[q];
	  tmid=rtb+(dt*=0.5);
	  if(RK4(t, tmid, &ftracer[q], GetVflowplusVsink)==1){
	    ftracer[q]=buffertracer[q];
	    FlagRK4[q]=1;
	    break;
	  }
	  fmid=(fdepth-ftracer[q].z);
	  if (fmid <= 0.0) rtb=tmid;
	  if (dt < (InFluxNumParams.TimeStep/16.0) || fmid == 0.0){
	    dt=tmid;	    
	    Flagfdepth[q]=1;
	    break;
	  }
	}
	ftime[q]=t+dt;
	break;
      }
    }
  }



  /* Free Velocities*/

  FreeMemoryVelocities(tau+4);

 /**********************************************
  * GRID 
  **********************************************/


  //initial grid

  vector<vectorXYZ> igrid;
  int igridDimx, igridDimy;
  int iKcompNormal=2*(InFluxNumParams.TimeStep>0)-1;
  
  if(MakeRegular2DGrid(InFluxNumParams.DomainBL, 
		       InFluxNumParams.InterSpacing, 
		       InFluxNumParams.DomainTR,
		       iKcompNormal,
		       &igrid,
		       &igridDimx,
		       &igridDimy)){//Grid construction 
    cout << "Regular igrid initial positions: [Fail]" << endl;
    return 1;
  }
  
  // max and min of lon and lat tracer final position

  vectorXYZ tracermax=ftracer[0];
  vectorXYZ tracermin=ftracer[0];

  for (unsigned int q=0; q<numtracers; q++) {
    if(FlagRK4[q]!=1 && Flagfdepth[q]!=0){

      if(ftracer[q].x>tracermax.x) tracermax.x=ftracer[q].x;
      if(ftracer[q].y>tracermax.y) tracermax.y=ftracer[q].y;
      if(ftracer[q].z>tracermax.z) tracermax.z=ftracer[q].z;
      
      if(ftracer[q].x<tracermin.x) tracermin.x=ftracer[q].x;
      if(ftracer[q].y<tracermin.y) tracermin.y=ftracer[q].y;
      if(ftracer[q].z<tracermin.z) tracermin.z=ftracer[q].z;

    }
    
  }


  vector<vectorXYZ> fgrid;
  int fgridDimx, fgridDimy;
  int fKcompNormal=2*(InFluxNumParams.TimeStep<0)-1;
  
  if(MakeRegular2DGrid(tracermin, 
		       InFluxNumParams.InterSpacing, 
		       tracermax,
		       fKcompNormal,
		       &fgrid,
		       &fgridDimx,
		       &fgridDimy)){//Grid construction 
    cout << "Regular fgrid initial positions: [Fail]" << endl;
    return 1;
  }

  // find index grid of initial position and final and add one to nini and nend in the corresponding grid box
  
  vector<double> nini(igrid.size(),0);
  vector<double> nend(fgrid.size(),0);
  vector<int> Flagnini(igrid.size(),0);
  vector<int> FlagnendFlag(fgrid.size(),0);
  
  vectorXYZ rini,rend;
  int iini,iend;
  int jini,jend;    
  int qini,qend;
  
  for (unsigned int q=0; q<numtracers; q++){    
    // locate initial cell
    rini=(itracer[q]-InFluxNumParams.DomainBL)/InFluxNumParams.InterSpacing;
    iini=floor(rini.x);
    jini=floor(rini.y);
    if(iini>=igridDimx || iini<0) continue;
    if(jini>=igridDimy || jini<0) continue;
    qini=iini+(jini*igridDimx);

    
    // locate arrival cell
    rend=(ftracer[q]-tracermin)/InFluxNumParams.InterSpacing;
    iend=floor(rend.x);
    jend=floor(rend.y);

    if(iend>=fgridDimx || iend<0) continue;
    if(jend>=fgridDimy || jend<0) continue;
    qend=iend+(jend*fgridDimx);

    if(FlagRK4[q]==1 && Flagfdepth[q]==0){
      Flagnini[qini]++;
      Flagnini[qend]++;
    }else{
      nini[qini]++;
      nend[qend]++;
    }
  }

  /**********************************************
   * WRITE RESULTS
   **********************************************/
  string rawfilename;
  size_t lastdot = SConfigurationFile.find_last_of(".");
  if(lastdot == string::npos){
    rawfilename=SConfigurationFile;
  } else {
    rawfilename=SConfigurationFile.substr(0,lastdot);
  }

  
  // IGRID CLUSTERING
  string VTKfgridfile=rawfilename+".vtk";  
  ofstream fgridfile(VTKfgridfile.c_str());
  fgridfile<<"# vtk DataFile Version 3.0"<<endl;
  fgridfile<<"fgrid aggregation 2D"<<endl; 
  fgridfile<<"ASCII"<<endl;
  fgridfile<<"DATASET STRUCTURED_POINTS"<< endl;
  fgridfile<<"DIMENSIONS "<< fgridDimx+1 <<" "<< fgridDimy+1 <<" "<< 1 <<endl;
  fgridfile<<"ORIGIN "<<fgrid[0].x<<" "<<fgrid[0].y<<" "<< fgrid[0].z <<endl;
  fgridfile<<"SPACING "<<InFluxNumParams.InterSpacing.x<<" "<<InFluxNumParams.InterSpacing.y<<" "<< InFluxNumParams.InterSpacing.z <<endl;
  fgridfile<<"CELL_DATA "<<(fgridDimx)*(fgridDimy)<<endl;
  fgridfile<<"SCALARS nend int"<<endl;
  fgridfile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<nend.size(); q++) {
      fgridfile<<nend[q]<<endl;
  }
  fgridfile.close();
  
  /*
  // STARTING POSITIONS
  string VTKinitfile="initialpositions.vtk";  
  ofstream initfile(VTKinitfile.c_str());
  
  initfile<<"# vtk DataFile Version 3.0"<<endl;
  initfile<<"Final grid"<<endl; 
  initfile<<"ASCII"<<endl;
  initfile<<"DATASET POLYDATA"<<endl;
  initfile<<"POINTS "<<itracer.size()<<" float"<<endl;
  for(unsigned int q=0; q<itracer.size(); q++){
    initfile<<itracer[q]<<endl;
  }
  initfile<<"VERTICES "<<itracer.size()<<" "<<itracer.size()*2<<endl;
  for(unsigned int q=0; q<itracer.size(); q++){
    initfile<<"1 "<<q<<endl;
  }
  initfile.close();

  // ENDING POSITIONS
  string VTKendfile="finalpositions.vtk";  
  ofstream endfile(VTKendfile.c_str());
  
  endfile<<"# vtk DataFile Version 3.0"<<endl;
  endfile<<"Final grid"<<endl; 
  endfile<<"ASCII"<<endl;
  endfile<<"DATASET POLYDATA"<<endl;
  endfile<<"POINTS "<<ftracer.size()<<" float"<<endl;
  for(unsigned int q=0; q<ftracer.size(); q++){
    endfile<<ftracer[q]<<endl;
  }
  endfile<<"VERTICES "<<ftracer.size()<<" "<<ftracer.size()*2<<endl;
  for(unsigned int q=0; q<ftracer.size(); q++){
    endfile<<"1 "<<q<<endl;
  }
  endfile.close();*/
  return 0;
}
