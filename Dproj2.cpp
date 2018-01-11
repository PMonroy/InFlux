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
  
  if(GetcmdlineParameters(argc, argv, &SConfigurationFile)) {//Get cmd line parameters 
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
#endif

  /********************************************************************************
   * GRID CONSTRUCTION
   ********************************************************************************/

  vector<vectorXYZ> grid;
  int GridDimx, GridDimy;
  int KcompNormal=2*(DprojParams.TimeStep>0)-1;
  
  if(MakeRegular2DGrid(DprojParams.DomainBL, 
		       DprojParams.InterSpacing, 
		       DprojParams.DomainTR,
		       KcompNormal,
		       &grid,
		       &GridDimx,
		       &GridDimy)){//Grid construction 
    cout << "Regular grid initial positions: [Fail]" << endl;
    return 1;
  }
  unsigned int numgridpoints=grid.size();

#ifdef DEBUG
  cout<<endl;
  cout << "Grid Initial Position: "<<endl; 
  cout << " Num. Tracers "<< grid.size() <<endl;
  cout << " Dim x "<< GridDimx << endl;
  cout << " Dim y "<< GridDimy << endl;
  cout << " depth " << grid[0].z<< endl;
#endif


  // Satellite grid construction  
  vector<vectorXYZ> satgrid;
  if(SatelliteGrid(DprojParams.DistSat, &grid, &satgrid)){
    cout << "Satellite grid of initial positions: [Fail]" << endl;
    return 1;
  }

  
  
  
  /********************************************************************************
   * SET UP TIME PARAMETERS
   ********************************************************************************/

  double StdVz;
  
  SetupVflow(SConfigurationFile);// SETUP VELOCITY MODEL: read parameters from file
  
  if(CalcStdVz(DprojParams.DepositionDate, &StdVz)){
    cout << " Calculation of Std Vz: [Fail]" << endl;
    return 1;
  }
  
#ifdef DEBUG
  cout << endl;
  cout << "SET UP TIME PARAMETERS: "<<endl; 
  cout << " Standard deviation of Vz at deposition time "<<StdVz<<endl;
#endif

  //int tau=1+int((DprojParams.DomainTR.z-DprojParams.DomainBL.z)/(DprojParams.Vsink-StdVz));
  //if(tau<0 ){
  //tau=360;
  //}

  int tau=3*int((DprojParams.DomainTR.z-DprojParams.DomainBL.z)/(DprojParams.Vsink));
  
  eqdate inidate;
  double tstart;
  double tend;
  int ascnd = DprojParams.TimeStep > 0;

  if(ascnd){
    unsigned int initime = DprojParams.DepositionDate.GetDay();
    inidate.SetDay(initime-2);
    tstart = 2.0;
    tend = (double) tau;
  }
  else{
    unsigned int initime = DprojParams.DepositionDate.GetDay()-tau;
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

if((LoadVLatLonGrid(DprojParams.DepositionDate))!=0){//Load velocity grid
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

  SetupLagrangianEngine(SConfigurationFile, DprojParams.Random);
  
  vector<vectorXYZ> tracer1,tracer2;
  vector<vectorXYZ> sattracer1,sattracer2;

  tracer1=grid;
  tracer2=grid;

  sattracer1=satgrid;
  sattracer2=satgrid;

   
  unsigned int numtracers=numgridpoints;

  vector<double> ftime(numtracers,0.0);
  vector<int> FlagRK4(numtracers,0);

  double fdepth=(DprojParams.TimeStep<0)?DprojParams.DomainTR.z:DprojParams.DomainBL.z;
  vector<int> Flagfdepth(numtracers,0);

  vector<vectorXYZ> Vtracer;
  vector<vectorXYZ> UnitNormal;
  vector<vectorXYZ> TangentX;
  vector<vectorXYZ> TangentY;
  
  Vtracer.reserve(numtracers);
  UnitNormal.reserve(numtracers);
  TangentX.reserve(numtracers);
  TangentY.reserve(numtracers);
  for(unsigned int q=0; q<numtracers; q++){
    Vtracer.push_back(vectorXYZ(0,0,0));
    UnitNormal.push_back(vectorXYZ(0,0,1));
    TangentX.push_back(vectorXYZ(DprojParams.DistSat,0,0));
    TangentY.push_back(vectorXYZ(0,DprojParams.DistSat,0));
  }
  

  vector<double> StretchingArea(numtracers,1.0);
  vector<double> StretchingLength(numtracers,1.0);
  vector<double> FactorDensity(numtracers,0);

#ifdef DEBUG
  cout << " fdepth " << fdepth<< endl;
#endif
  
  /****************************************************
   * TRACER LOOP
   ****************************************************/

  double t,period;
  
  omp_set_num_threads(20);
#pragma omp parallel for default(shared) private(t,period)// Parallelizing the code for computing trajectories
  for (unsigned int q=0; q<numtracers; q++) {
    period=0.0;
    for(t=tstart; ascnd==1?(t<tend):(t>=tend); t+=DprojParams.TimeStep) {
      
      //COMPUTE NEW POSITION
      tracer1[q]=tracer2[q];// it stores previous position in tracer1[q]
      sattracer1[4*q]=sattracer2[4*q];// it stores previous position in sattracer1[q]
      sattracer1[4*q+1]=sattracer2[4*q+1];// it stores previous position in sattracer1[q]
      sattracer1[4*q+2]=sattracer2[4*q+2];// it stores previous position in sattracer1[q]
      sattracer1[4*q+3]=sattracer2[4*q+3];// it stores previous position in sattracer1[q]
 
      if(RK4(t, DprojParams.TimeStep, &tracer2[q], GetVflowplusVsink)==1 ||
	 RK4(t, DprojParams.TimeStep, &sattracer2[4*q], GetVflowplusVsink)==1 ||
	 RK4(t, DprojParams.TimeStep, &sattracer2[4*q+1], GetVflowplusVsink)==1 ||
	 RK4(t, DprojParams.TimeStep, &sattracer2[4*q+2], GetVflowplusVsink)==1 ||
	 RK4(t, DprojParams.TimeStep, &sattracer2[4*q+3], GetVflowplusVsink)==1){
	ftime[q]=t;
	FlagRK4[q]=1;
	break;
      }          
      ftime[q]=t+DprojParams.TimeStep;

      double f=(fdepth-tracer1[q].z);
      double fmid=(fdepth-tracer2[q].z);
      double discriminant=f*fmid;   

      double NormCross;
      vectorXYZ ScaleFactor,deltaX,deltaY;
      
      if(discriminant<=0.0) {
	double rtb,dt;
	rtb=0.0;
	dt=DprojParams.TimeStep;
	double tmid;
	for (int j=0;j<20;j++) {
	  
	  if(Flagfdepth[q]==1) break;
	  tracer2[q]=tracer1[q];
	  sattracer2[4*q]=sattracer1[4*q];
	  sattracer2[4*q+1]=sattracer1[4*q+1];
	  sattracer2[4*q+2]=sattracer1[4*q+2];
	  sattracer2[4*q+3]=sattracer1[4*q+3];
	  tmid=rtb+(dt*=0.5);
	  if(RK4(t, tmid, &tracer2[q], GetVflowplusVsink)==1 ||
	     RK4(t, tmid, &sattracer2[4*q], GetVflowplusVsink)==1 ||
	     RK4(t, tmid, &sattracer2[4*q+1], GetVflowplusVsink)==1 ||
	     RK4(t, tmid, &sattracer2[4*q+2], GetVflowplusVsink)==1 ||
	     RK4(t, tmid, &sattracer2[4*q+3], GetVflowplusVsink)==1){
	    FlagRK4[q]=1;
	    break;
	  }
	  fmid=(fdepth-tracer2[q].z);
	  if (fmid <= 0.0) rtb=tmid;
	  if (dt < (DprojParams.TimeStep/16.0) || fmid == 0.0){
	    ftime[q]=t+tmid;
	    Flagfdepth[q]=1;
	    // Computing variables IF THE PARTILCE REACH THE DEPTH
	    if(GetVflowplusVsink(ftime[q],tracer2[q],&Vtracer[q])==1){
	      FlagRK4[q]=1;
	      break;
	    }      
	  }
	}
      }
      period+=DprojParams.TimeStep;
      if(period>=DprojParams.Period || Flagfdepth[q]==1){
	period=0.0;
	//COMPUTE ALL THE STUFS      
	TangentX[q].x=rearth*cos(rads*tracer2[q].y)*((sattracer2[4*q].x-sattracer2[4*q+2].x)*(rads/2.0)); 
	TangentX[q].y=rearth*((sattracer2[4*q].y-sattracer2[4*q+2].y)*(rads/2.0)); 
	TangentX[q].z=(sattracer2[4*q].z-sattracer2[4*q+2].z)/2.0;
	
	TangentY[q].x=rearth*cos(rads*tracer2[q].y)*((sattracer2[4*q+1].x-sattracer2[4*q+3].x)*(rads/2.0)); 
	TangentY[q].y=rearth*((sattracer2[4*q+1].y-sattracer2[4*q+3].y)*(rads/2.0)); 
	TangentY[q].z=(sattracer2[4*q+1].z-sattracer2[4*q+3].z)/2.0; 
      
	UnitNormal[q]=cross(TangentX[q],TangentY[q]);
	NormCross=sqrt(scalar(UnitNormal[q],UnitNormal[q]));
	UnitNormal[q]=(1.0/NormCross)*UnitNormal[q];
	
	StretchingArea[q]*=(NormCross/(DprojParams.DistSat*DprojParams.DistSat));
	double normX=sqrt(scalar(TangentX[q],TangentX[q]));
	double normY=sqrt(scalar(TangentY[q],TangentY[q]));
	vectorXYZ UnitTangentX=(1.0/normX)*TangentX[q];
	vectorXYZ UnitTangentY=(1.0/normY)*TangentY[q];
	vectorXYZ UnitTangentMax;

	if(normX>=normY){
	  StretchingLength[q]*=(normX/DprojParams.DistSat);
	  UnitTangentMax=UnitTangentX;
	}else{
	  StretchingLength[q]*=(normY/DprojParams.DistSat);
	  UnitTangentMax=UnitTangentY;
	}
	if(Flagfdepth[q]==1){
	  double normV=sqrt(scalar(Vtracer[q],Vtracer[q]));
	  vectorXYZ UnitV=(1.0/normV)*Vtracer[q];
	  vectorXYZ nv=cross(UnitNormal[q],UnitV);
	  double cosphi=scalar(UnitTangentMax,nv);
	  cosphi=0.0;
	  //double alpha=Vtracer[q].z/scalar(nv,Vtracer[q]);
	  double alpha=(Vtracer[q].z/normV)/((scalar(UnitNormal[q],Vtracer[q])/normV)*(scalar(UnitNormal[q],Vtracer[q])/normV)-1.0);
	  double beta;
	  beta=((StretchingArea[q]*cosphi)/StretchingLength[q]);
	  beta=beta*beta;
	  beta=beta+(StretchingLength[q]*StretchingLength[q])*(1.0-(cosphi*cosphi));
	  beta=sqrt(beta);
	  FactorDensity[q]=fabs(alpha)*beta;
	  break;
	}

	//ACTUALIZE VARIABLES

	if(normX>normY){
	  TangentX[q]=(DprojParams.DistSat/normX)*TangentX[q];
	  TangentY[q]=cross(UnitNormal[q],TangentX[q]);
	}else {
	  TangentX[q]=(DprojParams.DistSat/normY)*TangentY[q];
	  TangentY[q]=cross(UnitNormal[q],TangentX[q]);
	}
	
	ScaleFactor.x=degrees/(rearth*cos(rads*tracer2[q].y));
	ScaleFactor.y=degrees/rearth;
	ScaleFactor.z=1.0;
	
	deltaX=ScaleFactor*TangentX[q];
	deltaY=ScaleFactor*TangentY[q];
	
	sattracer2[4*q]=tracer2[q]+deltaX;
	sattracer2[4*q+1]=tracer2[q]+deltaY;
	sattracer2[4*q+2]=tracer2[q]-deltaX;
	sattracer2[4*q+3]=tracer2[q]-deltaY;
      }
    }
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
  for(unsigned int q=0; q<tracer2.size(); q++){
    posfile<<tracer2[q]<<" ";
    posfile<<Flagfdepth[q]<<" ";
    posfile<<FlagRK4[q]<<" ";
    posfile<<(isnan(FactorDensity[q])?-1.0:FactorDensity[q])<<endl;
  }
  posfile.close();
  
  
  // DATA FILE

  string datafilename=rawfilename+".dat";  
  ofstream datafile(datafilename.c_str());
  for(unsigned int q=0; q<tracer2.size(); q++){
    datafile<<grid[q]<<" ";
    datafile<<tracer2[q]<<" ";
    datafile<<Flagfdepth[q]<<" ";
    datafile<<FlagRK4[q]<<" ";
    datafile<<ftime[q]<<" ";
    datafile<<Vtracer[q]<<" ";
    datafile<<UnitNormal[q]<<" ";
    datafile<<(isnan(FactorDensity[q])?-1.0:FactorDensity[q])<<endl;
  }
  datafile.close();

  // VTK FILE
  string VTKfilename=rawfilename+".vtk";  
  ofstream offile(VTKfilename.c_str());
  
  offile<<"# vtk DataFile Version 3.0"<<endl;
  offile<<"Final grid"<<endl; 
  offile<<"ASCII"<<endl;
  offile<<"DATASET POLYDATA"<<endl;
  offile<<"POINTS "<<tracer2.size()<<" float"<<endl;
  for(unsigned int q=0; q<tracer2.size(); q++){
    offile<<tracer2[q].x<<" "<<tracer2[q].y<<" "<< DprojParams.DomainBL.z <<endl;
  }
  offile<<"VERTICES "<<tracer2.size()<<" "<<tracer2.size()*2<<endl;
  for(unsigned int q=0; q<tracer2.size(); q++){
    offile<<"1 "<<q<<endl;
  }
  offile<<"POINT_DATA "<<tracer2.size()<<endl;
  
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

  offile<<"SCALARS EndTime float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<ftime.size(); q++) {
      offile<<ftime[q]<<endl;
  }

  offile<<"SCALARS EndDepth float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<ftime.size(); q++) {
      offile<<tracer2[q].z<<endl;
  }

  offile<<"VECTORS Velocity float"<<endl;
  for(unsigned int q=0; q<Vtracer.size(); q++){
    offile<<Vtracer[q]<<endl;
  }

  offile<<"VECTORS UnitNormal float"<<endl;
  for(unsigned int q=0; q<UnitNormal.size(); q++){
    offile<<UnitNormal[q]<<endl;
  }
  
  offile<<"SCALARS FactorDensity float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<FactorDensity.size(); q++){
    if(isnan(FactorDensity[q]) || isinf(FactorDensity[q])){
	offile<<0.0<<endl;
    }else{
	offile<<FactorDensity[q]<<endl;
    }  
  }
  
  offile.close();

  
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
