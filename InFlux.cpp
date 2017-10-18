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
struct InFluxParameters {

  const vectorXYZ DomainBL;
  const vectorXYZ DomainTR;
  const vectorXYZ InterSpacing;
  const eqdate DepositionDate;
  const double TimeStep;
  const double Vsink;
  
  // Define a constructor that will load stuff from a configuration file.
  InFluxParameters(const string & InFluxParamsFileName)
  :DomainBL(getVectorXYZParam(InFluxParamsFileName, "DomainBL"))
  ,DomainTR(getVectorXYZParam(InFluxParamsFileName, "DomainTR"))
  ,InterSpacing(getVectorXYZParam(InFluxParamsFileName, "InterSpacing"))
  ,DepositionDate(getEqDateParam(InFluxParamsFileName, "DepositionDate"))
  ,TimeStep(getDoubleParam(InFluxParamsFileName, "TimeStep"))
  ,Vsink(getDoubleParam(InFluxParamsFileName, "Vsink"))
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
  cout << " Vsink "<<InFluxParams.Vsink<<endl;
#endif

  /********************************************************************************
   * GRID CONSTRUCTION
   ********************************************************************************/

  vector<vectorXYZ> grid;
  int GridDimx, GridDimy;
  int KcompNormal=2*(InFluxParams.TimeStep>0)-1;
  
  if(MakeRegular2DGrid(InFluxParams.DomainBL, 
		       InFluxParams.InterSpacing, 
		       InFluxParams.DomainTR,
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
  
  /********************************************************************************
   * SET UP TIME PARAMETERS
   ********************************************************************************/

  double StdVz;
  
  SetupVflow(SConfigurationFile);// SETUP VELOCITY MODEL: read parameters from file
  
  if(CalcStdVz(InFluxParams.DepositionDate, &StdVz)){
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

  int tau=2*int((InFluxParams.DomainTR.z-InFluxParams.DomainBL.z)/(InFluxParams.Vsink));
  
  eqdate inidate;
  double tstart;
  double tend;
  int ascnd = InFluxParams.TimeStep > 0;

  if(ascnd){
    unsigned int initime = InFluxParams.DepositionDate.GetDay();
    inidate.SetDay(initime-2);
    tstart = 2.0;
    tend = (double) tau;
  }
  else{
    unsigned int initime = InFluxParams.DepositionDate.GetDay()-tau;
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

  if((LoadVLatLonGrid(InFluxParams.DepositionDate))!=0){//Load velocity grid
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

  SetupLagrangianEngine(SConfigurationFile);
  
  vector<vectorXYZ> tracer1,tracer2;

  tracer1=grid;
  tracer2=grid;
  
   
  unsigned int numtracers=numgridpoints;
  vector<vectorXYZ> V1(numtracers),V2(numtracers);
  vector<vectorXYZ> dVdt1(numtracers),dVdt2(numtracers);
  vector<vectorXYZ> gradVz1(numtracers),gradVz2(numtracers);
  vector<vectorXYZ> gradh1(numtracers),gradh2(numtracers);
  
  vector<double> IntDivergence(numtracers),IntGradVz(numtracers), IntTimeTilt(numtracers); 
  vector<double> ftime(numtracers,0.0);
  vector<int> FlagRK4(numtracers,0);

  double fdepth=(InFluxParams.TimeStep<0)?InFluxParams.DomainTR.z:InFluxParams.DomainBL.z;
  vector<int> Flagfdepth(numtracers,0);
  

#ifdef DEBUG
  cout << " fdepth " << fdepth<< endl;
#endif
  

  /****************************************************
   * compute initial magnitudes
   ****************************************************/
  for (unsigned int q=0; q<numtracers; q++) {
       if(GetVflowplusVsink(tstart,tracer2[q],&V2[q])==1 
       || GradientVz(tstart,tracer2[q],&gradVz2[q])==1
       || TimeDerivativeV(tstart,tracer2[q], &dVdt2[q])==1){
      ftime[q]=double(tstart)+InFluxParams.TimeStep;
      FlagRK4[q]=1;      
      continue;
    }
    IntDivergence[q]=0.0;
    IntGradVz[q]=0.0;
    IntTimeTilt[q]=0.0;

    gradh1[q].x=0.0;
    gradh1[q].y=0.0;
    gradh1[q].z=0.0;
  }

  /****************************************************
   * TRACER LOOP
   ****************************************************/

  double t;
  
  omp_set_num_threads(5);
#pragma omp parallel for default(shared) private(t)// Parallelizing the code for computing trajectories
  for (unsigned int q=0; q<numtracers; q++) {
  /****************************************************
   * TIME LOOP
   ****************************************************/
  for(t=tstart; ascnd==1?(t<tend):(t>=tend); t+=InFluxParams.TimeStep) {
      
      //COMPUTE NEW POSITION
      tracer1[q]=tracer2[q];// it stores previous position in tracer1[q]
      if(RK4(t, InFluxParams.TimeStep, &tracer2[q], GetVflowplusVsink)==1){
	ftime[q]=double(t-tstart);
	FlagRK4[q]=1;
	break;
      }

      ftime[q]=t+InFluxParams.TimeStep;

      if(tracer2[q].z<=fdepth){
	if((tracer1[q].z-fdepth)<(fdepth-tracer2[q].z)){
	  tracer2[q]=tracer1[q];
	  ftime[q]-=InFluxParams.TimeStep;
	}
	Flagfdepth[q]=1;
        break;
      }

      //Store the previous value
      V1[q]=V2[q];
      gradVz1[q]=gradVz2[q];
      dVdt1[q]=dVdt2[q];
 
      // Computing variables
      if(GetVflowplusVsink(t+InFluxParams.TimeStep,tracer2[q],&V2[q])==1
	 || GradientVz(t+InFluxParams.TimeStep,tracer2[q],&gradVz2[q])==1
	 || TimeDerivativeV(t+InFluxParams.TimeStep,tracer2[q], &dVdt2[q])==1){
	ftime[q]=double(t-tstart)+InFluxParams.TimeStep;
	FlagRK4[q]=1;	
	break;
      }
      
      
      gradh1[q]=gradh2[q];//Store previous value
      
      gradh2[q].x=gradh1[q].x+((gradVz1[q].x+gradVz2[q].x)*(InFluxParams.TimeStep/2.0)); 
      gradh2[q].y=gradh1[q].y+((gradVz1[q].y+gradVz2[q].y)*(InFluxParams.TimeStep/2.0)); 
      gradh2[q].z=0.0;
      
      vectorXYZ Vnorm1=(1.0/V1[q].z)*V1[q];
      vectorXYZ Vnorm2=(1.0/V2[q].z)*V2[q];
           
      IntDivergence[q]+=(-gradVz1[q].z-gradVz2[q].z)*(InFluxParams.TimeStep/2.0);
      IntGradVz[q]+=(gradVz1[q].x*Vnorm1.x+gradVz1[q].y*Vnorm1.y+gradVz2[q].x*Vnorm2.x+gradVz2[q].y*Vnorm2.y)*(InFluxParams.TimeStep/2.0);
      

      vectorXYZ dVnorm1dt;
      vectorXYZ dVnorm2dt;

      dVnorm1dt.x=(1.0/V1[q].z)*dVdt1[q].x-(V1[q].x/(V1[q].z*V1[q].z))*dVdt1[q].z;
      dVnorm1dt.y=(1.0/V1[q].z)*dVdt1[q].y-(V1[q].y/(V1[q].z*V1[q].z))*dVdt1[q].z;
      dVnorm1dt.z=0.0;
      
      dVnorm2dt.x=(1.0/V2[q].z)*dVdt2[q].x-(V2[q].x/(V2[q].z*V2[q].z))*dVdt2[q].z;
      dVnorm2dt.y=(1.0/V2[q].z)*dVdt2[q].y-(V2[q].y/(V2[q].z*V2[q].z))*dVdt2[q].z;
      dVnorm2dt.z=0.0;
      
      double div1=scalar(gradh1[q],dVnorm1dt)/(scalar(gradh1[q],Vnorm1)-1.0);
      double div2=scalar(gradh2[q],dVnorm2dt)/(scalar(gradh2[q],Vnorm2)-1.0);
      
      IntTimeTilt[q]+=(div1+div2)*(InFluxParams.TimeStep/2.0);

  }
}



  /* Free Velocities*/

  FreeMemoryVelocities(tau+4);

  /**********************************************
   * WRITE RESULTS
   **********************************************/
  string VTKfilename="startgrid.vtk";  
  ofstream ofile(VTKfilename.c_str());
  
  ofile<<"# vtk DataFile Version 3.0"<<endl;
  ofile<<"Start grid"<<endl; 
  ofile<<"ASCII"<<endl;
  ofile<<"DATASET STRUCTURED_POINTS"<<endl;
  ofile<<"DIMENSIONS "<<GridDimx<<" "<<GridDimy<<" "<<1<<endl;
  ofile<<"ORIGIN "<<InFluxParams.DomainBL.x<<" "<<InFluxParams.DomainBL.y<<" "<<InFluxParams.DomainBL.z <<endl;
  ofile<<"SPACING "<<InFluxParams.InterSpacing.x<<" "<<InFluxParams.InterSpacing.y<<" "<<InFluxParams.InterSpacing.z <<endl;
  ofile<<"POINT_DATA "<<ftime.size()<<endl;
  ofile<<"SCALARS FinalTime float 1"<<endl;
  ofile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<ftime.size(); q++){
      ofile<<ftime[q]<<endl;
  }
  ofile<<"SCALARS FinalDepth float 1"<<endl;
  ofile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<tracer2.size(); q++){
      ofile<<tracer2[q].z<<endl;
  }
  ofile<<"VECTORS gradVz2 float"<<endl;
  //offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<tracer2.size(); q++){
    ofile<<tracer2[q]<<endl;
  }
  ofile<<"SCALARS intdivergence float 1"<<endl;
  ofile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<IntDivergence.size(); q++){
    ofile<<IntDivergence[q]<<endl;
  }

  ofile<<"SCALARS IntTimeTilt float 1"<<endl;
  ofile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<IntTimeTilt.size(); q++){
        ofile<<IntTimeTilt[q]<<endl;
  }

  /*ofile<<"VECTORS vel_vectors float"<<endl;
  //ofile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<V2.size(); q++){
    ofile<<V2[q]<<endl;
  }

  ofile<<"VECTORS dvdtnorm float"<<endl;
  //ofile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<dVdt2.size(); q++){
    
    vectorXYZ dVnorm2dt;
    
    dVnorm2dt.x=(1.0/V2[q].z)*dVdt2[q].x-(V2[q].x/(V2[q].z*V2[q].z))*dVdt2[q].z;
    dVnorm2dt.y=(1.0/V2[q].z)*dVdt2[q].y-(V2[q].y/(V2[q].z*V2[q].z))*dVdt2[q].z;
    dVnorm2dt.z=0.0;
    
    ofile<<dVnorm2dt<<endl;
  }
  ofile<<"VECTORS dvdt float"<<endl;
  //ofile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<dVdt2.size(); q++){
    ofile<<dVdt2[q]<<endl;
  }

  ofile<<"VECTORS gradh2 float"<<endl;
  //ofile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<gradh2.size(); q++){
    ofile<<gradh2[q]<<endl;
  }

  ofile<<"SCALARS IntGradVz float 1"<<endl;
  ofile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<IntGradVz.size(); q++){
    ofile<<IntGradVz[q]<<endl;
  }
  
  ofile<<"SCALARS IntTimeTilt float 1"<<endl;
  ofile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<IntTimeTilt.size(); q++){
    double intbuffer=exp(-IntTimeTilt[q]);
    if(isfinite(intbuffer)&& fabs(intbuffer)<1e2)
      ofile<< intbuffer <<endl;
    else
      ofile<< 1.0 <<endl;
  }
  
  ofile<<"SCALARS IntTotal float 1"<<endl;
  ofile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<IntDivergence.size(); q++){
    ofile<<exp(-(IntDivergence[q]+IntGradVz[q]+IntTimeTilt[q]))<<endl;
    }*/
  ofile<<"SCALARS flagDepth int 1"<<endl;
  ofile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<Flagfdepth.size(); q++) {
      ofile<<Flagfdepth[q]<<endl;
  }  
  ofile<<"SCALARS flagRK4 int 1"<<endl;
  ofile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<FlagRK4.size(); q++) {
      ofile<<FlagRK4[q]<<endl;
      }
  
  ofile.close();


  // FINAL GRID
  string VTKffilename="finalgrid.vtk";  
  ofstream offile(VTKffilename.c_str());
  
  offile<<"# vtk DataFile Version 3.0"<<endl;
  offile<<"Final grid"<<endl; 
  offile<<"ASCII"<<endl;
  offile<<"DATASET POLYDATA"<<endl;
  offile<<"POINTS "<<tracer2.size()<<" float"<<endl;
  for(unsigned int q=0; q<tracer2.size(); q++){
    offile<<tracer2[q].x<<" "<<tracer2[q].y<<" "<< fdepth <<endl;
  }
  offile<<"VERTICES "<<tracer2.size()<<" "<<tracer2.size()*2<<endl;
  for(unsigned int q=0; q<tracer2.size(); q++){
    offile<<"1 "<<q<<endl;
  }
  offile<<"POINT_DATA "<<ftime.size()<<endl;
  
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

  offile<<"SCALARS FinalTime float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<ftime.size(); q++){
      offile<<ftime[q]<<endl;
  }
  offile<<"SCALARS FinalDepth float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<tracer2.size(); q++){
      offile<<tracer2[q].z<<endl;
  }
  
  offile<<"SCALARS intdivergence float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<IntDivergence.size(); q++){
    offile<<IntDivergence[q]<<endl;
  }

  offile<<"SCALARS IntTimeTilt float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<IntTimeTilt.size(); q++){
        offile<<IntTimeTilt[q]<<endl;
  }

  offile<<"SCALARS IntGradVz float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<IntGradVz.size(); q++){
    offile<<IntGradVz[q]<<endl;
  }

  offile<<"VECTORS V2 float"<<endl;
  for(unsigned int q=0; q<V2.size(); q++){
    offile<<V2[q]<<endl;
  }

  offile<<"VECTORS gradVz2 float"<<endl;
  for(unsigned int q=0; q<gradVz2.size(); q++){
    offile<<gradVz2[q]<<endl;
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
