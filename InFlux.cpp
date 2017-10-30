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
  vector<vectorXYZ> Vinitial(numtracers),Vfinal(numtracers);
  vector<vectorXYZ> dVdt1(numtracers),dVdt2(numtracers);
  vector<vectorXYZ> gradVz1(numtracers),gradVz2(numtracers);
  vector<vectorXYZ> gradh1(numtracers),gradh2(numtracers);
  
  vector<double> IntDivergence(numtracers),IntGradVz(numtracers), IntTimeTilt(numtracers),IntTimeTilt2(numtracers); 
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
    IntTimeTilt2[q]=0.0;

    gradh1[q].x=0.0;
    gradh1[q].y=0.0;
    gradh1[q].z=0.0;
  }

  Vinitial=V2;
  
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
	ftime[q]=t;
	FlagRK4[q]=1;
	break;
      }
      


      // CHECK IF IT ARRIVE DEPTH 
      double f=(fdepth-tracer1[q].z);
      double fmid=(fdepth-tracer2[q].z);
      double discriminant=f*fmid;   
      double dt;

      dt=InFluxParams.TimeStep;      
      if(discriminant<=0.0) {	
	double rtb;
	rtb=0.0;
	double tmid;
	for (int j=0;j<20;j++) {	  
	  if(Flagfdepth[q]==1) break;
	  tracer2[q]=tracer1[q];
	  tmid=rtb+(dt*=0.5);
	  if(RK4(t, tmid, &tracer2[q], GetVflowplusVsink)==1){
	    FlagRK4[q]=1;
	    break;
	  }
	  fmid=(fdepth-tracer2[q].z);
	  if (fmid <= 0.0) rtb=tmid;
	  if (dt < (InFluxParams.TimeStep/16.0) || fmid == 0.0){
	    ftime[q]=t+tmid;
	    Flagfdepth[q]=1;
	  }
	}
      }
      ftime[q]=t+dt;
      
      //Store the previous value
      V1[q]=V2[q];
      gradVz1[q]=gradVz2[q];
      dVdt1[q]=dVdt2[q];
      
      // Computing variables
      if(GetVflowplusVsink(t+dt,tracer2[q],&V2[q])==1
	 || GradientVz(t+dt,tracer2[q],&gradVz2[q])==1
	 || TimeDerivativeV(t+dt,tracer2[q], &dVdt2[q])==1){
	ftime[q]=t+dt;
	FlagRK4[q]=1;	
	break;
      }
      
      
      gradh1[q]=gradh2[q];//Store previous value
      gradh2[q].x=gradh1[q].x+((gradVz1[q].x+gradVz2[q].x)*(dt/2.0)); 
      gradh2[q].y=gradh1[q].y+((gradVz1[q].y+gradVz2[q].y)*(dt/2.0)); 
      gradh2[q].z=0.0;
	
      vectorXYZ Vnorm1=(1.0/V1[q].z)*V1[q];
      vectorXYZ Vnorm2=(1.0/V2[q].z)*V2[q];
      
      IntDivergence[q]+=(gradVz1[q].z+gradVz2[q].z)*(dt/2.0);
      IntGradVz[q]+=(gradVz1[q].x*Vnorm1.x+gradVz1[q].y*Vnorm1.y+gradVz2[q].x*Vnorm2.x+gradVz2[q].y*Vnorm2.y)*(dt/2.0);
      
      vectorXYZ dVnorm1dt;
      vectorXYZ dVnorm2dt;
      
      dVnorm1dt.x=(1.0/V1[q].z)*dVdt1[q].x-(V1[q].x/(V1[q].z*V1[q].z))*dVdt1[q].z;
      dVnorm1dt.y=(1.0/V1[q].z)*dVdt1[q].y-(V1[q].y/(V1[q].z*V1[q].z))*dVdt1[q].z;
      dVnorm1dt.z=0.0;
      
      dVnorm2dt.x=(1.0/V2[q].z)*dVdt2[q].x-(V2[q].x/(V2[q].z*V2[q].z))*dVdt2[q].z;
      dVnorm2dt.y=(1.0/V2[q].z)*dVdt2[q].y-(V2[q].y/(V2[q].z*V2[q].z))*dVdt2[q].z;
      dVnorm2dt.z=0.0;
      
      double div1=scalar(gradh1[q],dVnorm1dt)/(1.0-scalar(gradh1[q],Vnorm1));
      double div2=scalar(gradh2[q],dVnorm2dt)/(1.0-scalar(gradh2[q],Vnorm2));
      
      IntTimeTilt[q]+=(div1+div2)*(dt/2.0);

      // NEWWWWW

      double numerator1;
      double numerator2;

      numerator1 = -gradh1[q].x*dVdt1[q].x-gradh1[q].y*dVdt1[q].y+dVdt1[q].z;
      numerator2 = -gradh2[q].x*dVdt2[q].x-gradh2[q].y*dVdt2[q].y+dVdt2[q].z;

      double denominator1;
      double denominator2;

      denominator1 = -gradh1[q].x*V1[q].x-gradh1[q].y*V1[q].y+V1[q].z;
      denominator2 = -gradh2[q].x*V2[q].x-gradh2[q].y*V2[q].y+V2[q].z;

      IntTimeTilt2[q]+=((numerator1/denominator1)+(numerator2/denominator2))*(dt/2.0);

      if(Flagfdepth[q]==1){
	Vfinal[q]=V2[q];
	break;
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
    double factordensity;
    factordensity = (Vfinal[q].z/Vinitial[q].z)/exp(IntTimeTilt2[q]);
    
    posfile<<(isnan(factordensity)?0.0:abs(factordensity))<<endl;
  }
  posfile.close();
    
  // FINAL GRID
  string VTKffilename=rawfilename+".vtk";  
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
    offile<<(isnan(IntTimeTilt[q])?0.0:IntTimeTilt[q])<<endl;
  }

  offile<<"SCALARS IntTimeTilt2 float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<IntTimeTilt2.size(); q++){
    offile<<(isnan(IntTimeTilt2[q])?0.0:IntTimeTilt2[q])<<endl;
  }

  offile<<"SCALARS FactorDensity float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<IntTimeTilt2.size(); q++){
    double factordensity;
    factordensity = (Vfinal[q].z/Vinitial[q].z)/exp(IntTimeTilt2[q]);
    
    offile<<(isnan(factordensity)?0.0:abs(factordensity))<<endl;
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
