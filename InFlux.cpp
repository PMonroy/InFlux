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
  cout << " Random "<<InFluxParams.Random<<endl;
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

  SetupLagrangianEngine(SConfigurationFile, InFluxParams.Random);
  
  vector<vectorXYZ> tracer1,tracer2;

  tracer1=grid;
  tracer2=grid;
     
  unsigned int numtracers=numgridpoints;
  vector<vectorXYZ> V1(numtracers),V2(numtracers);
  vector<vectorXYZ> Vinitial(numtracers),Vfinal(numtracers);
  vector<vectorXYZ> nu1(numtracers),nu2(numtracers);
  
  vector<double> logmu1(numtracers),logmu2(numtracers);
  vector<double> mu1(numtracers),mu2(numtracers);
  vector<vectorXYZ> UnitNormal(numtracers);
  

  vector<vectorXYZ> dVdt1(numtracers),dVdt2(numtracers);
  vector<vectorXYZ> gradVz1(numtracers),gradVz2(numtracers);
  
  vector<vectorXYZ> alpha1(numtracers),alpha2(numtracers);
  vector<vectorXYZ> beta1(numtracers),beta2(numtracers);
  vector<double> lambda1(numtracers),lambda2(numtracers);
  
  vector<double> logIntDiv(numtracers),logIntGradVz(numtracers), logIntH(numtracers);
  
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
    logIntDiv[q]=0.0;
    logIntGradVz[q]=0.0;
    logIntH[q]=0.0;

    nu2[q]=V2[q]/V2[q].z;

    UnitNormal[q].x=0.0;
    UnitNormal[q].y=0.0;
    UnitNormal[q].z=0.0;
    
    alpha2[q].x=dVdt2[q].x-nu2[q].x*dVdt2[q].z;
    alpha2[q].y=dVdt2[q].y-nu2[q].y*dVdt2[q].z;
    alpha2[q].z=0.0;

    logmu2[q]=0.0;
    mu2[q]=exp(logmu2[q]);

    beta2[q].x=0.0;
    beta2[q].y=0.0;
    beta2[q].z=mu2[q];
    
    lambda2[q]=scalar(alpha2[q],beta2[q])/scalar(V2[q],beta2[q]);


 
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
	if(Flagfdepth[q]==1) cout<<"SORPRESA1!!"<<endl;
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
      nu1[q]=nu2[q];
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
      nu2[q]=V2[q]/V2[q].z;
      
      //compute mu
      logmu1[q]=logmu2[q];
      mu1[q]=mu2[q];
      
      logmu2[q]=logmu1[q]-((gradVz1[q].z+gradVz2[q].z)*(dt/2.0));
      mu2[q]=exp(logmu2[q]);
      
      //compute beta(=gradh)
      beta1[q]=beta2[q];//Store previous value
      
      beta2[q].x=beta1[q].x-((mu1[q]*gradVz1[q].x+mu2[q]*gradVz2[q].x)*(dt/2.0)); 
      beta2[q].y=beta1[q].y-((mu1[q]*gradVz1[q].y+mu2[q]*gradVz2[q].y)*(dt/2.0)); 
      beta2[q].z=1.0;
      
      // compute alpha(=dvdt normalized with vz)
      
      alpha1[q]=alpha2[q];
      
      alpha2[q].x=dVdt2[q].x-nu2[q].x*dVdt2[q].z;
      alpha2[q].y=dVdt2[q].y-nu2[q].y*dVdt2[q].z;
      alpha2[q].z=0.0;

      // compute lambda

      lambda1[q]=lambda2[q];
      
      lambda2[q]=scalar(alpha2[q],factor(1.0/mu2[q],1.0/mu2[q],1.0)*beta2[q])/scalar(V2[q],factor(1.0/mu2[q],1.0/mu2[q],1.0)*beta2[q]);


      // compute logInts
      
      logIntGradVz[q]+=(scalar(gradVz1[q],nu1[q])+scalar(gradVz2[q],nu2[q]))*(dt/2.0);
      
      logIntDiv[q]+=(gradVz1[q].z+gradVz2[q].z)*(dt/2.0);
      
      logIntH[q]+=(lambda1[q]+lambda2[q])*(dt/2.0);
      
      if(Flagfdepth[q]==1){
	//Unit normal vector

	UnitNormal[q].z=1.0/(sqrt(1+scalar(factor(1.0/mu2[q],1.0/mu2[q],0.0)*beta2[q],factor(1.0/mu2[q],1.0/mu2[q],0.0)*beta2[q])));
	UnitNormal[q].x=(beta2[q].x/mu2[q])/UnitNormal[q].z;
	UnitNormal[q].y=(beta2[q].y/mu2[q])/UnitNormal[q].z;	
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

  string rawfilename;
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
    factordensity = exp(logIntGradVz[q]+logIntH[q]);
    if(isinf(factordensity) || isnan(factordensity)){
      posfile<<0.0<<" ";;
    } else {
      posfile<<factordensity<<" ";;
    }
    posfile<<UnitNormal[q]<<" ";
    posfile<<Vfinal[q]<<" ";
    posfile<<mu2[q]<<endl;
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
  
  offile<<"SCALARS logIntGradVz float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<logIntGradVz.size(); q++){
    if(isinf(logIntGradVz[q]) || isnan(logIntGradVz[q])){
      offile<<0.0<<endl;
    } else {
      offile<<logIntGradVz[q]<<endl;
    }
  }

  offile<<"SCALARS logIntDiv float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<logIntDiv.size(); q++){
  if(isinf(logIntDiv[q]) || isnan(logIntDiv[q])){
      offile<<0.0<<endl;
    } else {
    offile<<exp(logIntDiv[q])<<endl;
    }
  }

  offile<<"SCALARS logIntH float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<logIntH.size(); q++){
    if(isinf(logIntH[q]) || isnan(logIntH[q])){
      offile<<0.0<<endl;
    } else {
      offile<<logIntH[q]<<endl;
    }
  }

  offile<<"SCALARS mu float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<mu2.size(); q++){
      offile<<mu2[q]<<endl;
  }

  offile<<"VECTORS UnitNormal float"<<endl;
  for(unsigned int q=0; q<UnitNormal.size(); q++){
    offile<<UnitNormal[q]<<endl;
  }

  offile<<"SCALARS FactorDensity float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<logIntH.size(); q++){
    double factordensity;
    factordensity = exp(logIntGradVz[q]+logIntH[q]);
    if(isinf(factordensity) || isnan(factordensity)){
      offile<<0.0<<endl;
    } else {
      offile<<factordensity<<endl;
    }
  }

  offile<<"SCALARS FactorDensity2 float 1"<<endl;
  offile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<logIntH.size(); q++){
    double factordensity;
    factordensity = exp(logIntGradVz[q]);
    if(isinf(factordensity) || isnan(factordensity)){
      offile<<0.0<<endl;
    } else {
      offile<<factordensity<<endl;
    }
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
