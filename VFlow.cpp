#include <iostream>
#include <iomanip>
#include <netcdfcpp.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

#include "VFlow.hpp"
#include "Constants.hpp"
#include "EqDate.hpp"
#include "ParseParameters.hpp"

static const int NC_ERR = 2;

//GLOBAL VARIABLES
unsigned int DimLon;
unsigned int DimLonU;
unsigned int DimLat;
unsigned int DimLatV;
unsigned int DimDepth;
unsigned int DimTime;

vector <double> lon;
vector <double> lat;

vectorXYZ**** vflow;
double**** depth;

//char velocitydir[] = "/scratch/pmonroy/";
string VModelDir;

// Degree resolution
double DegreeResolution;


/* Structures: Esta estructura tiene que ser interna a velocity.cpp */
struct  vectorIJK {
    int i;
    int j;
    int k;
};
//FUNCTIONS PROTOTYPES
int GetIndices(unsigned long time, vectorXYZ point, vectorIJK *IntIndex, vectorXYZ *DecIndex);

//INLINE FUNCTIONS 
inline int LocateIndex(double x, double start, double step){
  return floor((x-start)/step);
}
inline int LocateIndex(double x, const vector <double> &xx){
  /* Given a vector xx[0...n-1], and given a value x, returns a value i such that x is between xx[i]
   * and xx[i+1]. xx must be monotonic, either increasing or decreasing. -1 or n is returned
   * to indicate that x is out of range.
   */

  if (x < xx.front()) 
    return -1;
  else if(x >= xx.back()) 
    return xx.size();
  else { 		
    unsigned long jl,jm,ju;
    jl=0;
    ju=xx.size()+1;
    int ascnd = (xx.back()>xx.front());
    while((ju-jl)>1) {
      jm =(ju+jl)>>1;
      if((x>=xx.at(jm-1))==ascnd) jl=jm;
      else ju=jm;
    }
    return int(jl-1);
  }
}

//FUNCTIONS
int SetupVflow(const string & vflowParamsFileName) {
  
  VModelDir=getStringParam(vflowParamsFileName, "VModelDir");
  DegreeResolution=getDoubleParam(vflowParamsFileName, "DegreeResolution");
  
  return 0;
}
int CalcStdVz(eqdate rdate, double *stdVz) {
  
  NcError err(NcError::verbose_nonfatal);

  // Open Netcdf file correspond to date rdate
  char ncfile[256];
  sprintf(ncfile, "%sextract_roms_avg_Y%1uM%1u.nc.1",VModelDir.c_str(), rdate.GetYear(),rdate.GetMonth());
  
  NcFile dataFile(ncfile, NcFile::ReadOnly);
  // Check to see if the file was opened.
  if(!dataFile.is_valid()){
    cout << "Error opening file:"<< ncfile <<endl; 
   return NC_ERR;
  }

  NcDim *NcDimLon;
  NcDim *NcDimLat;
  NcDim *NcDimDepth;

  if(!(NcDimLon = dataFile.get_dim("xi_rho"))
     || !(NcDimLat = dataFile.get_dim("eta_rho"))
     || !(NcDimDepth = dataFile.get_dim("s_rho"))){
    cout << "Error getting NcDimensions"<<endl;
    return NC_ERR;
  }

  vector <double> Vz(NcDimLon->size()*NcDimLat->size()*NcDimDepth->size());
  NcVar *NcVarVz;
  
  //Read Vz
  if (!(NcVarVz = dataFile.get_var("w"))
      || !NcVarVz->set_cur(rdate.GetMday()-1, 0, 0, 0)
      || !NcVarVz->get(&Vz[0], 1, NcDimDepth->size(), NcDimLat->size(), NcDimLon->size())){
    cout << "Error reading Vz variable"<<endl;
    return NC_ERR;
  }
  dataFile.close();


  // Compute Vz standard deviation

  double SumVz=0.0, SumVz2=0.0;
  int count=0;
  
  for(unsigned int q=0; q<Vz.size(); q++){
    // Vz[q]==0 corresponds to land condition
    if(Vz[q]!=0.0){
      count++;
      SumVz+=Vz[q];
      SumVz2+=(Vz[q]*Vz[q]);
    }
  }
  
  *stdVz=secondsday*sqrt((SumVz2/double(count))-(SumVz/double(count))*(SumVz/double(count)));
  
  return 0;
}
int LoadVLatLonGrid(eqdate rdate) {

  NcError err(NcError::verbose_nonfatal);

  // Open Netcdf file correspond to date rdate
  char ncfile[256];
  sprintf(ncfile, "%sextract_roms_avg_Y%1uM%1u.nc.1",VModelDir.c_str(), rdate.GetYear(),rdate.GetMonth());
  
  NcFile dataFile(ncfile, NcFile::ReadOnly);
  // Check to see if the file was opened.
  if(!dataFile.is_valid()){
    cout << "Error opening file:"<< ncfile <<endl; 
   return NC_ERR;
  }

  NcDim *NcDimLon;
  NcDim *NcDimLonU;
  NcDim *NcDimLat;
  NcDim *NcDimLatV;
  NcDim *NcDimDepth;
  NcDim *NcDimTime;

  if(!(NcDimLon = dataFile.get_dim("xi_rho"))
     || !(NcDimLonU = dataFile.get_dim("xi_u"))
     || !(NcDimLat = dataFile.get_dim("eta_rho"))
     || !(NcDimLatV = dataFile.get_dim("eta_v"))
     || !(NcDimDepth = dataFile.get_dim("s_rho"))
     || !(NcDimTime = dataFile.get_dim("time"))){
    cout << "Error getting NcDimensions"<<endl;
    return NC_ERR;
  }

  DimLon = NcDimLon->size();
  DimLonU = NcDimLonU->size();
  DimLat = NcDimLat->size();
  DimLatV = NcDimLatV->size();
  DimDepth = NcDimDepth->size();
  DimTime = NcDimTime->size();

  lon.resize(DimLon);
  lat.resize(DimLat);

  // Get latitude 
  NcVar *NcVarLon;
  if(!(NcVarLon = dataFile.get_var("lon_rho"))
     || !NcVarLon->set_cur(0,0)
     || !NcVarLon->get(&lon[0], 1, DimLon)){
    cout << "Error reading longitude variable"<<endl;
    return NC_ERR;
  }

  // Get longitude 
  NcVar *NcVarLat;
  if(!(NcVarLat = dataFile.get_var("lat_rho"))
     || !NcVarLat->set_cur(0,0) 
     || !NcVarLat->get(&lat[0], DimLat, 1)){
    cout << "Error reading latitude variable"<<endl;
    return NC_ERR;
  }

  return 0;
}
int LoadVelocities(eqdate startdate, int ntau) {

  unsigned int t,i,j,k;
  int ft_center,fk_center,fj_center,q_center;
  int fk_vmon,fj_vmon,q_vmon;
  int fj_umon,q_umon;

  char ncfile[256];

  NcError err(NcError::verbose_nonfatal);

  // Buffer variables

  vector <double> Hbuffer(DimTime*DimLon*DimLat*DimDepth);
  vector <double> Ubuffer(DimTime*DimLonU*DimLat*DimDepth);
  vector <double> Vbuffer(DimTime*DimLon*DimLatV*DimDepth);
  vector <double> Wbuffer(DimTime*DimLon*DimLat*DimDepth);

  NcVar *NcVarDepth, *NcVarU, *NcVarV, *NcVarW;
  unsigned int time, final_time;
  eqdate date;
  unsigned int ndays, sumndays=0;

  time = startdate.GetDay();
  date.SetDay(time);
  final_time = startdate.GetDay() + ntau;


  /* Allocate memory */
  vflow = new vectorXYZ ***[ntau];
  depth = new double ***[ntau];

  for(t=0; t<(unsigned int) ntau; t++){
    vflow[t] = new vectorXYZ **[DimLon];
    depth[t] = new double **[DimLon];
    for(i=0; i<DimLon; i++){
      vflow[t][i] = new vectorXYZ *[DimLat];
      depth[t][i] = new double *[DimLat];
      for(j=0; j<DimLat; j++){
	vflow[t][i][j] = new vectorXYZ [DimDepth];
	depth[t][i][j] = new double [DimDepth];
      }
    }
  }



  sumndays=0;
  while(time < final_time)
    {

      ndays = DimTime - (date.GetMday()-1);
      if((time + ndays) > final_time)
	ndays = final_time - time;
      
      sprintf(ncfile, "%sextract_roms_avg_Y%1uM%1u.nc.1",VModelDir.c_str(), date.GetYear(), date.GetMonth());
#ifdef DEBUG
      cout << "  NetCDF file " <<ncfile <<" start day="<< date.GetMday() <<" number of days=" << ndays;
#endif
      
      NcFile dataFile(ncfile, NcFile::ReadOnly);

      // Check to see if the file was opened.
      if(!dataFile.is_valid())
	return NC_ERR;
 
      //Read depth
      if (!(NcVarDepth = dataFile.get_var("depth"))
	  || !NcVarDepth->set_cur(date.GetMday()-1, 0, 0, 0)
	  || !NcVarDepth->get(&Hbuffer[0], ndays, DimDepth, DimLat, DimLon))
	return NC_ERR;

      //Read u
      if (!(NcVarU = dataFile.get_var("u"))
	  || !NcVarU->set_cur(date.GetMday()-1, 0, 0, 0)
	  || !NcVarU->get(&Ubuffer[0], ndays, DimDepth, DimLat, DimLonU))
	return NC_ERR;

      //Read v
      if (!(NcVarV = dataFile.get_var("v"))
	  || !NcVarV->set_cur(date.GetMday()-1, 0, 0, 0)
	  || !NcVarV->get(&Vbuffer[0], ndays, DimDepth, DimLatV, DimLon))
	return NC_ERR;

      //Read w
      if (!(NcVarW = dataFile.get_var("w"))
	  || !NcVarW->set_cur(date.GetMday()-1, 0, 0, 0)
	  || !NcVarW->get(&Wbuffer[0], ndays, DimDepth, DimLat, DimLon))
	return NC_ERR;

      dataFile.close();

      /* copy values and averages it */

      for(t=0; t<ndays; t++) {
	ft_center=t*DimDepth;
	
	for(k=0; k<DimDepth; k++){
	  fk_center=DimLat*(ft_center+k);
	  fk_vmon=DimLatV*(ft_center+k);
	    
	  for(j=0; j<DimLat; j++){
	    fj_center=DimLon*(fk_center+j);
	    fj_vmon=DimLon*(fk_vmon+j);
	    fj_umon=DimLonU*(fk_center+j);
	    for(i=0; i<DimLon; i++){
	      q_center=fj_center+i;
	      q_umon=fj_umon+i;
	      q_vmon=fj_vmon+i;
	      depth[t+sumndays][i][j][k]=Hbuffer[q_center];
	      vflow[t+sumndays][i][j][k].z=Wbuffer[q_center];
	      vflow[t+sumndays][i][j][k].x=(Ubuffer[q_umon-(i==(DimLon-1))]+Ubuffer[q_umon-(i!=0)])/2.0;
	      vflow[t+sumndays][i][j][k].y=(Vbuffer[q_vmon-(j==(DimLat-1))*DimLon]+Vbuffer[q_vmon-(j!=0)*DimLon])/2.0;

	      vflow[t+sumndays][i][j][k]=secondsday*vflow[t+sumndays][i][j][k];
	      
	    }
	  }
	}
      }

      time += ndays;
      date.SetDay(time);
      sumndays += ndays;
#ifdef DEBUG
      cout << " [OK]" <<endl;
#endif
    }

  return 0;
}

void FreeMemoryVelocities(int ntau) {
  int t,i,j;
  
  for(t=0; t<ntau; t++){
    for(i = 0; i < int(DimLon); i++){
      for(j = 0; j < int(DimLat); j++){
	delete[] depth[t][i][j];
	delete[] vflow[t][i][j];
      }
      delete[] depth[t][i];
      delete[] vflow[t][i]; 
    }
    delete[] depth[t];
    delete[] vflow[t]; 
  }
  delete[] depth;
  delete[] vflow;   
}

int GetIndices(unsigned long time, vectorXYZ point, vectorIJK *IntIndex, vectorXYZ *DecIndex) {

  /* point.x -> longitude
   * point.y -> latitude
   * point.z -> depth
   */

  int i,j,kbuffer, kmin, kmax;

  /* Locate index longitude*/
  IntIndex->i=LocateIndex(point.x, lon);

  if(IntIndex->i < 0 || IntIndex->i >= int(DimLon-1))
      return 1;

  DecIndex->x=(lon[IntIndex->i+1]-point.x)/(lon[IntIndex->i+1]-lon[IntIndex->i]);

  /* Locate index latitude*/
  IntIndex->j = LocateIndex(point.y, lat);

  if(IntIndex->j < 0 || IntIndex->j >= int(DimLat-1))
      return 1;
  // Transform to a more suitable coordinates
  
  double PointLat=log(fabs((1.0/cos(rads*point.y))+tan(rads*point.y)));
  double Lat0=log(fabs((1.0/cos(rads*lat[IntIndex->j]))+tan(rads*lat[IntIndex->j])));
  double Lat1=log(fabs((1.0/cos(rads*lat[IntIndex->j+1]))+tan(rads*lat[IntIndex->j+1])));
 
  DecIndex->y=(Lat1-PointLat)/(Lat1-Lat0);
  
  /* Locate index depth */
  kmax=0;
  kmin=DimDepth-1;

  int flagdepth=0;
  
  for(i=IntIndex->i; i<IntIndex->i+2; i++){
    for(j=IntIndex->j; j<IntIndex->j+2; j++)
      {
	vector<double> dpt;
	dpt.reserve (DimDepth);
	for (unsigned int k=0; k<DimDepth; k++){
	  dpt.push_back(depth[time][i][j][k]);
	}
	kbuffer = LocateIndex(point.z, dpt);
	if(kbuffer < 0 || kbuffer >= int(DimDepth-1))
	  flagdepth++;
	
	if(kbuffer<kmin && kbuffer>0) kmin=kbuffer;
	if(kbuffer>kmax && kbuffer<int(DimDepth)) kmax=kbuffer;
      }
  }

  if(flagdepth==4)
    return 1;
  
  //Calculation of hu and hl using kmin

  double hu=0.0, hl=0.0, phi;
  
  for(i=0; i<2; i++){
    for(j=0; j<2; j++){
      phi=((i>0)?DecIndex->x:(1-DecIndex->x));
      phi*=((j>0)?DecIndex->y:(1-DecIndex->y));
      hl+=depth[time][IntIndex->i+i][IntIndex->j+j][kmin]*phi;
      hu+=depth[time][IntIndex->i+i][IntIndex->j+j][kmin+1]*phi;
    }
  }
  
  
  if(kmin==kmax){ 
    IntIndex->k=kmin;
    DecIndex->z=(hu-point.z)/(hu-hl);    
  }
  else{
    while(!(point.z>=hl && point.z<hu)){
      if(point.z<hl){
	kmin=kmin-1;
	if(kmin<0||kmin>=int(DimDepth-1))
	  return 1;
	hu=hl; // and new calc of hl using kmin
	hl=0.0;
	for(i=0; i<2; i++){
	  for(j=0; j<2; j++){
	    phi=((i>0)?DecIndex->x:(1-DecIndex->x));
	    phi*=((j>0)?DecIndex->y:(1-DecIndex->y));
	    hl+=depth[time][IntIndex->i+i][IntIndex->j+j][kmin]*phi;
	  }
	}
      }
      if(point.z>=hu){
	kmin=kmin+1;
	if(kmin<0 || kmin>=int(DimDepth-1))
	  return 1;
	hl=hu; // and new calc of hu using kmin
	hu=0.0;
	for(i=0; i<2; i++){
	  for(j=0; j<2; j++){
	    phi=((i>0)?DecIndex->x:(1-DecIndex->x));
	    phi*=((j>0)?DecIndex->y:(1-DecIndex->y));
	    hu+=depth[time][IntIndex->i+i][IntIndex->j+j][kmin+1]*phi;
	  }
	}
      }
    }
    IntIndex->k=kmin;
    DecIndex->z=(hu-point.z)/(hu-hl);    
  }

  return 0;
}
int GetVelocity(double t,vectorXYZ point, vectorXYZ *vint) {

  vectorIJK IntIndex[2];
  vectorXYZ DecIndex[2];
  vectorXYZ vbuffer;
  unsigned long time;
  double alpha;
  int i,j,k,l;
  double phi;
  
  time=(unsigned long)t;
  alpha=((double)(time+1))-t;  
  
  if(GetIndices(time, point, &IntIndex[0], &DecIndex[0])==1)
    return 1;

  if(GetIndices(time+1, point, &IntIndex[1], &DecIndex[1])==1)  //increment 1 day,depend on vflow 
    return 1;

  vbuffer.x=0.0;
  vbuffer.y=0.0;
  vbuffer.z=0.0;

  for(l=0; l<2; l++){
    for(i=0; i<2; i++){
      for(j=0; j<2; j++){
	for(k=0; k<2; k++){
	  phi=((l>0)?alpha:(1-alpha));
	  phi*=((i>0)?DecIndex[l].x:(1-DecIndex[l].x));
	  phi*=((j>0)?DecIndex[l].y:(1-DecIndex[l].y));
	  phi*=((k>0)?DecIndex[l].z:(1-DecIndex[l].z));
	  vbuffer+=(vflow[time+l][IntIndex[l].i+i][IntIndex[l].j+j][IntIndex[l].k+k]*phi);
	}
      }
    }
  }
  *vint=vbuffer;

  return 0;
}
int GradientVz(double t,vectorXYZ point, vectorXYZ *grad){

  vectorIJK IntIndex[2];
  vectorXYZ DecIndex[2];
  vectorXYZ GradTilt;
  vectorXYZ GradRight;

  double  TanThetax,InvCosThetax,TanThetay,InvCosThetay;
  
  unsigned long time;

  double alpha;

  int i,di,i0,i1;
  int j,dj,j0,j1;
  int k,dk,k0,k1;
  int l,dl;

  double h1,h2;
  double deltalon,deltalat;
    
  double phit,philon,philat,phidepth;
  
  /* Calculates the vectors vb[15] and points ptmb[15] */
  time = (unsigned long) t;
  alpha=((double)(time+1))-t;  
  
  if(GetIndices(time, point, &IntIndex[0], &DecIndex[0])==1)
    return 1;

  if(GetIndices(time, point, &IntIndex[1], &DecIndex[1])==1)
    return 1;


  grad->x=0.0;
  grad->y=0.0;
  grad->z=0.0;
  
  for(dl=0; dl<2; dl++){
    l=time+dl;
    for(di=0; di<2; di++){
      i=IntIndex[dl].i+di;
      for(dj=0; dj<2; dj++){
	j=IntIndex[dl].j+dj;
	for(dk=0; dk<2; dk++){
	  k=IntIndex[dl].k+dk;

	  // Depth partial derivative:
	  k0=k-1;
	  if(k0<0||k0>int(DimDepth-1))
	    return 1;
	  
	  k1=k+1;	  
	  if(k1<0||k1>int(DimDepth-1))
	    return 1;

	  h1=depth[l][i][j][k]-depth[l][i][j][k0];
	  h2=depth[l][i][j][k1]-depth[l][i][j][k];
		  
	  GradTilt.z=(vflow[l][i][j][k1].z*(h1/h2)-vflow[l][i][j][k0].z*(h2/h1))/(h1+h2);
	  GradTilt.z=vflow[l][i][j][k].z*((1.0/h1)-(1.0/h2));

	  GradRight.z=GradTilt.z;
	  
	  // longitude partial derivative:
	  i0=i-1;
	  if(i0<0||i0>int(DimLon-1))
	    return 1;
	  
	  i1=i+1;
	  if(i1<0||i1>int(DimLon-1))
	    return 1;

	  deltalon=rearth*cos(rads*lat[j])*(rads*(lon[i1]-lon[i0]));	  
	  GradTilt.x=(vflow[l][i1][j][k].z-vflow[l][i0][j][k].z)/deltalon;

	  TanThetax=(depth[l][i1][j][k]-depth[l][i0][j][k])/deltalon;
	  InvCosThetax=sqrt(1.0+(TanThetax*TanThetax));

	  GradRight.x=(InvCosThetax*GradTilt.x)-(TanThetax*GradTilt.z);
	  
	  // latitude partial derivative:
	  j0=j-1;
	  if(j0<0||j0>int(DimLat-1))
	    return 1;
	  
	  j1=j+1;
	  if(j1<0||j1>int(DimLat-1))
	    return 1;

	  deltalat=rearth*rads*(lat[j1]-lat[j0]);
	  GradTilt.y=(vflow[l][i][j1][k].z-vflow[l][i][j0][k].z)/deltalat;

	  TanThetay=(depth[l][i][j1][k]-depth[l][i][j0][k])/deltalat;
	  InvCosThetay=sqrt(1.0+(TanThetay*TanThetay));

	  GradRight.y=(InvCosThetay*GradTilt.y)-(TanThetay*GradTilt.z);

	  // INTERPOL

	  phit=((dl>0)?alpha:(1-alpha));
	  philon=((di>0)?DecIndex[dl].x:(1-DecIndex[dl].x));
	  philat=((dj>0)?DecIndex[dl].y:(1-DecIndex[dl].y));
	  phidepth=((dk>0)?DecIndex[dl].z:(1-DecIndex[dl].z));
	  *grad+=(GradRight*(phit*philon*philat*phidepth));
	  
	}
      }
    }
  }
  
  return 0;
}
int TimeDerivativeV(double t,vectorXYZ point, vectorXYZ *dvdt){

  vectorIJK IntIndex[2];
  vectorXYZ DecIndex[2];
  vectorXYZ dvdtbuffer;

   unsigned long time;

  double alpha;

  int i,j,k,l;
  int di,dj,dk,dl;
  double phit,philon,philat,phidepth;
  
  
  /* Calculates the vectors vb[15] and points ptmb[15] */
  time = (unsigned long) t;
  alpha=((double)(time+1))-t;  
  
  if(GetIndices(time, point, &IntIndex[0], &DecIndex[0])==1)
    return 1;

  if(GetIndices(time, point, &IntIndex[1], &DecIndex[1])==1)
    return 1;


  dvdt->x=0.0;
  dvdt->y=0.0;
  dvdt->z=0.0;
  
  
  for(dl=0; dl<2; dl++){
    l=time+dl;
    for(di=0; di<2; di++){
      i=IntIndex[dl].i+di;
      for(dj=0; dj<2; dj++){
	j=IntIndex[dl].j+dj;
	for(dk=0; dk<2; dk++){
	  k=IntIndex[dl].k+dk;

	  dvdtbuffer=0.5*(vflow[l+1][i][j][k]-vflow[l-1][i][j][k]);
	  
	  // INTERPOL
	  phit=((dl>0)?alpha:(1-alpha));
	  philon=((di>0)?DecIndex[dl].x:(1-DecIndex[dl].x));
	  philat=((dj>0)?DecIndex[dl].y:(1-DecIndex[dl].y));
	  phidepth=((dk>0)?DecIndex[dl].z:(1-DecIndex[dl].z));
	  *dvdt+=(dvdtbuffer*(phit*philon*philat*phidepth));
	}
      }
    }
  }
  
  return 0;
}
