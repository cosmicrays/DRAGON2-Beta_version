/**
 * @file source.cc
 * @author Luca Maccione
 * @email luca.maccione@lmu.de
 * @brief In this file all the classes related to the model of the galaxy are implemented.
 */

#include "input.h"
#include "sources.h"
#include "grid.h"
#include "constants.h"
#include "geometry.h"

#include <fstream>
#include <iostream> 
#include <numeric>
#include <cstdlib> 
using namespace std;

// SN SOURCES

double TAstrophysicalSource::SourceDistribution(double x, double y, double z) {
  
  double radius = sqrt(x*x+y*y);
  double theta    = atan2(y,x) + M_PI;
  
  double Ferriere_distribution;
  double spiral_beta = 3.53;
  double spiral_z_h = 1.00;
  double robs = in->robs;
  

  //cout <<"$$$ SourceDistribution " << " " << SNR_model << " " << x << " " << y << " " << z << endl;

  switch (SNR_model) {
    
  case Lorimer: 
    // parametrization based on PSR catalogue: Lorimer et al., Mon. Not. R. Astron. Soc. 372, 777 (2006)
    if (radius >= 15.0) return 0.0;
    else return pow ( radius/robs , 1.9 ) * exp ( -5.00*(radius-robs)/robs - fabs(z)/0.2  );
    break;
    
  case Ferriere: 
    // parametrization based on PSR catalogue + disk stars: K. Ferriere, Rev.Mod.Phys. 73, 1031-1066 (2001)
    if (radius > 3.7) return (50.0*(0.79*exp(-pow(z/0.212,2.))+0.21*exp(-pow(z/0.636,2.)))*exp(-((radius)*(radius)-robs*robs)/(6.8*6.8)) + 7.3*exp(-(radius-robs)/4.5-fabs(z)/0.325));//)*pow(sin(2.*3.14*radius/16.+0.),2.);
    else return (177.5*(0.79*exp(-pow(z/0.212,2.))+0.21*exp(-pow(z/0.636,2.)))*exp(-pow((radius-3.7)/2.1, 2.)) + 7.3*exp(-(radius-robs)/4.5-fabs(z)/0.325));//*pow(sin(2.*3.14*radius/16.+0.),2.);
    break;
    
  case CaseBhattacharya: 
    // parametrization based on SNR catalogue: Case and  Bhattacharya, A&A Supplement, v.120, p.437-440 (1996)
    if (radius >= 15.0) return 0.0;
    else return pow ( radius/robs , 1.69 ) * exp ( -3.33*(radius-robs)/robs - fabs(z)/0.2  );
    break;
    
  case GiguereKaspi :  
    // Ref: Faucher-Giguere & Kaspi ApJ 643, 332 (2006)
    if (radius >= 15.0) return 0.0;
    else
      return pow(((radius+0.55)/(robs+0.55)),1.64)*exp(-4.0*((radius-robs)/(robs+0.55)))*exp(-fabs(z)/0.1);
    break;
    
  case PointSource : 
    // Point source at GC -- MW130827: You can specify the position now.
    cout <<"$$$ SourceDistribution-PointSource " << " " << SNR_model << endl;
    cout <<" PointSource " << pointsrc_x << " " << pointsrc_y << " " << pointsrc_z << endl;
    if ( fabs(x-pointsrc_x)<.05 && fabs(y-pointsrc_y)<.05 && fabs(z-pointsrc_z)<.05 ) return 1.;
    else return 0.;
    break;
    
  case Ring : 
    if ( radius >= ringmin && radius <= ringmax && fabs(z) <= 0.2 ) return 1.;
    else return 0.;
    break;
    
  case Rings :
    if (radius > 3.7)
      Ferriere_distribution = (50.0*(0.79*exp(-pow(z/0.212,2.))+0.21*exp(-pow(z/0.636,2.)))*exp(-((radius)*(radius)-robs*robs)/(6.8*6.8)) + 7.3*exp(-(radius-robs)/4.5-fabs(z)/0.325));
    else
      Ferriere_distribution = (177.5*(0.79*exp(-pow(z/0.212,2.))+0.21*exp(-pow(z/0.636,2.)))*exp(-pow((radius-3.7)/2.1, 2.)) + 7.3*exp(-(radius-robs)/4.5-fabs(z)/0.325));
    
    if (rings_period == 0) {
      if (fabs(y) < .1 && fabs(z) == 0. && fabs(x) < 10.) {
	//cout << "No rings" << endl;
	if (in->feedback >0) cout << "Ferriere distrib at radius = " << radius << " and z = " << z << " ---> " <<  Ferriere_distribution << endl;
      }
      //cout << endl;
      return Ferriere_distribution;
    }
    else return (Ferriere_distribution * pow(sin( 2.*3.14*radius/rings_period + rings_phase ),2.));
    break;
            
  case Galprop_:   
    // To be done : put here the model by GALPROP
    if (radius >= 15.0){
      return 0;
    }
    else {
      return pow ( radius/robs , 1.25 ) * exp ( -3.56*(radius-robs)/robs - fabs(z)/0.2  );
    }
    break;
    
  case Const :
    return 1.;
    break;
    
  case BlasiSmooth :
    // Add ref!
    return (1/(robs*robs)) * pow((radius/robs),2.) * exp(-spiral_beta*(radius-robs)/robs) * exp(-fabs(z)/spiral_z_h);
    break;
    
    //MW130828: Disable regular sources if I only want to look at one Extra Component:
  case OnlyExtra:
    return 0;
    break;
    
  default :
    return -1;
  }
}

TDMSource::TDMSource(TGrid* Coord_, Input* in_) : TSource() {
   
  in = in_;
  isDmsource = true;
  Coord = Coord_;
   
  vector<double> x = Coord->GetX();
  vector<double> y;
  if (Coord->GetType() == "3D") y = Coord->GetY();
  vector<double> z = Coord->GetZ();
   
  dimx = x.size();
  dimy = (Coord->GetType() == "3D") ? Coord->GetDimY() : 1;
  dimz = z.size();
   
  switch(in->dmprof) {
         
  case ISO :
    Alpha = 2.0;
    Beta  = 2.0;
    Gamma = 0.0;
    Rs    = 3.5;
    rc    = 0.0;
    rhoc  = 0.0;
    break;
         
  case NFW :
    Alpha = 1.0;
    Beta  = 3.0;
    Gamma = 1.0;
    Rs    = 20.;
    rc    = 0.1; // Check!
    rhoc  = Rs;  // Check!
    break;
         
  case Kra :
    Alpha = 2.0;
    Beta  = 3.0;
    Gamma = 0.4;
    Rs    = 10.;
    rc    = 0.1; // Check!
    rhoc  = Rs;  // Check!
    break;
         
  case Moore :
    Alpha = 1.5;
    Beta  = 3.0;
    Gamma = 1.5;
    Rs    = 28.;
    rc    = 0.1; // Check!
    rhoc  = Rs;  // Check!
    break;
         
  case Einasto :
    Alpha = -1.0;
    Beta  = 0.0;
    Gamma = -1.0;
    Rs    = 20.0;
    rc    = 0.0;
    rhoc  = Rs;
    break;
         
  default :
    break;
  }
   
  double rhos1 = in->rhos/in->mx;
  if (in->DMr == Decay) rhos1 /= in->taudec;
   
  for (int i = 0; i < dimx; i++) {
    double DeltaX = Coord->GetDeltaX(i);
    for (int j = 0; j < dimy; j++) {
      double DeltaY = (Coord->GetType() == "3D") ? Coord->GetDeltaY(j) : 0;
      double r = 0;
      if (y.size()) r = sqrt(x[i]*x[i]+y[j]*y[j]);
      else r = x[i];
      double DeltaZ;
      for (vector<double>::iterator zeta = z.begin(); zeta != z.end(); ++zeta)
	{
	  if(zeta==z.begin())      DeltaZ = *(zeta+1) - *(zeta);
	  else if(zeta==z.end()-1) DeltaZ = *(zeta) - *(zeta-1);
	  else                     DeltaZ = 0.5 * ( *(zeta+1) - *(zeta-1) );
            
	  double value = in->sigmav/(1.0+(in->DMr==Annihilation))*pow(rhos1*DM_profile_av(r,*zeta, DeltaZ, DeltaX),1+(in->DMr==Annihilation))*pow(kpc,3.)*Myr;
	  //ANNIHILATION CASE:
	  // [sigma v]     -->                    cm^3/s
	  // [(rhos1)^2]   --> (GeV/cm^3/GeV)^2 = 1/cm^6 
	  // Myr: conversion factor -->           s/Myr
	  // kpc^3: conversion factor -->           (cm/kpc)^3
	  // [value]       --> 1 / (kpc^3 Myr)  

	  // Source term: value [kpc^-3 Myr^-1] * spectrum [GeV^-1]	

	  //MW130725: disable DM ~ Spiral Arms for now
	  /*
	    if(Coord->GetType()=="3D")
	    {
	    value*= max( min( pow(Utility::Spiral_Arm_Density(x[i],y[j],(*zeta),in->SA_type), in->SA_DMsource), in->SA_cut_DMsource), 1./in->SA_cut_DMsource );
	    value*= pow( in->LB_DMsource, Coord->IsInLocalBubble(x[i],y[j],(*zeta)) );
	    }
	  */
	  source.push_back(value);
	}
    }
  }
   
  return ;
}

double TDMSource::DM_profile(double radius, double zeta) {
  double robs = in->robs;
  double rad = sqrt(radius*radius+zeta*zeta);
  double rat = rad/Rs;
  if (Gamma != 0 && rad < rc) rat = rc/Rs;
   
  double prof = 0.0;
  prof = 1.0/pow(rat, Gamma)/pow(1.0+ pow(rat,Alpha),(Beta-Gamma)/Alpha);
  double profsolar = 0.0;
  profsolar = 1.0/pow(robs/Rs, Gamma)/pow(1.0+ pow(robs/Rs,Alpha),(Beta-Gamma)/Alpha);
   
  if (Alpha == -1.0) { // Einasto profile!
    double alph = 0.17;
    prof = exp(-2.0/alph*(pow(rat,alph)-1));
    profsolar = exp(-2.0/alph*(pow(robs/Rs,alph)-1));
  }
   
  return prof/profsolar;
}

double TDMSource::DM_profile_av(double radius, double zeta, double Deltaz, double Deltar) {
  double rad = sqrt(radius*radius+zeta*zeta);
  double rat = rad/Rs;
   
  if (Gamma != 0 && rad < rc) rat = rc/Rs;
   
  double DM_profile_av_=0.0;
  int nuse=0;
   
  for (double zz=zeta-Deltaz/2.; zz<=zeta+Deltaz/2.; zz+=dzzGal)
    for (double rr=radius-Deltar/2.; rr<=radius+Deltar/2.; rr+=Deltar/10.)
      {
	if (rr<0.) continue;
	DM_profile_av_+=DM_profile(rr,zz);
	nuse++;
      }
  return DM_profile_av_/nuse;
}


TAstrophysicalSource::TAstrophysicalSource(TGrid* Coord_, Input* in_, TGeometry* geom, SNRType SNR_model_) : TSource() {
  
  in = in_;
  isDmsource = false;
  SNR_model = SNR_model_;
  Coord = Coord_;
  ringmin = in->ringmin;
  ringmax = in->ringmax;
  rings_period = in->rings_period;
  rings_phase  = in->rings_phase;
  if(SNR_model == Rings && in->feedback >0)
    cout << "***  Rings period " << rings_period << " and phase " << rings_phase << " ***" << endl;
  
  //MW130827: implementing placeable point sources

  cout << "%% AstrophysicalSource " << SNR_model << endl;

  if(SNR_model == PointSource) {

    cout <<" PointSource " << pointsrc_x << " " << pointsrc_y << " " << pointsrc_z << endl;
    if (Coord->GetType() == "2D") {

      double pointsrc_r = sqrt(pointsrc_x*pointsrc_x + pointsrc_y*pointsrc_y);
      //find low bins
      int ix,iz;
      ix = iz = 0;
      while(Coord->GetR()[ix + 1] <= pointsrc_r){ix++;}
      while(Coord->GetZ()[iz + 1] <= pointsrc_z){iz++;}
      
      pointsrc_r = Coord->GetR()[ix];
      pointsrc_z = Coord->GetZ()[iz];
      
      if (in->feedback >0)
	cout << "### 2D-grid  Point Source at (" << pointsrc_r << " , " << pointsrc_z << ") ###" << " " << ix << " " << iz << endl;

    } else {
      //find low bins
      int ix,iy,iz;
      ix = iy = iz = 0;
      while(Coord->GetX()[ix + 1] <= in->pointsrc_x){ix++;}
      while(Coord->GetY()[iy + 1] <= in->pointsrc_y){iy++;}
      while(Coord->GetZ()[iz + 1] <= in->pointsrc_z){iz++;}
      
      pointsrc_x = Coord->GetX()[ix];
      pointsrc_y = Coord->GetY()[iy];
      pointsrc_z = Coord->GetZ()[iz];
      
      if (in->feedback >0)
	cout << "***  Point Source at (" << pointsrc_x << " , " << pointsrc_y << " , " << pointsrc_z << ") ***" << endl;
   }

 }
   
  if (Coord->GetType() == "2D") {
      
    vector<double> r = Coord->GetR();
    vector<double> z = Coord->GetZ();
      
    dimx = r.size();
    dimz = z.size();
      
    for (unsigned int i = 0; i < dimx; ++i) {
      for (unsigned int j = 0; j < dimz; ++j) {
	double radial_distribution;
	double spiral_distribution;
	double spiral_arms;
	double spiral_zeta_h = 1.;
	double spiral_beta = 3.53;
	double sigma;
	double position_of_ring_0;
	double position_of_ring_1;
	double position_of_ring_2;
	double position_of_ring_3;
            
	source.push_back(SourceDistribution(r[i],0,z[j]));
      }
    }

    //Source integral over x,z
    double integral = accumulate(source.begin(),source.end(),0.);
    integral *= Coord->GetDeltaX(0);
    integral *= Coord->GetDeltaZ(0);
    integral *= 6.2831853;
    
    if (in->feedback >0) 
      cout << "Source has been integrated over space. I = " << integral << endl;
  }
  else {
    vector<double> x = Coord->GetX();
    vector<double> y = Coord->GetY();
    vector<double> z = Coord->GetZ();
    
    dimx = x.size();
    dimy = y.size();
    dimz = z.size();
    
    for (unsigned int i = 0; i < dimx; ++i) {
      for (unsigned int k = 0; k < dimy; ++k) {
	for (unsigned int j = 0; j < dimz; ++j) {
	  
	  //MW 130326
	  double value = SourceDistribution(x[i],y[k],z[j]);
	  value*= max( min( pow(geom->GetPattern(i,k,j), in->SA_source), in->SA_cut_source), 1./in->SA_cut_source );
	  if(Coord->IsInLocalBubble(x[i],y[k],z[j])) value*= in->LB_source;
	  
	  source.push_back(value);
	}
      }
    }
    
    //Source integral over x,y,z
    double integral = accumulate(source.begin(),source.end(),0.);
    integral *= Coord->GetDeltaX(0);
    integral *= Coord->GetDeltaY(0);
    integral *= Coord->GetDeltaZ(0);
    
    if (in->feedback >0) 
      cout << "Source has been integrated over space. I = " << integral << endl;
  }
  
  return ;
}
