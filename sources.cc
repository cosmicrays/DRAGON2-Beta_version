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
  double theta = atan2(y,x) + M_PI;
  
  double Ferriere_distribution;
  double spiral_beta = 3.53;
  double spiral_z_h = 1.00;
  double robs = in->robs;

  double Ferriere_classic;
  double rate_Ferr_typeI;
  double A_Ferr;
  double expo_sum_Ferr;
  double R0_Ferr;

  double R1_Yus;
  double a_Yus;
  double b_Yus;
  double z0_Yus;




  switch (SNR_model) {


  case Ferriere:
    // parametrization based on PSR catalogue + disk stars: K. Ferriere, Rev.Mod.Phys. 73, 1031-1066 (2001)

  	rate_Ferr_typeI = 7.3 * exp( - (radius-robs)/4.5 - fabs(z)/0.325 );
  	expo_sum_Ferr = 0.79 * exp(-pow(z/0.212,2.)) + 0.21 * exp(-pow(z/0.636,2.));

  	if (radius > 3.7)	{
    	A_Ferr = 50.0;
    	Ferriere_classic = A_Ferr * expo_sum_Ferr * exp( - pow( (radius-3.7)/2.1 , 2.) ) + rate_Ferr_typeI;
    					}
    else 	{
    	A_Ferr = 177.5;
    	Ferriere_classic = A_Ferr * expo_sum_Ferr * exp( - pow( (radius-robs)/6.8 , 2.) ) + rate_Ferr_typeI;
    		}
    return Ferriere_classic;
    break;
    
  

  // The following 4 cases are different parametrizations of the same general source-distribution function
  // the parameters are (R1_Yus, a_Yus, b_Yus, z0_Yus)
  // from Yusifov & Kucuk: A&A 459 p.545-553 (2004). DRAGON technical paper eq.(C.9)

  case Lorimer: 
    // parametrization based on PSR catalogue: Lorimer et al., Mon. Not. R. Astron. Soc. 372, 777 (2006)

  	R1_Yus = 0.;
  	a_Yus = 1.9;
  	b_Yus = 5.00;
  	z0_Yus = 0.2;
    if (radius >= 15.0) return 0.0;
    else return pow( (radius+R1_Yus) / (robs+R1_Yus) , a_Yus ) * exp( -b_Yus * ( (radius-robs)/(robs+R1_Yus) ) - fabs(z)/z0_Yus );
    break;
    
  case CaseBhattacharya: 
    // parametrization based on SNR catalogue: Case and  Bhattacharya, A&A Supplement, v.120, p.437-440 (1996)

  	R1_Yus = 0.;
  	a_Yus = 1.69;
  	b_Yus = 3.33;
  	z0_Yus = 0.2;
    if (radius >= 15.0) return 0.0;
    else return pow( (radius+R1_Yus) / (robs+R1_Yus) , a_Yus ) * exp( -b_Yus * ( (radius-robs)/(robs+R1_Yus) ) - fabs(z)/z0_Yus );
    break;
    
  case GiguereKaspi :  
    // parametrization based on PSR catalogue: Faucher-Giguere & Kaspi ApJ 643, 332 (2006)

  	R1_Yus = 0.55;
  	a_Yus = 1.64;
  	b_Yus = 4.00;
  	z0_Yus = 0.1;
    if (radius >= 15.0) return 0.0;
    else return pow( (radius+R1_Yus) / (robs+R1_Yus) , a_Yus ) * exp( -b_Yus * ( (radius-robs)/(robs+R1_Yus) ) - fabs(z)/z0_Yus );
    break;

  case StrongMoskalenko1998:  
  	// ad hoc parametrization to correct the "gradient problem": Strong & Moskalenko: ApJ 509:212-228 (1998), eq.(6) - Table 2

  	R1_Yus = 0.;
  	a_Yus = 0.5;
  	b_Yus = 1.00;
  	z0_Yus = 0.2;
    if (radius >= 20.0) return 0.0;
    else return pow( (radius+R1_Yus) / (robs+R1_Yus) , a_Yus ) * exp( -b_Yus * ( (radius-robs)/(robs+R1_Yus) ) - fabs(z)/z0_Yus );
    break;
            
            
  case BlasiAmato :
    // Blasi & Amato: JCAP01(2012)011, eq.(4.1) [radial profile] + eq.(4.3) [latitude profile]
    return 1./(robs*robs) * pow( radius/robs , 2. ) * exp( -spiral_beta*(radius-robs)/robs ) * exp( -fabs(z)/spiral_z_h );
    break;
    
  // customizable TEST cases
  case Const :
    return 1.;
    break;
    
  case OnlyExtra:
  // Disable regular sources if I only want to look at one Extra Component:
    return 0.0;
    break;
    
  default:
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

	/*switch(in->dmprof) {

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
	}*/

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
	double radius_spherical = sqrt(radius*radius+zeta*zeta);
	
	double profile;
	double profile_evaluated_at_sun;
	double profile_normalized;
	
	double Rs;
	
	double alpha = 0.17;
	
	switch(in->dmprof) {
	
		case ISO :
			
			Rs = 4.38; // kpc. From 1012.4515
			profile = 1./( 1. + pow(radius_spherical/Rs, 2.) );		
			profile_evaluated_at_sun = 1./( 1. + pow(robs/Rs, 2.) );				
			break;
			
		case NFW:
			
			Rs = 24.42; // kpc. From 1012.4515		
			profile = (Rs/radius_spherical)/pow((1. + radius_spherical/Rs), 2.);		
			profile_evaluated_at_sun = (Rs/robs)/pow((1. + robs/Rs), 2.);				
			break;
			
		
		case Einasto:
			
			profile = (-2./alpha)*( pow(radius_spherical/Rs, alpha) - 1.);
			profile_evaluated_at_sun = (-2./alpha)*( pow(robs/Rs, alpha) - 1.);
			break;
			
		default:
			
			profile = 1.;
			profile_evaluated_at_sun = 1.;
			break;
			
	}
	
	profile_normalized = profile/profile_evaluated_at_sun;
	return profile_normalized;
	
	/*double robs = in->robs;
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

	return prof/profsolar;*/
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


	cout << "%% AstrophysicalSource " << SNR_model << endl;


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
