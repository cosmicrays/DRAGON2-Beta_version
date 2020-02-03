/**
 * @file gas.cc
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief File where gas classes are implemented.
 */
#include "gas.h"
#include "grid.h"
#include "galaxy.h" // mw130620
#include "geometry.h"
#include "constants.h"
#include "input.h"
#include <fstream> //sk

using namespace std;

TH2Gas::TH2Gas(TGrid* Coord, Input* in, TGeometry* geom) : TGas() {
   
  vector<double> x;
  vector<double> y;
  vector<double> z;
   
  if (Coord->GetType() == "2D") {
    x = Coord->GetR();
    z = Coord->GetZ();
      
    dimr = x.size();
    dimx = dimr;
    dimy = 1;
    dimz = z.size();
      
    for (int k = 0; k < dimx; k++) {
      double r = x[k];
         
      for (int l = 0; l < dimz; l++) {
	double DeltaZ = Coord->GetDeltaZ(l);
            
	double height_H2 = 0;
	switch(in->gas_model) {
	case BronfFerr :
	  height_H2 = (r > 2.5) ? 0.059 : 0.015;
	  density.push_back(2.0* (150.*exp(-pow((r-0.05)/.2, 2.))+0.49*exp(-pow((r-5.3)/2.2, 2.))+.24*exp(-pow((r-8.1)/1.5, 2.)))*exp(-M_LN2*pow(z[l]/height_H2,2.)));
	  break;
                  
	case NS :
	  density.push_back(2.0*.94*(11.2*exp(-(r)*(r)/.874)+.83*exp(-pow((r-4)/3.2, 2.)))*exp(-M_LN2*pow(z[l]/(1.06/1000.*(10.8*exp(0.28*(r))+42.78)), 2.)));
	  break;
                  
	case Galprop :
	  density.push_back(2.0*TGas::nH2_av(r, z[l], DeltaZ, dzzGal, in));
	  break;
                  
	default:
	  density.push_back(-1);
	}
      }
    }
      
  }
  else {
    x = Coord->GetX();
    y = Coord->GetY();
    z = Coord->GetZ();
      
    dimx = x.size();
    dimy = y.size();
    dimz = z.size();
      
      
    for (int k = 0; k < dimx; k++) {
      for (int i = 0; i < dimy; i++) {
            
	double r = sqrt(x[k]*x[k]+y[i]*y[i]);
            
	for (int l = 0; l < dimz; l++) {
	  double DeltaZ = Coord->GetDeltaZ(l);
               
	  double height_H2 = 0;
	  switch(in->gas_model) {
	  case BronfFerr :
	    height_H2 = (r > 2.5) ? 0.059 : 0.015;
	    density.push_back(2.0* (150.*exp(-pow((r-0.05)/.2, 2.))+0.49*exp(-pow((r-5.3)/2.2, 2.))+.24*exp(-pow((r-8.1)/1.5, 2.)))*exp(-M_LN2*pow(z[l]/height_H2,2.)));
	    break;
                     
	  case NS :
	    density.push_back(2.0*.94*(11.2*exp(-(r)*(r)/.874)+.83*exp(-pow((r-4)/3.2, 2.)))*exp(-M_LN2*pow(z[l]/(1.06/1000.*(10.8*exp(0.28*(r))+42.78)), 2.)));
	    break;
                     
	  case Galprop :
	    density.push_back(2.0*TGas::nH2_av(r, z[l], DeltaZ, dzzGal, in));
                     
	    break;
                     
	  default:
	    density.push_back(-1);
	  }
	  //mw, 130326
	  double val = density.back();
	  val *= pow( in->LB_gasH2 * in->LB_gastotal, Coord->IsInLocalBubble(x[k],y[i],z[l]) );
	  density.pop_back();
	  density.push_back(val);
	}
      }
    }
    if (in->feedback >0) cout << "[MW] Total content of H2 gas before Spiral arms: " << GetTotalContent() << endl;
    if(in->SA_type!="None") geom->ApplySpiral(density,in->SA_gasH2, in->SA_cut_gasH2);
      
  }
   
  if (in->feedback >0) cout << "[MW] Total content of H2 gas: " << GetTotalContent() << endl;
   
}

THIGas::THIGas(TGrid* Coord, Input* in, TGeometry* geom) : TGas() {
   
  vector<double> x;
  vector<double> y;
  vector<double> z;
   
  if (Coord->GetType() == "2D") {
    x = Coord->GetR();
    z = Coord->GetZ();
      
    dimr = x.size();
    dimx = dimr;
    dimy = 1;
    dimz = z.size();
    for (int k = 0; k < dimx; k++) {
         
      double r = x[k];
         
      for (int l = 0; l < dimz; l++) {
	double DeltaZ = Coord->GetDeltaZ(l);
            
	double altezza_H1 = 0;
	double densita_H1 = 0;
            
	switch(in->gas_model) {
	case BronfFerr :
	  if (r < 2.5) altezza_H1 = 0.045;
	  else if (r < 8.5) altezza_H1 = 0.115;
	  else altezza_H1 = 0.115*exp((r-8.5)/6.7);
	  if (r < 3.) densita_H1 = 0.57*exp((r-3.)/0.4)+8.*exp(-pow(r/0.2, 2.));
	  else if (r < 13.) densita_H1 = 0.57;
	  else densita_H1 = 0.57*exp(-(r-13.)/4.);
	  density.push_back( densita_H1*exp(-M_LN2*pow(z[l]/altezza_H1, 2.)) );
	  break;
                  
	case NS :
	  densita_H1 =  .94*(.6*exp(-r/2.4) + .24*exp(-pow((r-9.5)/4.8, 2.)));
	  altezza_H1 = 1.06*(116.3 + 19.3*(r) +  4.1*(r)*(r) -0.05*(r)*(r)*(r))/1000.;
	  density.push_back(densita_H1*exp(-M_LN2*pow(z[l]/altezza_H1, 2.)));
	  break;
                  
	case Galprop :
	  density.push_back(TGas::nHI_av(r,z[l], DeltaZ, dzzGal));
	  break;
                  
	default:
	  density.push_back(-1);
	}
      }
    }
      
      
  }
  else {
    x = Coord->GetX();
    y = Coord->GetY();
    z = Coord->GetZ();
      
    dimx = x.size();
    dimy = y.size();
    dimz = z.size();
      
      
      
    for (int k = 0; k < dimx; k++) {
      for (int i = 0; i < dimy; i++) {
            
	double r = sqrt(x[k]*x[k]+y[i]*y[i]);
            
	for (int l = 0; l < dimz; l++) {
	  double DeltaZ = Coord->GetDeltaZ(l);
               
	  double altezza_H1 = 0;
	  double densita_H1 = 0;
               
	  switch(in->gas_model) {
	  case BronfFerr :
	    if (r < 2.5) altezza_H1 = 0.045;
	    else if (r < 8.5) altezza_H1 = 0.115;
	    else altezza_H1 = 0.115*exp((r-8.5)/6.7);
	    if (r < 3.) densita_H1 = 0.57*exp((r-3.)/0.4)+8.*exp(-pow(r/0.2, 2.));
	    else if (r < 13.) densita_H1 = 0.57;
	    else densita_H1 = 0.57*exp(-(r-13.)/4.);
	    density.push_back( densita_H1*exp(-M_LN2*pow(z[l]/altezza_H1, 2.)) );
	    break;
                     
	  case NS :
	    densita_H1 =  .94*(.6*exp(-r/2.4) + .24*exp(-pow((r-9.5)/4.8, 2.)));
	    altezza_H1 = 1.06*(116.3 + 19.3*(r) +  4.1*(r)*(r) -0.05*(r)*(r)*(r))/1000.;
	    density.push_back(densita_H1*exp(-M_LN2*pow(z[l]/altezza_H1, 2.)));
	    break;
                     
	  case Galprop :
	    density.push_back(TGas::nHI_av(r,z[l], DeltaZ, dzzGal));
                     
	    break;
                     
	  default:
	    density.push_back(-1);
	  }
	  //mw, 130326
	  double val = density.back();
	  val *= pow( in->LB_gasHI * in->LB_gastotal, Coord->IsInLocalBubble(x[k],y[i],z[l]) );
	  density.pop_back();
	  density.push_back(val);
	}
      }
    }
      
    if (in->feedback >0) cout << "[MW] Total content of HI gas before Spiral arms: " << GetTotalContent() << endl;
    if(in->SA_type!="None") geom->ApplySpiral(density,in->SA_gasHI,in->SA_cut_gasHI);
      
  }
  if (in->feedback >0) cout << "[MW] Total content of HI gas: " << GetTotalContent() << endl;
   
}

THIIGas::THIIGas(TGrid* Coord, Input* in, TGeometry* geom) : TGas() {
   
  vector<double> x;
  vector<double> y;
  vector<double> z;
   
  if (Coord->GetType() == "2D") {
    x = Coord->GetR();
    z = Coord->GetZ();
      
    dimr = x.size();
    dimx = dimr;
    dimy = 1;
    dimz = z.size();
      
    for (int k = 0; k < dimx; k++) {
      double r = x[k];
         
      for (int l = 0; l < dimz; l++) {
	double DeltaZ = Coord->GetDeltaZ(l);
            
	switch(in->gas_model) {
	case BronfFerr :
	case NS :
	  density.push_back(0);
	  break;
                  
	case Galprop :
	  density.push_back(TGas::nHII_av(r, z[l], DeltaZ, dzzGal));
	  //cout << density.back() << endl;
	  break;
                  
	default:
	  density.push_back(-1);
	}
      }
    }
      
      
  }
  else {
    x = Coord->GetX();
    y = Coord->GetY();
    z = Coord->GetZ();
      
    dimx = x.size();
    dimy = y.size();
    dimz = z.size();
      
      
    for (int k = 0; k < dimx; k++) {
      for (int i = 0; i < dimy; i++) {
            
	double r = sqrt(x[k]*x[k]+y[i]*y[i]);
            
	for (int l = 0; l < dimz; l++) {
	  double DeltaZ = Coord->GetDeltaZ(l);
               
	  switch(in->gas_model) {
	  case BronfFerr :
	  case NS :
	    density.push_back(0);
	    break;
                     
	  case Galprop :
	    density.push_back(TGas::nHII_av(r, z[l], DeltaZ, dzzGal));
                     
	    break;
                     
	  default:
	    density.push_back(-1);
	  }
	  //mw, 130326
	  double val = density.back();
	  val *= pow( in->LB_gasHII * in->LB_gastotal, Coord->IsInLocalBubble(x[k],y[i],z[l]) );
	  density.pop_back();
	  density.push_back(val);
	}
      }
    }
      
    if (in->feedback >0) cout << "[MW] Total content of HII gas before Spiral arms: " << GetTotalContent() << endl;
    if(in->SA_type!="None") geom->ApplySpiral(density,in->SA_gasHII, in->SA_cut_gasHII);
      
  }
   
  if (in->feedback >0) cout << "[MW] Total content of HII gas: " << GetTotalContent() << endl;
   
}

TGas::TGas(TGrid* Coord, Input* in) {
  //Construct TotalGas from other components!
   
  vector<double> x;
  vector<double> y;
  vector<double> z;
  double value = 0.0;
   
  if (Coord->GetType() == "2D")
    {
      x = Coord->GetR();
      z = Coord->GetZ();
      
      dimr = x.size();
      dimx = dimr;
      dimy = 1;
      dimz = z.size();
      
      
      for (int k = 0; k < dimx; k++)
	for (int l = 0; l < dimz; l++)
	  density.push_back(0);
    }
  else
    {
      x = Coord->GetX();
      y = Coord->GetY();
      z = Coord->GetZ();
      
      dimx = x.size();
      dimy = y.size();
      dimz = z.size();
      
      for (int k = 0; k < dimx; k++)
	for (int i = 0; i < dimy; i++)
	  for (int l = 0; l < dimz; l++)
	    density.push_back(0);
    }
   
  return;
}


// Galprop:
//
// GAS DENSITY MODELS
double TGas::nH2_Gal(double r, double z, Input* in) {
   
  int i;
  double nH2_ = 0.0, fR,fZ0,fZh;                                              // [B88]/Table 3
  double R[18] ={ 0.00, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75,       // kpc, col.1
		  6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75,10.25};
  double Y[18] ={ 0.00,  1.5,  3.3,  5.8,  5.5,  8.4,  9.0,  9.6,  8.6,       // CO, K km s^-1
		  9.1,  7.9,  9.2,  7.7,  5.0,  3.6,  4.8,  1.7,  0.0}; // (col.4)
  double Z0[18]={0.039,0.039,0.036,0.000,-.008,0.001,-.010,-.001,-.004,       // kpc, col.7
		 -.019,-.022,-.014,-.009,-.004,0.013,-.004,-.020,-.020};
  double Zh[18]={0.077,0.077,0.080,0.061,0.065,0.071,0.072,0.082,0.083,       // kpc, col.10
		 0.073,0.063,0.058,0.072,0.080,0.066,0.023,0.147,0.147};
   
  double H2toCO = X_CO(r, in);
  H2toCO *= 1.e20;
   
  if(r > R[17]) return nH2_;
  for (i=0; i<17; i++)  if(R[i] <= r && r <= R[i+1])  break;
   
  fR =  Y[i] + ( Y[i+1] - Y[i])/(R[i+1] - R[i])*(r - R[i]);
  fZ0= Z0[i] + (Z0[i+1] -Z0[i])/(R[i+1] - R[i])*(r - R[i]);
  fZh= Zh[i] + (Zh[i+1] -Zh[i])/(R[i+1] - R[i])*(r - R[i]);
  nH2_ =  fR * exp( -log(2.)*pow( (z-fZ0)/fZh, 2 ) )  *H2toCO/kpc;
  return nH2_< 0. ? 0.: nH2_;
}

double TGas::nHI_Gal(double Rkpc, double Zkpc)
{
  int i;                                                             // Table 1 [GB76]
  double R[30] ={ 0.0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0,  // kpc, col.1
		  6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,10.0,10.5,11.0,
		  11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,15.5,16.0},
    Y[30] ={ .10, .13, .14, .16, .19, .25, .30, .33, .32, .31,  // nHI, cm^-3
	     .30, .37, .38, .36, .32, .29, .38, .40, .25, .23,  // (col.3)
	     .32, .36, .32, .25, .16, .10, .09, .08, .06, .00};
    double fR, fZ,fZ1=0.,fZ2=0., R1,R2=R[29], Y1,Y2=Y[29];
    double nGB =0.33, nDL =0.57;       // cm^-3, disk density @ 4-8 kpc; [GB76], [DL90]
    double A1=0.395,     z1=0.212/2.,  // cm^-3, kpc; Z-distribution parameters from [DL90]
      A2=0.107,     z2=0.530/2.,
      B =0.064,     zh=0.403;
   
    for (i=0; i<29; i++)  if(R[i] <= Rkpc && Rkpc <= R[i+1])  break;
   
    R1 = (R[i]+R[i+1])/2;   Y1 = Y[i];
    if(Rkpc < R1)
      {
	if(i> 0)    { R2 = (R[i-1]+R[i])/2;   Y2 = Y[i-1]; }
	else        { R2 = R[0];              Y2 = Y[0];   }
      }
    else  if(i<28) { R2 = (R[i+1]+R[i+2])/2; Y2 = Y[i+1]; }
   
    fR = Y1 +(Y2 -Y1)/(R2 -R1)*(Rkpc -R1);                             // interpolation in R
   
    R2 = (R[28] +R[29]) /2;
    if(Rkpc > R2) fR = Y[28]*exp(-(Rkpc-R2)/3);                        // extrapolation in R
   
    // calculation of Z-dependence
    if(Rkpc <10.)                                                      // [DL90]
      fZ1 =A1*exp(-log(2.)*pow(Zkpc/z1,2))+A2*exp(-log(2.)*pow(Zkpc/z2,2))+B*exp(-fabs(Zkpc)/zh);
    if(Rkpc > 8.) fZ2=nDL*exp(-pow(Zkpc /(0.0523*exp(0.11*Rkpc)), 2)); // [C86] IMOS20010220
   
    if(Rkpc <= 8.) fZ = fZ1;
    else
      {   if(Rkpc >=10.) fZ = fZ2;
	else fZ = fZ1 +(fZ2 -fZ1)/2.*(Rkpc -8.);                       // interp. [DL90] & [C86]
      }
    return fZ *fR/nGB;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double TGas::nHII_Gal(double Rkpc,double Zkpc)
{
  double fne1=0.025, H1=1.00, A1=20.0;
  double fne2=0.200, H2=0.15, A2= 2.0;
  double R2=4.0;
  double ne1 = fne1 * exp(-fabs(Zkpc)/H1) * exp (-pow( Rkpc    /A1, 2));
  double ne2 = fne2 * exp(-fabs(Zkpc)/H2) * exp (-pow((Rkpc-R2)/A2, 2));
  return ne1+ne2;
}

double TGas::nH2_av(double R, double z,double dz,double dzz, Input* in)
{
  //   double R=sqrt(x*x + y*y);
  double nH2_av_=0.0;
  int nuse=0;
   
  for (double zz=z-dz/2.; zz<=z+dz/2.; zz+=dzz)
    {
      nH2_av_+=nH2_Gal(R,zz, in);
      nuse++;
    }
  return nH2_av_/nuse;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double TGas::nHI_av(double R, double z,double dz,double dzz)
{
   
  //   double R=sqrt(x*x + y*y);
  double nHI_av_=0.0;
  int nuse=0;
   
  for (double zz=z-dz/2.; zz<=z+dz/2.; zz+=dzz)
    {
      nHI_av_+=nHI_Gal(R,zz);
      nuse++;
    }
  return nHI_av_/nuse;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double TGas::nHII_av(double R, double z,double dz,double dzz)
{
  //   double R=sqrt(x*x + y*y);
  double nHII_av_=0.0;
  int nuse=0;
   
  for(double zz=z-dz/2.; zz<=z+dz/2.; zz+=dzz)
    {
      nHII_av_+=nHII_Gal(R,zz);
      nuse++;
    }
  return nHII_av_/nuse;
}


//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double TGas::X_CO(double r, Input* in) {
   
   
  switch (in->xco_mode) {
  case(SM96): // Strong&Mattox 1996
    return 1.9;
    break;
  case (galprop_2004): //Strong2004
    if (r < 3.5)
      return 0.4;
    else if (r < 5.5)
      return 0.6;
    else if (r < 7.5)
      return 0.8;
    else if (r < 9.5)
      return 1.5;
    else return 10.;
    break;
  case (galprop_2010): //compatible with Fermi-LAT gamma profile?
    if (r < 2.0)
      return 0.4;
    else if (r < 5.5)
      return 1.5;
    else if (r < 7.5)
      return 1.5;
    else if (r < 10.0)
      return 0.75;
    else if (r < 12.0)
      return 5.;
    else if (r < 17.0)
      return 12;
    else if (r < 19.0)
      return 60;
    else if (r < 30.0)
      return 200;
    else
      return 200;
    break;
  case (constant):
    return in->xco_constant;
    break;
  case (dragon):
    if (r < 2.0)
      return in->xco_inner;
    else
      return in->xco_outer;
    break;
  default:
    return 1.;
  }
}
