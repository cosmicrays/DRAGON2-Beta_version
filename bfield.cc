//
//  bfield.cc
//  DRAGON
//
//  Created by Luca Maccione on 04/06/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//
//  MW: merged farrar.cc and pshirkovfield.cc into this one.


#include <iostream>

#include "bfield.h"
#include "grid.h"
#include "constants.h"
#include "geometry.h"

#include <iostream>
#include <vector>
#include <math.h>

#include <fstream>

using namespace std;

TBField::TBField(TGrid* Coord) {

  //cout << "This is the generic magnetic field constructor " << endl;
    
    //eta = 0.0;
    
    if (Coord->GetType() == "3D") {
        x = Coord->GetX();
        y = Coord->GetY();
        dimx = x.size();
        dimy = y.size();
    }
    else {
        x = Coord->GetR();
        dimx = x.size();
        dimy = 1;
    }
    
    z = Coord->GetZ();
    dimz = z.size();

    
    MagField = std::vector<double>(dimx*dimy*dimz,0.0);
    //cout << "Created the MagField structure " << endl;
    
    if (Coord->GetType() == "3D") {
	cout << "initializing Breg and Breg_versors" << endl;
        Breg[0] = std::vector<double>(dimx*dimy*dimz,0.0);
        Breg[1] = std::vector<double>(dimx*dimy*dimz,0.0);
        Breg[2] = std::vector<double>(dimx*dimy*dimz,0.0);
        Breg_versors[0] = std::vector<double>(dimx*dimy*dimz,0.0);
        Breg_versors[1] = std::vector<double>(dimx*dimy*dimz,0.0);
        Breg_versors[2] = std::vector<double>(dimx*dimy*dimz,0.0);
    }
    
    //double Bdisk,
    //Bhalo;
    
    //const double RC = 5.0,
    //B0 = 2.0,
    //Z0 = 1.0,
    //R0 = 10.0;
    
    //const double B0H = 4., R0H = 8., Z0H = 1.3; // MEAN
    //const double B0H = 3., R0H = 6., Z0H = 1.;  // MIN
    //const double B0H = 12., R0H = 15., Z0H = 2.5; // MAX
    
    
    /*
     for (int i = 0; i < dimx; i++) {
     for (int l = 0; l < dimy; l++) {
     double radius = 0.0;
     if (y.size()) radius = sqrt(x[i]*x[i]+y[l]*y[l]);
     else radius = x[i];
     for (vector<double>::iterator zeta = z.begin(); zeta != z.end(); ++zeta) {
     
     double Z1H = (fabs(*zeta) < Z0H) ? 0.2 : 0.4;
     
     Bdisk = B0 * ( (radius < RC) ? exp(-fabs(*zeta)/Z0) : exp(-(radius-robs)/R0-fabs(*zeta)/Z0) );
     Bhalo = B0H / (1.+pow((fabs(*zeta)-Z0H)/Z1H,2)) * radius/R0H * exp(-1.-radius/R0H);
     MagField.push_back(Bdisk+Bhalo);
     //MagField.push_back(Bh*exp(-((radius)-robs)/rB)*exp(-fabs((*zeta))/zr));
     }
     }
     }
     */
}

std::vector<double> TBField::GetBregVersors(double x, double y, double z) {
    
    vector<double> B = GetField(x,y,z);
    //cout << "GetBregVersors " << endl;
    //cout << B[0] << " " << B[1] << " " << B[2] << endl;
    double Bmod = sqrt( pow(B[0],2) + pow(B[1],2) + pow(B[2],2) );
    //cout << "Bmod = " << Bmod << endl;	
    B[0] /= Bmod;
    B[1] /= Bmod;
    B[2] /= Bmod;
    return B;
}

double TBField::GetEnDensity(int i /**< linearized position index */) { return pow(MagField[i],2)/8.0/Pi*erg_to_eV; }
double TBField::GetEnDensity(int ir /**< radial position */, int iz /**< vertical position */) { return pow(MagField[index(ir,iz)],2)/8.0/Pi*erg_to_eV; }
double TBField::GetEnDensity(int ix /**< radial position */, int iy, int iz /**< vertical position */) { return pow(MagField[index(ix,iy,iz)],2)/8.0/Pi*erg_to_eV; }

TUniformField::TUniformField(double B, TGrid* Coord, TGeometry* geom) : TBField(Coord) {
    int counter = 0;
    for (int i = 0; i < dimx; i++) {
        for (int l = 0; l < dimy; l++) {
            for (int k = 0; k < dimz; k++) {
                MagField[counter] = B;

                if(Coord->GetType() == "3D")
                {
                    double spiral_factor_magfield = max( min( pow(geom->GetPattern(i,l,k), Coord->in_SA_MagField), Coord->in_SA_cut_MagField), 1./Coord->in_SA_cut_MagField );
                    MagField[counter] *= spiral_factor_magfield;
    
                    MagField[counter] *= pow( Coord->in_LB_MagField, Coord->IsInLocalBubble(x[i],y[l],z[k]) );
                }
                
                counter++;
            }
        }
    }
}


TPshirkovField::TPshirkovField(double B0_, double RC_, double Z0_, double R0_, double B0H_, double R0H_, double Z0H_, double B0turb_, double rscale_turb_, double zscale_turb_, double robs_, TGrid* Coord, TGeometry* geom) :
TBField(Coord),
B0(B0_),
RC(RC_),
Z0(Z0_),
R0(R0_),
B0H(B0H_),
R0H(R0H_),
Z0H(Z0H_),
B0turb(B0turb_),
rscale_turb(rscale_turb_),
zscale_turb(zscale_turb_),
robs(robs_)
{
  //cout << "Welcome to the Pshirkovfield constructor" << endl;
  //cout << "Sun is at R = " << robs << endl;
  //cout << "B0 = " << B0 << endl;
  //cout << "B0H = " << B0H << endl;	
    
    //eta = 0.0;
    //cout << B0H <<  " " << B0turb << endl;
	
    vector<double> x;
    vector<double> y;
    
    if (Coord->GetType() == "3D") {
        x = Coord->GetX();
        y = Coord->GetY();
    }
    else {
        x = Coord->GetR();
    }
    
    vector<double> z = Coord->GetZ();
    
    //    dimz = z.size();
    
    double Bdisk,
    Bhalo;
    
    //const double RC = 5.0,
    //B0 = 2.0,
    //Z0 = 1.0,
    //R0 = 10.0;
    
    //const double B0H = 4., R0H = 8., Z0H = 1.3; // MEAN
    //const double B0H = 3., R0H = 6., Z0H = 1.;  // MIN
    //const double B0H = 12., R0H = 15., Z0H = 2.5; // MAX
    int counter = 0;
    for (int i = 0; i < dimx; i++) {
        for (int l = 0; l < dimy; l++) {
            for (int k = 0; k < dimz; k++) {
                vector<double> B = (Coord->GetType() == "3D") ? GetField(x[i],y[l],z[k]) : GetField(x[i],0,z[k]);

		if (Coord->GetType() == "3D") {
                
	                Breg[0][counter] = B[0];
	                Breg[1][counter] = B[1];
	                Breg[2][counter] = B[2];
	                double Breg_mod = sqrt( pow(B[0],2) + pow(B[1],2) +pow(B[2],2) );
	                Breg_versors[0][counter] = Breg[0][counter]/Breg_mod;
	                Breg_versors[1][counter] = Breg[1][counter]/Breg_mod;
	                Breg_versors[2][counter] = Breg[2][counter]/Breg_mod;
        	        //MagField.push_back(Bh*exp(-((radius)-robs)/rB)*exp(-fabs((*zeta))/zr));
                	MagField[counter] = sqrt(pow(Breg_mod,2)+pow(B[3],2));
		}
		else {

                	MagField[counter] = sqrt(pow(B[0],2)+pow(B[1],2)+pow(B[2],2)+pow(B[3],2));

		}

                if(Coord->GetType() == "3D")
                { 
                    double spiral_factor_magfield = max( min( pow(geom->GetPattern(i,l,k), Coord->in_SA_MagField), Coord->in_SA_cut_MagField), 1./Coord->in_SA_cut_MagField );
                    MagField[counter] *= spiral_factor_magfield;

                    MagField[counter] *= pow( Coord->in_LB_MagField, Coord->IsInLocalBubble(x[i],y[l],z[k]) );
                }

                counter++;
            }
        }
    }
}

vector<double> TPshirkovField::GetField(double x, double y, double z) {
    
    double Z1H = (fabs(z) < Z0H) ? 0.2 : 0.4;
    double radius = sqrt(x*x+y*y);
    
    double Bdisk = B0 * ( (radius < RC) ? exp(-fabs(z)/Z0) : exp(-(radius-robs)/R0-fabs(z)/Z0) );
    double Bhalo = B0H / (1.+pow((fabs(z)-Z0H)/Z1H,2)) * radius/R0H * exp(1.-radius/R0H);
    double Breg_mod = Bdisk+Bhalo;
    double Brand = B0turb*exp(-(radius-robs)/rscale_turb)*exp(-fabs((z)/zscale_turb));
    double phi = atan2(y,x)+M_PI/2.0;
    
    vector<double> Bmixed;
    Bmixed.push_back(Breg_mod*cos(phi));
    Bmixed.push_back(Breg_mod*sin(phi));
    Bmixed.push_back(0);
    Bmixed.push_back(Brand);

    return Bmixed;
}


TFarrarField::TFarrarField(double betaFarrar, TGrid* Coord, TGeometry* geom) : TBField(Coord),
bring(0.1),
hdisk(0.40),
wdisk(0.27),
Bn(1.4),
Bs(-1.1),
rn(9.22),
rs(16.7),
wh(0.2),
z0(5.3),
BX(4.6),
Theta0X(49.0*M_PI/180.0),
rcX(4.8),
rX(2.9),
p(11.5*M_PI/180.0)
{
    
    bj = std::vector<double>(8,0);
    bj[0] = 0.1;
    bj[1] = 3.0;
    bj[2] = -0.9;
    bj[3] = -0.8;
    bj[4] = -2.0;
    bj[5] = -4.2;
    bj[6] = 0.0;
    bj[7] = 0.0;
    
   
   fj = std::vector<double>(8,0);
   fj[0] = 0.130;
   fj[1] = 0.165;
   fj[2] = 0.094;
   fj[3] = 0.122;
   fj[4] = 0.13;
   fj[5] = 0.118;
   fj[6] = 0.084;
   fj[7] = 0.156;
   
   
   for (int i = 0; i < 7; i++) bj[7] -= (fj[i]*bj[i]/fj[7]);
   
   
   rj = std::vector<double>(8,0);
   rj[0] = 5.1;
   rj[1] = 6.3;
   rj[2] = 7.1;
   rj[3] = 8.3;
   rj[4] = 9.8;
   rj[5] = 11.4;
   rj[6] = 12.7;
   rj[7] = 15.5;
   
    vector<double> x;
    vector<double> y;
   
    if (Coord->GetType() == "3D") {
        x = Coord->GetX();
        y = Coord->GetY();
        //dimx = x.size();
        //dimy = y.size();
    }
    else {
        x = Coord->GetR();
        //dimx = x.size();
        //dimy = 1;
    }
    
    vector<double> z = Coord->GetZ();
    
    int counter = 0;
    if (Coord->GetType() == "3D") {
        x = Coord->GetX();
        y = Coord->GetY();
        
        for (int i = 0; i < dimx; i++) {
            for (int j = 0; j < dimy; j++) {
                for (int k = 0; k < dimz; k++) {
                    
                    vector<double> Breg1 = GetField(x[i],y[j],z[k]);
                    
                    Breg[0][counter] = Breg1[0];
                    Breg[1][counter] = Breg1[1];
                    Breg[2][counter] = Breg1[2];
                    
                    double Breg_mod = pow(Breg1[0],2)+pow(Breg1[1],2)+pow(Breg1[2],2);
                    
                    Breg_versors[0][counter] = Breg1[0]/sqrt(Breg_mod);
                    Breg_versors[1][counter] = Breg1[1]/sqrt(Breg_mod);
                    Breg_versors[2][counter] = Breg1[2]/sqrt(Breg_mod);
                    
                    MagField[counter] = sqrt((1.0+betaFarrar)*Breg_mod);

                    if(Coord->GetType() == "3D")
                    { 
                        double spiral_factor_magfield = max( min( pow(geom->GetPattern(i,j,k), Coord->in_SA_MagField), Coord->in_SA_cut_MagField), 1./Coord->in_SA_cut_MagField );
                        MagField[counter] *= spiral_factor_magfield;

                        MagField[counter] *= pow( Coord->in_LB_MagField, Coord->IsInLocalBubble(x[i],y[j],z[k]) );
                    }
                
                    counter++;
                }
            }
        }        
    }
    else { // Need to use averages
        x = Coord->GetR();
        for (int i = 0; i < dimx; i++) {
            for (int k = 0; k < dimz; k++) {
                vector<double> Bav = GetAvField(x[i],z[k]);
                double Bav2 = pow(Bav[0],2)+pow(Bav[1],2)+pow(Bav[2],2);
                
                MagField[counter] = sqrt((1.0+betaFarrar)*Bav2);
                counter++;
                
            }
        }
    }
    
}

vector<double> TFarrarField::GetField(double x, double y, double z) {
    
   vector<double> Bret(3,0);
   
   double r = sqrt(x*x+y*y);
   //if (r > 20 || sqrt(r*r+z*z) < 1) return Bret;
   
   double Bdisk=0;
   double Bhalo=0;
   double Bxx=0;
   
   double rp=0;
   double ThetaX=0;
   double phi = atan2(y,x);
   double sinphi = sin(phi);
   double cosphi = cos(phi);
   double Ldisk = 1.0 / (1.0 + exp(-2.0*(fabs(z)-hdisk)/wdisk) );
   
   
   //if (r > 3.0) {
      if (r < 5.0) {
         Bret[0] = -bring*sinphi*(1.0-Ldisk);
         Bret[1] = bring*cosphi*(1.0-Ldisk);
      }
      else {
         
         double tanfactor = 1.0/tan(M_PI/2.-p);
         double rxneg = r*exp(-(phi-M_PI)*tanfactor);
         if (rxneg > rj[7]) rxneg = r*exp(-(phi+M_PI)*tanfactor);
         if (rxneg > rj[7]) rxneg = r*exp(-(phi+3.0*M_PI)*tanfactor);
         for (int loopcounter = 7; loopcounter>=0; loopcounter--) { if (rxneg < rj[loopcounter]) Bdisk = bj[loopcounter]; }
         
         Bdisk *= (5.0/r);
         Bret[0] = Bdisk*sin(p-phi)*( 1.0 - Ldisk );
         Bret[1] = Bdisk*cos(p-phi)*( 1.0 - Ldisk );
      }
   //}
   
   double Lhalo = 0;
   
   if (z>=0) {
      Bhalo = Bn;
      Lhalo = 1.0 / (1.0 + exp(-2.0*(r-rn)/wh) );
   }
   else {
      Bhalo = Bs;
      Lhalo = 1.0 / (1.0 + exp(-2.0*(r-rs)/wh) );
   }
   
   Bhalo *= (exp(-fabs(z)/z0)*Ldisk*(1.0-Lhalo));
   Bret[0] -= (Bhalo*sinphi);
   Bret[1] += (Bhalo*cosphi);
   
   
   double ztheta0x = fabs(z)/tan(Theta0X);
   double rpcX = rcX+ztheta0x;
   
   if (r<rpcX){ // interior region, with varying elevation angle
      rp   = r*rcX/rpcX;
      ThetaX = (z!=0) ? atan(fabs(z)/(r-rp)) : M_PI/2.;
      Bxx = BX*exp(-rp/rX)*pow(rcX/rpcX,2);
   }
   else {       // exterior region with constant elevation angle
      rp = r-ztheta0x;
      ThetaX = Theta0X;
      Bxx = BX*exp(-rp/rX)*(rp/r);
   }
   
   
   if (z<0) {
      Bret[0] -= (Bxx*cosphi*cos(ThetaX));
      Bret[1] -= (Bxx*sinphi*cos(ThetaX));
   }
   else {
      Bret[0] += (Bxx*cosphi*cos(ThetaX));
      Bret[1] += (Bxx*sinphi*cos(ThetaX));
   }
   Bret[2] += (Bxx*sin(ThetaX));
   
   Bret[0] *= 1.e-6;
   Bret[1] *= 1.e-6;
   Bret[2] *= 1.e-6;
   
   return Bret;
   
}

vector<double> TFarrarField::GetFieldCyl(double r, double phi, double z) {
   
   return TFarrarField::GetField(r*cos(phi), r*sin(phi), z);

}

vector<double> TFarrarField::GetAvField(double r, double z) {
    vector<double> Bav(3,0);
    if (sqrt(r*r+z*z) < 1 || r > 20) return Bav;
    
    const int Nav = 100;
    
    for (int i = 0; i < Nav; i++) {
        double phi = double(i)/double(Nav-1)*2.0*M_PI;
        vector<double> Breg1 = GetFieldCyl(r,phi,z);
        for (int j = 0; j < 3; j++) Bav[j] += Breg1[j]/double(Nav);
    }
    
    return Bav;
}


ToyModelField::ToyModelField(double Bx_, double By_, double Bz_, double B0turb_, TGrid* Coord, TGeometry* geom) :
  TBField(Coord),
Bx(Bx_),
By(By_),
Bz(Bz_),
B0turb(B0turb_)
{
    
    vector<double> x;
    vector<double> y;
    
    if (Coord->GetType() == "3D") {
        x = Coord->GetX();
        y = Coord->GetY();
    }
    else {
        x = Coord->GetR();
    }
    
    vector<double> z = Coord->GetZ();
	//cout << " dimx = " << dimx;
	//cout << " dimy = " << dimy;
    
    int counter = 0;
    for (int i = 0; i < dimx; i++) {
        for (int l = 0; l < dimy; l++) {
            for (int k = 0; k < dimz; k++) {
                vector<double> B = (Coord->GetType() == "3D") ? GetField(x[i],y[l],z[k]) : GetField(x[i],0,z[k]);
                
                Breg[0][counter] = B[0];
                Breg[1][counter] = B[1];
                Breg[2][counter] = B[2];
                double Breg_mod = sqrt( pow(B[0],2) + pow(B[1],2) +pow(B[2],2) );
                Breg_versors[0][counter] = Breg[0][counter]/Breg_mod;
                Breg_versors[1][counter] = Breg[1][counter]/Breg_mod;
                Breg_versors[2][counter] = Breg[2][counter]/Breg_mod;
                
                //MagField.push_back(Bh*exp(-((radius)-robs)/rB)*exp(-fabs((*zeta))/zr));
                MagField[counter] = sqrt(pow(Breg_mod,2)+pow(B[3],2)) * pow( Coord->in_LB_MagField, Coord->IsInLocalBubble(x[i],y[l],z[k]) );

                double spiral_factor_magfield = max( min( pow(geom->GetPattern(i,l,k), Coord->in_SA_MagField), Coord->in_SA_cut_MagField), 1./Coord->in_SA_cut_MagField );
                MagField[counter] *= spiral_factor_magfield;
                    
                counter++;
            }
        }
    }
}

vector<double> ToyModelField::GetField(double x, double y, double z) {
    
    vector<double> Bmixed;
    //Bmixed.push_back(B0reg);
    Bmixed.push_back(Bx);//
    Bmixed.push_back(By);
    Bmixed.push_back(Bz);
    Bmixed.push_back(B0turb);
    
    return Bmixed;
}


TSimpleField::TSimpleField(double b0, double rsun, TGrid* Coord, TGeometry* geom) : TBField(Coord) {
 
    r_scale=10.;
    z_scale=2.;
    B0=b0;
    robs=rsun;

    vector<double> x;
    vector<double> y;
    vector<double> z = Coord->GetZ();

    int counter = 0;

  if (Coord->GetType() == "3D") {
        x = Coord->GetX();
        y = Coord->GetY();
        
        for (int i = 0; i < dimx; i++) {
            for (int j = 0; j < dimy; j++) {
                for (int k = 0; k < dimz; k++) {
                    
                    vector<double> Breg1 = GetField(x[i],y[j],z[k]);
                    
                    Breg[0][counter] = Breg1[0];
                    Breg[1][counter] = Breg1[1];
                    Breg[2][counter] = Breg1[2];
                    
                    double Breg_mod = pow(Breg1[0],2)+pow(Breg1[1],2)+pow(Breg1[2],2);
                    
                    Breg_versors[0][counter] = Breg1[0]/sqrt(Breg_mod);
                    Breg_versors[1][counter] = Breg1[1]/sqrt(Breg_mod);
                    Breg_versors[2][counter] = Breg1[2]/sqrt(Breg_mod);
                    
                    MagField[counter] = sqrt(Breg_mod) * pow( Coord->in_LB_MagField, Coord->IsInLocalBubble(x[i],y[j],z[k]) );

                    double spiral_factor_magfield = max( min( pow(geom->GetPattern(i,j,k), Coord->in_SA_MagField), Coord->in_SA_cut_MagField), 1./Coord->in_SA_cut_MagField );
                    MagField[counter] *= spiral_factor_magfield;

                    counter++;
                    
                }
            }
        }        
    }
    else { 

ofstream datafile;
datafile.open("bfield_new.dat");

        x = Coord->GetR();
        for (int i = 0; i < dimx; i++) {
            for (int k = 0; k < dimz; k++) {
              
                MagField[counter] = B0*exp(-1.*(x[i]-robs)/r_scale)*exp(-1.*fabs(z[k])/z_scale);
datafile << x[i] << " " << z[k] << " " << MagField[counter] << endl;
                counter++;                
            }
        }
    }
}


vector<double> TSimpleField::GetField(double x, double y, double z) {

    vector<double> Breg1(3,0);
    
    double r = sqrt(x*x+y*y);
    if (fabs(r) < 1e-3) r = 1e-3;
    
    double theta = atan(fabs(z)/r);
    double phi = atan2(y,x)+M_PI/2.0;

    double B_simple = B0*exp(-1.*(r-robs)/r_scale)*exp(-1.*fabs(z)/z_scale);

//     cout << "[MW-DEBUG BFIELD] " << x << " " << y << " " << z << " | " << B0 << " " << r_scale << " " << z_scale << " " << robs << " " << phi << " " << theta << " | " << B_simple << endl;
    
    Breg1[0] += B_simple*cos(phi)*cos(theta);
    Breg1[1] += B_simple*sin(phi)*cos(theta);
    Breg1[2] += B_simple*sin(theta);

    return Breg1;
    
}
