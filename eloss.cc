/**
 * @file eloss.cc
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief In this file all the classes related to energy losses are implemented.
 */

#include "galaxy.h"
#include "grid.h"
#include "gas.h"
#include "input.h"
#include "eloss.h"
#include "constants.h"
#include "fitsio.h"
#include "bfield.h"
#include "geometry.h"

#include <fstream>
#include <iostream>

using namespace std;

TEnergyLoss::TEnergyLoss(TGrid* coord, Input *in) { 
    
    inp = in;
    
    dimz = coord->GetDimZ();
    dimE = coord->GetDimE();
    
    if (coord->GetType() == "3D") {
        dimy = coord->GetDimY();
        dimx = coord->GetDimX();
    }
    else {
        dimx = coord->GetDimR();
        dimy = 1;
    }
    if (in->feedback >0) cout << "creating dpdt vector... ";
    dpdt = vector<double>(dimx*dimz*dimE*dimy, 0.0);
    if (in->feedback >0) cout << "done. " << dimx*dimz*dimE*dimy << " elements." << endl;
}

TIonizationLoss::TIonizationLoss(TGrid* coord, vector<TGas*> gas, TGas* totalgas, Input* in, double A, double Z) : TEnergyLoss(coord,in) {
    
    const double factor = 1e-9*Myr;
    
    if (A != 0) { // hadrons
        
        vector<double> bet = coord->GetBeta();
        vector<double> gamma = coord->GetGamma();
        
        const double MA = A*mp*1.e3;  // nucleus mass MeV
        
        for (int k = 0; k < dimE; ++k) {
            double gammak = gamma[k];
            double betak = bet[k];
            double qmax = 2.*Mele*(gammak*gammak-1.)/(1. +2.*gammak*Mele/MA);
            double bh   = log(2.*Mele*(gammak*gammak-1.)*qmax /(EH*EH))  -2.*betak*betak;
            double bhe  = log(2.*Mele*(gammak*gammak-1.)*qmax /(EHe*EHe))-2.*betak*betak;  
            
            for (int i = 0; i < dimx; ++i) {
                for (int l = 0; l < dimy; ++l) {
                    for (int j = 0; j < dimz; ++j) {
                        int ind = coord->indexD(i,l,j);
                        double nh = totalgas->GetGas(ind) - gas[2]->GetGas(ind);
                        dpdt[index(i,l,j,k)] =  2.*PIR02MC2C*fabs(Z)*fabs(Z)/betak*nh*(bh +He_abundance*bhe) * factor / betak;
                    }
                }
            }
        }
    }
    else { // leptons
        
        vector<double> bet = coord->GetBetaEl();
        vector<double> gamma = coord->GetGammaEl();
        
        // IONIZATION LOSSES in the neutral H and He (Ginzburg 1979, p.360)
        for (int k = 0; k < dimE; ++k) {
            double gammak = gamma[k];
            double betk = bet[k];
            for (int i = 0; i < dimx; ++i) {
                for (int l = 0; l < dimy; ++l) {
                    for (int j = 0; j < dimz; ++j) {
                        int ind = coord->indexD(i,l,j);
                        double nh = totalgas->GetGas(ind) - gas[2]->GetGas(ind);
                        dpdt[index(i,l,j,k)] = 2.*PIR02MC2C/betk*(  nh*(1.0+ZHe*He_abundance)*( log(gammak-1.)- M_LN2 +1./8.
                                                                                        +2.*log(gammak*betk*Mele) )-nh*2.*log(EiH)-ZHe*nh*He_abundance*2.*log(EiHe)  ) * factor / betk;
// cout << "[MW-DEBUG ION] (A=" << A << ") " << k << " " << nh << " " << dpdt[index(i,l,j,k)] << endl;
                    }
                }
            }
        }
    }
}


TCoulombLoss::TCoulombLoss(TGrid* coord, vector<TGas*> gas, Input* in, double A, double Z) : TEnergyLoss(coord,in) {
    
    const double factor = 1e-9*Myr;
    if (A != 0) {
        
        vector<double> bet = coord->GetBeta();
        vector<double> gamma = coord->GetGamma();
        
        const double MA = A*mp*1.e3;  // nucleus mass MeV
        
        for (int k = 0; k < dimE; ++k) {
            double betk = bet[k];
            double gammak = gamma[k];
            double bet3 = pow(betk,3);
            double we = bet3/(xm + bet3);
            
            for (int i = 0; i < dimx; ++i) {
                for (int l = 0; l < dimy; ++l) {
                    for (int j = 0; j < dimz; ++j) {
                        int ind = coord->indexD(i,l,j);
                        
                        double coullog = -log(4.*PIR0H2C2*gas.back()->GetGas(ind) *(MA+2.*gammak*Mele)
                                              /(4.*gammak*gammak*bet3*betk*Mele*Mele*MA)) /2.;
                        dpdt[index(i,l,j,k)] = 4.*PIR02MC2C*fabs(Z)*fabs(Z)*gas.back()->GetGas(ind)/betk *coullog *we * factor/betk;
                        
                    }
                }
            }
        }
    }
    else {
        
        vector<double> bet = coord->GetBetaEl();
        vector<double> gamma = coord->GetGammaEl();
        
        // Coulomb energy losses in the cold H plasma limit (Ginzburg 1979, p.361)

        for (int k = 0; k < dimE; ++k) {
            double gammak = gamma[k];
            double betk = bet[k];
            
            for (int i = 0; i < dimx; ++i) {
                for (int l = 0; l < dimy; ++l) {
                    for (int j = 0; j < dimz; ++j) {
                        double coullog2 = 0.;
                        double nhi = gas.back()->GetGas(i,l,j);
                        if(nhi > 0.) 
                            coullog2=log(gammak*Mele/nhi*Mele/(4.*Pi*Rele*H2PiC*H2PiC))-3./4.;
                        dpdt[index(i,l,j,k)] = 2.*PIR02MC2C *nhi/betk *coullog2 * factor/betk;
// cout << "[MW-DEBUG COULOMB] (A=" << A << ") " << k << " " << coullog2 << " " << nhi << " " << dpdt[index(i,l,j,k)] << endl;
                    }
                }
            }
        }
    }
}

TBremsstrahlungLoss::TBremsstrahlungLoss(TGrid* coord, vector<TGas*> gas, TGas* totalgas, Input* in) : TEnergyLoss(coord,in) {
    
    const double factor = 1e-9*Myr;
    
    // Bremsstrahlung energy losses in neutral gas (Ginzburg 1979, p.386,409)
    double brem1 = 0.0;
    
    vector<double> bet = coord->GetBetaEl();
    vector<double> gamma = coord->GetGammaEl();
    
    for (int i = 0; i < dimx; ++i) {
        for (int l = 0; l < dimy; ++l) {
            for (int j = 0; j < dimz; ++j) {
                
                int ind = coord->indexD(i,l,j);
                double nhi = gas.back()->GetGas(ind);
                double nh = totalgas->GetGas(ind) - nhi;
                double nhe_cm3 = He_abundance*nh; // 1/cm^3, He number densit
                
                for (int k = 0; k < dimE; ++k) {
                    
                    double gammak = gamma[k];
                    
                    if(gammak < gam1) brem1 = gammak*Mele *4.*AFR02MC2C/Mele
                        *(2.*nh+ZHe*(ZHe+1.)*nhe_cm3)*(log(2.*gammak)-1./3.);
                    else 
                        if(gammak > gam2) brem1 = 1.e6 *gammak*Mele*C*(nh*MH/TH +nhe_cm3*MHe/THe);
                        else    // linear interpolation provides max 10% error to an exact value
                            brem1 = gam1*4.*AFR02MC2C
                            *(2.*nh+ZHe*(ZHe+1.)*nhe_cm3)*(log(2.*gam1)-1./3.)
                            *(gam2-gammak)/(gam2-gam1) +(gammak-gam1)/(gam2-gam1) 
                            *1.e6 *gam2*Mele *C*(nh*MH/TH +nhe_cm3*MHe/THe);
                    dpdt[index(i,l,j,k)] = (brem1 + gammak*Mele*4.*2.*AFR02MC2C/Mele*nhi*(log(2.*gammak)-1./3.) ) * factor/bet[k];
// cout << "[MW-DEBUG BREMS] " << k << " " << brem1 << " " << gammak << " " << dpdt[index(i,l,j,k)] << endl;

                    // Bremsstrahlung energy losses in hydrogen plasma (Ginzburg 1979, p.408)
                    //dpdt[index(i,j,k)] += gammak*Mele*4.*2.*AFR02MC2C/Mele*nhi*(log(2.*gammak)-1./3.) * 1.e-9*Myr/bet[k];
                }
            }
        }
    }
}

TSynchrotronLoss::TSynchrotronLoss(TGrid* coord, TBField* B, Input* in) : TEnergyLoss(coord,in) {
    //cout << "test synch " << endl;
    vector<double> bet = coord->GetBetaEl();
    vector<double> gamma = coord->GetGammaEl();

    vector<double> x_vec = coord->GetX();
    vector<double> y_vec = coord->GetY();
    vector<double> z_vec = coord->GetZ();
    vector<double> E_vec = coord->GetEk();		

    int dimx = coord->GetDimX();	 
    int dimy = coord->GetDimY();	 
    int dimz = coord->GetDimZ();	 

//MW130622
    unsigned int ixsun = 0;
    unsigned int iysun = 0;
    unsigned int izsun = 0;
    
    while(x_vec[++ixsun] <= inp->xobs){}
    ixsun--;
    
    if(coord->GetType()=="3D")
    {
        while(y_vec[++iysun] <= in->yobs){}
        iysun--;
    }
    
    while(z_vec[++izsun] <= in->zobs){}
    izsun--;

    //cout << endl;
    //cout << " *** Synch losses *** " << endl;
    //cout << ixsun << " " << iysun << " " <<  izsun << endl;

    //cout << "test synch" << endl;
    for (int ix = 0; ix < dimx; ++ix) {
        for (int iy = 0; iy < dimy; ++iy) {
            for (int iz = 0; iz < dimz; ++iz) {
                double bendens = B->GetEnDensity(ix,iy,iz);
		//cout << "Bendens = " << bendens << endl;
                for (int ip = 0; ip < dimE; ++ip) {
                    dpdt[index(ix,iy,iz,ip)] = 32./9.*PIR02C*bendens*(gamma[ip]*gamma[ip]-1.)* 1.e-9*Myr/bet[ip];
	            //cout << "[MW-DEBUG BF] " << ip << " " << bendens << " " << gamma[ip] << " " << dpdt[index(ix,iy,iz,ip)] << endl;
		    //if (ix == ixsun && iy == izsun && iz == iysun && ip%5 == 0)	
	                  //cout << x_vec[ix] << " " << y_vec[iy] << " " << z_vec[iz] << "; E= " << E_vec[ip] << "; dpdt= " << dpdt[index(ix,iy,iz,ip)] << " GeV/Myr " << endl;	
                }
            }
        }
    }
}

TICSLoss::TICSLoss(TGrid* coord, TISRF* isrf, Input* in) : TEnergyLoss(coord,in) {

    vector<double> x_vec = coord->GetX();
    vector<double> y_vec = coord->GetY();
    vector<double> z_vec = coord->GetZ();
    vector<double> E_vec = coord->GetEk();		

    int dimx = coord->GetDimX();	 
    int dimy = coord->GetDimY();	 
    int dimz = coord->GetDimZ();	 

//MW130622
    unsigned int ixsun = 0;
    unsigned int iysun = 0;
    unsigned int izsun = 0;
    while(x_vec[++ixsun] <= in->xobs){}
    ixsun--;
    
    if(coord->GetType()=="3D")
    {
        while(y_vec[++iysun] <= in->yobs){}
        iysun--;
    }
    
    while(z_vec[++izsun] <= in->zobs){}
    izsun--;

    //cout << endl;
    //cout << " *** IC losses *** " << endl;
    //cout << ixsun << " " << iysun << " " <<  izsun << endl;
    
    const double factor = eV_to_erg / h_planck * isrf->GetDnu()*Myr*1e-3;
    vector<double> nuarr = isrf->GetNuArray();
    
    //cout << " Starting loop over x,y,z,E to compute IC losses " << endl;
    for (int i = 0; i < dimx; ++i) {
        for (int l = 0; l < dimy; ++l) {
            for (int j = 0; j < dimz; ++j) {
                for (int k = 0; k < dimE; ++k) {
                    int ind = index(i,l,j,k);
                    // Need to integrate over frequency.
                    for (unsigned int nu_ind = 0; nu_ind < nuarr.size(); ++nu_ind) {
                        dpdt[ind] += factor*isrf->GetISRF(i, l, j, nu_ind)/nuarr[nu_ind]*isrf->GetElossCompton(nu_ind,k);
                        //if (i == ixsun && j == izsun && l == iysun)	
	                  //cout << x_vec[i] << " " << y_vec[l] << " " << z_vec[j] << "; E= " << E_vec[k] << "; dpdt= " << dpdt[index(i,j,j,k)] << " GeV/Myr " << endl;	
                    }
                }
            }
        }
    }
}


TISRF::TISRF(TGrid* Coord, string filename, TGeometry* geom, Input* in) {
    
    int status = 0;
    
    fitsfile *fptr = NULL;
    if ( fits_open_file(&fptr, filename.c_str(), READONLY, &status) ) fits_report_error(stderr, status);
    
    int NAXIS, NAXIS1, NAXIS2, NAXIS3, NAXIS4;
    float CRVAL1,CRVAL2,CRVAL3,CDELT1,CDELT2,CDELT3;
    char comment[100];
    
    if( fits_read_key(fptr,TINT,"NAXIS" ,&NAXIS ,comment,&status) ) fits_report_error(stderr, status);
    if( fits_read_key(fptr,TINT,"NAXIS1",&NAXIS1,comment,&status) ) fits_report_error(stderr, status);
    if( fits_read_key(fptr,TINT,"NAXIS2",&NAXIS2,comment,&status) ) fits_report_error(stderr, status);
    if( fits_read_key(fptr,TINT,"NAXIS3",&NAXIS3,comment,&status) ) fits_report_error(stderr, status);
    if( fits_read_key(fptr,TINT,"NAXIS4",&NAXIS4,comment,&status) ) fits_report_error(stderr, status);
    
    if( fits_read_key(fptr,TFLOAT,"CRVAL1",&CRVAL1,comment,&status) ) fits_report_error(stderr, status);
    if( fits_read_key(fptr,TFLOAT,"CRVAL2",&CRVAL2,comment,&status) ) fits_report_error(stderr, status);
    if( fits_read_key(fptr,TFLOAT,"CRVAL3",&CRVAL3,comment,&status) ) fits_report_error(stderr, status);
    if( fits_read_key(fptr,TFLOAT,"CDELT1",&CDELT1,comment,&status) ) fits_report_error(stderr, status);
    if( fits_read_key(fptr,TFLOAT,"CDELT2",&CDELT2,comment,&status) ) fits_report_error(stderr, status);
    if( fits_read_key(fptr,TFLOAT,"CDELT3",&CDELT3,comment,&status) ) fits_report_error(stderr, status);
    
    dimr_isrf = NAXIS1;
    dimz_isrf = NAXIS2;
    dimnu = NAXIS3;
    ncomp = NAXIS4;
    
    long nelements=dimr_isrf*dimz_isrf*dimnu*ncomp;
    long felement=1;
    float *isrf_in = new float[nelements]();
    float nulval=0;
    int anynul;
    
    if( fits_read_img(fptr,TFLOAT,felement,nelements,&nulval,isrf_in,&anynul,&status) ) fits_report_error(stderr, status);
    
    if (fits_close_file(fptr, &status)) fits_report_error(stderr, status);
    
    nu_array.assign(dimnu,0.0);

    // microns -> cm; nu=c/lambda
    for(int inu = 0; inu < dimnu; inu++) 
		nu_array[dimnu-1-inu]=C/(pow(10.,1.*CRVAL3+inu*CDELT3)*1.0e-4);  // nu is in Hz

    //cout << "nu vector for ISRF" << endl;	
    //for(int inu = 0; inu < dimnu; inu++) 
	         //cout << nu_array[inu] << " ";		
    //cout << endl;	    

    vector<double> x;
    vector<double> y;
    
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
    
    vector<double> z = Coord->GetZ();
    
    dimz = z.size();  
    
    ISRField.assign(dimnu*dimx*dimy*dimz, 0.0);

    vector<double> ISRF_r_profile(dimx*dimnu,0.0);
    vector<double> ISRF_z_profile(dimz*dimnu,0.0);

    if (in->ISRF_model == GalpropISRF) {

      for(int i = 0; i < ncomp; i++) {
        for(int ix = 0; ix < dimx; ix++) {    
            for(int iy = 0; iy < dimy; iy++) {    
                double r = 0;
                if (y.size()) r = sqrt(x[ix]*x[ix]+y[iy]*y[iy]);
                else r = x[ix];
                for(int iz = 0; iz < dimz; iz++) {
                    int irr=(int)((r -CRVAL1) /CDELT1+0.5);
                    int izz=(int)((fabs(z[iz])-CRVAL2) /CDELT2+0.5);
                    if(irr>NAXIS1-2) irr=NAXIS1-2;
                    if(izz>NAXIS2-2) izz=NAXIS2-2;
                    float rr=CRVAL1+irr*CDELT1;
                    float zz=CRVAL2+izz*CDELT2;
                    
                    for(int inu = 0; inu < dimnu; inu++) {
                        float v1=isrf_in[isrf_index(irr  ,izz  ,inu,i)];
                        float v2=isrf_in[isrf_index(irr+1,izz  ,inu,i)];
                        float v3=isrf_in[isrf_index(irr  ,izz+1,inu,i)];
                        float v4=isrf_in[isrf_index(irr+1,izz+1,inu,i)];
                        float v5=v1+(v2-v1)*(r-rr)/CDELT1;
                        float v6=v3+(v4-v3)*(r-rr)/CDELT1;
                        float value=v5+(v6-v5)*(fabs(z[iz])-zz)/CDELT2;

			//value *= 0.5;

                        if(value<0.0) value=0.0;

			if (i==0 || i==1)
				value *= 1.0;

			if (Coord->GetType() == "3D") { 

				double spiral_factor_isrf = 1.;
			
	     			if(i==0) //MW130625: modeling starlight component after spiral arms
	     			{
				       spiral_factor_isrf = max( min( pow(geom->GetPattern(ix,iy,iz), Coord->in_SA_ISRFStar), Coord->in_SA_cut_ISRF), 1./Coord->in_SA_cut_ISRF );
				}
				else if(i==1) //MW130625: modeling dust component after spiral arms
	     			{
				       spiral_factor_isrf = max( min( pow(geom->GetPattern(ix,iy,iz), Coord->in_SA_ISRFDust), Coord->in_SA_cut_ISRF), 1./Coord->in_SA_cut_ISRF );
     				}

		        	value *= spiral_factor_isrf;
			
			        //MW 130822: for Local Bubble, there is only one parameter, but still, don't change the CMB component!
			        if(i!=2) value *= pow( Coord->in_LB_ISRF, Coord->IsInLocalBubble(x[ix],y[iy],z[iz]) );
			}

                        // reverse scale from wavelength to frequency
                        ISRField[index(dimnu-1-inu, ix, iy, iz)] += value;
                      
                    }  //  inu
                }  //  iz
            }  //  iy
          }  //  ix
	}

        unsigned int irsun = (unsigned int) ((in->robs-x.front())/(x.back()-x.front())*(double)(dimx-1));
	unsigned int izsun = (unsigned int) ((in->zobs-z.front())/(z.back()-z.front())*(double)(dimz-1));

        for(int inu = 0; inu < dimnu; inu++) {
            for(int ix = 0; ix < dimx; ix++) {    
		ISRF_r_profile[inu*dimx+ix] = ISRField[index(inu, ix, dimy/2, izsun)] ;
		if (ix == irsun) {
			for(int iz = 0; iz < dimz; iz++) {
			    ISRF_z_profile[inu*dimz+iz] = ISRField[index(inu, irsun, dimy/2, iz)];
		         }
		}
	    }
	}		

        /*cout << "*** Galprop ISRF nu array *** " << endl;
        for (int inu=0; inu < dimnu; inu++) 
		cout << nu_array[inu] << " " ;

	cout << endl << "*** Galprop IRSF at Sun position ***" << endl;
        for (int inu=0; inu < dimnu; inu++) 
		cout << ISRField[index(inu,irsun,dimy/2,izsun)] << " " ;
	
	cout << endl << "*** Galprop IRSF profile along R ***" << endl;
        for (int inu=0; inu < dimnu; inu=inu+6) {
		cout << "* ------------------------- " << endl;
	        cout << "* frequency = " << nu_array[inu] << " Hz"<< endl;
		cout << "* ------------------------- " << endl;
	        cout << "* " ;
                for(int ix = 0; ix < dimx; ix=ix+4) {    
		  cout << " x = " << x[ix] << " --> " << ISRF_r_profile[inu*dimx+ix] << " ||";	
		}
		cout << endl;
	}

	cout << endl << "*** Galprop IRSF profile along z ***" << endl;
        for (int inu=0; inu < dimnu; inu=inu+6) {
		cout << "* ------------------------- " << endl;
	        cout << "* frequency = " << nu_array[inu] << " Hz"<< endl;
		cout << "* ------------------------- " << endl;
	        cout << "* " ;
                for(int iz = 0; iz < dimz; iz=iz+2) {    
		  cout << " z = " << z[iz] << " --> " << ISRF_z_profile[inu*dimz+iz] << " ||";	
		}
		cout << endl;
	}*/

    } 
    else if (in->ISRF_model == UniformISRF) {
	//implemented by D.Gaggero -- september 2013
	//WARNING: NO SPIRAL PATTERN IS APPLIED IN THIS MODE. Fix this!

	int n_isrf = 6;
	vector<double> t0isrf;
 	vector<double> factor;	
	t0isrf.assign(n_isrf,0.);
	factor.assign(n_isrf,0.);

	factor[0] = 1.;
        t0isrf[0]=2.726;
	
	factor[1] = 4.5e-5;
        t0isrf[1]=33.07;
	
	factor[2] = 1.2e-9;
        t0isrf[2]=313.32;
	
	factor[3] = 7.03e-13;
      	t0isrf[3]=3249.3; 
	
	factor[4] = 3.39e-14;
	t0isrf[4]=6150.4; 
	
	factor[5] = 8.67e-17;
      	t0isrf[5]=23209.0; 

	double hPlanck = 4.135638e-15;  // eV/Hz
	double kB      = 8.6173324e-5;  // eV/K 
 	// C is in cm/s from constants.h  

        vector<double> ISRF_vector;
	ISRF_vector.assign(dimnu,0.);

	cout << "*** ISRF vector from Delahaye et al. 2010 *** " << endl;
        for(int inu = 0; inu < dimnu; inu++) {
                // ISRF is in eV/(cm3 Hz)
		for (int i_isrf = 0; i_isrf < n_isrf; i_isrf++) {
  		  double u_nu = factor[i_isrf] * (8.* M_PI * hPlanck)/(C*C*C) * pow(nu_array[inu],3.) * ( 1./ ( exp(hPlanck*nu_array[inu]/(kB*t0isrf[i_isrf])) - 1. ) );
                  ISRF_vector[inu] += nu_array[inu] * u_nu;
		}
  		cout << ISRF_vector[inu] << " ";
        }
        cout << endl;


        for(int ix = 0; ix < dimx; ix++) {    
          for(int iy = 0; iy < dimy; iy++) {    
            for(int iz = 0; iz < dimz; iz++) {
              for(int inu = 0; inu < dimnu; inu++) {
                 ISRField[index(inu, ix, iy, iz)] = ISRF_vector[inu];
	      }	
	    }
	  }
	}

    }	
    
    // microns -> cm; nu=c/lambda
    //for(int inu = 0; inu < dimnu; inu++) nu_array[dimnu-1-inu]=C/(pow(10.,1.*CRVAL3+inu*CDELT3)*1.0e-4); --> already did before!
    
    delete[] isrf_in;
    
    Dnu = log(nu_array[1]/nu_array[0]);
    
    // Init Loss Compton
    int dimE = Coord->GetDimE();
    
    const double convfact = (h_planck * erg_to_eV * 1.0e-6 / Mele);
    
    vector<double> nu_me(nu_array);
    for (int i = 0; i < dimnu; i++) nu_me[i] *= convfact;
    
    vector<double> gamma_el = Coord->GetGammaEl();
    vector<double> beta_el = Coord->GetBetaEl();
    e_loss_Compton = vector<double>(dimE*dimnu,0.0);
    for (int i = 0; i < dimE; i++) {
        long double gam = gamma_el[i];
        long double bet = beta_el[i];
        long double gbar = gam*(1.0+bet);
        long double gambet = gam*bet;
        for (int k = 0; k < dimnu; k++) {
            long double num = nu_me[k];
            long double value = 0.0;
            if (gam*pow(num*gam,3)< 4e-5) 
                value = 
                Mele*0.5*Pi*Rele*Rele*C/pow(gam,2)/bet/pow(num,2)*(
                                                                   56.0/25.0*pow(num,5)*gambet*gam*(7.0*pow(gambet,4)
                                                                                                    +5.0*gam*gam*(-5.0+7.0*gam*gam)
                                                                                                    +5.0*gambet*gambet*(-5.0+14.0*gam*gam))
                                                                   -16.0/45.0*pow(num,4)*gambet*gam*gam*(-48.0 +63.0*gam*gam
                                                                                                         +bet*bet*(-16.0+63.0*gam*gam))
                                                                   +16.0/9.0*pow(num,3)*gambet*gam*(-3.0+(3.0+bet*bet)*gam*gam)     
                                                                   );
            else 
                value = Mele*0.5*Pi*Rele*Rele*C/pow(gam,2)/bet/pow(num,2)
                * ( gam * (f1(num*gbar)-f1(num/gbar)) -
                   num * (f2(num*gbar)-f2(num/gbar))
                   );
            e_loss_Compton[i*dimnu+k] = (value < 0.0) ? 0.0 : value;
        }
    }
    
    return ;
}
