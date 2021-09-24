/**
 * @file    dragon.cc
 * @author  Daniele Gaggero, Carmelo Evoli, Pedro De La Torre Luque
 * @email   daniele.gaggero@gmail.com
 * @brief   Implementation of the DRAGON class. See the .h file
 */

#include "crevolutor.h"

#include "grid.h"
#include "gas.h"
#include "galaxy.h"
#include "xsec.h"
#include "utilities.h"
#include "input.h"
#include "sources.h"
#include "diffusion.h"
#include "config.h"

#include <gsl/gsl_integration.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <string>
#include <sstream>
#include <fstream>

using namespace std;

TCREvolutor::TCREvolutor(Galaxy* gal1) {
   
   gal = gal1;
   in = gal1->GetInput();
   
   coord = gal->GetCoordinates();
   dimr = coord->GetDimR();
   dimz = coord->GetDimZ();
   dimE = coord->GetDimE();
   
   
   riac1    = vector<double>(dimr*dimz*dimE,0);
   riac2    = vector<double>(dimr*dimz*dimE,0);
   riac3    = vector<double>(dimr*dimz*dimE,0);
   Pdotup   = vector<double>(dimr*dimz*dimE,0);
   Pdotdown = vector<double>(dimr*dimz*dimE,0);
}

//2D engine -------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------   
void TCREvolutor::Run(vector<double>& N, vector<double>& N_previous, TInelasticCrossSection *xsec_object, const vector<double>& dpdt, const vector<double>& SecSource_, const vector<double>& spectrum, double A, double Z, double lifetime, int daughter, bool SecEl, int K_electron, bool isDM, bool isextra) {
//-----------------------------------------------------------------------------------------------------------------------------------   

  if (in->feedback >0) cout << "This is the 2D engine! " << endl;

   int uid = 1000*Z+A;
  
   double injfactor = CalcInjFactor(K_electron, SecEl, isDM, A, Z);
   
   if (in->feedback >0) cout << endl << "*** We are propagating this particle: *** " << endl;
   if (in->feedback >0) cout << "A = " << A << " Z = " << Z << " SecEl = " << SecEl << " K_electron = " << K_electron << " injfactor = " << injfactor << endl;
   if (A==0 && Z==-1 && SecEl==0 && in->feedback >0)
      cout << "Primary electrons" << endl;
   if (A==0 && Z==-1 && SecEl==1 && in->feedback >0)
      cout << "Secondary electrons" << endl;
   if (in->feedback >0) cout << "Starting propagation... " << endl;
   
   vector<double>::const_iterator itsecsource = max_element(SecSource_.begin(), SecSource_.end());
   if (injfactor == 0 && (*itsecsource) == 0 && in->feedback >0) {
      cout << "No primary nor secondary source. Returning." << endl;
      return ;
   }

   vector<double> xsec;
   map<pair<int,int>, vector<double> > xsec_extended;

   /*if (in->spallationxsec == Fluka)
	xsec_extended = xsec_object->GetXSec_extended();
   else*/
        xsec = xsec_object->GetXSec();
   
   cout << "inelastic xsec loaded" << endl;

   TGas*            totalgas = gal->GetTotalGas();
   map<int, double> gas_abundances = gal->GetGasAbundances();

   TSource* source = GetSourceTerm(isDM, isextra);
   TConvectionVelocity* vC = gal->GetVC();
   
   TDiffusionCoefficient* dperp = NULL;
   TReaccelerationCoefficient* dpp = NULL;
   vector<double> totmomentum;
   
   if (A == 0) {
      dperp = gal->GetDiffCoeffEl();
      dpp = gal->GetReaccCoeffEl();
      totmomentum = coord->GetMomentumEl();
   }
   else {
      dperp = gal->GetDiffCoeff();
      dpp = gal->GetReaccCoeff();
      totmomentum = coord->GetMomentum();
   }
   
   // Following quantities depend on rigidity, but are computed independently of the nucleus, to save memory and time. Here the necessary additional nucleus-dependent factors are restored.
   
   int bin=0;
   for (int k = 0; k < dimE; k++) if(totmomentum[k]<dperp->Getrho_b()) bin=k;
   
   double dperpfactor[dimr][dimE];
   double dppfactor[dimr][dimE];

   for (int i = 0; i < dimr; i++) {
     for (int k = 0; k < dimE; k++){
       if (A==0)    dperpfactor[i][k] = 1.0;
     
       else{
	 
	 if(k<=bin){
	   if (in->VariableDelta == true) dperpfactor[i][k] = pow(A/fabs(Z), in->delta_B + in->delta_A*in->robs);
	   else  	   dperpfactor[i][k] = pow(A/fabs(Z),dperp->GetDelta());
	 }
	 
	 else{
	   if (in->VariableDelta == true)     dperpfactor[i][k] = pow(A/fabs(Z),in->delta_B + in->delta_A*in->robs);
	   else   dperpfactor[i][k] = pow(A/fabs(Z),dperp->GetDelta_h());
	 }
	 
       }
       
       dppfactor[i][k]=1.0/dperpfactor[i][k];
     
     }
   }
   
   
   vector<double> gamma(coord->GetGamma());
   
   /* Some initialization */
   
   int ind = -1;
   int indspat = -1;
   int indspat_up=-1;//for nonequidistant grid, Iris 02/25/13
   int indspat_down=-1;//for nonequidistant grid, Iris 02/25/13
   double ds_up=0;
   double ds_down=0;
   double ds_double=0;
   
   int i = 0;
   int j = 0;
   int k = 0;
   int ip = 0;
   
   const double dvcdz = (vC==NULL) ? 0.0 : vC->GetDvdz()/3.0;
   
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,indspat,ind) schedule(dynamic) num_threads(OMP_NUM_THREADS)
#endif
   for (i = 0; i < dimr; i++) {
      for (j = 0; j < dimz; j++) {
         indspat = coord->index(i,j);
         double dpp1 = 0.0;
         
         for (k = 1; k < dimE-1; k++) {
	   if (dpp) dpp1 = dppfactor[i][k]*dpp->GetReaccelerationCoefficient(indspat);
            
            double momentumup = totmomentum[k+1];
            double momentumfix = totmomentum[k];
            double momentumdown = totmomentum[k-1];
            double upfix = momentumup-momentumfix;
            double updown = momentumup-momentumdown;
            double fixdown = momentumfix-momentumdown;
            
            ind = indspat*dimE+k;
            Pdotup[ind] = dpdt[ind+1]/(upfix);
            Pdotdown[ind] = dpdt[ind]/(upfix);
            if (A != 0) {
               Pdotup[ind] /= A;
               Pdotdown[ind] /= A;
            }
            if (vC) {
               Pdotup[ind] += momentumup*dvcdz/(upfix);
               Pdotdown[ind] += momentumfix*dvcdz/(upfix);
            }
            
            if (dpp) {
               double dppfix = dpp1*dpp->GetSpectrum(k);
               double dppdown = dpp1*dpp->GetSpectrum(k-1);
               riac1[ind] = 2.0*dppfix/(updown)/(upfix);
               riac2[ind] = ( - (dppfix-dppdown)/(fixdown*fixdown) +
                             2.0*dppfix*(1.0/(updown)*(1.0/(upfix) + 1.0/(fixdown)) +
                                         1.0/momentumfix/(fixdown))
                             );
               riac3[ind] = ( - (dppfix-dppdown)/(fixdown) +
                             2.0*(dppfix/(updown) +
                                  dppdown/momentumdown))/(fixdown);
            }
         }
         ind = indspat*dimE;
         Pdotup[ind] = dpdt[ind+1]/(totmomentum[1]-totmomentum[0]);
         Pdotdown[ind] = dpdt[ind]/(totmomentum[1]-totmomentum[0]);
         if (A != 0) {
            Pdotup[ind] /= A;
            Pdotdown[ind] /= A;
         }
         
         if (vC) {
            Pdotup[ind] += totmomentum[1]*dvcdz/(totmomentum[1]-totmomentum[0]);
            Pdotdown[ind] += totmomentum[0]*dvcdz/(totmomentum[1]-totmomentum[0]);
         }
         
         // Boundary conditions
         Pdotup[ind+dimE-1] = Pdotup[ind+dimE-2];
         Pdotdown[ind+dimE-1] = Pdotdown[ind+dimE-2];
         if (dpp) {
            riac1[ind] = riac1[ind+1];
            riac2[ind] = riac2[ind+1];
            riac3[ind] = riac3[ind+1];
            riac1[ind+dimE-1] = riac1[ind+dimE-2];
            riac2[ind+dimE-1] = riac2[ind+dimE-2];
            riac3[ind+dimE-1] = riac3[ind+dimE-2];
         }
      }
   }
   
   double dt = in->dtmax;
   
   double dtbar = 0.0;
   
   double halfdtbar = 0.0;
   double halfdt = 0.0;

   double halfdt_dperp_factor[dimr][dimE];
      for (int ii = 0; ii < dimr; ii++) {
	for(int e=0;e<dimE;e++) halfdt_dperp_factor[ii][e]=0.; }
   
   const double decay = (daughter!=0);
   int Niter = 0;
   
   double value = 0.0;
   
   long int count_temp = 0 ;
   
   while (dt > in->dtmin) {
      dtbar = dt/p;
      
      count_temp++;

      cout << "dt = " << dt << endl;  
      
      halfdtbar = 0.5*dtbar;
      halfdt = 0.5*dt;
      for (int ii = 0; ii < dimr; ii++) {
	for(int e=0;e<dimE;e++) halfdt_dperp_factor[ii][e] = halfdt*dperpfactor[ii][e];}
      
      for (Niter = 0; Niter < in->Nrept; ++Niter) {
         
         /*************************************************/
         // Here propagation in Z direction starts.
	 // ATTENTION -- there was a bug here before -- all indexes have to be declared as private (including ip) for parallelization!! Otherwise there are runtime errors on some systems
#ifdef _OPENMP
#pragma omp parallel default(shared) private(i,j,k,ip,indspat,ind,value) num_threads(OMP_NUM_THREADS)
#endif
         {
            vector<double> Rzz(dimz, 0.0);
            vector<double> dzz(dimz, 0.0);
            vector<double> uodzz(dimz, 0.0);
            vector<double> lodzz(dimz, 0.0);
            vector<double> yy(dimz, 0.0);
            
            
            vector<double> Rrr(dimr, 0.0);
            vector<double> drr(dimr, 0.0);
            vector<double> uodrr(dimr, 0.0);
            vector<double> lodrr(dimr, 0.0);
            vector<double> xx(dimr, 0.0);
            
            vector<double> de(dimE, 0.0);
            vector<double> ee(dimE, 0.0);
            vector<double> odeu(dimE, 0.0);
            vector<double> oded(dimE, 0.0);
            vector<double> Re(dimE, 0.0);
            
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif

            //for (ip = 0; ip < dimE; ip++) 
	    //cout << "\n\n These are my lines, ip = " << ip << " while xsec[ip] = " << xsec[ip] << "A = " << A << " Z = " << Z << "\n\n" <<endl;			


            for (ip = 0; ip < dimE; ip++) {
               
               //double halfdtbarxseck = halfdtbar*xsec[ip];
	       double halfdtbar_xsec;
	       //if (in->spallationxsec != Fluka)
	          halfdtbar_xsec = halfdtbar*xsec[ip];
		  //cout << "\n\n These are my lines, ip = " << ip << " while xsec[ip] = " << xsec[ip] << "\n\n" <<endl;			

	       //else {
		//  cout << "ip = " << ip <<endl;			
		//  halfdtbar_xsec = halfdtbar*( xsec_extended[make_pair(uid,1001)][ip]*gas_abundances[1001] +		
		//			xsec_extended[make_pair(uid,2004)][ip]*gas_abundances[2004]);
// + xsec_extended[make_pair(uid,6012)][ip]*gas_abundances[6012]);
		//}
               double halfdtbarlifetimegammak = halfdtbar*(decay)/(lifetime*gamma[ip]);
               double dtbarinjfactorspeck = dtbar*injfactor*spectrum[ip];
               double sp = dperp->GetSpectrum(ip);
               
               for (j = 0; j < dimr-1; j++) {
                  for (i = 0; i < dimz; i++) {
                     
                     //MW130705: consistency to 3D case
                     indspat = coord->index(j,i);
                     ind = indspat*dimE + ip;
                     
                     double CNalphaz1 = dperp->GetCNdiff_alpha1_z(ind);
                     double CNalphaz2 = dperp->GetCNdiff_alpha2_z(ind);
                     double CNalphaz3 = dperp->GetCNdiff_alpha3_z(ind);
                     
                     double vCi = 0.0; // vC(i)
                     double vCi1 = 0.0; // vC(i+1)
                     double vC1i = 0.0; // vC(i-1)
                     
                     if (vC)
                     {
                        vC1i = vC->GetCNconv_alpha1_z(indspat);
                        vCi  = vC->GetCNconv_alpha2_z(indspat);
                        vCi1 = vC->GetCNconv_alpha3_z(indspat);
                     }
                     
		     uodzz[i] = -CNalphaz3*halfdt_dperp_factor[j][ip] - halfdt*vCi1;
                     lodzz[i] = -CNalphaz1*halfdt_dperp_factor[j][ip] - halfdt*vC1i;

		     //totalgas->GetGas(indspat)*halfdtbarxseck
		     double gas_xsec = halfdtbar_xsec * totalgas->GetGas(indspat);	
		     dzz[i] = 1. + CNalphaz2*halfdt_dperp_factor[j][ip] + gas_xsec + halfdtbarlifetimegammak + halfdt*vCi;
                     
                     //         cout << "[MW-DEBUG] " << j << " " << i << " | " << CNalphaz1 << " " << CNalphaz2 << " " << CNalphaz3 << " " << vC1i << " " << vCi << " " << vCi1 << " | " << dzz[i] << endl;
                     
                     Rzz[i] = N[ind] * (2. - dzz[i]) + source->GetSource(indspat)*dtbarinjfactorspeck + dtbar*SecSource_[ind];
                     if (i < dimz-1) Rzz[i] -= N[ind+dimE] * uodzz[i] ;
                     if (i > 0) Rzz[i] -= N[ind-dimE] * lodzz[i];
                     
                  } // for i
                  
                  // Compute N
		  //cout << "test" << endl;
//#pragma omp critical
//{	
//		  if ((j==0) && (ip==0))
//		     cout<<"dt = "<<dt<<"; iteration = "<<Niter<<"; --> DEBUG: "<<lodzz[0]<<" "<<dzz[0]<<" "<<uodzz[0]<<" "<<Rzz[0]<<" "<<yy[0]<<"; dimz = "<<dimz<<endl;
//}
		  //exit(-1);
                  //Utility::solve_tridag(&(lodzz[0]), &(dzz[0]), &(uodzz[0]), &(Rzz[0]), &(yy[0]), dimz);
		  Utility::solve_tridag(lodzz, dzz, uodzz, Rzz, yy, dimz);
                  for (i = dimz-2; i > 0; --i) {
                     value = yy[i];
                     N[index(j,i,ip)] = (value > 0) ? value : 0.0;
                  }
               } // for k
            } // for j
            
            /***********************************************/
            /*********************************************/
            
            // Here propagation in Momentum direction starts: Reacceleration and/or energy losses
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for (i = 0; i < dimr-1; i++) {
               for (j = 1; j < dimz-1; j++) {
                  
                  indspat = coord->index(i,j);
                  //double totgas = totalgas->GetGas(indspat);
		  double gas_tot = totalgas->GetGas(indspat);
                  double dtbarprimsource = dtbar*injfactor*source->GetSource(indspat);
                  ind = indspat*dimE;
                  
                  for (ip = 0; ip < dimE; ip++) {
                     
                     int ind1 = ind + ip;

		     double gas_xsec;
	       	     //if (in->spallationxsec != Fluka)
	               gas_xsec = gas_tot * (halfdtbar*xsec[ip]);
	             //else		
		     //  gas_xsec = gas_tot * (halfdtbar*( xsec_extended[make_pair(uid,1001)][ip]*gas_abundances[1001] +		
		     //			xsec_extended[make_pair(uid,2004)][ip]*gas_abundances[2004])); 
//+xsec_extended[make_pair(uid,6012)][ip]*gas_abundances[6012]) ); 	
                     
                     double devect = 1.0
                     + gas_xsec + halfdtbar*(decay)/(lifetime*gamma[ip]) 
                     + halfdt*(riac2[ind1] + Pdotdown[ind1]);
                     double odeuvect = halfdt*(riac1[ind1] + Pdotup[ind1]);
                     de[ip] = devect;
                     odeu[ip] = -odeuvect;
                     oded[ip] = -halfdt*riac3[ind1];
                     
                     // cout << "[MW-DEBUG] " << ip << " " << i << " " << j << " | " << Pdotup[ind1] << " " << Pdotdown[ind1] << " " << riac1[ind1] << " " << riac2[ind1] << " " << riac3[ind1] << " " << totgas*xsec[ip] << " " << (decay)/(lifetime*gamma[ip]) << endl;
                     
                     Re[ip] = dtbarprimsource*spectrum[ip] + dtbar*SecSource_[ind1] + (2.0-devect)*N[ind1];
                     if (ip > 0) Re[ip] -= oded[ip]*N[ind1-1];
                     if (ip < dimE-1) Re[ip] += odeuvect*N[ind1+1];
                  }
                  
                  //Utility::solve_tridag(&(oded[0]), &(de[0]), &(odeu[0]), &(Re[0]), &(ee[0]), dimE);
		  Utility::solve_tridag(oded, de, odeu, Re, ee, dimE);
                  
                  for (ip = 0; ip < dimE; ip++) {
                     value = ee[ip];
                     N[ind+ip] = (value > 0) ? value : 0.0;
                  }
               }
            }
            
            /********************************************/
            
            // Here propagation in R direction starts.
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for (ip = 0; ip < dimE; ip++) {
               
               //double halfdtbarxseck = halfdtbar*xsec[ip];
	       double halfdtbar_xsec;
	       //if (in->spallationxsec != Fluka)
	          halfdtbar_xsec = halfdtbar*xsec[ip];
	       //else		
	//	  halfdtbar_xsec = halfdtbar*( xsec_extended[make_pair(uid,1001)][ip]*gas_abundances[1001] +		
	//				xsec_extended[make_pair(uid,2004)][ip]*gas_abundances[2004]);
// +	xsec_extended[make_pair(uid,6012)][ip]*gas_abundances[6012]);

               double halfdtbarlifetimegammak = halfdtbar*(decay)/(lifetime*gamma[ip]);
               double dtbarinjfactorspeck = dtbar*injfactor*spectrum[ip];
               double sp = dperp->GetSpectrum(ip);
               
               for (j = 1; j < dimz-1; j++) {
                  for (i = 0; i < dimr; i++) {
                     
                     //MW130705: consistency to 3D case
                     indspat = coord->index(i,j);
                     ind = indspat*dimE + ip;
                     
                     double CNalphar1 = dperp->GetCNdiff_alpha1_r(ind);
                     double CNalphar2 = dperp->GetCNdiff_alpha2_r(ind);
                     double CNalphar3 = dperp->GetCNdiff_alpha3_r(ind);
		     //double halfdtdperpfactorphi = halfdt_dperp_factor[i][ip]*dperp->GetPhi(indspat)*sp;  // Variable phi not used!

		     double gas_xsec = halfdtbar_xsec * totalgas->GetGas(indspat);	                     
		     drr[i] = 1.0 + CNalphar2*halfdt_dperp_factor[i][ip] + gas_xsec + halfdtbarlifetimegammak;
                  
                     if(i==0){
			uodrr[i] = -CNalphar3*halfdt_dperp_factor[i][ip];
                        uodrr[i] *= 2; //Symmetry condition at R = 0
                     }
		     else uodrr[i] = -CNalphar3*halfdt_dperp_factor[i][ip] // -halfdtdperpfactorphi;   // Variable phi not used!
                     
                     
		     lodrr[i] = -CNalphar1*halfdt_dperp_factor[i][ip]; // +halfdtdperpfactorphi;  // Variable phi not used!
                  
                     Rrr[i] = N[ind] * (2.0-drr[i]) + source->GetSource(indspat)*dtbarinjfactorspeck + dtbar*SecSource_[ind];
                     if (i < dimr-1) Rrr[i] -= N[ind+dimE*dimz]*uodrr[i];
                     if (i > 0)      Rrr[i] -= N[ind-dimE*dimz]*lodrr[i];
                     
                     // 				        cout << "[MW-DEBUG] " << ip << " " << i << " " << j << " | " << CNalphar1 << " " << CNalphar2 << " " << CNalphar3 << " " << drr[i] << " " << lodrr[i] << " " << uodrr[i] << " | " << Rrr[i] << endl;
                     
                  }
                  
                  // Calculate N
                  //Utility::solve_tridag(&(lodrr[0]), &(drr[0]), &(uodrr[0]), &(Rrr[0]), &(xx[0]), dimr);
		  Utility::solve_tridag(lodrr, drr, uodrr, Rrr, xx, dimr);
                  for (int i = dimr-2; i >= 0; --i) {
                     value = xx[i];
                     N[index(i,j,ip)] = (value > 0) ? value : 0.0;
                  }
               }
            }
         } // #pragma omp parallel
      } // for Niter
      
      dt *= in->dtfactor;
      
   } // while (dt > dtmin)
   
   return ;
}

TCREvolutorADI::TCREvolutorADI(Galaxy* gal1) {
   
   gal = gal1;
   in = gal1->GetInput();
   
   tolerance = in->tol;
   
   coord = gal->GetCoordinates();
   dimr = coord->GetDimR();
   dimz = coord->GetDimZ();
   dimE = coord->GetDimE();
   
   N1       = vector<double>(dimr*dimz*dimE,0);
   N2       = vector<double>(dimr*dimz*dimE,0);
   riac1    = vector<double>(dimr*dimz*dimE,0);
   riac2    = vector<double>(dimr*dimz*dimE,0);
   riac3    = vector<double>(dimr*dimz*dimE,0);
   Pdotup   = vector<double>(dimr*dimz*dimE,0);
   Pdotdown = vector<double>(dimr*dimz*dimE,0);
   
}

double TCREvolutorADI::FindMax(TDiffusionCoefficient* dperp) {
   
   vector<double> spatial = dperp->GetDiffusionCoefficient();
   vector<double> spectrum = dperp->GetSpectrum();
   vector<double>::iterator it = max_element(spatial.begin(), spatial.end());
   vector<double>::iterator itsp = max_element(spectrum.begin(), spectrum.end());
   return (*itsp)*(*it);
}

double TCREvolutorADI::FindTRiacc(TReaccelerationCoefficient* dpp, vector<double>& totmomentum) {
   double result = 1e20;
   double aux = 0.0;
   
   for (int i = 1; i < dimE; i++) {
      aux = pow(totmomentum[i]-totmomentum[i-1],2.0)/dpp->GetSpectrum(i);
      result = min(aux, result);
   }
   
   vector<double> spat = dpp->GetReaccelerationCoefficient();
   vector<double>::iterator it = max_element(spat.begin(), spat.end());
   
   return result/(*it);
}

double TCREvolutorADI::FindTLoss(const vector<double>& dpdt, vector<double>& totmomentum) {
   
   double result = 1e20;
   double aux = 0.0;
   for (int k = 0; k < dimE; k++) {
      for (int i = 0; i < dimr; i++) {
         for (int j = 0; j < dimz; j++) {
            aux = totmomentum[k]/dpdt[index(i,j,k)];
            result = min(aux, result);
         }
      }
   }
   
   return result;
}

double TCREvolutorADI::FindTInt(TGas* gas, const vector<double>& xsec) {
   
   vector<double> gasvec = gas->GetGas();
   vector<double>::iterator itgas = max_element(gasvec.begin(), gasvec.end());
   vector<double>::const_iterator itxsec = max_element(xsec.begin(), xsec.end());
   return 1.0/((*itgas)*(*itxsec));
}

//modified
void TCREvolutorADI::Run(vector<double>& N, vector<double>& N_previous,  TInelasticCrossSection *xsec_object, const vector<double>& dpdt, const vector<double>& SecSource_, const vector<double>& spectrum, double A, double Z, double lifetime, int daughter, bool SecEl, int K_electron, bool isDM, bool isextra) {
   
   // //MW130706: Disable the whole ADI stuff because we never use it and it relies on Phi,Psi which I got rid of.
   //
   //     vector<double>::const_iterator itsecsource = max_element(SecSource_.begin(), SecSource_.end());
   //     if (gal->GetSourceAbundance(int(A+1000*Z)) == 0 && (*itsecsource) == 0) return ;
   //
   //     TDiffusionCoefficient* dperp = NULL;
   //     TReaccelerationCoefficient* dpp = NULL;
   //
   //     if (A == 0) {
   //         dperp = gal->GetDiffCoeffEl();
   //         dpp = gal->GetReaccCoeffEl();
   //     }
   //     else {
   //         dperp = gal->GetDiffCoeff();
   //         dpp = gal->GetReaccCoeff();
   //     }
   //
   //     // Following quantities depend on rigidity, but are computed independently of the nucleus, to save memory and time. Here the necessary additional nucleus-dependent factors are restored.
   //
   //     const double dperpfactor = (A==0) ? 1.0 : pow(A/fabs(Z),dperp->GetDelta());
   //     const double dppfactor = 1.0/dperpfactor; // = A*A/dperpfactor/A/A
   //
   //     /*
   //      #ifdef HAVE_DS
   //      const double injfactor = (SecEl) ? 0 : 1;
   //      #else
   //      const double injfactor = (SecEl) ? 0 : (A==0) ? gal->GetSourceAbundance(int(A+1000*Z)) : gal->GetSourceAbundance(int(A+1000*Z))/A;
   //      #endif
   //      */
   //     double injfactor = CalcInjFactor(K_electron, SecEl, isDM, A, Z);
   //     /*
   //      double injfactor = 0.0;
   //      if (SecEl) injfactor = 0.0;
   //      else if (A > 0) injfactor = gal->GetSourceAbundance(int(A+1000*Z))/A;
   //      else injfactor = 1.0;
   //      */
   //
   //     TGas* totalgas = gal->GetTotalGas();
   //     TSource* source = GetSourceTerm(isDM, isextra);
   //
   //     //    TSource* source = (isDM) ? gal->GetDMSource() : gal->GetSource();
   //     TConvectionVelocity* vC = gal->GetVC();
   //
   //     // TEST!!!
   //     cout << A << " " << Z << " " << SecEl << " " << injfactor << endl;
   //     //if (K_electron > 0) injfactor = 0.; // the nucleus that already attached an electron only has a secondary source term!!
   //
   //     vector<double> totmomentum;
   //     if (A == 0) totmomentum = coord->GetMomentumEl();
   //     else totmomentum = coord->GetMomentum();
   //
   //     vector<double> gamma(coord->GetGamma());
   //
   //     /* Some initialization */
   //
   //     int ind = -1;
   //     int indspat = -1;
   //
   //     for (int i = 0; i < dimr; i++) {
   //         for (int j = 0; j < dimz; j++) {
   //             indspat = coord->index(i,j);
   //             for (int k = 1; k < dimE-1; k++) {
   //                 double momentumup = totmomentum[k+1];
   //                 double momentumfix = totmomentum[k];
   //                 double momentumdown = totmomentum[k-1];
   //                 double upfix = momentumup-momentumfix;
   //                 double updown = momentumup-momentumdown;
   //                 double fixdown = momentumfix-momentumdown;
   //
   //                 ind = indspat*dimE+k;
   //                 Pdotup[ind] = dpdt[ind+1]/(upfix);
   //                 Pdotdown[ind] = dpdt[ind]/(upfix);
   //                 if (A != 0) {
   //                     Pdotup[ind] /= A;
   //                     Pdotdown[ind] /= A;
   //                 }
   //                 if (vC) {
   //                     Pdotup[ind] += momentumup/3.0*vC->GetDvdz()/(upfix);
   //                     Pdotdown[ind] += momentumfix/3.0*vC->GetDvdz()/(upfix);
   //                 }
   //
   //                 if (dpp) {
   //                     double dppfix = dppfactor*dpp->GetReaccelerationCoefficient(indspat, k);
   //                     double dppdown = dppfactor*dpp->GetReaccelerationCoefficient(indspat, k-1);
   //                     riac1[ind] = 2.0*dppfix/(updown)/(upfix);
   //                     riac2[ind] = ( - (dppfix-dppdown)/pow(fixdown, 2.0) +
   //                                   2.0*dppfix/(updown)*(1.0/(upfix) + 1.0/(fixdown)) +
   //                                   2.0*dppfix/momentumfix/(fixdown)
   //                                   );
   //                     riac3[ind] = ( - (dppfix-dppdown)/pow(fixdown, 2.0) +
   //                                   2.0*dppfix/(updown)/(fixdown) +
   //                                   2.0*dppdown/momentumdown/(fixdown)
   //                                   );
   //                 }
   //             }
   //             ind = indspat*dimE;
   //             Pdotup[ind] = dpdt[indspat+1]/(totmomentum[1]-totmomentum[0]);
   //             Pdotdown[ind] = dpdt[ind]/(totmomentum[1]-totmomentum[0]);
   //             if (A != 0) {
   //                 Pdotup[ind] /= A;
   //                 Pdotdown[ind] /= A;
   //             }
   //
   //             if (vC) {
   //                 Pdotup[ind] += totmomentum[1]/3.0*vC->GetDvdz()/(totmomentum[1]-totmomentum[0]);
   //                 Pdotdown[ind] += totmomentum[0]/3.0*vC->GetDvdz()/(totmomentum[1]-totmomentum[0]);
   //             }
   //
   //             // Boundary conditions
   //             Pdotup[ind+dimE-1] = Pdotup[ind+dimE-2];
   //             Pdotdown[ind+dimE-1] = Pdotdown[ind+dimE-2];
   //             if (dpp) {
   //                 riac1[ind] = riac1[ind+1];
   //                 riac2[ind] = riac2[ind+1];
   //                 riac3[ind] = riac3[ind+1];
   //                 riac1[ind+dimE-1] = riac1[ind+dimE-2];
   //                 riac2[ind+dimE-1] = riac2[ind+dimE-2];
   //                 riac3[ind+dimE-1] = riac3[ind+dimE-2];
   //             }
   //         }
   //     }
   //
   //     double decay = (daughter!=0);
   //     //  const double maxDperp = dperpfactor*FindMax(dperp);
   //     const double minTRiacc = (dpp) ? 1.0/dppfactor*FindTRiacc(dpp, totmomentum) : 1e20;
   //     const double minTLoss = (A == 0) ? FindTLoss(dpdt, totmomentum) : A*FindTLoss(dpdt, totmomentum);
   //     const double minTInt = (A == 0) ? 1e20 : FindTInt(totalgas, xsec);
   //     const double minTDec = (daughter != 0) ? lifetime : 1e20;
   //
   //     const double ADIfactor = 1.0/p;
   //
   //     bool conv = false;
   //
   //     double vCi = 0.0; // vC(i)
   //     double vCi1 = 0.0; // vC(i+1)
   //     double vC1i = 0.0; // vC(i-1)
   //     double sp = 0.0;
   //
   //     int counter = 0;
   //
   //     vector<double> Rzz(dimz-2, 0.0);
   //     vector<double> dzz(dimz-2, 0.0);
   //     vector<double> uodzz(dimz-2, 0.0);
   //     vector<double> lodzz(dimz-2, 0.0);
   //     vector<double> yy(dimz-2, 0.0);
   //
   //
   //     vector<double> Rrr(dimr-1, 0.0);
   //     vector<double> drr(dimr-1, 0.0);
   //     vector<double> uodrr(dimr-1, 0.0);
   //     vector<double> lodrr(dimr-1, 0.0);
   //     vector<double> xx(dimr-1, 0.0);
   //
   //     vector<double> de(dimE-2, 0.0);
   //     vector<double> ee(dimE-2, 0.0);
   //     vector<double> odeu(dimE-2, 0.0);
   //     vector<double> oded(dimE-2, 0.0);
   //     vector<double> Re(dimE-2, 0.0);
   //
   //
   //     while (!conv) {
   //         double error = 0.0;
   //
   //         /*************************************************/
   //         // Here propagation in Z direction starts.
   //         for (int k = 1; k < dimE-1; k++) {
   //
   //
   //
   //             for (int j = 0; j < dimr-1; j++) {
   //             double dr = coord->GetDeltaR_central(j);
   //             double coeffz = 1/dr/dr;
   //                 for (int i = 1; i < dimz-1; i++) {
   //                 double dz = coord->GetDeltaZ_central(i);
   //                 double coeffr = 1/dz/dz;
   //
   // //MW130620: dtADI --> in->alpha and multiply dr*dz on it later on. - have to move this inside of the loop
   //
   //                     double dtbar = ADIfactor * in->alpha * dr * dz;//min(dtADI, min(minTRiacc, min(minTLoss, min(minTInt,minTDec) ) ) );
   //                     double halfdtbar = 0.5*dtbar;
   //
   //                     double halfdtbarxseck = halfdtbar*xsec[k];
   //                     double halfdtbarlifetimegammak = halfdtbar*(decay)/(lifetime*gamma[k]);
   //                     double dtbarinjfactorspeck = dtbar*injfactor*spectrum[k];
   //                     sp = dperpfactor*dperp->GetSpectrum(k);
   //
   //
   //                     indspat = coord->index(j,i);
   //                     ind = indspat*dimE+k;
   //                     double Nind = N[ind];
   //                     double psi = dperp->GetPsi(indspat);
   //
   //                     double diffcoeff = dperp->GetDiffusionCoefficient(indspat)*sp;
   //                     double coeffrdiffcoeff = coeffr*diffcoeff;
   //
   //                     if (vC) {
   //                         if (i > (dimz-1)/2) {
   //                             vCi1 = 0.0;
   //                             vC1i = vC->GetVC(j,i-2)/dz;
   //                             vCi = vC->GetVC(j,i)/dz;
   //                         }
   //                         else if (i < (dimz-1)/2) {
   //                             vC1i = 0.0;
   //                             vCi1 = -vC->GetVC(j,i+1)/dz;
   //                             vCi = -vC->GetVC(j,i)/dz;
   //                         }
   //                         else {
   //                             vCi = 0.0;
   //                             vCi1 = 0.0;
   //                             vC1i = 0.0;
   //                         }
   //                     }
   //
   //                     double dzzcoeff = dtbar*coeffrdiffcoeff + (totalgas->GetGas(indspat)*halfdtbarxseck + halfdtbarlifetimegammak) + halfdtbar*vCi;
   //                     dzz[i-1] = 1.0 + dzzcoeff;
   //
   //                     if (i < dimz-2) uodzz[i-1] = -halfdtbar*(coeffrdiffcoeff + psi*sp) - halfdtbar*vCi1;
   //                     if (i > 1) lodzz[i-1] = -halfdtbar*(coeffrdiffcoeff - dperp->GetPsi(indspat-1)*sp) - halfdtbar*vC1i;
   //
   //                     Rzz[i-1] = Nind*(1.0 - dzzcoeff)
   //                     + N[ind+dimE]*(halfdtbar*(coeffrdiffcoeff + psi*sp) + vCi1)
   //                     + N[ind-dimE]*(halfdtbar*(coeffrdiffcoeff - psi*sp) + vC1i)
   //                     + source->GetSource(indspat)*dtbarinjfactorspeck
   //                     + dtbar*SecSource_[ind]
   //                     ;
   //
   //                     if (j == 0) Rzz[i-1] += 2.0*dtbar*coeffz*diffcoeff*(N[ind+dimz*dimE] - Nind);
   //                     else Rzz[i-1] += dtbar*(
   //                                             (coeffz*diffcoeff+dperp->GetPhi(indspat)*sp)*N[ind+dimz*dimE]
   //                                             - 2.0*coeffz*diffcoeff*Nind
   //                                             + (coeffz*diffcoeff-dperp->GetPhi(indspat)*sp)*N[ind-dimz*dimE]);
   //
   //                     Rzz[i-1] += dtbar*(
   //                                        (riac1[ind]+Pdotup[ind])*N[ind+1]
   //                                        - (riac2[ind]+Pdotdown[ind])*Nind
   //                                        + riac3[ind]*N[ind-1]);
   //                 } // dimz-1
   //
   //                 // Calcolare N
   //                 //Utility::solve_tridag(lodzz, dzz, uodzz, Rzz, yy);
   //                 Utility::solve_tridag(&(lodzz[0]), &(dzz[0]), &(uodzz[0]), &(Rzz[0]), &(yy[0]), dimz);
   //                 for (int i = dimz-2; i > 0; --i) N1[index(j,i,k)] = yy[i-1];
   //             } // dimr-1
   //         } // dimE
   //
   //         /***********************************************/
   //         /*********************************************/
   //
   //         // Here propagation in Momentum direction starts: Reacceleration and/or energy losses
   //         for (int i = 0; i < dimr-1; i++) {
   //         double dr = coord->GetDeltaR_central(i);
   //         double coeffz = 1/dr/dr;
   //             for (int j = 1; j < dimz-1; j++) {
   //             double dz = coord->GetDeltaZ_central(j);
   //             double coeffr = 1/dz/dz;
   //
   // //MW130620: dtADI --> in->alpha and multiply dr*dz on it later on. - have to move this inside of the loop
   //                 double dtbar = ADIfactor * in->alpha * dr * dz;//min(dtADI, min(minTRiacc, min(minTLoss, min(minTInt,minTDec) ) ) );
   //                 double halfdtbar = 0.5*dtbar;
   //
   //                 indspat = coord->index(i,j);
   //                 double totgas = totalgas->GetGas(indspat);
   //                 ind = indspat*dimE;
   //
   //                 for (int k = 1; k < dimE-1; k++) {
   //                     int ind1 = ind + k;
   //                     double N1ind1 = N1[ind1];
   //                     sp = dperpfactor*dperp->GetSpectrum(k);
   //
   //                     double dtbarprimsource = dtbar*injfactor*source->GetSource(indspat);
   //                     double defact = halfdtbar*(riac2[ind1] + Pdotdown[ind1] + totgas*xsec[k] + (decay)/(lifetime*gamma[k]));
   //
   //                     de[k-1] = 1.0 + defact;
   //
   //                     if (k < dimE-2) odeu[k-1] = - halfdtbar*(riac1[ind1] + Pdotup[ind1]);
   //                     if (k > 1) oded[k-1] = - halfdtbar*riac3[ind1-1];
   //
   //
   //                     Re[k-1] = N1ind1*(1.0 - defact)
   //                     + N1[ind1+1]*halfdtbar*(riac1[ind1] + Pdotup[ind1])
   //                     + N1[ind1-1]*halfdtbar*riac3[ind1]
   //                     + dtbarprimsource*spectrum[k] // primary source
   //                     + dtbar*SecSource_[ind1] // secondary source
   //                     ;
   //
   //                     double diffcoeff = dperp->GetDiffusionCoefficient(indspat)*sp;
   //
   //                     if (i == 0) Re[k-1] += dtbar*coeffz*diffcoeff*(2.0*N1[ind1+dimz*dimE] - 2.0*N1ind1);
   //                     else Re[k-1] += dtbar*(
   //                                            (coeffz*diffcoeff + dperp->GetPhi(indspat)*sp)*N1[ind1+dimz*dimE]
   //                                            - 2.0*coeffz*diffcoeff*N1ind1
   //                                            + (coeffz*diffcoeff - dperp->GetPhi(indspat)*sp)*N1[ind1-dimz*dimE]
   //                                            );
   //
   //                     Re[k-1] += dtbar*(
   //                                       (coeffr*diffcoeff + dperp->GetPsi(indspat)*sp)*N1[ind1+dimE]
   //                                       - 2.0*coeffr*diffcoeff*N1ind1
   //                                       + (coeffr*diffcoeff - dperp->GetPsi(indspat)*sp)*N1[ind1-dimE]
   //                                       );
   //                     if (vC) {
   //                         if (j > (dimz-1)/2) Re[k-1] -= dtbar/dz*(vC->GetVC(i,j)*N1ind1 - vC->GetVC(i,j-1)*N1[ind1-dimE]);
   //                         else if (j < (dimz-1)/2) Re[k-1] += dtbar/dz*(vC->GetVC(i,j)*N1ind1 - vC->GetVC(i,j+1)*N1[ind1+dimE]);
   //                     }
   //                 } // dimE-1
   //
   //                 //Utility::solve_tridag(oded, de, odeu, Re, ee);
   //                 Utility::solve_tridag(&(oded[0]), &(de[0]), &(odeu[0]), &(Re[0]), &(ee[0]), dimE);
   //                 for (int k = 1; k < dimE-1; k++) N2[index(i,j,k)] = ee[k-1];
   //
   //             } // dimz-1
   //         } // dimr-1
   //
   //
   //         /********************************************/
   //
   //         // Here propagation in R direction starts.
   //         for (int k = 1; k < dimE-1; k++) {
   //
   //             for (int j = 1; j < dimz-1; j++) {
   //             double dz = coord->GetDeltaZ_central(j);
   //             double coeffr = 1/dz/dz;
   //                 for (int i = 0; i < dimr-1; i++) {
   //                 double dr = coord->GetDeltaR_central(i);
   //                 double coeffz = 1/dr/dr;
   //
   // //MW130620: have to move this (has been dtADI) inside of the loop
   //                 double dtbar = ADIfactor * in->alpha * dr * dz;//min(dtADI, min(minTRiacc, min(minTLoss, min(minTInt,minTDec) ) ) );
   //                 double halfdtbar = 0.5*dtbar;
   //                 double halfdtbarxseck = halfdtbar*xsec[k];
   //                 double halfdtbarlifetimegammak = halfdtbar*(decay)/(lifetime*gamma[k]);
   //                 double dtbarinjfactorspeck = dtbar*injfactor*spectrum[k];
   //                 sp = dperpfactor*dperp->GetSpectrum(k);
   //
   //                     indspat = coord->index(i,j);
   //                     ind = indspat*dimE+k;
   //
   //                     double N2ind = N2[ind];
   //
   //                     double diffcoeff = dperp->GetDiffusionCoefficient(indspat)*sp;
   //                     double coeffzdiffcoeff = coeffz*diffcoeff;
   //
   //                     double drrcoeff = dtbar*coeffzdiffcoeff + (totalgas->GetGas(indspat)*halfdtbarxseck + halfdtbarlifetimegammak);
   //                     drr[i] = 1.0 + drrcoeff;
   //
   //                     double phi = dperp->GetPhi(indspat);
   //
   //                     if (i == 0) uodrr[i] = -dtbar*coeffzdiffcoeff; // i == 0, Symmetry condition at R = 0
   //                     else uodrr[i] = -halfdtbar*(coeffzdiffcoeff + phi*sp);
   //                     if (i > 0) lodrr[i] = -halfdtbar*(coeffzdiffcoeff - phi*sp);
   //
   //                     Rrr[i] = N2ind*(1.0 - drrcoeff)
   //                     + N2[ind+dimz*dimE]*((i==0)*dtbar*coeffzdiffcoeff + (i>0)*halfdtbar*(coeffzdiffcoeff + phi*sp))
   //                     + source->GetSource(indspat)*dtbarinjfactorspeck
   //                     + dtbar*SecSource_[ind]
   //                     ;
   //
   //                     if (i > 0) Rrr[i] += N2[ind-dimz*dimE]*halfdtbar*(coeffzdiffcoeff - phi*sp);
   //                     Rrr[i] += dtbar*(
   //                                      (coeffr*diffcoeff + dperp->GetPsi(indspat)*sp)*N2[ind+dimE]
   //                                      - 2.0*coeffr*diffcoeff*N2ind
   //                                      + (coeffr*diffcoeff - dperp->GetPsi(indspat)*sp)*N2[ind-dimE]
   //                                      );
   //
   //                     Rrr[i] += dtbar*(
   //                                      (riac1[ind] + Pdotup[ind])*N2[ind+1]
   //                                      - (riac2[ind] + Pdotdown[ind])*N2ind
   //                                      + riac3[ind]*N2[ind-1]
   //                                      );
   //
   //                     if (vC) {
   //                         if (j > (dimz-1)/2) Rrr[i] -= dtbar/dz*(vC->GetVC(i,j)*N2ind - vC->GetVC(i,j-1)*N2[ind-dimE]);
   //                         else if (j < (dimz-1)/2) Rrr[i] += dtbar/dz*(vC->GetVC(i,j)*N2ind - vC->GetVC(i,j+1)*N2[ind+dimE]);
   //                     }
   //                 } // dimr-1
   //
   //                 // Calculate N
   //                 //                Utility::solve_tridag(lodrr, drr, uodrr, Rrr, xx);
   //                 Utility::solve_tridag(&(lodrr[0]), &(drr[0]), &(uodrr[0]), &(Rrr[0]), &(xx[0]), dimr);
   //
   //                 for (int i = dimr-2; i >= 0; --i) {
   //                     error = max(fabs(1.0 - xx[i]/N[index(i,j,k)]), error);
   //                     N[index(i,j,k)] = xx[i];
   //                 }
   //
   //            } // dimz-1
   //        } // dimE-1
   //
   //        if (counter%100 == 0) cout << error << endl;
   //        counter++;
   //        if (error < tolerance) conv = true;
   //    } // while (!conv)
   //
   //    cout << "Steps: " << counter << endl;
   
   return ;
}

// MODIFIED 17-09-2012

//**********************************************************************************************************************************************************
//**********************************************************************************************************************************************************
//**********************************************************************************************************************************************************

TCREvolutor3D::TCREvolutor3D(Galaxy* gal1)  {
   
   gal = gal1;
   in = gal1->GetInput();
   
   coord = gal->GetCoordinates();
   
   dimx = coord->GetDimX();
   dimy = coord->GetDimY();
   dimz = coord->GetDimZ();
   dimE = coord->GetDimE();
   
   riac1    = vector<double>(dimx*dimy*dimz*dimE,0);
   riac2    = vector<double>(dimx*dimy*dimz*dimE,0);
   riac3    = vector<double>(dimx*dimy*dimz*dimE,0);
   Pdotup   = vector<double>(dimx*dimy*dimz*dimE,0);
   Pdotdown = vector<double>(dimx*dimy*dimz*dimE,0);
   
   
}

//**********************************************************************************************************************************************************
//**********************************************************************************************************************************************************
//**********************************************************************************************************************************************************

// implemented by Daniele Gaggero (daniele.gaggero@sissa.it)   -- september 2012

// last modified by Daniele Gaggero (daniele.gaggero@sissa.it) -- november 2012 -- tested point source propagation against analytical solution -- tested against 2D code in realistic setup

// anisotropic diffusion implemented by Luca Maccione and Daniele Gaggero -- november 2012 -- !!!under testing (dec 2012)!!!

// moving sources added by D.Gaggero -- 17/01/2013

// changed structure to improve performance - M. Weinreuter 05/07/2013

void TCREvolutor3D::Run(vector<double>& N, vector<double>& N_previous, TInelasticCrossSection *xsec_object, const vector<double>& dpdt, const vector<double>& SecSource_, const vector<double>& spectrum, double A, double Z, double lifetime, int daughter, bool SecEl, int K_electron, bool isDM, bool isExtra) {


   // This is for the horizon(E) plot -- just remove that after the plot is done
   if (A == 0 && Z == -1) {
      TDiffusionCoefficient* dperp = gal->GetDiffCoeffEl();
      vector<double> totmomentum = coord->GetMomentumEl();
      vector<double> dperp_vec = dperp->GetSpectrum();
      vector<double> x = (gal->GetCoordinates()->GetX());
      vector<double> y = (gal->GetCoordinates()->GetY());
      vector<double> z = (gal->GetCoordinates()->GetZ());
      double deltax = ( (x.back() - x.front()) / (dimx-1) );
      double deltay = ( (y.back() - y.front()) / (dimy-1) );
      double deltaz = ( (z.back() - z.front()) / (dimz-1) );
      unsigned int ixsun = (unsigned int) ((in->xobs-x.front())/(2.0*x.back())*(double)(dimx-1));
      unsigned int iysun = (unsigned int) ((in->yobs-y.front())/(2.0*y.back())*(double)(dimy-1));
      unsigned int izsun = (unsigned int) ((in->zobs-z.front())/(2.0*z.back())*(double)(dimz-1));
      long int indspat = coord->indexD(ixsun,iysun,izsun);
      cout << ixsun << " " << iysun << " " << izsun << endl;
      cout << "momentum vector in GeV *** " << endl;	
      for (int ip = 0; ip < totmomentum.size(); ip++)
	cout << totmomentum[ip] << ", " ;
      cout << endl;
      cout << "eloss vector at Sun position in GeV/Myr *** " << endl;	
      for (int ip = 0; ip < totmomentum.size(); ip++) {
	long int ind = indspat*dimE+ip;
	cout << dpdt[ip] << ", " ;
      }
      cout << endl;
      cout << "Diff coefficient at Sun position in kpc^2/Myr *** " << endl;	
      for (int ip = 0; ip < totmomentum.size(); ip++) {
	long int ind = indspat*dimE+ip;
	cout << dperp_vec[ip] << ", " ;
      }
      cout << endl;
      cout << "Horizon at Sun position in kpc *** " << endl;	
      for (int ip = 0; ip < totmomentum.size(); ip++) {
	long int ind = indspat*dimE+ip;
	double horizon = sqrt(dperp_vec[ip] * totmomentum[ip]/dpdt[ip]);	
	cout << horizon << ", " ;
      }
      cout << endl;
   }
 
   // initialization stuff
  
  if (in->feedback >0) cout << "We are starting the 3D ISOTROPIC engine. Diffusion coefficient is a scalar, propagation is computed in 3+1D grid (x,y,z,p)." <<endl;

   cout << "injection spectrum" << endl;
   for (int i=0; i<spectrum.size(); i++)
	cout << spectrum[i] << " " ;
   cout << endl;

   vector<double> xsec;
   map<pair<int,int>, vector<double> > xsec_extended;

   // if (in->spallationxsec == Fluka) {
   // 	xsec_extended = xsec_object->GetXSec_extended();
   // 	cout << "Fluka cross sections don't work in 3D mode! " << endl;
   // 	return;
   // }	
   // else
        xsec = xsec_object->GetXSec();


   cout << "inelastic xsec loaded" << endl;

   //cout << "eloss term" << endl;
   //for (int i=0; i<dpdt.size(); i++)
	//if (i%1000==0) cout << dpdt[i] << " " ;
   //cout << endl;
   
   double injfactor = CalcInjFactor(K_electron, SecEl, isDM, A, Z);
   
   if (in->feedback >0) cout << endl << "*** We are propagating this particle: *** " << endl;
   if (in->feedback >0) cout << "A = " << A << " Z = " << Z << " SecEl = " << SecEl << " K_electron = " << K_electron << " injfactor = " << injfactor << endl;
   if (A==0 && Z==-1 && SecEl==0 && in->feedback >0)
      cout << "Primary electrons" << endl;
   if (A==0 && Z==-1 && SecEl==1 && in->feedback >0) 
      cout << "Secondary electrons" << endl;
   if (A==0 && Z==-1 && SecEl==0 && isExtra && in->feedback >0) 
       cout << "Extra component leptons" << endl;
   if (A==0 && isDM)
       cout << "DM leptons" << endl;
   if (A==1 && isDM)
       cout << "DM antiprotons" << endl;
   if (in->feedback >0) cout << "Starting propagation... " << endl;
   
   vector<double>::const_iterator itsecsource = max_element(SecSource_.begin(), SecSource_.end());
   if (injfactor == 0 && (*itsecsource) == 0 && in->feedback >0) {
      cout << "No primary nor secondary source. Returning." << endl;
      return ;
   }
   
   TGas* totalgas = gal->GetTotalGas();
   TSource* source = GetSourceTerm(isDM, isExtra);
   TConvectionVelocity* vC = gal->GetVC(); //careful -- is it 3D-ready??
   const double dvcdz = (vC==NULL) ? 0.0 : vC->GetDvdz()/3.0;
   
   TDiffusionCoefficient* dperp = NULL;
   
   TReaccelerationCoefficient* dpp = NULL;
   
   vector<double> totmomentum;
   
   if (A == 0) {
      dperp = gal->GetDiffCoeffEl();
      dpp = gal->GetReaccCoeffEl();
      totmomentum = coord->GetMomentumEl();
   }
   else {
      dperp = gal->GetDiffCoeff();
      dpp = gal->GetReaccCoeff();
      totmomentum = coord->GetMomentum();
   }
   
   // Following quantities depend on rigidity, but are computed independently of the nucleus, to save memory and time. Here the necessary additional nucleus-dependent factors are restored.
   
	int bin=0;
	for (int ip = 0; ip < dimE; ip++) if(totmomentum[ip]<dperp->Getrho_b()) bin=ip;
   
   double dperpfactor[dimE];
   double dppfactor[dimE];
   
	for (int ip = 0; ip < dimE; ip++){
		if(ip<=bin) dperpfactor[ip] = (A==0) ? 1.0 : pow(A/fabs(Z),dperp->GetDelta());
		else dperpfactor[ip] = (A==0) ? 1.0 : pow(A/fabs(Z),dperp->GetDelta_h());
		dppfactor[ip]=1.0/dperpfactor[ip];
	}
   
   vector<double> gamma(coord->GetGamma());
   
   // main loops
   
   int ind = -1;
   int indspat = -1;
   int indspat_up=-1;//for nonequidistant grid, Iris 02/25/13
   int indspat_down=-1;//for nonequidistant grid, Iris 02/25/13
   double ds_up=0;
   double ds_down=0;
   double ds_double=0;
   int i =  0;
   int j =  0;
   int k =  0;
   int ip = 0;
   
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,ip,indspat,ind) schedule(dynamic) num_threads(NUMTHREADS)
#endif
   for (i = 0; i < dimx; i++) {
      for (j = 0; j < dimy; j++) {
         for (k = 0; k < dimz; k++) {
            indspat = coord->indexD(i,j,k);
            double dpp1 = 0.0;
            
            for (ip = 1; ip < dimE-1; ip++) {
               if (dpp) dpp1 = dppfactor[ip]*dpp->GetReaccelerationCoefficient(indspat);
               
               double momentumup = totmomentum[ip+1];
               double momentumfix = totmomentum[ip];
               double momentumdown = totmomentum[ip-1];
               double upfix = momentumup-momentumfix;
               double updown = momentumup-momentumdown;
               double fixdown = momentumfix-momentumdown;
               
               ind = indspat*dimE+ip;
               Pdotup[ind] = dpdt[ind+1]/(upfix);
               Pdotdown[ind] = dpdt[ind]/(upfix);
               if (A != 0) {
                  Pdotup[ind] /= A;
                  Pdotdown[ind] /= A;
               }
               if (vC) {
                  Pdotup[ind] += momentumup*dvcdz/(upfix);
                  Pdotdown[ind] += momentumfix*dvcdz/(upfix);
               }
               
               if (dpp) {
                  double dppfix = dpp1*dpp->GetSpectrum(ip);
                  double dppdown = dpp1*dpp->GetSpectrum(ip-1);
                  riac1[ind] = 2.0*dppfix/(updown)/(upfix);
                  riac2[ind] = ( - (dppfix-dppdown)/(fixdown*fixdown) +
                                2.0*dppfix*(1.0/(updown)*(1.0/(upfix) + 1.0/(fixdown)) +
                                            1.0/momentumfix/(fixdown))
                                );
                  riac3[ind] = ( - (dppfix-dppdown)/(fixdown) +
                                2.0*(dppfix/(updown) +
                                     dppdown/momentumdown))/(fixdown);
                  
                  // if(A==0) cout << "[MW-DEBUG REACC] " << i << "," << j << "," << k << " " << ip << " dppfactor " << dppfactor[ip] << " " << " reacccoeff " << dpp->GetReaccelerationCoefficient(indspat) << " ... so dpp1 " << dpp1 << " | " << totmomentum[ip] << " " << momentumup << " " << momentumdown << " | " << dpp->GetSpectrum(ip) << endl;
                  
               }
            }
            ind = indspat*dimE;
            Pdotup[ind] = dpdt[ind+1]/(totmomentum[1]-totmomentum[0]);
            Pdotdown[ind] = dpdt[ind]/(totmomentum[1]-totmomentum[0]);
            if (A != 0) {
               Pdotup[ind] /= A;
               Pdotdown[ind] /= A;
            }
            
            if (vC) {
               Pdotup[ind] += totmomentum[1]*dvcdz/(totmomentum[1]-totmomentum[0]);
               Pdotdown[ind] += totmomentum[0]*dvcdz/(totmomentum[1]-totmomentum[0]);
            }
            
            // Boundary conditions
            Pdotup[ind+dimE-1] = Pdotup[ind+dimE-2];
            Pdotdown[ind+dimE-1] = Pdotdown[ind+dimE-2];
            if (dpp) {
               riac1[ind] = riac1[ind+1];
               riac2[ind] = riac2[ind+1];
               riac3[ind] = riac3[ind+1];
               riac1[ind+dimE-1] = riac1[ind+dimE-2];
               riac2[ind+dimE-1] = riac2[ind+dimE-2];
               riac3[ind+dimE-1] = riac3[ind+dimE-2];
            }
            
         }
      }
   }
   
   double dt = in->dtmax;
   
   double dtbar = 0.0;
   
   double halfdtbar = 0.0;
   double halfdt = 0.0;
   
   double halfdt_dperp_factor[dimE];
	for(int e=0;e<dimE;e++) halfdt_dperp_factor[e]=0.;
   
   const double decay = (daughter!=0);
   int Niter = 0;
   
   double value = 0.0;
   
   long counter = -1;

   //double total_injected_energy_so_far = 0.;
   //double total_injected_energy_NOW = 0.;
   //double total_energy_NOW = 0.;
   
   //**********************
   while (dt > in->dtmin) {
      //**********************
      
      counter++;

      cout << counter << " dt =  " << dt << endl;
      
      dtbar = dt/4;  //check p!!!
      
      halfdtbar = 0.5*dtbar;
      halfdt    = 0.5*dt;
      for(int e=0;e<dimE;e++) halfdt_dperp_factor[e] = halfdt*dperpfactor[e];
      
      //*******************************************
      for (Niter = 0; Niter < in->Nrept; ++Niter) {
         //*******************************************

 
#ifdef _OPENMP
#pragma omp parallel default(shared) private(i,j,k,ip,indspat,ind) num_threads(NUMTHREADS)
#endif
         {
            
            //cout << "Niter " << Niter << endl;
            
            vector<double> Rxx(dimx,0.0);
            vector<double> dxx(dimx,0.0);
            vector<double> uodxx(dimx,0.0);
            vector<double> lodxx(dimx,0.0);
            vector<double> xx(dimx,0.0);
            
            vector<double> Ryy(dimy,0.0);
            vector<double> dyy(dimy,0.0);
            vector<double> uodyy(dimy,0.0);
            vector<double> lodyy(dimy,0.0);
            vector<double> yy(dimy,0.0);
            
            vector<double> Rzz(dimz,0.0);
            vector<double> dzz(dimz,0.0);
            vector<double> uodzz(dimz,0.0);
            vector<double> lodzz(dimz,0.0);
            vector<double> zz(dimz,0.0);
            
            vector<double> de(dimE, 0.0);
            vector<double> ee(dimE, 0.0);
            vector<double> odeu(dimE, 0.0);
            vector<double> oded(dimE, 0.0);
            vector<double> Re(dimE, 0.0);
            
            
            // **************************
            // Propagation in x direction
            // **************************
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for  (ip = 0; ip < dimE; ip++) {
               
               double halfdtbar_xsec_ip             = halfdtbar*xsec[ip];
               double halfdtbar_lifetime_gamma_ip   = halfdtbar*(decay)/(lifetime*gamma[ip]);
               double dtbar_injfactor_spec_ip       =     dtbar*injfactor*spectrum[ip];

	       for  (k = 1; k < dimz-1; k++) {
                  for (j = 1; j < dimy-1; j++) {
                     
                     for (i = 0; i < dimx;  i++) {

                        //MW130624
                        indspat = coord->indexD(i,j,k);
                        ind = indspat*dimE+ip;
 
                        double CNalphax1 = dperp->GetCNdiff_alpha1_x(ind);
                        double CNalphax2 = dperp->GetCNdiff_alpha2_x(ind);
                        double CNalphax3 = dperp->GetCNdiff_alpha3_x(ind);
                        
                        dxx[i] = 1. + CNalphax2*halfdt_dperp_factor[ip]+  totalgas->GetGas(indspat)*halfdtbar_xsec_ip  +  halfdtbar_lifetime_gamma_ip ;//CHECK!!! IG
                        uodxx[i] = -CNalphax3*halfdt_dperp_factor[ip];
                        lodxx[i] = -CNalphax1*halfdt_dperp_factor[ip];
                        
                        Rxx[i] = N[index(i,j,k,ip)]*(2.0-dxx[i]) + source->GetSource(indspat)*dtbar_injfactor_spec_ip + SecSource_[ind]*dtbar;
                        
                        if (i < dimx-1) Rxx[i] -= N[index(i+1,j,k,ip)] * uodxx[i] ;
                        if (i > 0)      Rxx[i] -= N[index(i-1,j,k,ip)] * lodxx[i] ;
                     }
                     
                     //Utility::solve_tridag( &(lodxx[0]), &(dxx[0]), &(uodxx[0]), &(Rxx[0]), &(xx[0]), dimx );
		     Utility::solve_tridag(lodxx, dxx, uodxx, Rxx, xx, dimx );
                     
                     for (int i = dimx-2; i >0; --i) {
                        value = xx[i];
                        N[index(i,j,k,ip)] = (value > 0) ? value : 0.0;
                     }
                  }
               }

            } // for ip
            
            // **************************
            // Propagation in y direction
            // **************************
            
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for  (ip = 0; ip < dimE; ip++) {
               
               double halfdtbar_xsec_ip            = halfdtbar*xsec[ip];
               double halfdtbar_lifetime_gamma_ip   = halfdtbar*(decay)/(lifetime*gamma[ip]);
               double dtbar_injfactor_spec_ip       =     dtbar*injfactor*spectrum[ip];
               
               for  (k = 1; k < dimz-1; k++) {
                  for (i = 1; i < dimx-1; i++) {
                     
                     for (j = 0; j < dimy;  j++) {
                        
                        //MW130624
                        indspat = coord->indexD(i,j,k);
                        ind = indspat*dimE+ip;
                        
                        double CNalphay1 = dperp->GetCNdiff_alpha1_y(ind);
                        double CNalphay2 = dperp->GetCNdiff_alpha2_y(ind);
                        double CNalphay3 = dperp->GetCNdiff_alpha3_y(ind);
                        
                        dyy[j] = 1. + CNalphay2*halfdt_dperp_factor[ip]+  totalgas->GetGas(indspat)*halfdtbar_xsec_ip  +  halfdtbar_lifetime_gamma_ip ;//CHECK!!! IG
                        uodyy[j] = -CNalphay3*halfdt_dperp_factor[ip];
                        lodyy[j] = -CNalphay1*halfdt_dperp_factor[ip];
                        
                        Ryy[j] = N[index(i,j,k,ip)]*(2.0-dyy[j]) + source->GetSource(indspat)*dtbar_injfactor_spec_ip + SecSource_[ind]*dtbar;
                        
                        if (j < dimy-1) Ryy[j] -= N[index(i,j+1,k,ip)] * uodyy[j] ;
                        if (j > 0)      Ryy[j] -= N[index(i,j-1,k,ip)] * lodyy[j] ;
                     }
                     
                     //Utility::solve_tridag( &(lodyy[0]), &(dyy[0]), &(uodyy[0]), &(Ryy[0]), &(yy[0]), dimy );
		     Utility::solve_tridag(lodyy, dyy, uodyy, Ryy, yy, dimy);
                     
                     for (int j = dimy-2; j >0; --j) {
                        value = yy[j];
                        N[index(i,j,k,ip)] = (value > 0) ? value : 0.0;
                     }
                  }
               }
               
            }// for ip
            
            // **************************
            // Propagation in z direction
            // **************************
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for  (ip = 0; ip < dimE; ip++) {
               
               double halfdtbar_xsec_ip            = halfdtbar*xsec[ip];
               double halfdtbar_lifetime_gamma_ip   = halfdtbar*(decay)/(lifetime*gamma[ip]);
               double dtbar_injfactor_spec_ip       =     dtbar*injfactor*spectrum[ip];
               
               for (i = 1; i < dimx-1;  i++) {
                  for (j = 1; j < dimy-1; j++) {
                     
                     for  (k = 0; k < dimz; k++) {
                        //MW130620
                        indspat = coord->indexD(i,j,k);
                        ind = indspat*dimE+ip;
                        
                        double CNalphaz1 = dperp->GetCNdiff_alpha1_z(ind);
                        double CNalphaz2 = dperp->GetCNdiff_alpha2_z(ind);
                        double CNalphaz3 = dperp->GetCNdiff_alpha3_z(ind);
                        
                        //moved to galaxy.cc MW130624
                        double vCk = 0.0; // vC(i)
                        double vCk1 = 0.0; // vC(i+1)
                        double vC1k = 0.0; // vC(i-1)
                        
                        if (vC)
                        {
                           vC1k = vC->GetCNconv_alpha1_z(indspat);
                           vCk  = vC->GetCNconv_alpha2_z(indspat);
                           vCk1 = vC->GetCNconv_alpha3_z(indspat);
                        }
                        
                        //                     if(counter==0 && i==28 && j==20) cout << "TEST A " << A << " ip " << ip << " k " << k << " | " << CNalphaz1 << " " << CNalphaz2 << " " << CNalphaz3 << " " << vC1k << " " << vCk << " " << vCk1 << " " << xsec[ip] << " " << decay << " " << lifetime << " " << gamma[ip] << " " << injfactor << " " << spectrum[ip] << endl;
                        
                        dzz[k] = 1. + CNalphaz2*halfdt_dperp_factor[ip]+  totalgas->GetGas(indspat)*halfdtbar_xsec_ip  +  halfdtbar_lifetime_gamma_ip + halfdt*vCk;//CHECK!!! IG
                        uodzz[k] = -CNalphaz3*halfdt_dperp_factor[ip]-halfdt*vCk1;
                        lodzz[k] = -CNalphaz1*halfdt_dperp_factor[ip]-halfdt*vC1k;
                        
                        //  				        cout << "[MW-DEBUG Z] " << i << " " << j << " " << k << " | " << CNalphaz1 << " " << CNalphaz2 << " " << CNalphaz3 << " " << vC1k << " " << vCk << " " << vCk1 << " | " << lodzz[k] << " " << dzz[k] << " "  << uodzz[k] << endl;
                        
                        Rzz[k] = N[index(i,j,k,ip)]*(2.0-dzz[k]) + source->GetSource(indspat)*dtbar_injfactor_spec_ip + SecSource_[ind]*dtbar;
                        
                        if (k < dimz-1) Rzz[k] -= N[index(i,j,k+1,ip)] * uodzz[k] ;
                        if (k > 0)      Rzz[k] -= N[index(i,j,k-1,ip)] * lodzz[k] ;
                     }
                     
                     //Utility::solve_tridag( &(lodzz[0]), &(dzz[0]), &(uodzz[0]), &(Rzz[0]), &(zz[0]), dimz );
		     Utility::solve_tridag(lodzz, dzz, uodzz, Rzz, zz, dimz);
                     
                     for (int k = dimz-2; k >0; --k) {
                        value = zz[k];
                        N[index(i,j,k,ip)] = (value > 0) ? value : 0.0;
                     }
                  }
               }
               
            }
            
            // **************************
            // Propagation in p direction
            // **************************
            
            
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for (i = 1; i < dimx-1;  i++) {
               for (j = 1; j < dimy-1; j++) {
                  for (k = 1; k < dimz-1; k++) {
                     
                     indspat = coord->indexD(i,j,k);
                     double totgas = totalgas->GetGas(indspat);
                     double dtbarprimsource = dtbar*injfactor*source->GetSource(indspat);
                     ind = indspat*dimE;
                     
                     for (ip = 0; ip < dimE; ip++) {
                        int ind1 = ind + ip;
                        
                        double devect = 1.0
                        + halfdtbar*(totgas*xsec[ip] + (decay)/(lifetime*gamma[ip]))
                        + halfdt*(riac2[ind1] + Pdotdown[ind1]);
                        double odeuvect = halfdt*(riac1[ind1] + Pdotup[ind1]);
                        de[ip] = devect;
                        odeu[ip] = -odeuvect;
                        oded[ip] = -halfdt*riac3[ind1];
                        
                        // cout << "[MW-DEBUG P] (" << A << "," << Z << ") " << i << " " << j << " " << k << " " << ip << " | " << Pdotdown[ind1] << " " << Pdotup[ind1] << " " << xsec[ip] << " " << gamma[ip] << " | " << oded[ip] << " " << de[ip] << " "  << odeu[ip] << endl;
                        
                        
                        Re[ip] = dtbarprimsource*spectrum[ip] + dtbar*SecSource_[ind1] + (2.0-devect)*N[ind1];
                        if (ip > 0) Re[ip] -= oded[ip]*N[ind1-1];
                        if (ip < dimE-1) Re[ip] += odeuvect*N[ind1+1];
                     }
                     //Utility::solve_tridag(&(oded[0]), &(de[0]), &(odeu[0]), &(Re[0]), &(ee[0]), dimE);
		     Utility::solve_tridag(oded, de, odeu, Re, ee, dimE);
                     
                     for (ip = 0; ip < dimE; ip++) {
                        value = ee[ip];
                        N[ind+ip] = (value > 0) ? value : 0.0;
                     }
                  }
               } // dimy-1
            } // dimr-1
            
            
         } //#pragma omp parallel
         
         /////////////////////////////////////////////////////////
         //      writing temporary fits files
         //
         //      use like this:
         //         <Nrept value=60/>  -- number of steps
         //         <Dtmax value=0.1/> -- timestep in Myr
         //
         //      Dtfactor and Dtmin don't matter for timestepstore,
         //
         /////////////////////////////////////////////////////////
         
         if (in->timestepstore && counter==0 && ((Z==1 && A==1) || (Z==-1 && A==0)) ) {
            
            string name;
            stringstream namestream;
            
            namestream << "!output/" << in->filename.c_str() << "_temp_" << 1000*Z + A;
            if(SecEl) namestream << "s";
            namestream << "_step" << Niter << ".fits.gz";
            
            name = namestream.str();
            
            fitsfile* output_ptr;
            
            int status = 0;
            
            const long naxis = 4;
            long size_axes[naxis] = {dimx,dimy,dimz,dimE};
            long nelements = size_axes[0]*size_axes[1]*size_axes[2]*size_axes[3];
            long fpixel = 1;
            int bitpix = FLOAT_IMG;
            
            cerr << "Writing temporary output file " << name << ", nelements = " << nelements << endl;
            
            if (fits_create_file(&output_ptr, name.c_str(), &status)) fits_report_error(stderr, status);
            /*            if(status==105) //file could not be created. Probably already exists.
             {
             cerr << " deleting " << name.c_str() << " and try once again...";
             if( remove(name.c_str()) ) cerr << " successfully removed. Now?";
             cerr << endl;
             
             if (fits_create_file(&output_ptr, name.c_str(), &status)) fits_report_error(stderr, status);
             }*/
            
            if (fits_create_img(output_ptr, bitpix, naxis, size_axes, &status)) fits_report_error(stderr, status);
            
            //             if (fits_write_key(output_ptr, TINT, (char*) "dimE", &dimE, NULL, &status)) fits_report_error(stderr, status);
            //             if (fits_write_key(output_ptr, TINT, (char*) "dimx", &dimx, NULL, &status)) fits_report_error(stderr, status);
            //             if (fits_write_key(output_ptr, TINT, (char*) "dimy", &dimy, NULL, &status)) fits_report_error(stderr, status);
            //             if (fits_write_key(output_ptr, TINT, (char*) "dimz", &dimz, NULL, &status)) fits_report_error(stderr, status);
            
            float time_myr = (float)Niter*dt;
            if (fits_write_key(output_ptr, TFLOAT, (char*) "Time_Myr", &time_myr, NULL, &status)) fits_report_error(stderr, status);
            
            float* array = new float[nelements]();
            if (fits_write_img(output_ptr, TFLOAT, fpixel, nelements, array, &status)) fits_report_error(stderr, status);
            
            // <-- header HDU ; now distribution HDU -->
            
            if (fits_create_img(output_ptr, bitpix, naxis, size_axes, &status)) fits_report_error(stderr, status);
            
            int I_Z = Z;
            int I_A = A;
            int I_S = (SecEl==1);
            
            if (fits_write_key(output_ptr, TINT,   (char*) "Z_",            &I_Z, NULL, &status)) fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TINT,   (char*) "A",             &I_A, NULL, &status)) fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TINT,   (char*) "Sec",           &I_S, NULL, &status)) fits_report_error(stderr, status);
            
            ofstream tempfile;
            
            stringstream tempnamestream;
            tempnamestream << "output/" << in->filename.c_str() << "_temp_" << 1000*Z + A;
            if(SecEl) tempnamestream << "s";
            tempnamestream << "_step" << Niter << ".plain";
            char tfn[100];
            tempnamestream >> tfn;
            
            tempfile.open(tfn);
            for(int ix=0;ix<dimx;ix++)
               for(int iy=0;iy<dimy;iy++)
                  for(int iz=0;iz<dimz;iz++)
                  {
                     for(int iE=0;iE<dimE;iE++)
                        if(N[index(ix,iy,iz,iE)]>0) tempfile << coord->GetX()[ix] << '\t' << coord->GetY()[iy] << '\t' <<  coord->GetZ()[iz] << '\t' <<  coord->GetEk()[iE] << '\t' << N[index(ix,iy,iz,iE)] << endl;
                     
                     //                 tempfile << endl;
                  }
            tempfile.close();
            
            int temp_counter = 0;
            for (int k = 0; k < dimz; ++k) {
               for (int j = 0; j < dimy; ++j) {
                  for (int i = 0; i < dimx; ++i) {
                     for (int ip = 0; ip < dimE; ++ip) {
                        array[temp_counter] = N[index(i,j,k,ip)];
                        temp_counter++;
                     }
                  }
               }
            }
            
            if (fits_write_img(output_ptr, TFLOAT, fpixel, nelements, array, &status)) fits_report_error(stderr, status);
            
            if (fits_close_file(output_ptr, &status)) fits_report_error(stderr, status);
            
            delete [] array;
            
         }
         /////////////////////////////////////////////////////////
         // <--    writing temporary fits files
         /////////////////////////////////////////////////////////
         
         
      }// for Niter
      
      dt *= in->dtfactor;
      
   }//while (dt > dtmin)
   
   return;
}

// ***********************************************
// To be used in the case of Anisotropic Diffusion
// ***********************************************
void TCREvolutor3D::Run(vector<double>& N, vector<double>& N_previous, TInelasticCrossSection *xsec_object, const vector<double>& dpdt, const vector<double>& SecSource_, const vector<double>& spectrum, TDiffusionCoefficient* dperp, TReaccelerationCoefficient* dpp, double A, double Z, double lifetime, int daughter, bool SecEl, int K_electron, bool isDM, bool isExtra, Galaxy* gal) {
   
   //initialization stuff

   cout << endl << "We are starting the 3D ***ANISOTROPIC*** engine!" << endl;
   
   //In TestMode the propagation in Momentum is switched off; no source term, the evolution of a Dirac delta initial condition is followed
   //Nrept is taken from xml file. Deltat decreases as indicated in the xml file
   if (gal->GetTestMode() == true && gal->IsSourceMoving() == false)  {
      cout << "************************************** " << endl;
      cout << "*** Test mode, source is NOT moving! *" << endl;
      cout << "************************************** " << endl;
   }
   
   //In TestMode with moving source, the propagation in Momentum is switched off; the evolution of a *moving dirac delta* source term is followed
   //Nrept is fixed to 1. Deltat does not change (***set to  DeltatMax***)
   if (gal->GetTestMode() == true && gal->IsSourceMoving() == true) {
      cout << "************************************** " << endl;
      cout << "*** Test mode, source is now moving! * " << endl;
      cout << "************************************** " << endl;
   }
   
   //If the appropriate flag is specified, the Dark Matter is propagated as a moving clump
   //Nrept is fixed to 1. Deltat does not change (***set to  clump_deltat***)
   if (isDM == true && gal->IsMovingClump() == true) {
      cout << "************************************** " << endl;
      cout << "*** Moving clump of Dark Matter!!! *** " << endl;
      cout << "************************************** " << endl;
   }
   
   vector<double> xsec;
   map<pair<int,int>, vector<double> > xsec_extended;

   if (in->spallationxsec == Fluka) {
	xsec_extended = xsec_object->GetXSec_extended();
	cout << "Fluka cross sections don't work in 3D mode! " << endl;
	return;
   }	
   else
        xsec = xsec_object->GetXSec();

   double injfactor = CalcInjFactor(K_electron, SecEl, isDM, A, Z);
   
   cout << endl << "*** We are propagating this particle: *** " << endl;
   cout << "A = " << A << " Z = " << Z << " SecEl = " << SecEl << " K_electron = " << K_electron << " injfactor = " << injfactor << endl;
   if (A==0 && Z==-1 && SecEl==0)
      cout << "Primary electrons" << endl;
   if (A==0 && Z==-1 && SecEl==1)
      cout << "Secondary electrons" << endl;
   cout << "Starting propagation... " << endl;
   
   vector<double>::const_iterator itsecsource = max_element(SecSource_.begin(), SecSource_.end());
   if (injfactor == 0 && (*itsecsource) == 0) {
      cout << "No primary nor secondary source. Returning." << endl;
      return ;
   }
   
   TGas* totalgas = gal->GetTotalGas();
   TSource* source = GetSourceTerm(isDM, isExtra);
   
   TConvectionVelocity* vC = gal->GetVC();
   const double dvcdz = (vC==NULL) ? 0.0 : vC->GetDvdz()/3.0;
   
   vector<double> totmomentum;
   
   if (A == 0) totmomentum = coord->GetMomentumEl();
   else  totmomentum = coord->GetMomentum();
   
   // Following quantities depend on rigidity, but are computed independently of the nucleus, to save memory and time. Here the necessary additional nucleus-dependent factors are restored.
   
   vector<double> gamma(coord->GetGamma());
   
   // main loops
   
   int ind = -1;
   int indspat = -1;
   
   int i =  0;
   int j =  0;
   int k =  0;
   int ip = 0;
   
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,ip,indspat,ind) schedule(dynamic) num_threads(NUMTHREADS)
#endif
   for (i = 0; i < dimx; i++) {
      for (j = 0; j < dimy; j++) {
         for (k = 0; k < dimz; k++) {
            indspat = coord->indexD(i,j,k);
            double dpp1 = 0.0;
            if (dpp) dpp1 = dpp->GetReaccelerationCoefficient(indspat);
            
            for (ip = 1; ip < dimE-1; ip++) {
               double momentumup = totmomentum[ip+1];
               double momentumfix = totmomentum[ip];
               double momentumdown = totmomentum[ip-1];
               double upfix = momentumup-momentumfix;
               double updown = momentumup-momentumdown;
               double fixdown = momentumfix-momentumdown;
               
               ind = indspat*dimE+ip;
               Pdotup[ind] = dpdt[ind+1]/(upfix);
               Pdotdown[ind] = dpdt[ind]/(upfix);
               if (A != 0) {
                  Pdotup[ind] /= A;
                  Pdotdown[ind] /= A;
               }
               if (vC) {
                  Pdotup[ind] += momentumup*dvcdz/(upfix);
                  Pdotdown[ind] += momentumfix*dvcdz/(upfix);
               }
               
               if (dpp) {
                  double dppfix = dpp1*dpp->GetSpectrum(ip);
                  double dppdown = dpp1*dpp->GetSpectrum(ip-1);
                  riac1[ind] = 2.0*dppfix/(updown)/(upfix);
                  riac2[ind] = ( - (dppfix-dppdown)/(fixdown*fixdown) +
                                2.0*dppfix*(1.0/(updown)*(1.0/(upfix) + 1.0/(fixdown)) +
                                            1.0/momentumfix/(fixdown))
                                );
                  riac3[ind] = ( - (dppfix-dppdown)/(fixdown) +
                                2.0*(dppfix/(updown) +
                                     dppdown/momentumdown))/(fixdown);
                  
                  // if(A==0) cout << "[MW-DEBUG REACC] (" << A << "," << Z << ") " << i << "," << j << "," << k << " " << ip << /*" dppfactor " << dppfactor[ip] <<*/ " " << " reacccoeff " << dpp->GetReaccelerationCoefficient(indspat) << " ... so dpp1 " << dpp1 << " | " << totmomentum[ip] << " " << momentumup << " " << momentumdown << " | " << dpp->GetSpectrum(ip) << endl;
               }
            }
            ind = indspat*dimE;
            Pdotup[ind] = dpdt[ind+1]/(totmomentum[1]-totmomentum[0]);
            Pdotdown[ind] = dpdt[ind]/(totmomentum[1]-totmomentum[0]);
            if (A != 0) {
               Pdotup[ind] /= A;
               Pdotdown[ind] /= A;
            }
            
            if (vC) {
               Pdotup[ind] += totmomentum[1]*dvcdz/(totmomentum[1]-totmomentum[0]);
               Pdotdown[ind] += totmomentum[0]*dvcdz/(totmomentum[1]-totmomentum[0]);
            }
            
            // Boundary conditions
            Pdotup[ind+dimE-1] = Pdotup[ind+dimE-2];
            Pdotdown[ind+dimE-1] = Pdotdown[ind+dimE-2];
            if (dpp) {
               riac1[ind] = riac1[ind+1];
               riac2[ind] = riac2[ind+1];
               riac3[ind] = riac3[ind+1];
               riac1[ind+dimE-1] = riac1[ind+dimE-2];
               riac2[ind+dimE-1] = riac2[ind+dimE-2];
               riac3[ind+dimE-1] = riac3[ind+dimE-2];
            }
         }
      }
   }
   
   
   double dt = in->dtmax;
   double dt_official = dt;
   
   if (isDM == true && gal->IsMovingClump() == true){
      dt = in->clump_deltat;
      dt_official = in->clump_deltat;
   }
   
   //cout << " dt = " << dt << endl;
   
   double dtbar = 0.0;
   
   double halfdtbar = 0.0;
   double halfdt = 0.0;
   
   //double halfdt_dperp_factor = 0.0;
   
   const double decay = (daughter!=0);
   int Niter = 0;
   
   double value = 0.0;
   
   long counter = 0;
   long counter_total = 0;
   
   //for Moving Source (added by D.Gaggero -- 17/01/2013)
   //MW130620: What are these deltas used for? -- for now, set them equidistantially
   vector<double> x = (gal->GetCoordinates()->GetX());
   vector<double> y = (gal->GetCoordinates()->GetY());
   vector<double> z = (gal->GetCoordinates()->GetZ());
   double deltax = ( (x.back() - x.front()) / (dimx-1) );
   double deltay = ( (y.back() - y.front()) / (dimy-1) );
   double deltaz = ( (z.back() - z.front()) / (dimz-1) );
   double delta_ = in->clump_size;
   //end of initializations for Moving Source

   double total_number_of_injected_particles_so_far = 0.;
 
   //**********************
   while (dt_official > in->dtmin || (isDM == true && gal->IsMovingClump() == true)) {
      //**********************
      
      //counter++;
      //if (counter%2==0)
      //cout << counter << " dt =  " << dt << endl;
      //
      cout << "...dt =  " << dt_official << " Myr... " << endl;


      if (isDM == true && gal->IsMovingClump() == true && gal->GetInput()->analytical_refinement == true)
      	dt = (in->clump_deltat)/2.; 
	//if the analytical refinement is ON, dt is split in two. During the first half the point source is spread analytically over the grid, then the numerical computation starts.	
      
      dtbar = dt/4;  //check p!!!
      
      if (gal->GetTestMode() == true)
         dtbar = dt/3;
      
      halfdtbar = 0.5*dtbar;
      halfdt    = 0.5*dt;
      //halfdt_dperp_factor = halfdt * dperp_factor;
      
      //*****************************************
      //for the moving source
      //(added by D.Gaggero -- 17/01/2013)
      double x_0   = gal->GetSourceX0();
      double y_0   = gal->GetSourceY0();
      double z_0   = gal->GetSourceZ0();
      double source_vx = gal->GetSourceVX(); // 10000.; // kpc/Myr
      double source_vy = gal->GetSourceVY(); // 0.; // kpc/Myr
      double source_vz = gal->GetSourceVZ(); // 0.; // kpc/Myr
      double x_now = x_0 + (source_vx * dt_official * counter_total);
      double y_now = y_0 + (source_vy * dt_official * counter_total);
      double z_now = z_0 + (source_vz * dt_official * counter_total);
      
      if (gal->GetTestMode() == true && gal->IsSourceMoving() == true) {
         cout << "t [Myr] = " << dt_official*counter_total << endl;
         cout << "source_vx = " << source_vx << endl;
         cout << "source_vy = " << source_vy << endl;
         cout << "source_vz = " << source_vz << endl;
         cout << "source x coordinate [kpc] = " << x_now << endl;
         cout << "source y coordinate [kpc] = " << y_now << endl;
         cout << "source z coordinate [kpc] = " << z_now << endl;
      }
      //*****************************************
      //*****************************************
      //for a moving DM clump
      //(added by D.Gaggero -- 11/02/2013)
      double clump_x_0   = gal->GetClumpX0();
      double clump_y_0   = gal->GetClumpY0();
      double clump_z_0   = gal->GetClumpZ0();
      double clump_vx = gal->GetClumpVX(); // 10000.; // kpc/Myr
      double clump_vy = gal->GetClumpVY(); // 0.; // kpc/Myr
      double clump_vz = gal->GetClumpVZ(); // 0.; // kpc/Myr
      double clump_x_now = clump_x_0 + (clump_vx * dt_official * counter_total);
      int ixclump =  (clump_x_now - x[0])/(deltax);
      double clump_y_now = clump_y_0 + (clump_vy * dt_official * counter_total);
      int iyclump =  (clump_y_now - y[0])/(deltay);
      double clump_z_now = clump_z_0 + (clump_vz * dt_official * counter_total);
      int izclump =  (clump_z_now - z[0])/(deltaz);
      if (isDM == true && gal->IsMovingClump() == true) {
  	 cout << "t [Myr] = " << dt_official*counter_total << endl;
	 cout << "delta values = " << deltax << " " << deltay << " " << deltaz << endl;
	 cout << "grid number of points = " << dimx << " " << dimy << " " << dimz << endl;
         cout << "grid extreme values = " << x[0] << " " << y[0] << " " << z[0] << endl;
         cout << "grid extreme values = " << x[dimx-1] << " " << y[dimy-1] << " " << z[dimz-1] << endl;
         cout << "clump_vx = " << clump_vx << endl;
         cout << "clump_vy = "  << clump_vy << endl;
         cout << "clump_vz = "  << clump_vz << endl;
         cout << "clump x coordinate [kpc] = " << clump_x_now << " -- index = " << ixclump << endl;
         cout << "clump y coordinate [kpc] = " << clump_y_now << " -- index = " << iyclump << endl;
         cout << "clump z coordinate [kpc] = " << clump_z_now << " -- index = " << izclump << endl;
         cout << "** DM spectrum **" << endl;
         for (int i_p = 0; i_p < dimE; i_p++) {
            cout << spectrum[i_p] << " ";
         }
         cout << endl;
      }
      //******************************************
      
      int Nrept = in->Nrept;
      
      if (isDM == true && gal->IsMovingClump() == true)
         Nrept = 1;
      
      if (gal->GetTestMode() == true && gal->IsSourceMoving() == true)
         Nrept = 1;

//------------------------------------------TOTAL NUMBER OF INJECTED PARTICLES in the leptonic case-------------------------------------------
/*
      double total_number_of_injected_particles = 0.;
      double energy_integral=0.;
      if (isDM == true && gal->IsMovingClump() == true && gal->GetInput()->analytical_refinement == false) {
	vector<double> Ek = coord->GetEk();	
        for (int ip = 0; ip < dimE; ip++) {
    	  for (int k = 0; k < dimz; k++) {
            for (int j = 0; j < dimy; j++) {
              for (int i = 0; i < dimx; i++) {
	          long int linearized_index = index(i,j,k,ip);
		  clump_source_at_current_time[linearized_index] = norm_DM_clump*(1./(pow(2*M_PI,1.5)*pow(delta_,3))) * exp( -((x[i]-clump_x_now)*(x[i]-clump_x_now))/(2*delta_*delta_) - ((y[j]-clump_y_now)*(y[j]-clump_y_now))/(2*delta_*delta_) - ((z[k]-clump_z_now)*(z[k]-clump_z_now))/(2*delta_*delta_) ); //EDIT HANI
		   total_number_of_injected_particles += Ek[ip]*coord->GetDeltaE()*spectrum[ip]*coord->GetDeltaZ_central(k)*coord->GetDeltaY_central(j)*coord->GetDeltaX_central(i)*(pow(kpc,3.)/1.e6)*clump_source_at_current_time[linearized_index]/(Clight/1.e2/4./M_PI)*dt_official; //pure number	
		   }
             	 }
	       }
	     }
       }
*/
//----------------------------------------------------------------------------------------------------------------------------------------------
//
//----------------------------------------------------------------------------------------------------------------------------------------------
	    // DG.05.10.2013 -- the moving clump source term is computed for each timestep
	    // If the proper flag about analytical refinement was specified in the xml, during the first half of the timestep the particles emitted from the clump are spread analytically over the grid
	            
	    vector<double> clump_source_at_current_time(dimx*dimy*dimz*dimE,0.0);
	    double conversion_factor = C * pow(kpc,-3)/4./Pi*1.e4;
	    double norm_DM_clump = conversion_factor * in->clump_norm *  pow(M_Sun_in_GeV,2)/pow(kpc,3)*(in->sigmav)/(2*pow(in->mx,2))*Myr;
	    // norm_DM_clump has correct unit -->  (m^-2 s^-1 sr^-1 kpc^3 Myr^-1)
	    // so Q = clump_source * spectrum has unit --> m^-2 s^-1 Myr^-1 sr^-1 GeV^-1 
	    // so N = final result has standard DRAGON units --> m^-2 s^-1 sr^-1 GeV^-1  //again, CHECK THIS!!!!!
	    //
	    double clump_corrective_factor = 1.;
            if (isDM == true && gal->IsMovingClump() == true) {
	      /*cout << "norm_DM_clump = " << norm_DM_clump << endl;
 	      double total_number_of_particles = 0.;
              for  (int ip_ = 0; ip_ < dimE; ip_++) {
	       vector<double> Ek = coord->GetEk();	
               for  (int k_ = 0;  k_ < dimz;  k_++) {
                  for (int j_ = 0;  j_ < dimy;  j_++) {
                     for (int i_ = 0; i_ < dimx;  i_++) {
			  long int linearized_index = index(i_,j_,k_,ip_);
                          total_number_of_particles += Ek[ip_]*coord->GetDeltaE()*coord->GetDeltaZ_central(k_)*coord->GetDeltaY_central(j_)*coord->GetDeltaX_central(i_)*(pow(kpc,3.)/1.e6)*N[linearized_index]/(Clight/1.e2/4./M_PI); //pure number	
                  }
                }
               }
              }
	      cout << "---Total number of propagated particles NOW = " << total_number_of_particles << endl;	 
	      */double total_number_of_injected_particles = 0.;
	      double total_number_of_injected_particles_in_theory = 0.;//integrating analytically the delta
	      if (!(in->stop_after_timestep>0 && counter_total>=in->stop_after_timestep)){
	       double energy_integral=0.;
	       vector<double> Ek = coord->GetEk();
	       for (int ip = 0; ip < dimE; ip++) 
		energy_integral += Ek[ip]*coord->GetDeltaE()*spectrum[ip];
	       total_number_of_injected_particles_in_theory = norm_DM_clump*(pow(kpc,3.)/1.e6)/(Clight/1.e2/4./M_PI)*dt_official*energy_integral; //pure number	
	      }
	      if (isDM == true && gal->IsMovingClump() == true && gal->GetInput()->analytical_refinement == false) {
		  vector<double> Ek = coord->GetEk();	
         	  for (int ip = 0; ip < dimE; ip++) {
    		    for (int k = 0; k < dimz; k++) {
         	      for (int j = 0; j < dimy; j++) {
             	        for (int i = 0; i < dimx; i++) {
			    long int linearized_index = index(i,j,k,ip);
			    clump_source_at_current_time[linearized_index] = norm_DM_clump*(1./(pow(2*M_PI,1.5)*pow(delta_,3))) * exp( -((x[i]-clump_x_now)*(x[i]-clump_x_now))/(2*delta_*delta_) - ((y[j]-clump_y_now)*(y[j]-clump_y_now))/(2*delta_*delta_) - ((z[k]-clump_z_now)*(z[k]-clump_z_now))/(2*delta_*delta_) ); //EDIT HANI
			    if (in->stop_after_timestep>0 && counter_total>=in->stop_after_timestep)
					clump_source_at_current_time[linearized_index] = 0.;
			    total_number_of_injected_particles += Ek[ip]*coord->GetDeltaE()*spectrum[ip]*coord->GetDeltaZ_central(k)*coord->GetDeltaY_central(j)*coord->GetDeltaX_central(i)*(pow(kpc,3.)/1.e6)*clump_source_at_current_time[linearized_index]/(Clight/1.e2/4./M_PI)*dt_official; //pure number	
		          }
             	        }
		      }
		    }
	      }
	      if (isDM == true && gal->IsMovingClump() == true && gal->GetInput()->analytical_refinement == true) {
	        vector<double> Ek = coord->GetEk();
         	  for (int ip = 0; ip < dimE; ip++) { //analytical solution from Aharonian.Atoyan.1995 --- CHECK BETTER!!!
		    double p_   = totmomentum[ip];
		    double p_0  = p_ + dt*dpdt[ip]; 
		    double Ek_  = sqrt( MeleGeV*MeleGeV + p_*p_ ) - MeleGeV;
		    double Ek_0 = sqrt( MeleGeV*MeleGeV + p_0*p_0 ) - MeleGeV;
		    int i_      = ( log(Ek_)  - log(coord->GetEkMin()) )/(coord->GetDeltaE());
		    int i_0     = ( log(Ek_0) - log(coord->GetEkMin()) )/(coord->GetDeltaE());
		    double diffusion_coefficient_ip = max( dperp->GetDperpTotal(ip), dperp->GetDparTotal(ip) );
		    double r_diff = 2. * sqrt( diffusion_coefficient_ip * (Ek_0-Ek_) / dpdt[ip] );
		    for (int k = 0; k < dimz; k++) {
         	      for (int j = 0; j < dimy; j++) {
             	        for (int i = 0; i < dimx; i++) {
			    double r = sqrt( pow(x[i]-clump_x_now,2.) + pow(y[j]-clump_y_now,2.) + pow(z[k]-clump_z_now,2.) );	
			    long int linearized_index = index(i,j,k,ip);
			    clump_source_at_current_time[linearized_index] = norm_DM_clump * dpdt[i_0]/(pow(M_PI,1.5)*dpdt[i_0]*pow(r_diff,3.)) * exp(-r*r/(r_diff*r_diff));
			    if (in->stop_after_timestep>0 && counter_total>=in->stop_after_timestep)
					clump_source_at_current_time[linearized_index] = 0.;
 			    
			    N[linearized_index]	+= clump_source_at_current_time[linearized_index]*dt;	
			    total_number_of_injected_particles += Ek[ip]*coord->GetDeltaE()*spectrum[ip]*coord->GetDeltaZ_central(k)*coord->GetDeltaY_central(j)*coord->GetDeltaX_central(i)*(pow(kpc,3.)/1.e6)*clump_source_at_current_time[linearized_index]/(Clight/1.e2/4./M_PI)*dt_official; //pure number -- referring to the whole delta_t
             	        }
		      }
		    }
		  }
	      }
	      cout << "Total number of particles injected at all energies during this timestep = " << total_number_of_injected_particles << endl;	
	      cout << "Should be " << total_number_of_injected_particles_in_theory << endl;
	      //
	      if (total_number_of_injected_particles_in_theory>0. && total_number_of_injected_particles>0.)
	       clump_corrective_factor =  1.;//= total_number_of_injected_particles_in_theory/total_number_of_injected_particles;		
	      cout << "Applying a corrective factor of " << clump_corrective_factor << endl;
	      total_number_of_injected_particles_so_far += total_number_of_injected_particles*clump_corrective_factor;	            
	      cout << "Total number of particles injected so far after correction = " << total_number_of_injected_particles_so_far << endl; 
	    }
//----------------------------------------------------------------------------------------------------------------------------------------- 
		
      //*******************************************
      for (Niter = 0; Niter < Nrept; ++Niter) {
         //*******************************************
         
         counter_total++;
         
         //for the convergence check; DG 01.09.2013
         if (Niter == Nrept-1 && dt_official*in->dtfactor <= in->dtmin && gal->IsMovingClump() == false) {
            //if this is the last iteration then save N as N_previous for the convergence check!
            cout << "This is the #" << Niter << " iteration for dt = " << dt_official << "; last computation" << endl;
            N_previous = N;
         }
         
#ifdef _OPENMP
#pragma omp parallel default(shared) private(i,j,k,ip,indspat,ind) num_threads(NUMTHREADS)
#endif
         {
            
            //cout << "Niter " << Niter << endl;
            
            vector<double> Rxx(dimx,0.0);
            vector<double> dxx(dimx,0.0);
            vector<double> uodxx(dimx,0.0);
            vector<double> lodxx(dimx,0.0);
            vector<double> xx(dimx,0.0);
            
            vector<double> Ryy(dimy,0.0);
            vector<double> dyy(dimy,0.0);
            vector<double> uodyy(dimy,0.0);
            vector<double> lodyy(dimy,0.0);
            vector<double> yy(dimy,0.0);
            
            vector<double> Rzz(dimz,0.0);
            vector<double> dzz(dimz,0.0);
            vector<double> uodzz(dimz,0.0);
            vector<double> lodzz(dimz,0.0);
            vector<double> zz(dimz,0.0);
            
            vector<double> de(dimE, 0.0);
            vector<double> ee(dimE, 0.0);
            vector<double> odeu(dimE, 0.0);
            vector<double> oded(dimE, 0.0);
            vector<double> Re(dimE, 0.0);

            // **************************
            // Propagation in x direction
            // **************************
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for  (ip = 0; ip < dimE; ip++) {
               
               double halfdtbar_xsec_ip             = halfdtbar*xsec[ip];
               double halfdtbar_lifetime_gamma_ip   = halfdtbar*(decay)/(lifetime*gamma[ip]);
               double dtbar_injfactor_spec_ip       =     dtbar*injfactor*spectrum[ip];
               //double sp = dperp->GetSpectrum(ip);
               
               for  (k = 1; k < dimz-1; k++) {
                  double dz = coord->GetDeltaZ_central(k);
                  
                  for (j = 1; j < dimy-1; j++) {
                     double dy = coord->GetDeltaY_central(j);
                     
                     for (i = 0; i < dimx;  i++) {
                        double dx = coord->GetDeltaX_central(i);
                        
                        indspat = coord->indexD(i,j,k);
                        ind = indspat*dimE+ip;
                        
                        double alphax = dperp->GetAlpha_XX(indspat, ip);
                        
                        /*if (k == dimz/2 && ip == dimE/2) {
                         cout << " i = " << i << " j = " << j << " k = " << k << " ip = " << ip << endl;
                         cout << " alpha_x = " << alphax << endl;
                         }*/
                        
                        double halfdt_coeffx_dp =  halfdt*alphax/dx/dx;
                        double phix = dperp->GetPhiX(indspat,ip);
                        double halfdt_dperp_factor_phix = halfdt*phix;
                        double dtbar_coeffx_dp_gas_lifetime  =  2.0*halfdt_coeffx_dp  +  totalgas->GetGas(indspat)*halfdtbar_xsec_ip  +  halfdtbar_lifetime_gamma_ip ;
                        
                        
                        if (isDM == true && gal->IsMovingClump() == true)
                           double dtbar_coeffx_dp_gas_lifetime  =  2.0*halfdt_coeffx_dp;//  +  totalgas->GetGas(indspat)*halfdtbar_xsec_ip  +  halfdtbar_lifetime_gamma_ip ;
                        
                        if (gal->GetTestMode() == true)
                           dtbar_coeffx_dp_gas_lifetime  =  2.0*halfdt_coeffx_dp;
                        
                        dxx[i] = 1. + dtbar_coeffx_dp_gas_lifetime;
                        uodxx[i] = -halfdt_coeffx_dp - halfdt_dperp_factor_phix;
                        lodxx[i] = -halfdt_coeffx_dp + halfdt_dperp_factor_phix;
                        
                        Rxx[i] = N[index(i,j,k,ip)]*(1.0 - dtbar_coeffx_dp_gas_lifetime) + source->GetSource(indspat)*dtbar_injfactor_spec_ip + SecSource_[ind]*dtbar;
                        
                        
                        
                        //**************************************************************
                        if (gal->GetTestMode() == true && gal->IsSourceMoving() == false)
                           //In TestMode with non-moving source: **no source term**, the evolution of a Dirac delta initial condition is followed
                           //**************************************************************
                           
                           Rxx[i] = N[index(i,j,k,ip)]*(1.0 - dtbar_coeffx_dp_gas_lifetime);// + source->GetSource(indspat)*dtbar_injfactor_spec_ip + SecSource_[ind]*dtbar;
                        
                        //**************************************************************
                        if (gal->GetTestMode() == true && gal->IsSourceMoving() == true) {
                           //In TestMode with moving source, the evolution of a *moving dirac delta* source term is followed
                           //**************************************************************
                           
                           double source_at_current_time  = (1./(sqrt(3.14)*delta_)) * exp( -((x[i]-x_now)*(x[i]-x_now))/(delta_*delta_) - ((y[j]-y_now)*(y[j]-y_now))/(delta_*delta_) - ((z[k]-z_now)*(z[k]-z_now))/(delta_*delta_) );
                           Rxx[i] = N[index(i,j,k,ip)]*(1.0 - dtbar_coeffx_dp_gas_lifetime) + source_at_current_time*dtbar;
                        }
                        
                        
                        //**************************************************************
                        if (isDM == true && gal->IsMovingClump() == true) {
                           // if the proper flag was specified, if it is a DM particle, it is considered as a DM moving clump!
                           //**************************************************************
                           
                           //double source_at_current_time  = (1./(sqrt(3.14)*delta_)) * exp( -((x[i]-clump_x_now)*(x[i]-clump_x_now))/(delta_*delta_) - ((y[j]-clump_y_now)*(y[j]-clump_y_now))/(delta_*delta_) - ((z[k]-clump_z_now)*(z[k]-clump_z_now))/(delta_*delta_) );
                           //if (counter*dt > 0.5)
                           //source_at_current_time = 0.;
                           Rxx[i] = N[index(i,j,k,ip)]*(1.0 - dtbar_coeffx_dp_gas_lifetime) + clump_source_at_current_time[index(i,j,k,ip)]*dtbar_injfactor_spec_ip*clump_corrective_factor;
                        }
                        
                        //**************************************************************
                        
                        if (i < dimx-1) {
                           Rxx[i] -= N[index(i+1,j,k,ip)] * uodxx[i] ;
                           Rxx[i] += dperp->GetAlpha_XY(indspat,ip)*halfdt/(4.0*dy*dx)*(N[index(i+1,j+1,k,ip)] - N[index(i+1,j-1,k,ip)]);
                           Rxx[i] += dperp->GetAlpha_XZ(indspat,ip)*halfdt/(4.0*dz*dx)*(N[index(i+1,j,k+1,ip)] - N[index(i+1,j,k-1,ip)]);
                        }
                        if (i > 0) {
                           Rxx[i] -= N[index(i-1,j,k,ip)] * lodxx[i] ;
                           Rxx[i] -= dperp->GetAlpha_XY(indspat,ip)*halfdt/(4.0*dy*dx)*(N[index(i-1,j+1,k,ip)] - N[index(i-1,j-1,k,ip)]);
                           Rxx[i] -= dperp->GetAlpha_XZ(indspat,ip)*halfdt/(4.0*dz*dx)*(N[index(i-1,j,k+1,ip)] - N[index(i-1,j,k-1,ip)]);
                        }
                        
                        // if(A==0) cout << "[MW-DEBUG X] (" << A << "," << Z << ") " << i << " " << j << " " << k << " | " /*<< CNalphax1 << " " << CNalphax2 << " " << CNalphax3 << " | "*/ << lodxx[i] << " " << dxx[i] << " "  << uodxx[i] << " " << Rxx[i] << " | " << gal->GetCoordinates()->GetDeltaX_central(i) << " " << coord->GetDeltaX_central(i) << "=" << dx << endl;
                        
                     }
                     
                     //Utility::solve_tridag( &(lodxx[0]), &(dxx[0]), &(uodxx[0]), &(Rxx[0]), &(xx[0]), dimx );
		     Utility::solve_tridag(lodxx, dxx, uodxx, Rxx, xx, dimx);
                     
                     for (int i = dimx-2; i >0; --i) {
                        value = xx[i];
                        N[index(i,j,k,ip)] = (value > 0) ? value : 0.0;
                     }
                  }
               }
               
            } // for ip
            
            
            // **************************
            // Propagation in y direction
            // **************************
            
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for  (ip = 0; ip < dimE; ip++) {
               
               
               double halfdtbar_xsec_ip            = halfdtbar*xsec[ip];
               double halfdtbar_lifetime_gamma_ip   = halfdtbar*(decay)/(lifetime*gamma[ip]);
               double dtbar_injfactor_spec_ip       =     dtbar*injfactor*spectrum[ip];
               //double sp = dperp->GetSpectrum(ip);
               
               
               for  (k = 1; k < dimz-1; k++) {
                  double dz = coord->GetDeltaZ_central(k);
                  for (i = 1; i < dimx-1; i++) {
                     double dx = coord->GetDeltaX_central(i);
                     
                     for (j = 0; j < dimy;  j++) {
                        double dy = coord->GetDeltaY_central(j);
                        
                        indspat = coord->indexD(i,j,k);
                        ind = indspat*dimE+ip;
                        
                        double alphay= dperp->GetAlpha_YY(indspat, ip);
                        
                        
                        /*if (k == dimz/2 && ip == dimE/2) {
                         cout << " i = " << i << " j = " << j << " k = " << k << " ip = " << ip << endl;
                         cout << " alpha_y = " << alphay << endl;
                         }*/
                        
                        double halfdt_coeffy_dp =        halfdt*alphay/dy/dy;
                        double phiy = dperp->GetPhiY(indspat,ip);
                        double halfdt_dperp_factor_phiy = halfdt*phiy;
                        double dtbar_coeffy_dp_gas_lifetime  =  2.0*halfdt_coeffy_dp  +  totalgas->GetGas(indspat)*halfdtbar_xsec_ip  +  halfdtbar_lifetime_gamma_ip ;
                        
                        
                        if (isDM == true && gal->IsMovingClump() == true)
                           double dtbar_coeffy_dp_gas_lifetime  =  2.0*halfdt_coeffy_dp;//  +  totalgas->GetGas(indspat)*halfdtbar_xsec_ip  +  halfdtbar_lifetime_gamma_ip ;
                        
                        if (gal->GetTestMode() == true)
                           dtbar_coeffy_dp_gas_lifetime  =  2.0*halfdt_coeffy_dp;
                        
                        dyy[j] = 1. + dtbar_coeffy_dp_gas_lifetime;
                        uodyy[j] = -halfdt_coeffy_dp - halfdt_dperp_factor_phiy;
                        lodyy[j] = -halfdt_coeffy_dp + halfdt_dperp_factor_phiy;
                        
                        Ryy[j] = N[index(i,j,k,ip)]*(1.0 - dtbar_coeffy_dp_gas_lifetime) + source->GetSource(indspat)*dtbar_injfactor_spec_ip  + SecSource_[ind]*dtbar;;
                        
                        
                        if (gal->GetTestMode() == true && gal->IsSourceMoving() == false)
                           Ryy[j] = N[index(i,j,k,ip)]*(1.0 - dtbar_coeffy_dp_gas_lifetime); // + source->GetSource(indspat)*dtbar_injfactor_spec_ip  + SecSource_[ind]*dtbar;;
                        
                        if (gal->GetTestMode() == true && gal->IsSourceMoving() == true) {
                           double source_at_current_time  = (1./(sqrt(3.14)*delta_)) * exp( -((x[i]-x_now)*(x[i]-x_now))/(delta_*delta_) - ((y[j]-y_now)*(y[j]-y_now))/(delta_*delta_) - ((z[k]-z_now)*(z[k]-z_now))/(delta_*delta_) );
                           Ryy[j] = N[index(i,j,k,ip)]*(1.0 - dtbar_coeffy_dp_gas_lifetime) + source_at_current_time*dtbar;
                        }
                        
                        if (isDM == true && gal->IsMovingClump() == true) {
                           //double source_at_current_time  = (1./(sqrt(3.14)*delta_)) * exp( -((x[i]-clump_x_now)*(x[i]-clump_x_now))/(delta_*delta_) - ((y[j]-clump_y_now)*(y[j]-clump_y_now))/(delta_*delta_) - ((z[k]-clump_z_now)*(z[k]-clump_z_now))/(delta_*delta_) );
                           //if (counter*dt > 0.5)
                           //source_at_current_time = 0;
                           Ryy[j] = N[index(i,j,k,ip)]*(1.0 - dtbar_coeffy_dp_gas_lifetime) + clump_source_at_current_time[index(i,j,k,ip)]*dtbar_injfactor_spec_ip*clump_corrective_factor;
                        }
                        
                        
                        if (j < dimy-1) {
                           Ryy[j] -= N[index(i,j+1,k,ip)] * uodyy[j] ;
                           Ryy[j] += dperp->GetAlpha_XY(indspat,ip)*halfdt/(4.0*dy*dx)*(N[index(i+1,j+1,k,ip)] - N[index(i-1,j+1,k,ip)]);
                           Ryy[j] += dperp->GetAlpha_YZ(indspat,ip)*halfdt/(4.0*dy*dz)*(N[index(i,j+1,k+1,ip)] - N[index(i,j+1,k-1,ip)]);
                        }
                        if (j > 0) {
                           Ryy[j] -= N[index(i,j-1,k,ip)] * lodyy[j] ;
                           Ryy[j] -= dperp->GetAlpha_XY(indspat,ip)*halfdt/(4.0*dy*dx)*(N[index(i+1,j-1,k,ip)] - N[index(i-1,j-1,k,ip)]);
                           Ryy[j] -= dperp->GetAlpha_YZ(indspat,ip)*halfdt/(4.0*dz*dy)*(N[index(i,j-1,k+1,ip)] - N[index(i,j-1,k-1,ip)]);
                        }
                        
                        // cout << "[MW-DEBUG Y] " << i << " " << j << " " << k << " | "/* << CNalphay1 << " " << CNalphay2 << " " << CNalphay3 << " | " */<< lodyy[j] << " " << dyy[j] << " "  << uodyy[j] << " " << Ryy[j] << " " << endl;
                        
                     }
                     
		     //Utility::solve_tridag( &(lodyy[0]), &(dyy[0]), &(uodyy[0]), &(Ryy[0]), &(yy[0]), dimy);	
                     Utility::solve_tridag(lodyy, dyy, uodyy, Ryy, yy, dimy);
                     
                     for (int j = dimy-2; j >0; --j) {
                        value = yy[j];
                        N[index(i,j,k,ip)] = (value > 0) ? value : 0.0;
                     }
                  }
               }
               
            }// for ip
            
            // **************************
            // Propagation in z direction
            // **************************
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for  (ip = 0; ip < dimE; ip++) {
               
               double halfdtbar_xsec_ip            = halfdtbar*xsec[ip];
               double halfdtbar_lifetime_gamma_ip   = halfdtbar*(decay)/(lifetime*gamma[ip]);
               double dtbar_injfactor_spec_ip       =     dtbar*injfactor*spectrum[ip];
               //double sp = dperp->GetSpectrum(ip);
               
               for (i = 1; i < dimx-1;  i++) {
                  double dx = coord->GetDeltaX_central(i);
                  
                  for (j = 1; j < dimy-1; j++) {
                     double dy = coord->GetDeltaY_central(j);
                     
                     for  (k = 0; k < dimz; k++) {
                        double dz = coord->GetDeltaZ_central(k);
                        
                        indspat = coord->indexD(i,j,k);
                        ind = indspat*dimE+ip;
                        
                        double alphaz = dperp->GetAlpha_ZZ(indspat, ip);
                        
                        double halfdt_coeffz_dp =        halfdt*alphaz/dz/dz;
                        double phiz = dperp->GetPhiZ(indspat, ip);
                        double halfdt_dperp_factor_phiz = halfdt*phiz;
                        double dtbar_coeffz_dp_gas_lifetime  =  2.0*halfdt_coeffz_dp +  totalgas->GetGas(indspat)*halfdtbar_xsec_ip  +  halfdtbar_lifetime_gamma_ip ;
                        
                        if (isDM == true && gal->IsMovingClump() == true)
                           double dtbar_coeffz_dp_gas_lifetime  =  2.0*halfdt_coeffz_dp;//  +  totalgas->GetGas(indspat)*halfdtbar_xsec_ip  +  halfdtbar_lifetime_gamma_ip ;
                        
                        if (gal->GetTestMode() == true)
                           dtbar_coeffz_dp_gas_lifetime = 2.0*halfdt_coeffz_dp;
                        
                        //MW130720: Didn't have any convection?
                        double vCk = 0.0; // vC(i)
                        double vCk1 = 0.0; // vC(i+1)
                        double vC1k = 0.0; // vC(i-1)
                        
                        if (vC)
                        {
                           vC1k = vC->GetCNconv_alpha1_z(indspat);
                           vCk  = vC->GetCNconv_alpha2_z(indspat);
                           vCk1 = vC->GetCNconv_alpha3_z(indspat);
                        }
                        
                        dzz[k] = 1. + dtbar_coeffz_dp_gas_lifetime  + halfdt*vCk;
                        uodzz[k] = -halfdt_coeffz_dp - halfdt_dperp_factor_phiz -halfdt*vCk1;
                        lodzz[k] = -halfdt_coeffz_dp + halfdt_dperp_factor_phiz -halfdt*vC1k;
                        
                        Rzz[k] = N[index(i,j,k,ip)]*(2.0 - dzz[k]) + source->GetSource(indspat)*dtbar_injfactor_spec_ip  + SecSource_[ind]*dtbar;
                        
                        if (gal->GetTestMode() == true && gal->IsSourceMoving() == false)
                           Rzz[k] = N[index(i,j,k,ip)]*(1.0 - dtbar_coeffz_dp_gas_lifetime);// + source->GetSource(indspat)*dtbar_injfactor_spec_ip  + SecSource_[ind]*dtbar;;
                        
                        if (gal->GetTestMode() == true && gal->IsSourceMoving() == true) {
                           double source_at_current_time  = (1./(sqrt(3.14)*delta_)) * exp( -((x[i]-x_now)*(x[i]-x_now))/(delta_*delta_) - ((y[j]-y_now)*(y[j]-y_now))/(delta_*delta_) - ((z[k]-z_now)*(z[k]-z_now))/(delta_*delta_) );
                           Rzz[k] = N[index(i,j,k,ip)]*(1.0 - dtbar_coeffz_dp_gas_lifetime) + source_at_current_time*dtbar;
                        }
                        
                        if (isDM == true && gal->IsMovingClump() == true) {
                           //double source_at_current_time  = (1./(sqrt(3.14)*delta_)) * exp( -((x[i]-clump_x_now)*(x[i]-clump_x_now))/(delta_*delta_) - ((y[j]-clump_y_now)*(y[j]-clump_y_now))/(delta_*delta_) - ((z[k]-clump_z_now)*(z[k]-clump_z_now))/(delta_*delta_) );
                           //if (counter*dt > 0.5)
                           //source_at_current_time = 0;
                           Rzz[k] = N[index(i,j,k,ip)]*(1.0 - dtbar_coeffz_dp_gas_lifetime) + clump_source_at_current_time[index(i,j,k,ip)]*dtbar_injfactor_spec_ip*clump_corrective_factor;
                        }
                        
                        if (k < dimz-1) {
                           Rzz[k] -= N[index(i,j,k+1,ip)] * uodzz[k] ;
                           Rzz[k] += dperp->GetAlpha_XZ(indspat,ip)*halfdt/(4.0*dz*dx)*(N[index(i+1,j,k+1,ip)] - N[index(i-1,j,k+1,ip)]);
                           Rzz[k] += dperp->GetAlpha_YZ(indspat,ip)*halfdt/(4.0*dy*dz)*(N[index(i,j+1,k+1,ip)] - N[index(i,j-1,k+1,ip)]);
                        }
                        if (k > 0) {
                           Rzz[k] -= N[index(i,j,k-1,ip)] * lodzz[k] ;
                           Rzz[k] -= dperp->GetAlpha_XZ(indspat,ip)*halfdt/(4.0*dz*dx)*(N[index(i+1,j,k-1,ip)] - N[index(i-1,j,k-1,ip)]);
                           Rzz[k] -= dperp->GetAlpha_YZ(indspat,ip)*halfdt/(4.0*dz*dy)*(N[index(i,j+1,k-1,ip)] - N[index(i,j-1,k-1,ip)]);
                        }
                        
                        // if(A==0) cout << "[MW-DEBUG Z] (" << A << "," << Z << ") " << i << " " << j << " " << k << " | " /*<< CNalphaz1 << " " << CNalphaz2 << " " << CNalphaz3 << " " */ << vC1k << " " << vCk << " " << vCk1 << " | " << lodzz[k] << " " << dzz[k] << " "  << uodzz[k] << " " << Rzz[k] << " | " << N[index(i,j,k,ip)] << " " << source->GetSource(indspat) << " " << dtbar_injfactor_spec_ip << " " << SecSource_[ind] << " " << dtbar << endl;
                        
                     }
                     
                     //Utility::solve_tridag( &(lodzz[0]), &(dzz[0]), &(uodzz[0]), &(Rzz[0]), &(zz[0]), dimz );
		     Utility::solve_tridag(lodzz, dzz, uodzz, Rzz, zz, dimz);                     

                     for (int k = dimz-2; k >0; --k) {
                        value = zz[k];
                        N[index(i,j,k,ip)] = (value > 0) ? value : 0.0;
                     }
                  }
               }
               
            }
            
            // **************************
            // Propagation in p direction
            // **************************
            
            //In TestMode the propagation in Momentum is switched off
            if (gal->GetTestMode() == false) {
               
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
               for (i = 1; i < dimx-1;  i++) {
                  for (j = 1; j < dimy-1; j++) {
                     for (k = 1; k < dimz-1; k++) {
                        
                        indspat = coord->indexD(i,j,k);
                        double totgas = totalgas->GetGas(indspat);
                        double dtbarprimsource = dtbar*injfactor*source->GetSource(indspat);
                        
                       
                        ind = indspat*dimE;
                        
                        for (ip = 0; ip < dimE; ip++) {

                           if (isDM == true && gal->IsMovingClump() == true) {
                            //double source_at_current_time  = (1./(sqrt(3.14)*delta_)) * exp( -((x[i]-clump_x_now)*(x[i]-clump_x_now))/(delta_*delta_) - ((y[j]-clump_y_now)*(y[j]-clump_y_now))/(delta_*delta_) - ((z[k]-clump_z_now)*(z[k]-clump_z_now))/(delta_*delta_) );
                            dtbarprimsource = dtbar*injfactor*clump_source_at_current_time[index(i,j,k,ip)]*clump_corrective_factor;
                           }
                           
                           int ind1 = ind + ip;
                           
                           double devect = 1.0
                           + halfdtbar*(totgas*xsec[ip] + (decay)/(lifetime*gamma[ip]))
                           + halfdt*(riac2[ind1] + Pdotdown[ind1]);
                           double odeuvect = halfdt*(riac1[ind1] + Pdotup[ind1]);
                           de[ip] = devect;
                           odeu[ip] = -odeuvect;
                           oded[ip] = -halfdt*riac3[ind1];
                           
                           // if(i==20 && j==20) cout << "[MW-DEBUG P] (" << A << "," << Z << ") " << i << " " << j << " " << k << " " << ip << " | " << Pdotdown[ind1] << " " << Pdotup[ind1] << " " << xsec[ip] << " " << gamma[ip] << " | " << oded[ip] << " " << de[ip] << " "  << odeu[ip] << endl;
                           
                           Re[ip] = dtbarprimsource*spectrum[ip] + dtbar*SecSource_[ind1] + (2.0-devect)*N[ind1];
                           
                           //Maybe redundant, check!!!
                           if (isDM == true && gal->IsMovingClump() == true)
                              Re[ip] = dtbarprimsource*spectrum[ip] + (2.0-devect)*N[ind1];
                           
                           if (ip > 0) Re[ip] -= oded[ip]*N[ind1-1];
                           if (ip < dimE-1) Re[ip] += odeuvect*N[ind1+1];
                        }
                        
                        //Utility::solve_tridag(&(oded[0]), &(de[0]), &(odeu[0]), &(Re[0]), &(ee[0]), dimE);
			Utility::solve_tridag(oded, de, odeu, Re, ee, dimE);
                        
                        for (ip = 0; ip < dimE; ip++) {
                           value = ee[ip];
                           N[ind+ip] = (value > 0) ? value : 0.0;
                        }
                     }
                  } // dimy-1
               } // dimr-1
               
            }
            
            
         } //#pragma omp parallel
         
         //********************************
	 if (isDM == true && gal->IsMovingClump() == true) {   
	    double total_number_of_particles = 0.;
            for  (int ip = 0; ip < dimE; ip++) {
	       vector<double> Ek = coord->GetEk();	
               for  (int k = 0;  k < dimz;  k++) {
                  for (int j = 0;  j < dimy;  j++) {
                     for (int i = 0; i < dimx;  i++) {
			  long int linearized_index = index(i,j,k,ip);
                          total_number_of_particles += Ek[ip]*coord->GetDeltaE()*coord->GetDeltaZ_central(k)*coord->GetDeltaY_central(j)*coord->GetDeltaX_central(i)*(pow(kpc,3.)/1.e6)*N[linearized_index]/(Clight/1.e2/4./M_PI); //pure number	
                  }    
                }
              }
            }
	    cout << "---Total number of propagated particles NOW = " << total_number_of_particles << endl;	 
	 }
	 //********************************  

         //for the convergence check; DG 01.09.2013
         if (Niter == Nrept-1 && dt_official*in->dtfactor <= in->dtmin && gal->IsMovingClump() == false) {
            //if we are at the last iteration, do the convergence check
            //N_previous and N contain the particle vector at the last-but-one and last iteration
            cout << endl << "*** Convergence check ***" << endl;
            for  (int ip = 0; ip < dimE; ip++) {
               cout << "***" << endl;
               cout << "Energy index -> " << ip << endl;
               double average_difference = 0.;
               double maximum_difference = 0.;
               long counter = 0;
               for  (int k = 0;  k < dimz;  k++) {
                  for (int j = 0;  j < dimy;  j++) {
                     for (int i = 0; i < dimx;  i++) {
                        long indspat = coord->indexD(i,j,k);
                        long ind = indspat*dimE+ip;
                        if (N[ind] != 0) {
                           counter++;
                           double difference = fabs((N_previous[ind] - N[ind])/N[ind]);
                           if (difference > maximum_difference)
                              maximum_difference = difference;
                           average_difference += difference;
                        }
                     }
                  }
               }
               if (counter > 0) average_difference /= counter;
               cout << "Nonzero matrix elements -> " << counter << endl;
               cout << "Average difference -> " << average_difference*100. << "% " << endl;
               cout << "Maximum difference -> " << maximum_difference*100. << "% " << endl;
            }
            cout << "*** End of convergence check ***" << endl;
         }
         
      }// for Niter
      
      bool StopHere = false;
      
      if ( (isDM == true && gal->IsMovingClump() == true) || (gal->GetTestMode() == true && gal->IsSourceMoving() == true) ) {
         
         //if the source is moving, the deltat is kept constant!!
         dt *= 1.; dt_official *= 1.;
         if (counter_total > 5000) StopHere = true;
      }
      else {
         dt *= in->dtfactor;
         dt_official *= in->dtfactor;
      }
      
      if (StopHere == true)
         break;
      
      //counter_total++;
      
      if ((gal->GetTestMode() == true && A==1 && Z==1 && !SecEl) || (A==0 && Z==-1 && isDM == true && gal->IsMovingClump()==true))  { // In Test mode, prints the primary proton density matrix in the center of the energy range every 2 timesteps
			
         if (counter_total == 1 || counter_total%20==0) {
            
            string name;
            stringstream namestream;
            
            cout << "****************************************************** " << endl;
            cout << "Writing temporary output file corresponding to " << dt_official*counter_total << " Myr" << endl;
            cout << "****************************************************** " << endl;
	    cout << "counter = " << counter_total << endl;	
            
            namestream << "output/temporary_output_aniso_" << counter_total << "_near_2.fits.gz";
            
            name = namestream.str();
            
            fitsfile* output_ptr;
            
            int status = 0;
            
            const long naxis = 4;
            long size_axes[naxis] = {dimE,dimx,dimy,dimz};
            long nelements = size_axes[0]*size_axes[1]*size_axes[2]*size_axes[3];
            cout << "nelements = " << nelements ;
            //long totelements = nelements * dimE;
            //cout << "totelements = " << totelements << endl;
            long fpixel = 1;
            int bitpix = FLOAT_IMG;
            
            //int ip_ = dimE/2;
            vector<double> energy_vec = coord->GetEk();
            cout << " Writing temp file " << endl; //<< " at energy = " << energy_vec[ip_] << endl;
            
            if (fits_create_file(&output_ptr, name.c_str(), &status)) fits_report_error(stderr, status);
            
            //***********************************************************************
            //create the first HDU, with metadata only -- like DRAGON standard output
            //***********************************************************************
            if (fits_create_img(output_ptr, bitpix, naxis, size_axes, &status)) fits_report_error(stderr, status);
            
            if (fits_write_key(output_ptr, TINT, (char*) "dimx",     &dimx,             NULL, &status))       fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TINT, (char*) "dimy",     &dimy,             NULL, &status))       fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TINT, (char*) "dimz",     &dimz,             NULL, &status))       fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TINT, (char*) "dimE",     &dimE,             NULL, &status))       fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TINT, (char*) "ixclump",     &ixclump,             NULL, &status)) fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TINT, (char*) "iyclump",     &iyclump,             NULL, &status)) fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TINT, (char*) "izclump",     &izclump,             NULL, &status)) fits_report_error(stderr, status);

	    unsigned int ixsun = (unsigned int) ((in->xobs-x.front())/(2.0*x.back())*(double)(dimx-1));
	    unsigned int iysun = (unsigned int) ((in->yobs-y.front())/(2.0*y.back())*(double)(dimy-1));
	    unsigned int izsun = (unsigned int) ((in->zobs-z.front())/(2.0*z.back())*(double)(dimz-1));
	    cout << "Sun position indexes = " << ixsun << " " << iysun << " " << izsun << endl;
            if (fits_write_key(output_ptr, TINT, (char*) "ixsun",     &ixsun,             NULL, &status)) fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TINT, (char*) "iysun",     &iysun,             NULL, &status)) fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TINT, (char*) "izsun",     &izsun,             NULL, &status)) fits_report_error(stderr, status);

            
            double xmax  = x.back() ;
            double xmin  = x.front();
            double zmax  = z.back();
            double zmin  = z.front();
            double Ekmin =  (gal->GetCoordinates()->GetEkMin());
            double Ekmax =  (gal->GetCoordinates()->GetEkMax());
            double DeltaE = (gal->GetCoordinates()->GetDeltaE());
            double Ekin_f = exp(DeltaE);
            
            if (fits_write_key(output_ptr, TDOUBLE, (char*) "xmin",     &xmin,             NULL, &status))   fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TDOUBLE, (char*) "xmax",     &xmax,             NULL, &status))   fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TDOUBLE, (char*) "zmin",     &zmin,             NULL, &status))   fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TDOUBLE, (char*) "zmax",     &zmax,             NULL, &status))   fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TDOUBLE, (char*) "Ekmin",    &Ekmin,             NULL, &status))  fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TDOUBLE, (char*) "Ekmax",    &Ekmax,             NULL, &status))  fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TDOUBLE, (char*) "Ekin_fac", &Ekin_f,             NULL, &status)) fits_report_error(stderr, status);
            
            float* array1 = new float[nelements];
            if (fits_write_img(output_ptr, TFLOAT, fpixel, nelements, array1, &status)) fits_report_error(stderr, status);
            
            //***********************************************************************
            //create data HDU
            //***********************************************************************
            if (fits_create_img(output_ptr, bitpix, naxis, size_axes, &status)) fits_report_error(stderr, status);
            
            int Z_pos = 1;
            int A_pos = 0;
            
            if (fits_write_key(output_ptr, TINT, (char*) "Z_",    &Z_pos,             NULL, &status))       fits_report_error(stderr, status);
            if (fits_write_key(output_ptr, TINT, (char*) "A",     &A_pos,             NULL, &status))       fits_report_error(stderr, status);
            
            int temp_counter = 0;
            float* array = new float[nelements];
            for (int k = 0;  k  < dimz; ++k) {
               for (int j = 0;  j < dimy; ++j) {
                  for (int i = 0; i < dimx; ++i) {
                     for (int ie = 0;  ie < dimE; ++ie) {
                        array[temp_counter] = N[index(i,j,k,ie)];
                        temp_counter++;
                     }
                  }
               }
            }
            
            if (fits_write_img(output_ptr, TFLOAT, fpixel, nelements, array, &status)) fits_report_error(stderr, status);
            
            if (fits_close_file(output_ptr, &status)) fits_report_error(stderr, status);
            
            delete [] array;
         }
      }
      
   }//while (dt > dtmin)
   
   return;
}



