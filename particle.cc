/**
 * @file particle.cc
 * @author Luca Maccione, Daniele Gaggero
 * @email luca.maccione@desy.de
 * @email daniele.gaggero@sissa.it
 * @brief Implementation of TParticle class. See the .h file
 */


#include "constants.h"
#include "geometry.h"
#include "particle.h"
#include "galaxy.h"
#include "eloss.h"
#include "xsec.h"
#include "crevolutor.h"
#include "grid.h"
#include "gas.h"
#include "input.h"
#include "spectrum.h"
#include "errorcode.h"
#include "diffusion.h"

#ifdef HAVE_ROOT
#include "TFile.h"
#include "TNtupleD.h"
#include "TString.h"
#endif

using namespace std;

//******************************************************************************************************************************************************
vector<double> TParticle::GetEnergyGrid() const {
   return _fGalaxy->GetCoordinates()->GetEk();
}

TParticle::TParticle(int A_ /**< Mass number */, int Z_ /**< Charge */, Galaxy* gal /**< A model for the galaxy */, Input* in_ /**< User input */, bool issec_, vector<TXSecBase*> xsecmodel, TNucleiList* l, int K_electron_, bool isDM_, bool isextra, bool isTPP_) :
A(A_),
Z(Z_),
uid(1000*Z_+A_),
dimz(gal->GetCoordinates()->GetDimZ()),
dimE(gal->GetCoordinates()->GetDimE()),
K_electron(K_electron_) {

  in = in_;

  if (in->feedback >0) cout << endl;
  if (in->feedback >0) cout << "********************************************" << endl;
  if (in->feedback >0) cout << "*** Welcome to the Particle constructor! ***" << endl;
  if (in->feedback >0) cout << "********************************************" << endl << endl;
  
  if (in->DiffT == Anisotropic) {
      _fDiff = new TDiffusionCoefficient3D(gal->GetCoordinates(), in, gal->GetSource(), gal->GetBField(), gal->GetGeometry(), A, Z, (A==0));
      
      vector<double> momentum = (A!=0) ? gal->GetCoordinates()->GetMomentum() : gal->GetCoordinates()->GetMomentumEl();
      _fDpp = new TReaccelerationCoefficient(momentum, _fDiff, gal->GetGeometry(), in);
   }
   else {
      _fDiff = NULL;
      _fDpp = NULL;
   }
   isDM = isDM_;
    if (in->feedback >0) cout << "IsDM? " << isDM << endl;
   //modified
   issec = issec_;
   isExtra = isextra;
    if (in->feedback >0) cout << "IsEXTRA? " << isExtra << endl;
   isTPP = isTPP_;
   
   _fGalaxy = gal;
   decmode = l->GetDecayMode(uid);
   lifetime = l->GetLifeTime(uid)/M_LN2;


   // possible decay modes ----- DG 24.10.2013 - SK 07.11.2013 ------------------------------------------------------   
    if (decmode == BP)  daughter = 1000*(Z-1)+A;
    else if (decmode == BM) {
        if ((Z == 14 && A == 32) || (Z == 18 && A ==42)) daughter = 1000*(Z+2)+A; //double beta- for these isotopes
        else daughter = 1000*(Z+1)+A;
	// 32Si --> 32P (14 yr) --> 32S || 42Ar(2b-)42Ca  100%  Ar-K -Ca
    }
    else if (decmode == EC) {     
        if (K_electron <= 0) {  
            lifetime = -1.; // If the nucleus is unstable for EC the "naked" version of these isotopes does not decay
            daughter = 0;
        }
        else { //modified 
            lifetime = l->GetLifeTime(uid)/M_LN2; // If the nucleus is unstable for EC the decay may happen only after it has attached an electron 
            daughter = 1000*(Z-1)+A;
        }
    }
    else if (decmode == ECBM) {     
        if (K_electron <= 0) {  
         lifetime = l->GetLifeTime_naked(uid)/M_LN2;; // Beta- decay is possible for the naked nucleus 
         daughter = 1000*(Z+1)+A;
        }
        else { //modified 
            lifetime = l->GetLifeTime(uid)/M_LN2; // If the nucleus is unstable for EC the decay may happen only after it has attached an electron 
            daughter = 1000*(Z-1)+A;
        }
    }
    else if (decmode == ECBP) { 
        if (K_electron <= 0) {  
         lifetime = l->GetLifeTime_naked(uid)/M_LN2;; //  Beta+ decay is possible for the naked nucleus 
         daughter = 1000*(Z-1)+A;
         if(Z == 28 && A == 56) daughter = 1000*(Z-2)+A;  //56Ni makes a double Beta+ decay
        }
        else { //modified 
            lifetime = l->GetLifeTime(uid)/M_LN2; // If the nucleus is unstable for EC the decay may happen only after it has attached an electron 
            daughter = 1000*(Z-1)+A;
        }
    }
    else daughter = 0;
   //-----------------------------------------------------------------------------------------------------------------

   if(!in->RDECAY){ //fk 130701
      lifetime = -1; // -1 means stable!
      daughter =0;
      if (in->feedback >0) cout << "Radioactive Decay is turned off so every nucleus is stable!" << endl;
   }
   
   if (in->feedback >0) cout << "Computing spectrum " << endl;
   
   if (issec) sp = new TSpectrum(gal->GetCoordinates()); // Secondary particles
   else {
     if (!isDM) {
       if (isExtra) sp = new TSpectrum(gal->GetCoordinates(), in, gal->GetInjSpectrum_rho(1000), gal->GetInjSpectrum_alpha(1000), in->cutoff_rig, true, true); // Standard EXTRA component particles
       else {
	 //sp = new TSpectrum(gal->GetCoordinates(), in, gal->GetInjSpectrum(uid), gal->GetInjLowSpectrum(uid), (A==0)); // Standard particles
	 sp = new TSpectrum(gal->GetCoordinates(), in, gal->GetInjSpectrum_rho(uid), gal->GetInjSpectrum_alpha(uid), in->cutoff_rig, (A==0), false); // Standard particles
       }
     }
     else {
       if (in->feedback >0) cout << "Particles coming from DM. Building spectrum..." << endl;
         // DM particles
         if (in->DMs == EWCorrections) sp = new TSpectrum(gal->GetCoordinates(), int(in->mx), in->dmmode, (A==0)*151 + (A==1)*154);
         else if (in->DMs == DarkSUSY || in->DMs == Delta) sp = new TSpectrum(gal->GetCoordinates(), in, (A==0)*151 + (A==1)*154);
         else if (in->DMs == SelfTable) {
	  if (in->feedback >0) cout << "Particles coming from DM. Getting datafile with inj spectrum specified by the user... " << endl;
	  if (A==0)  //Electrons and Positrons from DM annihilation/decay	
		 sp = new TSpectrum(gal->GetCoordinates(), in->MySelfTableDMel, 0);
	  if (A==1 && Z==-1) //Antiprotons from DM annihilation/decay	
		 sp = new TSpectrum(gal->GetCoordinates(), in->MySelfTableDMap, 0);
	  if (A==2 && Z==-1)	
		 sp = new TSpectrum(gal->GetCoordinates(), in->MySelfTableDMdbar, 0);
	 }
      }
   }
   
   eloss.push_back(new TIonizationLoss(gal->GetCoordinates(), gal->GetGas(), gal->GetTotalGas(), in, A, Z));
   eloss.push_back(new TCoulombLoss(gal->GetCoordinates(), gal->GetGas(), in, A, Z));
   
   if (in->feedback >0) cout << "[debug] elosses - Bremsstrahlung" << endl;
   if (A == 0) {
      eloss.push_back(new TBremsstrahlungLoss(gal->GetCoordinates(), gal->GetGas(), gal->GetTotalGas(), in));
      if (in->feedback >0) cout << "[debug] elosses - Synchrotron" << endl;
      eloss.push_back(new TSynchrotronLoss(gal->GetCoordinates(), gal->GetBField(), in));
      if (in->feedback >0) cout << "[debug] elosses - Inverse Compton" << endl;
      eloss.push_back(new TICSLoss(gal->GetCoordinates(), gal->GetISRF(), in));
   }
   
   if (in->feedback >0) cout << "[debug] elosses - Inelastic Cross Section" << endl;
   
   _fInXSec = new TInelasticCrossSection(gal->GetCoordinates(), in, uid, K_electron, xsecmodel); //modified
   
   if (in->feedback >0) cout << "End of constructor " << endl;
}
/**< The constructor normally used in DRAGON */

//modified
TParticle::TParticle(TParticle& part) :
A(part.A),
Z(part.Z),
uid(part.uid),
dimz(part.dimz),
dimE(part.dimE),
daughter(part.daughter),
lifetime(part.lifetime),
isExtra(part.isExtra),
isTPP(part.isTPP),
isDM(part.isDM),
decmode(part.decmode) {
   _fGalaxy = part._fGalaxy;
   _fDiff = part._fDiff;
   _fDpp = part._fDpp;
   _fInXSec = part._fInXSec;
   eloss = part.eloss;
   density = part.density;
   sp = part.sp;
   issec = part.issec;
   K_electron = part.K_electron; //modified
   in = part.in;
}

TParticle::~TParticle() {
   if (sp) delete sp;
   if (_fInXSec) delete _fInXSec;
   for (vector<TEnergyLoss*>::iterator i = eloss.begin(); i != eloss.end(); ++i) {
      if (*i) delete *i;
   }
   eloss.clear();
   density.clear();
}

//******************************************************************************************************************************************************

void TParticle::Evolve(vector<TParticle*> part, vector<TCREvolutorBasis*> crev_, TSpallationNetwork* spnet, vector<TXSecBase*> xsecmodel, bool isSecondIteration) {

   if (isSecondIteration)
     if (in->feedback >0) cout << "This is the SECOND ITERATION " << endl;
   
   TEnergyLoss el(*(eloss.front()));
   for (vector<TEnergyLoss*>::iterator i = eloss.begin()+1; i != eloss.end(); ++i) el += (*(*i));
   
   if(!in->ELOSS) el.Set_dpdt_Zero(); //fk 130701
   
   //modified
   if (in->feedback >0) cout << "We are building the secondary source term " << endl;
   const vector<double> secsource = ComputeSecondarySource(part, spnet, xsecmodel);
   
   if (in->DiffT == Isotropic) {
      for (vector<TCREvolutorBasis*>::iterator crev = crev_.begin(); crev != crev_.end(); ++crev) (*crev)->Run(density, density_previous,  _fInXSec, el.GetDpdt(), secsource, sp->GetSpectrum(), double(A), double(Z), lifetime, daughter, bool(issec), K_electron, isDM, isExtra); //modified
   }
   else {
      for (vector<TCREvolutorBasis*>::iterator crev = crev_.begin(); crev != crev_.end(); ++crev) (*crev)->Run(density, density_previous, _fInXSec, el.GetDpdt(), secsource, sp->GetSpectrum(), _fDiff, _fDpp, double(A), double(Z), lifetime, daughter, bool(issec), K_electron, isDM, isExtra, _fGalaxy); //modified
   }
   
   /*
   //
   // THIS PART IS ALREADY IN THE DESCRUCTOR and creates problems with the second iteration
   // 
   for (vector<TEnergyLoss*>::iterator i = eloss.begin(); i != eloss.end(); ++i) {
      if (*i) delete *i;
   }
   eloss.clear();
   
   if (_fDiff) delete _fDiff;
   if (_fDpp) delete _fDpp;*/
   
   return ;
}

//******************************************************************************************************************************************************


/**< The constructor normally used in DRAGON */
TParticle2D::TParticle2D(int A_ /**< Mass number */, int Z_ /**< Charge */, Galaxy* gal /**< A model for the galaxy */, Input* in /**< User input */, bool uid_prec, vector<TXSecBase*> xsecmodel, TNucleiList* l, int K_electron_, bool isDM, bool isextra,bool isTPP) :
TParticle(A_, Z_, gal, in, uid_prec, xsecmodel, l, K_electron_, isDM, isextra,isTPP),
dimr(gal->GetCoordinates()->GetDimR()) {

  density = vector<double>(dimr*dimz*dimE, 0.0);
}
/**< The constructor normally used in DRAGON */

//******************************************************************************************************************************************************

double TParticle2D::FindNormalization(const double& normrig, const double& normval) {
   
   if (!(uid == 1001 || uid == -1000)) {
      cerr << "Asking for normalization but no protons nor primary electrons found" << endl;
      exit(NOPROTELNORM);
   }
   
   vector<double> r = _fGalaxy->GetCoordinates()->GetR();
   unsigned int irsun = (unsigned int) ((in->robs-r.front())/(r.back()-r.front())*(double)(dimr-1));
   
   vector<double> z = _fGalaxy->GetCoordinates()->GetZ();
   vector<double> Ek = _fGalaxy->GetCoordinates()->GetEk();
   
   unsigned int izsun = (unsigned int) ((in->zobs-z.front())/(z.back()-z.front())*(double)(dimz-1));
 
   cout << "irsun = " << irsun << " izsun = " << izsun << endl;

   for (int ind = 0; ind < Ek.size(); ind++) 
     cout << "spectrum at irsun, izsun = " << density[TParticle::index(irsun,izsun,ind)] << endl;
  
   double r1 = (r[irsun+1]-in->robs)/(r[irsun+1]-r[irsun]);
   double r2 = (in->robs-r[irsun])/(r[irsun+1]-r[irsun]);
   
   int ilow = int(log(normrig/Ek[0])/log(Ek[1]/Ek[0]));

   cout << "energy index = " << ilow << endl;
   
   double valuelow = log10( (density[TParticle::index(irsun,izsun,ilow)]*r1+density[TParticle::index(irsun+1,izsun,ilow)]*r2) );
   double valuehigh = log10( (density[TParticle::index(irsun,izsun,ilow+1)]*r1+density[TParticle::index(irsun+1,izsun,ilow+1)]*r2) );

   cout << "valuelow and valuehigh " << valuelow << " " << valuehigh << endl;

   double e1 = log(Ek[ilow+1]/normrig)/log(Ek[ilow+1]/Ek[ilow]);
   double e2 = log(normrig/Ek[ilow])/log(Ek[ilow+1]/Ek[ilow]);
   
   double value = pow(10,valuelow*e1 + valuehigh*e2);
   double factor = normval/value;
   //cout << value << " " << inp->spect_norm << endl;
   //cout << factor << " " << factorint << " comparison of interpolation." << endl;
   
   //cout << "Normalization required in the xml: " << inp->spect_norm << endl;
   if (in->feedback >0) cout << "Reference value in arbitrary units  = " << value << "; Normalization factor = " << factor << endl << endl;
   
   return factor;
}

//******************************************************************************************************************************************************
/*
double TParticle2D::FindNormalization(TParticle* electrons, const double& normel) {
   //to be used for the extra component
   //
   if (!(uid == -1000)) {
      cerr << "Asking for normalization of extra component but no primary electrons found" << endl;
      exit(NOPROTELNORM);
   }
   
   cout << "2D normalization" << endl;
   
   vector<double> r = _fGalaxy->GetCoordinates()->GetR();
   unsigned int irsun = (unsigned int) ((robs-r.front())/(r.back()-r.front())*(double)(dimr-1));
   
   vector<double> z = _fGalaxy->GetCoordinates()->GetZ();
   vector<double> Ek = _fGalaxy->GetCoordinates()->GetEk();
   
   unsigned int izsun = (unsigned int) ((zobs-z.front())/(z.back()-z.front())*(double)(dimz-1));
   
   double r1 = (r[irsun+1]-robs)/(r[irsun+1]-r[irsun]);
   double r2 = (robs-r[irsun])/(r[irsun+1]-r[irsun]);
   
   int ilow = int(log(inp->sp_ref_rig_el_extra/Ek[0])/log(Ek[1]/Ek[0]));
   
   double valuelow = log10( (density[TParticle::index(irsun,izsun,ilow)]*r1+density[TParticle::index(irsun+1,izsun,ilow)]*r2) );
   double valuehigh = log10( (density[TParticle::index(irsun,izsun,ilow+1)]*r1+density[TParticle::index(irsun+1,izsun,ilow+1)]*r2) );
   
   double valuelow_std = 0.0;
   double valuehigh_std = 0.0;
   if (electrons) {
      valuelow_std = log10( (electrons->GetDensity(irsun,izsun,ilow)*r1+electrons->GetDensity(irsun+1,izsun,ilow)*r2) );
      valuehigh_std = log10( (electrons->GetDensity(irsun,izsun,ilow+1)*r1+electrons->GetDensity(irsun+1,izsun,ilow+1)*r2) );
   }
   
   double e1 = log(Ek[ilow+1]/inp->sp_ref_rig_el_extra)/log(Ek[ilow+1]/Ek[ilow]);
   double e2 = log(inp->sp_ref_rig_el_extra/Ek[ilow])/log(Ek[ilow+1]/Ek[ilow]);
   
   double value = pow(10,valuelow*e1 + valuehigh*e2);
   double value_std = (electrons) ? pow(10,valuelow_std*e1 + valuehigh_std*e2) : 0;
   double factor = (inp->spect_norm_el_extra - normel*value_std)/value;
   //cout << value << " " << inp->spect_norm << endl;
   //cout << factor << " " << factorint << " comparison of interpolation." << endl;
   cout << "Normalization required in the xml: " << inp->spect_norm_el_extra << endl;
   cout << "Reference value in arbitrary units  = " << value << "; Normalization factor = " << factor<< endl << endl;
   return factor;
}
*/
//******************************************************************************************************************************************************

vector<double> TParticle2D::ComputeSecondarySource(vector<TParticle*> part, TSpallationNetwork* spnet, vector<TXSecBase*> xsecmodel) {
   /**< Here also inverse beta-decay (K-electron capture) must be taken into account */
   
   vector<double> result(density.size(), 0.0);
   // If it is the only particle to be propagated, or if it is primary electrons, or if it is DM, do not add secondary contribution
   if (part.size() == 1 || (uid == -1000 && !issec) || isDM) return result;
   TGas* totalgas = _fGalaxy->GetTotalGas();
   TGrid* coord = _fGalaxy->GetCoordinates();
#undef DEBUG
#ifdef DEBUG
   vector<double> Ek = coord->GetEk();
   
#endif
   
   if (issec && fabs(uid) != 1000) { // Secondary protons or tertiary antiprotons or tertiary antideuterons: they are next to their primary brothers...

     vector<double> Ek(coord->GetEk());  
     
     if (uid == 1001 && issec){
	 cout << "Secondary protons!" << endl;

	 if ( in->spallationxsec == Fluka){//Fluka model
       
	   //pair<int,int> coupleppr(1001,1001);  // Secondary protons, from prim. protons
	   //pair<int,int> couplepHe(2004,1001);  // Secondary protons, from Helium	
	   //pair<int,int> coupleapap(-999, -999);  // Tertiary antiprotons, from sec. antiprotons
	   
	   //cout << "*** FLUKA Spallation cross section: ***" << endl;
	   
	   vector<TParticle*>::iterator iprim = part.end()-2; // Primary species
	   
	   vector<double> beta(coord->GetBeta());
           
	   for (int k = 0; k < dimr; ++k) {
	     for (int l = 0; l < dimz; ++l) {
	       int indspat = coord->index(k,l);
	       double gasdensity = totalgas->GetGas(indspat);
	       
	       for (int i = 0; i < dimE; ++i)  {
		 double betaigasdensity = beta[i]*gasdensity;
		 int ind = indspat*dimE+i;
		 vector<double> spall_spectrum( spnet->GetXSecApEl(uid, uid, i) );
		 // if (k==0 && l==0)
		 // 	for (int ie=0; ie<spall_spectrum.size(); ie++)
		 // 	  cout << i << " " << spall_spectrum[ie] << endl;
		 //if (k==0 && l==0) cout << "E= " <<  Ek[i] << " spall xsec from p = " << spall_spectrum[i] << endl;
		 for (int ii = i; ii < dimE; ++ii) result[ind] += (*iprim)->GetDensity(indspat*dimE+ii) * betaigasdensity * spall_spectrum[ii];
		 //std::cout << "\nsecondary protons, Fluka " << result[ind] << std::endl;		 
	       }
	     }
	   }
	 }
	 
	 else {  
	   cout << "DRAGON CROSEC Secondary protons" << endl;
	   vector<TParticle*>::iterator iprim = part.end()-2; // Primary species	   
	   vector<double> beta(coord->GetBeta());
	   
	   vector<double> spall_spectrum( spnet->GetXSec(uid, uid) );
	   
	   for (int k = 0; k < dimr; ++k) {
	     for (int l = 0; l < dimz; ++l) {
	       int indspat = coord->index(k,l);
	       double gasdensity = totalgas->GetGas(indspat);
	       
	       for (int i = 0; i < dimE; ++i)  {
		 double betaigasdensity = beta[i]*gasdensity;
		 int ind = indspat*dimE+i;
		 //if (k==0 && l==0) cout << "E= " <<  Ek[i] << " spall xsec from p = " << spall_spectrum[i] << endl;
		 for (int ii = i; ii < dimE; ++ii) result[ind] += (*iprim)->GetDensity(indspat*dimE+ii) * betaigasdensity * spall_spectrum[ii];
		 //std::cout << "\nsecondary protons, no Fluka" << result[ind] << std::endl;
	       }
	     }
	   }
	 }
	 /*       vector<double> beta(coord->GetBeta());
	 //loop over the parent particles. Only p and He are relevant
	 for (vector<TParticle*>::iterator ipart = part.begin(); ipart != part.end()-1; ++ipart) {
	 if ((*ipart)->GetUid() > 2004 || (*ipart)->GetUid() < 0) continue;
	 for (int i = 0; i < dimE; i++) {
	 vector<double> spall_spectrum( spnet->GetXSecApEl((*ipart)->GetUid(), uid, i) );
	 if (spall_spectrum.size() == dimE) {
	 for (int k = 0; k < dimr; k++) {
	 for (int l = 0; l < dimz; l++) {
	 int indspat = coord->index(k,l);
	 double gasdensity = totalgas->GetGas(indspat);
	 double betaigasdensity = beta[i]*gasdensity;
	 int ind = indspat*dimE+i;
	 // ENERGY INTEGRAL
	 for (int ii=0; ii < dimE; ii++) {
	 //result[ind] += (*ipart)->GetDensity(indspat*dimE+ii) * gasdensity * spall_spectrum[ii];
	 result[ind] += (*ipart)->GetDensity(indspat*dimE+ii) * betaigasdensity * spall_spectrum[ii];
	 if ((*ipart)->GetUid() == 1001 && uid==1001 && k==0 && l==0) 
	 cout << " E_parent= " <<  Ek[i] << " E_daughter= " <<  Ek[ii] << " spall xsec from p = " << spall_spectrum[ii] << endl;
	 
	 }
	 }
	 }
	 }
	 }
	 }*/
	 
	 //cout << "*** End of spallation cross section: ***" << endl;

     }//End Sec protons
   
     else if (uid == -999 && issec){
       cout << "Tertiary antiprotons" << endl;
       
       if ( in->apy == FlukaAp){//Fluka model
	 //pair<int,int> coupleapap(-999, -999);  // Tertiary antiprotons, from sec. antiprotons
	 cout << "*** FLUKA Tertiary Ap Spallation cross section: ***" << endl;
	 
	 vector<TParticle*>::iterator iprim = part.end()-2; // Primary species
	 
	 vector<double> beta(coord->GetBeta());
         
	 for (int k = 0; k < dimr; ++k) {
	   for (int l = 0; l < dimz; ++l) {
	     int indspat = coord->index(k,l);
	     double gasdensity = totalgas->GetGas(indspat);
	     
	     for (int i = 0; i < dimE; ++i)  {
	       double betaigasdensity = beta[i]*gasdensity;
	       int ind = indspat*dimE+i;
	       vector<double> spall_spectrum( spnet->GetXSecApEl(uid, uid, i) );
	       for (int ii = i; ii < dimE; ++ii) result[ind] += (*iprim)->GetDensity(indspat*dimE+ii) * betaigasdensity * spall_spectrum[ii];
	       //std::cout << "\ntertiary antiprotons, Fluka " << result[ind] << std::endl;
	       
	     }
	   }
	 }
       }
	 
       else {  
	 cout << "*** DRAGON approach for tertiary Ap cross sections! " << endl;
	 int a=0;
	 vector<TParticle*>::iterator iprim = part.end()-2; // Primary species
	 vector<double> beta(coord->GetBeta());
	 
	 vector<double> spall_spectrum( spnet->GetXSec(uid, uid) );
	 
	 for (int k = 0; k < dimr; ++k) {
	   for (int l = 0; l < dimz; ++l) {
	     int indspat = coord->index(k,l);
	     double gasdensity = totalgas->GetGas(indspat);
	     
	     for (int i = 0; i < dimE; ++i)  {
	       double betaigasdensity = beta[i]*gasdensity;
	       int ind = indspat*dimE+i;
	       //if (k==0 && l==0) cout << "E= " <<  Ek[i] << " spall xsec from p = " << spall_spectrum[i] << endl;
	       for (int ii = i; ii < dimE; ++ii){ result[ind] += (*iprim)->GetDensity(indspat*dimE+ii) * betaigasdensity * spall_spectrum[ii];}
	       //if (result[ind] > 0){
	       //  a = 1;
	       //  std::cout << "\ntertiary antiprotons, no Fluka " << spall_spectrum[ii] << " " << (*iprim)->GetDensity(indspat*dimE+ii) << std::endl;}}
	       //if (a == 1) 	       exit(-1);
	     }
	   }
	 }
       }
       
     }//End Tertiary Ap
     
   } // if (issec && uid != -1000)
   else if (uid > 1000) { // nuclei spallation
     
     //modified
     vector<double> gamma(coord->GetGamma());
     vector<double> beta(coord->GetBeta());
     vector<double> energy(coord->GetEk());
     
     //modified
     //if K_capture is present, the source term due to electron attachment/stripping is computed
     if (K_electron > 0) {
         
         double factor = 1.e-27*Clight;
         
         TParticle* naked_nucleus;
         if (part.size() > 1) naked_nucleus = part[part.size()-2];
         if (in->feedback >0) cout << "Computing source term for dressed nucleus " << endl;
         //cout << "source term coming from nucleus with A =  " << naked_nucleus->GetA() << "; Z = " << naked_nucleus->GetZ() << endl;
         
         for (int k = 0; k < dimr; ++k) {
            for (int l = 0; l < dimz; ++l) {
               int indspat = coord->index(k,l);
               double Gasdensity = totalgas->GetGas(indspat);
               for (int i = 0; i < dimE; ++i) {
                  int ind = indspat*dimE+i;
                  double attach_H=0., attach_He=0., strip_H=0., strip_He=0.;
                  xsecmodel[0]->Kcapture_cs(energy[i],naked_nucleus->GetZ(),1,&attach_H,&strip_H);
                  xsecmodel[0]->Kcapture_cs(energy[i],naked_nucleus->GetZ(),2,&attach_He,&strip_He);
                  result[ind] = Gasdensity * (attach_H + He_abundance*attach_He) * factor * beta[i] * naked_nucleus->GetDensity(ind);
               }
            }
         }
         
      }
    
      //if K_capture is present and K_electron == 1, the source term ends here: it only contains contribution from the corresponding naked nucleus
      
      else { // if K_capture is present and K_electron == 0, the source term must contain also the spallation and decay contribution from heavier nuclei
         
         for (vector<TParticle*>::iterator ipart = part.begin(); ipart != part.end()-1; ++ipart) { // // //
	   
	   if ((*ipart)->GetDaughter() == uid && in->feedback >0)
	     cout << "---> Nucleus " << ((*ipart))->GetUid() << " decays into current particle " << uid << endl; 	
	   
            vector<double> spall_spectrum( spnet->GetXSec( (*ipart)->GetUid(), uid ) );
            if (spall_spectrum.size() == dimE || (*ipart)->GetDaughter() == uid) {
               
               double Afactor = double((*ipart)->GetA())/double(A); // To convert between kinetic energy per nucleon (which is approx. constant) to momentum
               for (int k = 0; k < dimr; ++k) {
                  for (int l = 0; l < dimz; ++l) {
                     int indspat = coord->index(k,l);
                     double Afactorgasdensity = Afactor*totalgas->GetGas(indspat);
                     
                     for (int i = 0; i < dimE; ++i) {
                        int ind = indspat*dimE+i;
                        if (spall_spectrum.size() == dimE) result[ind] += Afactorgasdensity * spall_spectrum[i] * (*ipart)->GetDensity(ind); // spallation
			  
                        if ((*ipart)->GetDaughter() == uid) result[ind] += (*ipart)->GetDensity(ind)/(*ipart)->GetLifetime()/gamma[i]; // decay
                     }
                  }
               }
            }
         }
         
      } // if (K_electron == 0)
   } // else if (uid > 1000)
   else { // Secondary antiprotons, electrons and positrons, antideuterons, positrons from TPP
    	
      
      if (!isTPP) {
	vector<double> beta(coord->GetBeta());
	
	// ofstream outfile("prova_Pohl.dat", ios::out);
	for (vector<TParticle*>::iterator ipart = part.begin(); ipart != part.end()-1; ++ipart) {
	  if ((*ipart)->GetUid() > 2004 || (*ipart)->GetUid() < 0) continue;
	  for (int i = 0; i < dimE; i++) {
	    vector<double> spall_spectrum( spnet->GetXSecApEl((*ipart)->GetUid(), uid, i) );
	    
	    if (spall_spectrum.size() == dimE) {
	      for (int k = 0; k < dimr; k++) {
		for (int l = 0; l < dimz; l++) {
		  int indspat = coord->index(k,l);
		  double betaigasdensity = totalgas->GetGas(indspat) * beta[i];
		  int ind = indspat*dimE+i;
		  // ENERGY INTEGRAL
		  for (int ii=0; ii < dimE; ii++) {	result[ind] += (*ipart)->GetDensity(indspat*dimE+ii) * betaigasdensity * spall_spectrum[ii];
#ifdef DEBUG
		    if ((*ipart)->GetUid() == 1001 && uid==1000) outfile << Ek[i] << " " << Ek[ii] << " " << spall_spectrum[ii] << endl;
#endif
		  }
		  //std::cout << result[ind] << " for nuclei " << (*ipart)->GetUid() << " to " << uid << std::endl;
		}
	      }
	    }
	    //else std::cout << spall_spectrum.size() << std::endl; exit(-1);
	  }
	}
	//            outfile.close();
      }
      else {	//TPP: Secondary positrons coming from TPP -- implemented by D.Gaggero -- november 2011
         //        if (uid == 1000 && issec) { //TPP: Secondary positrons coming from TPP
         if (in->feedback >0) cout << "Building source term for TPP positrons..." << endl;
         
         vector<double> energy(coord->GetEk());
         vector<double> R_vec(coord->GetR());
         vector<double> z_vec(coord->GetZ());
         double DlogE = coord->GetDeltaE();
         
         TISRF* Gal_ISRF = _fGalaxy->GetISRF();
         vector <double> nu_vector = Gal_ISRF->GetNuArray();
         double DlogNu = Gal_ISRF->GetDnu();
         int dimNu = nu_vector.size();
         double hPlanck = 4.135667e-15; // eV/Hz
         vector <double> ISRF_vector = Gal_ISRF->GetISRF(); // nu U_nu; units: Hz eV cm^-3 Hz^-1
         
         for (vector<TParticle*>::iterator ipart=part.begin(); ipart!=part.end()-1; ++ipart) {
            
            if ((*ipart)->GetUid() == -1000) { //TPP positrons only come from (primary and secondary) electrons
               
               if (in->feedback >0) cout << "Calculating TPP cross section..." << endl;
               
               vector<double> TPP_cross_section(spnet->GetXSecTPP(nu_vector));
               
               if (in->feedback >0) cout << "TPP cross section calculated" << endl;
               
               vector<double> parent_particle_density (dimE*dimr*dimz);
               for (int k = 0; k < dimr; k++) {
                  for (int l = 0; l < dimz; l++) {
                     for (int ii = 0; ii < dimE; ii++) {
                        int indspat = coord->index(k,l);
                        parent_particle_density[indspat*dimE+ii] = (*ipart)->GetDensity(indspat*dimE+ii);
                     }
                  }
               }
               
               for (int i = 0; i < dimE; i++) {
                  
                  for (int k = 0; k < dimr; k++) {
                     for (int l = 0; l < dimz; l++) {
                        
                        int indspat = coord->index(k,l);
                        int ind = indspat*dimE+i;
                        
                        for (int inu = 0; inu<dimNu; inu++) { // integral over frequency_vector of ISRF
                           
                           double photon_density = ISRF_vector[(inu*dimr+k)*dimz+l] / ((hPlanck * nu_vector[inu])*nu_vector[inu]);  // rho_nu; units: cm^-3 Hz^-1
                           
                           double energy_integral = 0.;
                           
                           for (int ii=0; ii < dimE; ii++)	{  // integral over energy of parent particle
                              
         		      double contrib = parent_particle_density[indspat*dimE+ii] * TPP_cross_section[(i*dimNu + inu)*dimE +ii] * DlogE;   // cm^-3 * cm^3/Myr DeltalogE:adimensional
                              
                              energy_integral += contrib;
                              
                           }
                           
                           double contrib_2 = photon_density * energy_integral * nu_vector[inu] * DlogNu; // cm^-3 * 1/Myr Deltalognu: adimensional
                           result[ind] += contrib_2;
                           
                        }
                     } //for l = 0,dimz
                  } // for k = 0,dimr
               } //for i=0,dimE
               
            }
            
            if (in->feedback >0) cout << "Building source term for TPP positrons succeeded." << endl;
         }
         
      }
    	
    	
   } // else
   
   return result;
}


//******************************************************************************************************************************************************

void TParticle2D::Print(fitsfile* output_ptr, double norm) { /**< Print all the information relevant to that nucleus: charge, mass, source abundance, injection spectrum and propagated density. */
   int status = 0;
   
   const long naxis = 3;
   long size_axes[naxis] = {dimE,dimr,dimz};
   long nelements = size_axes[0]*size_axes[1]*size_axes[2];
   long fpixel = 1;
   int bitpix = FLOAT_IMG;
   if (fits_create_img(output_ptr, bitpix, naxis, size_axes, &status)) fits_report_error(stderr, status);
   
   double sab = _fGalaxy->GetSourceAbundance(uid);
    
   if (fits_write_key(output_ptr, TINT,    (char*) "Z_",      &Z,                 NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "A",       &A,                 NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "Sec",     &issec,             NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "DM",      &isDM,              NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "EXTRA",   &isExtra,           NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "TPP",     &isTPP,             NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "S_Ab",    &sab,               NULL, &status)) fits_report_error(stderr, status);
   int counter = 0;
   float* array = new float[nelements];
   double weight = 1;
   if (A != 0) weight = double(A);
   for (int l = 0; l < dimz; ++l) {
      for (int k = 0; k < dimr; ++k) {
         for (int j = 0; j < dimE; ++j) {
            array[counter] = float(weight*density[TParticle::index(k,l,j)]*norm);
            counter++;
         }
      }
   }
   
   if (fits_write_img(output_ptr, TFLOAT, fpixel, nelements, array, &status)) fits_report_error(stderr, status);
   
   delete [] array;
   
   return ;
}

//******************************************************************************************************************************************************

void TParticle2D::PrintSpectrum(fitsfile* output_ptr, double norm) { /**< Print all the information relevant to that nucleus: charge, mass, source abundance, injection spectrum and propagated spectrum at Sun position. */
   int status = 0;
   const long naxis = 1;
   long size_axes[naxis] = {dimE};
   long nelements = size_axes[0];
   long fpixel = 1;
   int bitpix = FLOAT_IMG;
   if (fits_create_img(output_ptr, bitpix, naxis, size_axes, &status)) fits_report_error(stderr, status);
   
   double sab = (!issec) ? _fGalaxy->GetSourceAbundance(uid) : 0.0;
   /*double injind_rho_0 = _fGalaxy->GetInjSpectrum_rho_0(uid);
   double injind_rho_1 = _fGalaxy->GetInjSpectrum_rho_1(uid);
   double injind_rho_2 = _fGalaxy->GetInjSpectrum_rho_2(uid);
   double injind_alpha_0 = _fGalaxy->GetInjSpectrum_alpha_0(uid);
   double injind_alpha_1 = _fGalaxy->GetInjSpectrum_alpha_1(uid);
   double injind_alpha_2 = _fGalaxy->GetInjSpectrum_alpha_2(uid);
   double injind_alpha_3 = _fGalaxy->GetInjSpectrum_alpha_3(uid);
   */
   if (fits_write_key(output_ptr, TINT,    (char*) "Z_",      &Z,                 NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "A",       &A,                 NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "Sec",     &issec,             NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "DM",     &isDM,             NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "EXTRA",     &isExtra,             NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "TPP",     &isTPP,             NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "S_Ab",    &sab,               NULL, &status)) fits_report_error(stderr, status);
   /*if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_rho_0", &injind_rho_0,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_rho_1", &injind_rho_1,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_rho_2", &injind_rho_2,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_alpha_0", &injind_alpha_0,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_alpha_1", &injind_alpha_1,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_alpha_2", &injind_alpha_2,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_alpha_3", &injind_alpha_3,            NULL, &status)) fits_report_error(stderr, status);
   */
   float* array = new float[nelements];
   double weight = 1;
   
   vector<double> r = _fGalaxy->GetCoordinates()->GetR();
   vector<double> z = _fGalaxy->GetCoordinates()->GetZ();
   vector<double> Ek =_fGalaxy->GetCoordinates()->GetEk();
   
   unsigned int irsun = (unsigned int) ((in->robs-r.front())/(r.back()-r.front())*(double)(dimr-1));
   unsigned int izsun = (unsigned int) ((in->zobs-z.front())/(2.0*z.back())*(double)(dimz-1));
   
   double r1 = (r[irsun+1]-in->robs)/(r[irsun+1]-r[irsun]);
   double r2 = (in->robs-r[irsun])/(r[irsun+1]-r[irsun]);
   if (A != 0) weight = double(A);

   if (isDM && in->feedback >0)
     cout << "**** particle originating from DARK MATTER *****" << endl;
   
   if (in->feedback >0) cout << "************************************************" << endl;
   if (in->feedback >0) cout << "***Writing spectrum at Solar System position... " << endl;
   if (in->feedback >0) cout << "************************************************" << endl;
   if (in->feedback >0) cout << "A = " <<  A << " Z = " << Z << " issec = " << issec << endl;
   if (in->feedback >0) cout << "************************************************" << endl;
   
   if (in->feedback >0) cout << "Energy -- Normalized flux " << endl;
   for (int j = 0; j < dimE; ++j){
     array[j] = float(weight*norm*( density[TParticle::index(irsun,izsun,j)]*r1 + density[TParticle::index(irsun+1,izsun,j)]*r2 ));
     if (in->feedback >0) cout <<  Ek[j] << " " <<  array[j] << endl;
   }
   if (in->feedback >0) cout << endl;
   
   if (fits_write_img(output_ptr, TFLOAT, fpixel, nelements, array, &status)) fits_report_error(stderr, status);
   delete [] array;
   return;
}

vector<double> TParticle2D::GetSpectrumAtSunPosition() {

  vector<double> result(dimE);
  vector<double> r =  _fGalaxy->GetCoordinates()->GetR();
  vector<double> z =  _fGalaxy->GetCoordinates()->GetZ();
  vector<double> Ek = _fGalaxy->GetCoordinates()->GetEk();
  
  unsigned int irsun = (unsigned int) ((in->robs-r.front())/(r.back()-r.front())*(double)(dimr-1));
  unsigned int izsun = (unsigned int) ((in->zobs-z.front())/(2.0*z.back())*(double)(dimz-1));
  
  double r1 = (r[irsun+1]-in->robs)/(r[irsun+1]-r[irsun]);
  double r2 = (in->robs-r[irsun])/(r[irsun+1]-r[irsun]);
  
  double weight = 1;
  if (A != 0) weight = double(A);
  
  for (int j = 0; j < dimE; ++j)
    result[j] = weight*(density[TParticle::index(irsun,izsun,j)]*r1 + density[TParticle::index(irsun+1,izsun,j)]*r2);
  return result;
}

double TParticle2D::GetFluxAtSunPosition(int j) {

   vector<double> r =  _fGalaxy->GetCoordinates()->GetX();
   vector<double> z =  _fGalaxy->GetCoordinates()->GetZ();
   vector<double> Ek = _fGalaxy->GetCoordinates()->GetEk();
   unsigned int irsun = (unsigned int) ((in->robs-r.front())/(r.back()-r.front())*(double)(dimr-1));
   unsigned int izsun = (unsigned int) ((in->zobs-z.front())/(2.0*z.back())*(double)(dimz-1));
   double r1 = (r[irsun+1]-in->robs)/(r[irsun+1]-r[irsun]);
   double r2 = (in->robs-r[irsun])/(r[irsun+1]-r[irsun]);
   double weight = 1;
   if (A != 0) weight = double(A);
   return density[TParticle::index(irsun,izsun,j)]*r1 + density[TParticle::index(irsun+1,izsun,j)]*r2;

}


//******************************************************************************************************************************************************
/*
std::vector <double> TParticle2D::GetSpectra(double norm,string run_id,int writeflag,double *eng) {
   
   
   int status = 0;
   const long naxis = 1;
   long size_axes[naxis] = {dimE};
   long nelements = size_axes[0];
   long fpixel = 1;
   int bitpix = FLOAT_IMG;
   
   double sab = (!issec) ? _fGalaxy->GetSourceAbundance(uid) : 0.0;
   
   double* array = new double[nelements];
   double weight = 1;
   
   vector<double> r = _fGalaxy->GetCoordinates()->GetR();
   unsigned int irsun = (unsigned int) ((robs-r.front())/(r.back()-r.front())*(double)(dimr-1));
   unsigned int izsun = (unsigned int) ((zobs-_fGalaxy->GetCoordinates()->GetZ().front())/(2.0*_fGalaxy->GetCoordinates()->GetZ().back())*(double)(dimz-1));  //original code
   
   
   double r1 = (r[irsun+1]-robs)/(r[irsun+1]-r[irsun]);
   double r2 = (robs-r[irsun])/(r[irsun+1]-r[irsun]);
   if (A != 0) weight = double(A);
   for (int j = 0; j < dimE; ++j) array[j] = double(weight*norm*( density[index(irsun,izsun,j)]*r1 + density[index(irsun+1,izsun,j)]*r2 ));
   
	//for (int j = 0; j < dimE; ++j)	cout << r1 << " " << r2 << " " << index(irsun,izsun,j) << " " << density[index(irsun,izsun,j)] << " " << weight << endl;
   
   if(writeflag){
      
      bool write = true;
      
      ofstream file;
      char buff_save[1000];
      
      if(Z==-1 && A==0 && issec==0) sprintf(buff_save,"./ASCII_spectra/%s/electron_prim_density.dat",run_id.c_str());
      else if(Z==-1 && A==0 &&issec==1)  sprintf(buff_save,"./ASCII_spectra/%s/electron_sec_density.dat",run_id.c_str());
      else if(Z==1 && A==1 && issec==0)  sprintf(buff_save,"./ASCII_spectra/%s/proton_prim_density.dat",run_id.c_str());
      else if(Z==1 && A==1 && issec==1)  sprintf(buff_save,"./ASCII_spectra/%s/proton_sec_density.dat",run_id.c_str());
      else write = false;
      
      if(write)
      {
         file.open(buff_save);
         for(int i=0;i<20;i++) for(int j=0;j<dimE;j++) file << i << " " << j << " " << eng[j] << " " << double(weight*norm*( density[index(i,izsun,j)]*r1 + density[index(i+1,izsun,j)]*r2 )) << endl;
         file.close();
      }
   }
   
   vector <double> spectrum;
   
   
   if ((Z==1 && A==1)&& (issec==0)) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);	// prim proton
   if ((Z==1 && A==1)&& (issec==1)) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);	// sec proton
   if ((Z==-1 && A==1) && (issec==0)) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);	// sec antiproton
   if ((Z==-1 && A==1) && (issec==1)) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);	// ter antiproton
   if ((Z==-1 && A==0) && (issec==0)) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);	//prim electron
   if ((Z==-1 && A==0) && (issec==1)) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);	// sec electron
   if (Z==1 && A==0) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			// positron
   if (Z==6 && A==12) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Carbon 12
   if (Z==6 && A==13) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Carbon 13
   if (Z==6 && A==14) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Carbon 14
   if (Z==5 && A==10) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Boron 10
   if (Z==5 && A==11) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Carbon 11
   if (Z==4 && A==7) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Be_7
   if (Z==4 && A==8) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Be_8
   if (Z==4 && A==9) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Be_9
   if (Z==4 && A==10) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Be_10
   if (Z==4 && A==11) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Be_11
   if (Z==13 && A==26) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);		//Al_26
   if (Z==13 && A==27) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);	 	//Al_27
   
   
   
   delete [] array;
   return spectrum;
}
*/
//******************************************************************************************************************************************************

//******************************************************************************************************************************************************
//******************************************************************************************************************************************************

TParticle3D::TParticle3D(int A_ /**< Mass number */, int Z_ /**< Charge */, Galaxy* gal /**< A model for the galaxy */, Input* in /**< User input */, bool uid_prec, vector<TXSecBase*> xsecmodel, TNucleiList* l, int K_electron_, bool isDM, bool isextra, bool isTPP) :
TParticle(A_, Z_, gal, in, uid_prec, xsecmodel, l, K_electron_, isDM, isextra, isTPP),
dimx(gal->GetCoordinates()->GetDimX()),
dimy(gal->GetCoordinates()->GetDimY()) {
   
   cout << "Welcome to the Particle 3D constructor! " << endl;
   
   int dimz = (gal->GetCoordinates()->GetDimZ());

   vector<double> x = (gal->GetCoordinates()->GetX());
   vector<double> y = (gal->GetCoordinates()->GetY());
   vector<double> z = (gal->GetCoordinates()->GetZ());
   
   //MW130620: What are these deltas used for? -- for now, set them equidistantially
   double deltax = ( (x.back() - x.front()) / (dimx-1) );
   double deltay = ( (y.back() - y.front()) / (dimy-1) );
   double deltaz = ( (z.back() - z.front()) / (dimz-1) );
   
   cout << "Initializing the density vector... " << endl;
   
   density = vector<double>(dimx*dimy*dimz*dimE, 0.0);
   density_previous = vector<double>(dimx*dimy*dimz*dimE, 0.0);
   
   if (gal->GetTestMode() == true) {
      
      cout << "Initializing to a gaussian... " << endl;
      
      //double delta_ = deltax * 2.;
      double delta_ = deltax/100.;
      cout << "Amplitude of initial gaussian = " << delta_ << endl;
      
      double init_value = 0.;
      
      double x1 = 0.;
      double x2 = 0.;
      double norm1 = 1.;
      double norm2 = 0.;
      //double x3 = -10;
      
      for (int ix=0; ix<dimx; ix++) {
         for (int iy=0; iy<dimy; iy++) {
            for (int iz=0; iz<dimz; iz++) {
               for (int ip=0; ip<dimE; ip++) {
                  
                  //				double x1 = -20.;
                  //				double x2 =  20.;
                  
                  if ( gal->IsSourceMoving() == false) {
                     
                     init_value = norm1 * (1./(sqrt(3.14)*delta_)) * exp( -((x[ix]-x1)*(x[ix]-x1))/(delta_*delta_) - (y[iy]*y[iy])/(delta_*delta_) - (z[iz]*z[iz])/(delta_*delta_) );
                     
                     
                     init_value += norm2 * (1./(sqrt(3.14)*delta_)) * exp( -((x[ix]-x2)*(x[ix]-x2))/(delta_*delta_) - (y[iy]*y[iy])/(delta_*delta_) - (z[iz]*z[iz])/(delta_*delta_) );
                     
                  }
                  else {
                     
                     init_value = 0.; // (1./(sqrt(3.14)*delta_)) * exp( -((x[ix]-x3)*(x[ix]-x3))/(delta_*delta_) - (y[iy]*y[iy])/(delta_*delta_) - (z[iz]*z[iz])/(delta_*delta_) );
                     
                  }
                  
                  // if (ip == 0 && iz == dimz/2 && init_value > 1.e-1)
                  //	cout << " x = " << x[ix] << " y = " << y[iy] << " z = " << z[iz] << " value = " << init_value  << endl;
                  
                  density[index(ix,iy,iz,ip)] = init_value;  //initialization ot a Gaussian centered on the origin in order to compare with the analytical solution
               }
            }
         }
      }
      
   }
   
}

//******************************************************************************************************************************************************

double TParticle3D::FindNormalization(const double& sp_ref_rig_norm, const double& spect_norm) {
   
   if (!(uid == 1001 || uid == -1000)) {
      cerr << "Asking for normalization but no protons nor primary electrons found" << endl;
      exit(NOPROTELNORM);
   }
   
   cout << "3D normalization " << endl;
   cout << "Observer position: " << in->xobs << " " << in->yobs << " " << in->zobs << endl;	
   
   // Matze: rewrite stuff for nonequidistant grid
   
   vector <double> x = _fGalaxy->GetCoordinates()->GetX();
   vector <double> y = _fGalaxy->GetCoordinates()->GetY();
   vector <double> z = _fGalaxy->GetCoordinates()->GetZ();
   vector <double> Ek = _fGalaxy->GetCoordinates()->GetEk();
   
   unsigned int ixsun = 0;
   while(x[++ixsun] <= in->xobs){}
   ixsun--;
   unsigned int iysun = 0;
   while(y[++iysun] <= in->yobs){}
   iysun--;
   unsigned int izsun = 0;
   while(z[++izsun] <= in->zobs){}
   izsun--;
   
   double rx1 = (x[ixsun+1]-in->xobs)/(x[ixsun+1]-x[ixsun]);
   double ry1 = (y[iysun+1]-in->yobs)/(y[iysun+1]-y[iysun]);
   double rz1 = (z[izsun+1]-in->zobs)/(z[izsun+1]-z[izsun]);
   double rx2 = 1-rx1;
   double ry2 = 1-ry1;
   double rz2 = 1-rz1;
   cout << rx1 << " = rx1; rx2 = " << rx2 << endl;
   cout << ry1 << " = ry1; ry2 = " << ry2 << endl;
   cout << rz1 << " = rz1; rz2 = " << rz2 << endl;
   
   //double sp_ref_rig_norm, spect_norm;
   //sp_ref_rig_norm = (uid == 1001) ? inp->sp_ref_rig_norm : inp->sp_ref_rig_el;
   //spect_norm      = (uid == 1001) ? inp->spect_norm      : inp->spect_norm_el;
   
   int ilow = int(log(sp_ref_rig_norm/Ek[0])/log(Ek[1]/Ek[0]));
   cout << "ilow, i.e. energy index  " << ilow << endl;
   cout << "Ek[ilow] = " << Ek[ilow] << endl;
   
   cout << "Sun coordinates: " <<  ixsun << " " << iysun << " " << izsun << endl;

   cout << "+++ Propagated spectrum at Sun +++" << endl;
   for (int i=0; i < Ek.size(); i++)
	cout << density[index(ixsun,  iysun,  izsun,  i)] <<  " " ;
   cout << "+++ End of propagated spectrum +++" << endl;	
   
   double valuelow = log10( (density[index(ixsun,  iysun,  izsun,  ilow)]*rx1*ry1
                             + density[index(ixsun,  iysun+1,izsun,  ilow)]*rx1*ry2
                             + density[index(ixsun+1,iysun,  izsun,  ilow)]*rx2*ry1
                             + density[index(ixsun+1,iysun+1,izsun,  ilow)]*rx2*ry2)*rz1
                           +(density[index(ixsun,  iysun,  izsun+1,ilow)]*rx1*ry1
                             + density[index(ixsun,  iysun+1,izsun+1,ilow)]*rx1*ry2
                             + density[index(ixsun+1,iysun,  izsun+1,ilow)]*rx2*ry1
                             + density[index(ixsun+1,iysun+1,izsun+1,ilow)]*rx2*ry2)*rz2);
   cout << "10^valuelow i.e. density at sun pos., ilow energy = " << pow(10,valuelow) << endl;
   double valuehigh = log10((density[index(ixsun,  iysun,  izsun  ,ilow+1)]*rx1*ry1
                             + density[index(ixsun,  iysun+1,izsun  ,ilow+1)]*rx1*ry2
                             + density[index(ixsun+1,iysun,  izsun  ,ilow+1)]*rx2*ry1
                             + density[index(ixsun+1,iysun+1,izsun  ,ilow+1)]*rx2*ry2)*rz1
                            +(density[index(ixsun,  iysun,  izsun+1,ilow+1)]*rx1*ry1
                              + density[index(ixsun,  iysun+1,izsun+1,ilow+1)]*rx1*ry2
                              + density[index(ixsun+1,iysun,  izsun+1,ilow+1)]*rx2*ry1
                              + density[index(ixsun+1,iysun+1,izsun+1,ilow+1)]*rx2*ry2)*rz2);
   cout << "10^valuehigh i.e. density at sun pos., ihigh energy = " << pow(10,valuehigh) << endl;
   
   double e1 = log(Ek[ilow+1]/sp_ref_rig_norm)/log(Ek[ilow+1]/Ek[ilow]);
   double e2 = log(sp_ref_rig_norm/Ek[ilow])/log(Ek[ilow+1]/Ek[ilow]);
   cout << e1 << " = e1; e2 = " << e2 << endl;

   double value = pow(10,valuelow*e1 + valuehigh*e2);
   cout << "10^(valuelow*e1 + valuehigh*e2) = " << value << endl;
   
   double factor = spect_norm/value;
   
   cout << "normalization required in the xml: " << spect_norm << endl;
   cout << "interpolated value at ixsun, ixsun+1, iysun, iysun+1, ilow, ilow+1: " << value << endl;
   cout << "normalization factor = norm/value: " << factor << endl;

   cout << "+++ Propagated NORMALIZED spectrum at Sun +++" << endl;
   for (int i=0; i < Ek.size(); i++)
	cout << density[index(ixsun,  iysun,  izsun,  i)]*factor <<  " " ;
   cout << "+++ End of propagated spectrum +++" << endl;
   int dimx = x.size();
   int dimy = y.size();
   int dimz = z.size();	
   cout << "+++ Propagated NORMALIZED spectrum at the Gal.center +++" << endl;
   for (int i=0; i < Ek.size(); i++)
	cout << density[index(dimx/2+1,  dimy/2+1,  izsun,  i)]*factor <<  " " ;
   cout << "+++ End of propagated spectrum +++" << endl;	
   cout << "+++ Propagated NORMALIZED spectrum at intermediate latitude +++" << endl;
   for (int i=0; i < Ek.size(); i++)
	cout << density[index(ixsun-1,  iysun-1,  (dimz+izsun)/2,  i)]*factor <<  " " ;
   cout << "+++ End of propagated spectrum +++" << endl;	   
   return factor;
}

/*
double TParticle3D::FindNormalization(TParticle* electrons, const double& normel) {
   //To be used for the extra component
   
   if (!(uid == -1000)) {
      cerr << "Asking for normalization of the extra component but no primary electrons found" << endl;
      exit(NOPROTELNORM);
   }
   
   // Matze: rewrite what needed to be revised,
   // 10/12/12: copied that from above
   // To Be Improved!! Now it is only a copy of the previous routine!!
   
   vector <double> x = _fGalaxy->GetCoordinates()->GetX();
   vector <double> y = _fGalaxy->GetCoordinates()->GetY();
   vector<double> z = _fGalaxy->GetCoordinates()->GetZ();
   vector<double> Ek = _fGalaxy->GetCoordinates()->GetEk();
   unsigned int ixsun = 0;
   unsigned int iysun = 0;
   unsigned int izsun = 0;
   while(x[++ixsun] <= xobs){}
   ixsun--;
   while(y[++iysun] <= yobs){}
   iysun--;
   while(z[++izsun] <= zobs){}
   izsun--;
   
   double rx1 = (x[ixsun+1]-xobs)/(x[ixsun+1]-x[ixsun]);
   double ry1 = (y[iysun+1]-yobs)/(y[iysun+1]-y[iysun]);
   double rz1 = (z[izsun+1]-zobs)/(z[izsun+1]-z[izsun]);
   double rx2 = 1-rx1;
   double ry2 = 1-ry1;
   double rz2 = 1-rz1;
   cout << rx1 << " = rx1; rx2 = " << rx2 << endl;
   cout << ry1 << " = ry1; ry2 = " << ry2 << endl;
   cout << rz1 << " = rz1; rz2 = " << rz2 << endl;
   
   double sp_ref_rig_norm, spect_norm;
   sp_ref_rig_norm = (uid == 1001) ? inp->sp_ref_rig_norm : inp->sp_ref_rig_el_extra;
   spect_norm      = (uid == 1001) ? inp->spect_norm      : inp->spect_norm_el_extra;
   
   int ilow = int(log(sp_ref_rig_norm/Ek[0])/log(Ek[1]/Ek[0]));
   cout << "ilow, i.e. energy index  " << ilow << endl;
   
   // To be checked
   cout << "sun coordinates (for electrons): " <<  ixsun << " " << iysun << " " << izsun << endl;
   double valuelow = log10( (density[index(ixsun,iysun,izsun,ilow)]*rx1*ry1
                             + density[index(ixsun,iysun+1,izsun,ilow)]*rx1*ry2
                             + density[index(ixsun+1,iysun,izsun,ilow)]*rx2*ry1
                             + density[index(ixsun+1,izsun+1,izsun,ilow)]*rx2*ry2)*rz1
                           + (density[index(ixsun,iysun,izsun+1,ilow)]*rx1*ry1
                              + density[index(ixsun,iysun+1,izsun+1,ilow)]*rx1*ry2
                              + density[index(ixsun+1,iysun,izsun+1,ilow)]*rx2*ry1
                              + density[index(ixsun+1,izsun+1,izsun+1,ilow)]*rx2*ry2)*rz2);
   cout << "valuelow = " << pow(10,valuelow) << endl;
   double valuehigh = log10( (density[index(ixsun,iysun,izsun,ilow+1)]*rx1*ry1
                              + density[index(ixsun,iysun+1,izsun,ilow+1)]*rx1*ry2
                              + density[index(ixsun+1,iysun,izsun,ilow+1)]*rx2*ry1
                              + density[index(ixsun+1,izsun+1,izsun,ilow+1)]*rx2*ry2)*rz1
                            + (density[index(ixsun,iysun,izsun+1,ilow+1)]*rx1*ry1
                               + density[index(ixsun,iysun+1,izsun+1,ilow+1)]*rx1*ry2
                               + density[index(ixsun+1,iysun,izsun+1,ilow+1)]*rx2*ry1
                               + density[index(ixsun+1,izsun+1,izsun+1,ilow+1)]*rx2*ry2)*rz2);
   cout << "valuehigh = " << pow(10,valuehigh) << endl;
   
   double e1 = log(Ek[ilow+1]/sp_ref_rig_norm)/log(Ek[ilow+1]/Ek[ilow]);
   double e2 = log(sp_ref_rig_norm/Ek[ilow])/log(Ek[ilow+1]/Ek[ilow]);
   
   cout << e1 << " = e1; e2 = " << e2 << endl;
   
   double value = pow(10,valuelow*e1 + valuehigh*e2);
   double factor = spect_norm/value;
   
   cout << "normalization required in the xml: " << spect_norm << endl;
   cout << "interpolated value at ixsun, ixsun+1, iysun, iysun+1, ilow, ilow+1: " << value << endl;
   cout << "normalization factor = norm/value: " << factor << endl << endl;
   
   return factor;
}
*/
//******************************************************************************************************************************************************

vector<double> TParticle3D::ComputeSecondarySource(vector<TParticle*> part, TSpallationNetwork* spnet, vector<TXSecBase*> xsecmodel) {
   /**< Here also inverse beta-decay (K-electron capture) must be taken into account */
   
   cout << "We are starting the 3D computation of the secondary source " << endl;
   
   vector<double> result(density.size(), 0.0);
   
   // If it is the only particle to be propagated, or if it is primary electrons, do not add secondary contribution
   if (part.size() == 1 || (uid == -1000 && !issec) || isDM || isExtra) return result;
   
   TGas* totalgas = _fGalaxy->GetTotalGas();
   TGrid* coord = _fGalaxy->GetCoordinates();
   
   if (issec && fabs(uid) != 1000) { // Secondary protons or tertiary antiprotons or tertiary antideuterons: they are next to their primary brothers...
 

     if (uid == 1001 && issec){
       if (in->spallationxsec != Fluka)   { 
	 
	 cout << "*** DRAGON procedure for secondary protons! " << endl;
	 
	 
	 vector<TParticle*>::iterator iprim = part.end()-2; // Primary species
	 
	 vector<double> beta(coord->GetBeta());
	 vector<double> spall_spectrum( spnet->GetXSec(uid, uid) );
	 
	 //cout << "Sec. proton or tert. antiproton/antideuteron: " << uid << " "<< issec << endl;
	 
	 for (int k = 0; k < dimx; ++k) {
	   for (int j = 0; j < dimy; ++j) {
	     for (int l = 0; l < dimz; ++l) {
	       int indspat = coord->indexD(k,j,l);
	       double gasdensity = totalgas->GetGas(indspat);
	       
	       
	       //if (l == dimz/2 && j==0)
	       //	cout << k << " " << j << "  " << gasdensity << endl;
	       
	       //if (l == dimz/2 && j==dimy/2)
	       //	cout << k << " " << j << "  " << gasdensity << endl;
               
	       for (int i = 0; i < dimE; ++i)  {
		 double betaigasdensity = beta[i]*gasdensity;
		 int ind = indspat*dimE+i;
		 for (int ii = i+1; ii < dimE; ++ii) {
		   result[ind] += (*iprim)->GetDensity(indspat*dimE+ii) * betaigasdensity * spall_spectrum[ii];
		 }
		 
	       }
	     }
	   }
	 }
       } 
       else { //Fluka model
	 
	 cout << "Fluka Secondary protons!" << endl;
	 
	 //pair<int,int> coupleppr(1001,1001);  // Secondary protons, from prim. protons
	 //pair<int,int> couplepHe(2004,1001);  // Secondary protons, from Helium	
	 //pair<int,int> coupleapap(-999, -999);  // Tertiary antiprotons, from sec. antiprotons
	 
	 //cout << "*** FLUKA Spallation cross section: ***" << endl;
	 
	 vector<TParticle*>::iterator iprim = part.end()-2; // Primary species
	 
	 vector<double> beta(coord->GetBeta());
	 //vector<double> spall_spectrum( spnet->GetXSec(uid, uid) );
	 
	 //cout << "*** Spallation cross section: ***" << endl;
	 
	 for (int k = 0; k < dimx; ++k) {
	   for (int j = 0; j < dimy; ++j) {
	     for (int l = 0; l < dimz; ++l) {
	       int indspat = coord->indexD(k,j,l);
	       double gasdensity = totalgas->GetGas(indspat);
	       
	       
	       //if (l == dimz/2 && j==0)
	       //	cout << k << " " << j << "  " << gasdensity << endl;
	       
	       //if (l == dimz/2 && j==dimy/2)
	       //	cout << k << " " << j << "  " << gasdensity << endl;
	       
	       for (int i = 0; i < dimE; ++i)  {
		 double betaigasdensity = beta[i]*gasdensity;
		 int ind = indspat*dimE+i;
		 
		 vector<double> spall_spectrum( spnet->GetXSecApEl(uid, uid, i) );
		 // if (k==0 && l==0)
		 // 	 for (int ie=0; ie<spall_spectrum.size(); ie++)
		 // 	   cout << i << " " << spall_spectrum[ie] << endl;
		 
		 for (int ii = i+1; ii < dimE; ++ii) {
		   result[ind] += (*iprim)->GetDensity(indspat*dimE+ii) * betaigasdensity * spall_spectrum[ii];
		 }
		 
	       }
	     }
	   }
	 }
       }
     }
     else if (uid == -999 && issec)
       cout << "Tertiary antiprotons" << endl;
     if (in->apy != FlukaAp)   { 
       
       cout << "*** DRAGON procedure for tertiary antiprotons! " << endl;
       
       vector<TParticle*>::iterator iprim = part.end()-2; // Primary species
       
       vector<double> beta(coord->GetBeta());
       vector<double> spall_spectrum( spnet->GetXSec(uid, uid) );
       
       //cout << "Sec. proton or tert. antiproton/antideuteron: " << uid << " "<< issec << endl;
       
       for (int k = 0; k < dimx; ++k) {
	 for (int j = 0; j < dimy; ++j) {
	   for (int l = 0; l < dimz; ++l) {
	     int indspat = coord->indexD(k,j,l);
	     double gasdensity = totalgas->GetGas(indspat);
	     
	     
	     //if (l == dimz/2 && j==0)
	     //	cout << k << " " << j << "  " << gasdensity << endl;
	     
	     //if (l == dimz/2 && j==dimy/2)
	     //	cout << k << " " << j << "  " << gasdensity << endl;
             
	     for (int i = 0; i < dimE; ++i)  {
	       double betaigasdensity = beta[i]*gasdensity;
	       int ind = indspat*dimE+i;
	       for (int ii = i+1; ii < dimE; ++ii) {
		 result[ind] += (*iprim)->GetDensity(indspat*dimE+ii) * betaigasdensity * spall_spectrum[ii];
	       }
	       
	     }
	   }
	 }
       }
     } else { //Fluka model
       
       cout << "Fluka tertiaty antiprotons!" << endl;
       
       //pair<int,int> coupleppr(1001,1001);  // Secondary protons, from prim. protons
       //pair<int,int> couplepHe(2004,1001);  // Secondary protons, from Helium	
       //pair<int,int> coupleapap(-999, -999);  // Tertiary antiprotons, from sec. antiprotons
       
       //cout << "*** FLUKA Spallation cross section: ***" << endl;
       
       vector<TParticle*>::iterator iprim = part.end()-2; // Primary species
       
       vector<double> beta(coord->GetBeta());
       
       for (int k = 0; k < dimx; ++k) {
	 for (int j = 0; j < dimy; ++j) {
	   for (int l = 0; l < dimz; ++l) {
	     int indspat = coord->indexD(k,j,l);
	     double gasdensity = totalgas->GetGas(indspat);
	     
	     
	     //if (l == dimz/2 && j==0)
	     //	cout << k << " " << j << "  " << gasdensity << endl;
	     
	     //if (l == dimz/2 && j==dimy/2)
	     //	cout << k << " " << j << "  " << gasdensity << endl;
             
	     for (int i = 0; i < dimE; ++i)  {
	       double betaigasdensity = beta[i]*gasdensity;
	       int ind = indspat*dimE+i;
	       
	       vector<double> spall_spectrum( spnet->GetXSecApEl(uid, uid, i) );
	       // if (k==0 && l==0)
	       // 	 for (int ie=0; ie<spall_spectrum.size(); ie++)
	       // 	   cout << i << " " << spall_spectrum[ie] << endl;
	       
	       for (int ii = i+1; ii < dimE; ++ii) {
		 result[ind] += (*iprim)->GetDensity(indspat*dimE+ii) * betaigasdensity * spall_spectrum[ii];
	       }
	       
	     }
	   }
	 }
       }
     }
     
   } // if (issec && uid != -1000)
   
   else if (uid > 1000) { // nuclei spallation
     
     vector<double> gamma(coord->GetGamma());
     vector<double> beta(coord->GetBeta());
     vector<double> energy(coord->GetEk());
     
     //modified
     //if K_capture is present, the source term due to electron attachment/stripping is computed
     if (K_electron > 0) {
       
       double factor = 1.e-27*Clight;
         
         TParticle* naked_nucleus;
         if (part.size() > 1) naked_nucleus = part[part.size()-2];
         cout << "Computing source term for dressed nucleus " << endl;
         //cout << "source term coming from nucleus with A =  " << naked_nucleus->GetA() << "; Z = " << naked_nucleus->GetZ() << endl;
         /*
          // To be redone
          for (int k = 0; k < dimr; ++k) {
          for (int l = 0; l < dimz; ++l) {
          int indspat = coord->index(k,l);
          double Gasdensity = totalgas->GetGas(indspat);
          for (int i = 0; i < dimE; ++i) {
          int ind = indspat*dimE+i;
          double attach_H=0., attach_He=0., strip_H=0., strip_He=0.;
          xsecmodel[0]->Kcapture_cs(energy[i],naked_nucleus->GetZ(),1,&attach_H,&strip_H);
          xsecmodel[0]->Kcapture_cs(energy[i],naked_nucleus->GetZ(),2,&attach_He,&strip_He);
          result[ind] = Gasdensity * (attach_H + He_abundance*attach_He) * factor * beta[i] * naked_nucleus->GetDensity(ind);
          }
          }
          }
          */
         for (int k = 0; k < dimx; ++k) {
            for (int j = 0; j < dimy; ++j) {
               for (int l = 0; l < dimz; ++l) {
                  int indspat = coord->indexD(k,j,l);
                  double Gasdensity = totalgas->GetGas(indspat);
                  
                  for (int i = 0; i < dimE; ++i)  {
		    //double betaigasdensity = beta[i]*gasdensity;
                     int ind = indspat*dimE+i;
                     
                     double attach_H=0., attach_He=0., strip_H=0., strip_He=0.;
                     xsecmodel[0]->Kcapture_cs(energy[i],naked_nucleus->GetZ(),1,&attach_H,&strip_H);
                     xsecmodel[0]->Kcapture_cs(energy[i],naked_nucleus->GetZ(),2,&attach_He,&strip_He);
                     result[ind] = Gasdensity * (attach_H + He_abundance*attach_He) * factor * beta[i] * naked_nucleus->GetDensity(ind);
                  }
               }
            }
         }
         
      }
      
      //modified
      //if K_capture is present and K_electron == 1, the source term ends here: it only contains contribution from the corresponding naked nucleus
      
      if (K_electron <= 0) { // if K_capture is present and K_electron == 0, the source term must contain also the spallation and decay contribution from heavier nuclei
         
         for (vector<TParticle*>::iterator ipart = part.begin(); ipart != part.end()-1; ++ipart) { // // //
            
            vector<double> spall_spectrum( spnet->GetXSec( (*ipart)->GetUid(), uid ) );
            
            if (spall_spectrum.size() == dimE || (*ipart)->GetDaughter() == uid) {
               
               double Afactor = double((*ipart)->GetA())/double(A); // To convert between kinetic energy per nucleon (which is approx. constant) to momentum
               
               for (int k = 0; k < dimx; ++k) {
                  for (int j = 0; j < dimy; ++j) {
                     for (int l = 0; l < dimz; ++l) {
                        int indspat = coord->indexD(k,j,l);
                        double Afactorgasdensity = Afactor*totalgas->GetGas(indspat);
                        for (int i = 0; i < dimE; ++i) {
                           int ind = indspat*dimE+i;
                           if (spall_spectrum.size() == dimE)  result[ind] += Afactorgasdensity * spall_spectrum[i] * (*ipart)->GetDensity(ind); // spallation
                           if ((*ipart)->GetDaughter() == uid) result[ind] += (*ipart)->GetDensity(ind)/(*ipart)->GetLifetime()/gamma[i]; // decay
                        }
                     }
                  }
               }
            }
         }
      } // if (K_electron == 0)
      
   } // else if (uid > 1000)
   else { // Secondary antiprotons, electrons and positrons, antideuterons
      
      for (vector<TParticle*>::iterator ipart = part.begin(); ipart != part.end()-1; ++ipart) {
         if ((*ipart)->GetUid() > 2004 || (*ipart)->GetUid() < 0) continue;
         for (int i = 0; i < dimE; i++) {
            vector<double> spall_spectrum( spnet->GetXSecApEl((*ipart)->GetUid(), uid, i) );
            if (spall_spectrum.size() == dimE) {
               for (int k = 0; k < dimx; k++) {
                  for (int j = 0; j < dimy; ++j) {
                     for (int l = 0; l < dimz; l++) {
                        int indspat = coord->indexD(k,j,l);
                        double gasdensity = totalgas->GetGas(indspat);
                        int ind = indspat*dimE+i;
                        // ENERGY INTEGRAL
                        for (int ii=i+1; ii < dimE; ii++)	result[ind] += (*ipart)->GetDensity(indspat*dimE+ii) * gasdensity * spall_spectrum[ii];
                     }
                  }
               }
            }
         }
      }
   } // else
   
   return result;
}


//******************************************************************************************************************************************************

void TParticle3D::Print(fitsfile* output_ptr, double norm) { /**< Print all the information relevant to that nucleus: charge, mass, source abundance, injection spectrum and propagated density. */
   int status = 0;
   
   // cout << "************************ " << endl;  
   // cout << "Particle3D print; norm = " << norm << endl;
   // cout << "A = " <<  A << " Z = " << Z << " issec = " << issec << endl;
   // cout << "************************ " << endl;  
   // cout << endl << endl;   

   const long naxis = 4;
   long size_axes[naxis] = {dimE,dimx,dimy,dimz};
   long nelements = size_axes[0]*size_axes[1]*size_axes[2]*size_axes[3];
   long fpixel = 1;
   int bitpix = FLOAT_IMG;
   if (fits_create_img(output_ptr, bitpix, naxis, size_axes, &status)) fits_report_error(stderr, status);
   
   double sab = _fGalaxy->GetSourceAbundance(uid);
   
   /*double injind_rho_0 = _fGalaxy->GetInjSpectrum_rho_0(uid);
   double injind_rho_1 = _fGalaxy->GetInjSpectrum_rho_1(uid);
   double injind_rho_2 = _fGalaxy->GetInjSpectrum_rho_2(uid);
   double injind_alpha_0 = _fGalaxy->GetInjSpectrum_alpha_0(uid);
   double injind_alpha_1 = _fGalaxy->GetInjSpectrum_alpha_1(uid);
   double injind_alpha_2 = _fGalaxy->GetInjSpectrum_alpha_2(uid);
   double injind_alpha_3 = _fGalaxy->GetInjSpectrum_alpha_3(uid);
   */
   if (fits_write_key(output_ptr, TINT,    (char*) "Z_",      &Z,                 NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "A",       &A,                 NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "Sec",     &issec,             NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "DM",     &isDM,             NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "EXTRA",     &isExtra,             NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "TPP",     &isTPP,             NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "S_Ab",    &sab,               NULL, &status)) fits_report_error(stderr, status);
   /*if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_rho_0", &injind_rho_0,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_rho_1", &injind_rho_1,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_rho_2", &injind_rho_2,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_alpha_0", &injind_alpha_0,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_alpha_1", &injind_alpha_1,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_alpha_2", &injind_alpha_2,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_alpha_3", &injind_alpha_3,            NULL, &status)) fits_report_error(stderr, status);
   */
   
   int counter = 0;
   float* array = new float[nelements];
   double weight = 1;
   if (A != 0) weight = double(A);
   
   double* M_number = new double[dimE]; //MW: measure particle count in Galaxy without Local Bubble
   for (int j = 0; j < dimE; ++j) M_number[j] = 0;

   vector <double> x = _fGalaxy->GetCoordinates()->GetX();
   vector <double> y = _fGalaxy->GetCoordinates()->GetY();
   vector <double> z = _fGalaxy->GetCoordinates()->GetZ();
   vector <double> Ek = _fGalaxy->GetCoordinates()->GetEk();

   unsigned int ixsun = 0;
   while(x[++ixsun] <= in->xobs){}
   ixsun--;
   unsigned int iysun = 0;
   while(y[++iysun] <= in->yobs){}
   iysun--;
   unsigned int izsun = 0;
   while(z[++izsun] <= in->zobs){}
   izsun--;
    
   map<int, vector<double> > N_near;
   map<int, double> anisotropy;
 
   for (int l = 0; l < dimz; ++l) {
      for (int k = 0; k < dimy; ++k) {
         for (int i = 0; i < dimx; ++i) {

            for (int j = 0; j < dimE; ++j) {
               //array[counter] = float(weight*density[index(i,k,l,j)]*norm*C/4.0/M_PI*pow(kpc,-3.0)*1e4);
               array[counter] = float(weight*density[index(i,k,l,j)]*norm);
               
               //                     if (not _fGalaxy->GetCoordinates()->IsInLocalBubble_Indexed(i,k,l)
               //                     and not _fGalaxy->GetCoordinates()->IsInLocalBubble_Indexed(i+1,k,l)
               //                     and not _fGalaxy->GetCoordinates()->IsInLocalBubble_Indexed(i-1,k,l)
               //                     and not _fGalaxy->GetCoordinates()->IsInLocalBubble_Indexed(i,k+1,l)
               //                     and not _fGalaxy->GetCoordinates()->IsInLocalBubble_Indexed(i,k-1,l)
               //                     and not _fGalaxy->GetCoordinates()->IsInLocalBubble_Indexed(i,k,l+1)
               //                     and not _fGalaxy->GetCoordinates()->IsInLocalBubble_Indexed(i,k,l-1)) M_number[j] += array[counter];//MW 130407
               
	       if ( (fabs(i-ixsun) < 2) && (fabs(k-iysun) < 2) && (fabs(l-izsun)< 1) ) 	
		N_near[j].push_back(density[index(i,k,l,j)]*norm*weight);

               counter++;
            }
         }
      }
   }
   
   for (int j = 0; j < dimE; ++j) {
     // cout << "energy: " << Ek[j] << "; anisotropy: ";
     double N_min_element = *min_element(N_near[j].begin(), N_near[j].end());	
     double N_max_element = *max_element(N_near[j].begin(), N_near[j].end());	
     anisotropy[j] = ( N_max_element - N_min_element ) / ( N_max_element + N_min_element );   
     // cout << anisotropy[j] << endl;
     //for (int i=0; i<N_near[j].size(); i++)
	//cout << N_near[j][i] << " " ;
     //cout << endl;
     //cout << N_min_element << endl;
     //cout << N_max_element << endl;

   }

   cout << endl;
   if (fits_write_img(output_ptr, TFLOAT, fpixel, nelements, array, &status)) fits_report_error(stderr, status);
   
   delete [] array;
   
   //MW130407 - maybe move that out of particle.cc?x
   //     cout << "UID " << Z*1000+A << endl;
   //     for (int j = 0; j < dimE; ++j) cout << "M_number for energy " << _fGalaxy->GetCoordinates()->GetEk()[j] << " : " << M_number[j] << endl;
   
   return ;
}

void TParticle3D::PrintSpectrum(fitsfile* output_ptr, double norm) { /**< Print all the information relevant to that nucleus: charge, mass, source abundance, injection spectrum and propagated spectrum at Sun position. */
   
   //     cout << "PrintSpectrum " << endl;
   //     cout << Z << "<--- Z || A --->  " << A << endl;
   
   int status = 0;
   const long naxis = 1;
   long size_axes[naxis] = {dimE};
   long nelements = size_axes[0];
   long fpixel = 1;
   int bitpix = FLOAT_IMG;
   if (fits_create_img(output_ptr, bitpix, naxis, size_axes, &status)) fits_report_error(stderr, status);
   
   double sab = (!issec) ? _fGalaxy->GetSourceAbundance(uid) : 0.0;
   
   /*double injind_rho_0 = _fGalaxy->GetInjSpectrum_rho_0(uid);
   double injind_rho_1 = _fGalaxy->GetInjSpectrum_rho_1(uid);
   double injind_rho_2 = _fGalaxy->GetInjSpectrum_rho_2(uid);
   double injind_alpha_0 = _fGalaxy->GetInjSpectrum_alpha_0(uid);
   double injind_alpha_1 = _fGalaxy->GetInjSpectrum_alpha_1(uid);
   double injind_alpha_2 = _fGalaxy->GetInjSpectrum_alpha_2(uid);
   double injind_alpha_3 = _fGalaxy->GetInjSpectrum_alpha_3(uid);
   */
   if (fits_write_key(output_ptr, TINT,    (char*) "Z_",      &Z,                 NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "A",       &A,                 NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "Sec",     &issec,             NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "DM",     &isDM,             NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "EXTRA",     &isExtra,             NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TINT,    (char*) "TPP",     &isTPP,             NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "S_Ab",    &sab,               NULL, &status)) fits_report_error(stderr, status);
   /*if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_rho_0", &injind_rho_0,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_rho_1", &injind_rho_1,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_rho_2", &injind_rho_2,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_alpha_0", &injind_alpha_0,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_alpha_1", &injind_alpha_1,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_alpha_2", &injind_alpha_2,            NULL, &status)) fits_report_error(stderr, status);
   if (fits_write_key(output_ptr, TDOUBLE, (char*) "SpecInd_alpha_3", &injind_alpha_3,            NULL, &status)) fits_report_error(stderr, status);
   */
   float* array = new float[nelements];
   double weight = 1;
   
   vector<double> r = _fGalaxy->GetCoordinates()->GetX();
   vector<double> z = _fGalaxy->GetCoordinates()->GetZ();
   vector<double> Ek = _fGalaxy->GetCoordinates()->GetEk();
   
   unsigned int irsun = (unsigned int) ((in->robs-r.front())/(r.back()-r.front())*(double)(dimx-1));
   
   //  unsigned int irsun = (unsigned int) (robs/Rmax*(double)(dimr-1));
   unsigned int izsun = (unsigned int) ((in->zobs-z.front())/(2.0*z.back())*(double)(dimz-1));
   
   double r1 = (r[irsun+1]-in->robs)/(r[irsun+1]-r[irsun]);
   double r2 = (in->robs-r[irsun])/(r[irsun+1]-r[irsun]);
   if (A != 0) weight = double(A);

   if (isDM && in->feedback >0)
   cout << "**** particle originating from DARK MATTER *****" << endl;

   cout << "************************************************" << endl;
   cout << "*** Writing spectrum at Solar System position..." << endl;
   cout << "************************************************" << endl;
   cout << "A = " <<  A << " Z = " << Z << " issec = " << issec << endl;
   cout << "************************************************" << endl;
   
   cout << "Energy -- Normalized flux " << endl;
   for (int j = 0; j < dimE; ++j)
   {
      array[j] = float(weight*norm*( density[index(irsun,(dimy-1)/2,izsun,j)]*r1 + density[index(irsun+1,(dimy-1)/2,izsun,j)]*r2 )/**C/4.0/M_PI*pow(kpc,-3.0)*1e4*/);
      cout << Ek[j] << "  " << array[j] << endl;
   }
   cout << endl;
   
   //MW130725: what for??
   if (Z == 5 || Z == 6)
   {
      cout << Z << " <--- Z | A ---> " << A << endl;
      for (int j = 0; j < dimE; ++j)
         cout << j << " " << (weight*norm*( density[index(irsun,(dimy-1)/2,izsun,j)]*r1 + density[index(irsun+1,(dimy-1)/2,izsun,j)]*r2 )) << endl;
      cout << endl;
   }
   
   if (fits_write_img(output_ptr, TFLOAT, fpixel, nelements, array, &status)) fits_report_error(stderr, status);
   delete [] array;
   return ;
}

vector<double> TParticle3D::GetSpectrumAtSunPosition() {

   vector<double> result;
   vector<double> r =  _fGalaxy->GetCoordinates()->GetX();
   vector<double> z =  _fGalaxy->GetCoordinates()->GetZ();
   vector<double> Ek = _fGalaxy->GetCoordinates()->GetEk();
   unsigned int irsun = (unsigned int) ((in->robs-r.front())/(r.back()-r.front())*(double)(dimx-1));
   unsigned int izsun = (unsigned int) ((in->zobs-z.front())/(2.0*z.back())*(double)(dimz-1));
   double r1 = (r[irsun+1]-in->robs)/(r[irsun+1]-r[irsun]);
   double r2 = (in->robs-r[irsun])/(r[irsun+1]-r[irsun]);
   double weight = 1;
   if (A != 0) weight = double(A);
   for (int j = 0; j < dimE; ++j)
      result[j] = density[index(irsun,(dimy-1)/2,izsun,j)]*r1 + density[index(irsun+1,(dimy-1)/2,izsun,j)]*r2;
   return result;

}

double TParticle3D::GetFluxAtSunPosition(int j) {

   vector<double> r =  _fGalaxy->GetCoordinates()->GetX();
   vector<double> z =  _fGalaxy->GetCoordinates()->GetZ();
   vector<double> Ek = _fGalaxy->GetCoordinates()->GetEk();
   unsigned int irsun = (unsigned int) ((in->robs-r.front())/(r.back()-r.front())*(double)(dimx-1));
   unsigned int izsun = (unsigned int) ((in->zobs-z.front())/(2.0*z.back())*(double)(dimz-1));
   double r1 = (r[irsun+1]-in->robs)/(r[irsun+1]-r[irsun]);
   double r2 = (in->robs-r[irsun])/(r[irsun+1]-r[irsun]);
   double weight = 1;
   if (A != 0) weight = double(A);
   return density[index(irsun,(dimy-1)/2,izsun,j)]*r1 + density[index(irsun+1,(dimy-1)/2,izsun,j)]*r2;

}



//******************************************************************************************************************************************************
/*
std::vector <double> TParticle3D::GetSpectra(double norm,string run_id,int writeflag,double *eng) {
   
   const long naxis = 1;
   long size_axes[naxis] = {dimE};
   long nelements = size_axes[0];
   
   float* array = new float[nelements];
   double weight = 1;
   
   vector<double> r = _fGalaxy->GetCoordinates()->GetX();
   unsigned int irsun = (unsigned int) ((robs-r.front())/(r.back()-r.front())*(double)(dimx-1));
   
   //  unsigned int irsun = (unsigned int) (robs/Rmax*(double)(dimr-1));
   unsigned int izsun = (unsigned int) ((zobs-_fGalaxy->GetCoordinates()->GetZ().front())/(2.0*_fGalaxy->GetCoordinates()->GetZ().back())*(double)(dimz-1));
   
   double r1 = (r[irsun+1]-robs)/(r[irsun+1]-r[irsun]);
   double r2 = (robs-r[irsun])/(r[irsun+1]-r[irsun]);
   if (A != 0) weight = double(A);
   for (int j = 0; j < dimE; ++j) array[j] = float(weight*norm*( density[index(irsun,(dimy-1)/2,izsun,j)]*r1 + density[index(irsun+1,(dimy-1)/2,izsun,j)]*r2 ));
  
	cout << "suns pos x: " << _fGalaxy->GetCoordinates()->GetX()[irsun] << " y: " <<   _fGalaxy->GetCoordinates()->GetX()[(dimy-1)/2] << " z: " << _fGalaxy->GetCoordinates()->GetZ()[izsun] << endl;
 
   if(writeflag){
      
      bool write = true;
      
      ofstream file;
      char buff_save[1000];
		
      if(Z==-1 && A==0 && issec==0) sprintf(buff_save,"./ASCII_spectra/%s/electron_prim_density.dat",run_id.c_str());
      else if(Z==-1 && A==0 &&issec==1)  sprintf(buff_save,"./ASCII_spectra/%s/electron_sec_density.dat",run_id.c_str());
      else if(Z==1 && A==1 && issec==0)  sprintf(buff_save,"./ASCII_spectra/%s/proton_prim_density.dat",run_id.c_str());
      else if(Z==1 && A==1 && issec==1)  sprintf(buff_save,"./ASCII_spectra/%s/proton_sec_density.dat",run_id.c_str());
      else write = false;
      
      if(write)
      {
         file.open(buff_save);
         for(int i=0;i<20;i++) for(int j=0;j<dimE;j++) file << i << " " << j << " " << eng[j] << " " << array[j] << endl;
         file.close();
      }
   }
   
   vector <double> spectrum;
   
   if ((Z==1 && A==1)&& (issec==0)) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);	// prim proton
   if ((Z==1 && A==1)&& (issec==1)) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);	// sec proton
   if ((Z==-1 && A==1) && (issec==0)) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);	// sec antiproton
   if ((Z==-1 && A==1) && (issec==1)) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);	// ter antiproton
   if ((Z==-1 && A==0) && (issec==0)) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);	//prim electron
   if ((Z==-1 && A==0) && (issec==1)) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);	// sec electron
   if (Z==1 && A==0) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			// positron
   if (Z==6 && A==12) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Carbon 12
   if (Z==6 && A==13) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Carbon 13
   if (Z==6 && A==14) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Carbon 14
   if (Z==5 && A==10) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Boron 10
   if (Z==5 && A==11) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Carbon 11
   if (Z==4 && A==7) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Be_7
   if (Z==4 && A==8) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Be_8
   if (Z==4 && A==9) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Be_9
   if (Z==4 && A==10) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Be_10
   if (Z==4 && A==11) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);			//Be_11
   if (Z==13 && A==26) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);		//Al_26
   if (Z==13 && A==27) for (int j = 0; j < dimE; ++j) spectrum.push_back(array[j]);	 	//Al_27
   
   delete [] array;
   return spectrum;
}
*/
//******************************************************************************************************************************************************

void TParticle::PrintELosses() {
   
   ofstream outfile("Electron_eloss.dat", ios::out);
   
   int dimEloss = eloss.size();
   vector<double> r = _fGalaxy->GetCoordinates()->GetR();
   vector<double> z = _fGalaxy->GetCoordinates()->GetZ();
   vector<double> Ekin = _fGalaxy->GetCoordinates()->GetEk();
   for (int iE = 0; iE < dimE; ++iE) {
      for (int i = 0; i < r.size(); ++i) {
         for (int j = 0; j < dimz; ++j) {
            outfile << Ekin[iE] << " " << r[i] << " " << z[j] << " ";
            for (int k = 0; k < dimEloss; k++) outfile << eloss[k]->GetDpdt(i,j,iE) << " ";
            outfile << endl;
         }
      }
   }
   
   outfile.close();
   exit(-1) ;
   return ;
}

