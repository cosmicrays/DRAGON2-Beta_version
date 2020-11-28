/**
 * @file    dragon.cc
 * @author  Daniele Gaggero, Carmelo Evoli, Pedro De La Torre Luque
 * @email   daniele.gaggero@gmail.com
 * @brief   Implementation of the DRAGON class. See the .h file
 */



#include "dragon.h"
#include "grid.h"
#include "nucleilist.h"
#include "galaxy.h"
#include "particle.h"
#include "xsec.h"
#include "crevolutor.h"
#include "input.h"
#include "utilities.h"
#include "sources.h"
#include "spectrum.h"
#include <iomanip>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

//#include "include/spline.C"
//#include "include/spline.H"
//#include "include/functions.C"
//#include "include/daten.C"

#include <sys/stat.h> //sk for mkdir
#include <sys/types.h>

using namespace std;

//*****************************************************************************************************************************************
//

DRAGON::DRAGON(Input* in_ /**< Wrapper to user input. */) {

	in = in_;

	status = 0;

	norm = 1.0;
	normel = 0.0;
	normextra = 0.0;
	normDM = 0.0;

	if (in->feedback >0) cout << "Getting nuclei list... " << endl;
	list = new TNucleiList(in);
	if (in->feedback >0) cout << " done" << endl;

	if (in->feedback >0) cout << "Initializing the Galaxy...";
	gal = new Galaxy(in, list);
	if (in->feedback >0) cout << "done" << endl;

	if (in->feedback >0) cout << "Initializing Algorithm...";
	if (in->gridtype == "2D") {
		if(in->OpSplit) alg.push_back( new TCREvolutor(gal) );
		if (in->ADI) alg.push_back( new TCREvolutorADI(gal) );
	}
	else {

		if (in->OpSplit) alg.push_back( new TCREvolutor3D(gal) );
		if (in->ADI) alg.push_back( new TCREvolutor3DADI(gal) );
	}
	if (in->feedback >0) cout << "done" << endl;

	if (in->fullstore) {
		string name = make_filename(in->filename.c_str(), ".fits.gz");
		if (fits_create_file(&output_ptr, name.c_str(), &status)) fits_report_error(stderr, status);
	} else output_ptr = NULL;

	if (in->partialstore) {
		string namesp = make_filename(in->filename.c_str(), "_spectrum.fits.gz");
		if (fits_create_file(&output_ptr_sp, namesp.c_str(), &status)) fits_report_error(stderr, status);

	} else {
		output_ptr_sp = NULL;
	}
}


//*****************************************************************************************************************************************

DRAGON::~DRAGON() {

	if (gal)
		delete gal;
	if (list)
		delete list;

	for (vector<TCREvolutorBasis*>::iterator i = alg.begin(); i != alg.end(); ++i)
		delete *i;
	alg.clear();

	for (vector<TParticle*>::iterator i = particles.begin(); i != particles.end(); ++i)
		delete *i;
	particles.clear();

	if (output_ptr) {
		if (fits_close_file(output_ptr, &status))
			fits_report_error(stderr, status);
	}
	if (output_ptr_sp) {
		if (fits_close_file(output_ptr_sp, &status))
			fits_report_error(stderr, status);
	}
}

//*****************************************************************************************************************************************

TParticle* DRAGON::FindParticle(int uid, bool sec, int isDM, int isextra) {
	for (vector<TParticle*>::reverse_iterator ripart = particles.rbegin(); ripart != particles.rend(); ++ripart) {
		if ( (*ripart)->GetUid() == uid && (*ripart)->GetIsSec() == sec && (*ripart)->IsDM() == isDM && (*ripart)->IsExtra() == isextra )
			return *ripart;
	}

	return NULL;
}


//*****************************************************************************************************************************************

vector<TParticle*> DRAGON::FindSpecies(int Z) {
	vector<TParticle*> result;
	for (vector<TParticle*>::reverse_iterator ripart = particles.rbegin(); ripart != particles.rend(); ++ripart) {
		if ( (*ripart)->GetZ() == Z )
			result.push_back(*ripart);
	}

	return result;
}

TParticle* DRAGON::CreateParticle(const string& grtype, int A_ /**< Mass number */, int Z_ /**< Charge */, Galaxy* gal /**< A model for the galaxy */, Input* in /**< User input */, bool issec, vector<TXSecBase*> xsecmodel, TNucleiList* l, int K_electron, bool isDM, bool isextra, bool isTPP_) {

	if (grtype == "2D") return new TParticle2D(A_, Z_, gal, in, issec, xsecmodel, l, K_electron, isDM, isextra, isTPP_);
	if (grtype == "3D") return new TParticle3D(A_, Z_, gal, in, issec, xsecmodel, l, K_electron, isDM, isextra, isTPP_);

	return NULL;

}


//*****************************************************************************************************************************************
//******************************************** The main method of DRAGON ******************************************************************
//*****************************************************************************************************************************************

void DRAGON::Run() {



	vector<TXSecBase*> xsecmodel;
	xsecmodel.push_back(new TGalpropXSec(gal->GetCoordinates(), in));

	//cout << "Galprop cross sections loaded" << endl;

	//if (in->spallationxsec == Webber03) xsecmodel.push_back(new TWebber03());
	//if (in->spallationxsec == Fluka) {
	//	cout << "Fluka cross sections " << endl;
	//    xsecmodel.push_back(new FlukaXSec(gal->GetCoordinates()));
	//}

	//cout << "Cross sections ok" << endl;

	vector<int> list_nuc = list->GetList();

	TSpallationNetwork* spallnet = new TSpallationNetwork(gal->GetCoordinates(), in, xsecmodel, list_nuc);

	cout << "Spallation network ok " << endl;

	const string grty = in->gridtype;

	//-- Loop over particle list!--------------------------------------------------------------------------------
	for (vector<int>::iterator inuc = list_nuc.begin(); inuc != list_nuc.end(); ++inuc) {

		int A = -1000;
		int Z = -1000;
		Utility::id_nuc(*inuc, A, Z); /**< Compute the ID of the nucleus. */

		cout << "- Starting with nucleus A = " << A << " Z = " << Z << endl;

		if (in->feedback >0) cout << endl << "-----------------------------------------" << endl;
		if (in->feedback >0) cout         << "- Starting with nucleus A = " << A << " Z = " << Z << endl;
		if (in->feedback >0) cout         << "- Starting propagation..." << endl;
		if (in->feedback >0) cout         << "-----------------------------------------" << endl;
		//    int prev_uid = (particles.size() > 0) ? particles.back()->GetUid() : 0;

		if (*inuc > 1000) { // nuclei

			DECMODE decay_mode = list->GetDecayMode(*inuc);
			double life = list->GetLifeTime(*inuc);

			if (in->feedback >0){
				if (decay_mode != STABLE) cout << "Lifetime = " << life*1e6 << " yr" << endl;
				else cout << "Stable nucleus" << endl;
				if (decay_mode == EC || decay_mode == ECBM || decay_mode == ECBP) cout << "Nucleus may attach an electron and decay via EC" << endl;
			}

			if ((decay_mode != EC && decay_mode != ECBM && decay_mode != ECBP) || ((life*1e6 < t_half_limit) && !(in->Kcapture))) { //modified

				//bool issec = (prev_uid == *inuc);
				//bool isDM = ( (prop_DMap && *inuc == -999) || (prop_DMel && *inuc == -1000) || (prop_DMdeuterons && *inuc == -998) );
				//modified

				particles.push_back(CreateParticle(grty, A, Z, gal, in, false, xsecmodel, list, -1, false, false)); /**< Creates the nucleus and add it to the particle vector. */
				//cout << "Evolve" << endl; 
				particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */

				if (*inuc == 1001) {
					particles.push_back(CreateParticle(grty, A, Z, gal, in, true, xsecmodel, list, -1, false, false)); /**< Creates the nucleus and add it to the particle vector. */
					particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
				}

			}
			else { /* if the decay mode is  K_capture and the lifetime is long enough */

				//     if (in->Kcapture) {
				if (in->feedback >0) cout << "Since the lifetime is long enough, the nucleus is propagated twice" << endl;
				if (in->feedback >0) cout << "Propagating naked nucleus" << endl;
				particles.push_back(CreateParticle(grty, A, Z, gal, in, false, xsecmodel, list, 0, false, false)); // naked nucleus. K_electron = 0
				particles.back()->Evolve(particles, alg, spallnet, xsecmodel);

				//	prev_uid = particles.back()->GetUid();

				if (in->feedback >0) cout << "Propagating nucleus with an attached electron" << endl;
				particles.push_back(CreateParticle(grty, A, Z, gal, in, false, xsecmodel, list, 1, false, false)); // K_electron = 1
				particles.back()->Evolve(particles, alg, spallnet, xsecmodel);
			}
		}
		else {

			if (*inuc == -999) {
				// Antiprotons

				if (in->prop_secap) {
					particles.push_back(CreateParticle(grty, A, Z, gal, in, false, xsecmodel, list, -1, false, false));
					particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */

					particles.push_back(CreateParticle(grty, A, Z, gal, in, true, xsecmodel, list, -1, false, false));
					particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */

				}
				if (in->prop_DMap) {
					particles.push_back(CreateParticle(grty, A, Z, gal, in, false, xsecmodel, list, -1, true, false));
					particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */

					particles.push_back(CreateParticle(grty, A, Z, gal, in, true, xsecmodel, list, -1, true, false));
					particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */

				}

			}
			else if (*inuc == -998) {
				// Antideuterons

				if (in->prop_secdeuteron) {
					particles.push_back(CreateParticle(grty, A, Z, gal, in, false, xsecmodel, list, -1, false, false));
					particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */

					particles.push_back(CreateParticle(grty, A, Z, gal, in, true, xsecmodel, list, -1, false, false));
					particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */

				}
				if (in->prop_DMdeuteron) {
					particles.push_back(CreateParticle(grty, A, Z, gal, in, false, xsecmodel, list, -1, true, false));
					particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */

					particles.push_back(CreateParticle(grty, A, Z, gal, in, true, xsecmodel, list, -1, true, false));
					particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */

				}
			}
			else {
				// Leptons

				if (in->prop_extracomp && *inuc == -1000) {
					particles.push_back(CreateParticle(grty, A, Z, gal, in, false, xsecmodel, list, -1, false, true));
					particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
				}
				if (in->prop_lep) {

					// Secondary electron or positron
					particles.push_back(CreateParticle(grty, A, Z, gal, in, true, xsecmodel, list, -1, false, false));
					particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */

					//TPP: positrons from TPP
					if (*inuc == 1000 && in->prop_TPP) {
						if (in->feedback >0) cout << "TPP positrons" << endl;
						particles.push_back(CreateParticle(grty, A, Z, gal, in, true, xsecmodel, list, -1, false,false,true));
						if (in->feedback >0) cout << "Evolving TPP positrons " << endl;
						particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
					}

					// Primary electrons
					if (*inuc == -1000) {
						particles.push_back(CreateParticle(grty, A, Z, gal, in, false, xsecmodel, list, -1, false, false));
						particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
					}
				}

				if (in->prop_DMel && *inuc == -1000 ) {
					if (in->feedback >0) cout << "DM electrons" << endl;
					particles.push_back(CreateParticle(grty, A, Z, gal, in, false, xsecmodel, list, -1, true, false));
					particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
				}

				if (in->prop_DMel && *inuc == 1000 ) {
					if (in->feedback >0) cout << "DM positrons" << endl;
					particles.push_back(CreateParticle(grty, A, Z, gal, in, false, xsecmodel, list, -1, true, false));
					particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
				}	

			}
		}
		if (in->feedback >0) cout << "Propagation done.\n" << endl;
	}
	//-- End of loop over particle list!-----------------------------------------------------------------

	//-- Second Run!--------DG24.10.2013-----------------------------------------------------------------
	//---------------------------------------------------------------------------------------------------
	//-- The second run is useful to refine the computation of very heavy nuclei. In this way the beta- decays
	//-- of nuclei with lower Z into nuclei with higher Z are taken into account properly!
	//-- The second run is useless for B/C protons and leptons; 
	//-- since it makes the code much slower so it is DISABLED in the default xml
	//-- In order to enable it, put the <DoubleRun> flag after the <NuclearChain> block!
	//---------------------------------------------------------------------------------------------------

	if (in->DoubleRun == true) {

		if (in->feedback >0) cout<<endl<<endl<<" ******* STARTING SECOND ITERATION ******* "<<endl<<endl;

		//int count_nuc = 0;

		for (vector<TParticle*>::iterator i_current_part = particles.begin(); i_current_part != particles.end()-1; ++i_current_part) { 

			int A = (*i_current_part)->GetA();
			int Z = (*i_current_part)->GetZ();
			int uid = (*i_current_part)->GetUid();	
			if (in->feedback >0){
				cout << endl << "--SECOND ITERATION--SECOND ITERATION-----" << endl;
				cout         << "-----------------------------------------" << endl;
				cout         << "- Starting with nucleus A = " << A << " Z = " << Z << endl;
				cout         << "- Starting propagation..." << endl;
				cout         << "-----------------------------------------" << endl;
				cout         << "--SECOND ITERATION--SECOND ITERATION-----" << endl;
			}
			(*i_current_part)->Evolve(particles, alg, spallnet, xsecmodel, true);

		} 
	}  
	//-- End of second loop over particle vector!-----------------------------------------
	//------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------

	/**
	 * Normalization part. If it is DM then just convert units, otherwise find protons and electrons and normalize to observed flux.
	 */

	//DG17.10.2013
	//	
	// BEFORE normalization DM particles are number densities per energy unit: 1/GeV/kpc^3
	// ------
	// applying kpc^-3 ---> 1/cm^3/GeV
	// applying C/4pi  ---> 1/cm^2/s/sr/GeV
	// applying 1.e4   ---> 1/m^2/s/sr/GeV
	// -----
	// AFTER normalization DM particles have DRAGON units for the flux: 1/m^2/s/sr/GeV

	normDM = C * pow(kpc,-3.) /4./Pi*1.e4;    //C is in cm/s
	if (in->propDM && in->feedback >0) cout << "normDM = " << normDM << endl << endl;
	//double normtest = 0;

	TParticle* protons = FindParticle(1001, false) ;
	if (protons) norm = protons->FindNormalization(in->sp_ref_rig_norm, in->spect_norm);
	if (in->feedback >0) cout << "check normalization! now norm proton = " << norm << endl;

	TParticle* electrons = FindParticle( -1000, false);
	if (electrons) normel = electrons->FindNormalization(in->sp_ref_rig_el, in->spect_norm_el);   
	if (in->feedback >0) cout << "check normalization! now norm electron = " << normel << endl;

	TParticle* electronsExtra = FindParticle( -1000, false, false, true);
	if (electronsExtra) normextra = electronsExtra->FindNormalization(in->sp_ref_rig_el_extra, in->spect_norm_el_extra);
	if (in->feedback >0) cout << "check normalization! now norm extracomponent = " << electronsExtra << endl;
	
	if (in->DontNormalize) {
		
		// DOES NOT NORMALIZE to local CR flux
		
		// if you want to inject a given energy budget in [particles / (GeV kpc^3 Myr)] the source term has to be provided in those units!!!

		double norm_factor = 10000. * (1./pow(kpc,3.)) * (C/(4.*M_PI));	// (kpc in cm)^3 * c/4pi measured in    cm / ( s sr)  
		
		// Q in DRAGON is expressed in                             --> particles / (GeV kpc^3 Myr) 

		// so the propagated N is a number density per unit energy --> particles / (GeV kpc^3)

		// with norm_factor:
		// multiplying by 1/kpc3 transforms N into particles / (GeV cm^3)
		// multiplying by C/4pi  transforms N into particles / (GeV cm^2 s sr)
		// multiplying by 10000. transforms N into particles / (GeV m^2 s sr)


		norm = norm_factor;
		normel = norm_factor;
		normextra = norm_factor;


	}	

	//cout << "Finding normalization for the extra component" << endl;
	//cout << in->sp_ref_rig_el_extra << " <- energy || norm -> " << in->spect_norm_el_extra << endl	;

	//DG14.01.2014	--- injected energy rate --- --- --- --- --- --- --- --- --- 
	//updated 16.01.2014

	vector<double> Ek;
	vector<double> spectrum_;
	vector<double> source_distribution;

	/*if (protons && electrons && electronsExtra) {
    for (int ind=0; ind<3; ind++) {

      double injection_rate = 0.; //we want it in GeV/Myr
      cout << ind << endl; 
      Ek = gal->GetCoordinates()->GetEk();
      //cout << Ek.size() << endl;	
      int dimx = gal->GetCoordinates()->GetDimX();	
      int dimy = gal->GetCoordinates()->GetDimY();	
      int dimz = gal->GetCoordinates()->GetDimZ();	
      int dimE = gal->GetCoordinates()->GetDimE();
      //
      //cout << "test" << endl;
      if (ind==0) {cout <<"computing energy for proton injection" <<endl; spectrum_ = protons->GetSpectrum()->GetSpectrum();}
      if (ind==1) {cout <<"computing energy for electron injection" <<endl; spectrum_ = electrons->GetSpectrum()->GetSpectrum();}
      if (ind==2) {cout <<"computing energy for extra comp. injection" <<endl; spectrum_ = electronsExtra->GetSpectrum()->GetSpectrum();}
      //
      if (ind==0 || ind==1) source_distribution = gal->GetSource()->GetSource(); 
      if (ind==2) source_distribution = gal->GetSourceExtra()->GetSource(); 
      //
      double kpc_m = 3.08568e19;
      //   
      double normalization;
      if (ind==0) normalization	 = norm / (C/(100.*4.*M_PI));	
      if (ind==1) normalization	 = normel / (C/(100.*4.*M_PI));	
      if (ind==2) normalization	 = normextra / (C/(100.*4.*M_PI));	

      // normextra takes N from arbitrary units to 1 / (m2 s sr GeV)
      // normalization takes N from arbitrary units to  1 / (m3 GeV) -- number density of particles per unit energy
      // if normalization is applied to Q, since [Q] = [N]/[t], Q gets the following units: 1 / (m3 GeV Myr) 	
      //
      cout << "********************************************************** " << endl;
      cout << "Normalization factor A for the integral of Q [source term] " << endl;
      cout << "********************************************************** " << endl;
      cout << "units: [A] = 1 / ( m^3 GeV Myr X(arb units of Q) ) " << endl;
      cout << "source term               -> Q_phyisical [1/(m^3 GeV Myr) ]  =  A * Q [X] " << endl;
      cout << "propagated number density -> N_physical  [1/(m^3 GeV) ]      =  A * N [X Myr]  " << endl;
      cout << "********************************************************** " << endl;
      cout << "if the I = integral (E[GeV] Q dx dy dz dE) is in [GeV*X m^3 GeV] " << endl;
      cout << "then A*I is in [GeV/Myr] " << endl;
      //
      if (ind==0) cout << " A_p = "     << normalization<<endl;//	 = norm / (C/(100.*4.*M_PI));	
      if (ind==1) cout << " A_e = "     << normalization<<endl;//	 = normel / (C/(100.*4.*M_PI));	
      if (ind==2) cout << " A_extra = " << normalization<<endl;//	 = normextra / (C/(100.*4.*M_PI));	
      //
      for (int ip = 0; ip < dimE; ip++) {
        //cout << ip << "  ";
	  for (int k = 0; k < dimz; k++) {
	  for (int j = 0; j < dimy; j++) {
	    for (int i = 0; i < dimx; i++) {
	      //
	      long int linearized_index = gal->GetCoordinates()->indexD(i,j,k,ip);
              long int linearized_SPATIAL_index = gal->GetCoordinates()->indexD(i,j,k);
	      double deltax = gal->GetCoordinates()->GetDeltaX_central(i) * kpc_m; //kpc -> m
	      double deltay = gal->GetCoordinates()->GetDeltaY_central(j) * kpc_m; //kpc -> m
	      double deltaz = gal->GetCoordinates()->GetDeltaZ_central(k) * kpc_m; //kpc -> m
	      //  
	      double source_term_arbitrary_units = spectrum_[ip] * source_distribution[linearized_SPATIAL_index];
	      double source_term_normalized = source_term_arbitrary_units * normalization; // in particles / (m3 GeV Myr)
	      // X --> X * 1/(m^3 GeV Myr X) = 1/(m^3 GeV Myr)	
	      double energy_source_term_normalized = source_term_arbitrary_units * normalization * Ek[ip] ; // in GeV / (m3 GeV Myr)
	      //
	      double num =  Ek[ip]*gal->GetCoordinates()->GetDeltaE() * deltax*deltay*deltaz * energy_source_term_normalized;
	      // num is in GeV * m3 * GeV / (m3 GeV Myr) --> GeV/Myr
	      injection_rate += num;
	      //   
	    }
	  }
	  }
        //cout << endl;
      }

      if (ind==0) cout << "+++ The injection rate for the protons is " << injection_rate << " GeV/Myr +++" << endl;
      if (ind==1) cout << "+++ The injection rate for the electrons is " << injection_rate << " GeV/Myr +++" << endl;
      if (ind==2) cout << "+++ The injection rate for the extra component is " << injection_rate << " GeV/Myr +++" << endl;
    }
  }*/

	Ek.clear();
	spectrum_.clear();
	source_distribution.clear();
	//norm = 1.;

	/**
	 * Print nuclear densities and/or spectra in output file.
	 * Lighter nuclei are printed first.
	 */
	if (in->feedback >0) cout << endl << endl;

	for (vector<TParticle*>::reverse_iterator ripart = particles.rbegin(); ripart != particles.rend(); ++ripart) {
		if ((*ripart)->IsDM()) {
			if (output_ptr) (*ripart)->Print(output_ptr, normDM);
			if (output_ptr_sp) (*ripart)->PrintSpectrum(output_ptr_sp, normDM);

		}
		else {
			if ((*ripart)->IsExtra()) {
				if (output_ptr) (*ripart)->Print(output_ptr, normextra);
				if (output_ptr_sp) (*ripart)->PrintSpectrum(output_ptr_sp, normextra);

			}
			else {
				if ( ( (*ripart)->GetUid() == -1000 && !(*ripart)->GetIsSec()) || (*ripart)->IsTPP() ) {
					// primary electrons and secondary positrons from TPP are normalized to electrons

					if (output_ptr) (*ripart)->Print(output_ptr, normel);
					if (output_ptr_sp) (*ripart)->PrintSpectrum(output_ptr_sp, normel);

				}
				else {
					if (output_ptr) (*ripart)->Print(output_ptr, norm);
					if (output_ptr_sp) (*ripart)->PrintSpectrum(output_ptr_sp, norm);
				}
			}
		}
	}

	/* ASCII output */
	//*****************************************************************************************************************************************

	vector<double> e_GeV = gal->GetCoordinates()->GetEk();
	char ascii_filename[1000];
	sprintf(ascii_filename,"output/%s.txt",in->filename.c_str());
	cout << "Writing ASCII output file: " << ascii_filename << endl;
	ofstream outfile(ascii_filename, ios::out);
	outfile << "Energy [GeV] ";
	for (vector<TParticle*>::reverse_iterator current_part = particles.rbegin(); current_part != particles.rend(); ++current_part) {
		if ( (*current_part)->IsDM() && (*current_part)->GetUid() == -1000 )
			outfile << "DM_e-        " ;
		if ( (*current_part)->IsDM() && (*current_part)->GetUid() == 1000 )
			outfile << "DM_e+        " ;
		if ( (*current_part)->IsDM() && (*current_part)->GetUid() == -999 )
			outfile << "DM_pbar      " ;
		if ( !((*current_part)->IsDM()) && (*current_part)->GetUid() == -1000 && !(*current_part)->GetIsSec() &&  !(*current_part)->IsExtra() )
			outfile << "Pri_e-       " ;
		if ( !((*current_part)->IsDM()) && (*current_part)->GetUid() == -1000 &&  (*current_part)->GetIsSec() &&  !(*current_part)->IsExtra() )
			outfile << "Sec_e-       " ;
		if ( !((*current_part)->IsDM()) && (*current_part)->GetUid() ==  1000 &&  (*current_part)->GetIsSec() &&  !(*current_part)->IsExtra() )
			outfile << "Sec_e+       " ;
		if ( !((*current_part)->IsDM()) && (*current_part)->GetUid() == -1000 &&  !(*current_part)->GetIsSec() &&   (*current_part)->IsExtra() )
			outfile << "Extra_e+     " ;
		if ( !((*current_part)->IsDM()) && (*current_part)->GetUid() == -999 && !(*current_part)->GetIsSec()  &&  !(*current_part)->IsExtra() )
			outfile << "Sec_pbar     " ;
		if ( !((*current_part)->IsDM()) && (*current_part)->GetUid() == -999 &&  (*current_part)->GetIsSec()  &&  !(*current_part)->IsExtra() )
			outfile << "Ter_pbar     " ;
		if ( !((*current_part)->IsDM()) && (*current_part)->GetUid() == 1001 && !(*current_part)->GetIsSec()  &&  !(*current_part)->IsExtra() )
			outfile << "Pri_p        " ;
		if ( !((*current_part)->IsDM()) && (*current_part)->GetUid() == 1001 &&  (*current_part)->GetIsSec()  &&  !(*current_part)->IsExtra() )
			outfile << "Sec_p        " ;
		if ( (*current_part)->GetUid() > 1001)
			outfile << "NUC_" << (*current_part)->GetUid() <<"     ";
	}
	outfile << endl;
	for (int ie=0; ie<e_GeV.size(); ie++) {
		//
		outfile << scientific << setprecision(6) << e_GeV[ie] << " ";
		//
		for (vector<TParticle*>::reverse_iterator current_part = particles.rbegin(); current_part != particles.rend(); ++current_part) {
			//
			double flux = (*current_part)->GetFluxAtSunPosition(ie);
			//
			if ((*current_part)->IsDM())
				flux *= normDM;
			else {
				if ((*current_part)->IsExtra())			
					flux *= normextra;
				else {
					if ( ( (*current_part)->GetUid() == -1000 && !(*current_part)->GetIsSec()) || (*current_part)->IsTPP() )
						flux *= normel;
					else
						flux *= norm;
				}
			}	 		
			outfile << flux << " ";
		}	
		outfile << endl;
	}
	outfile.close(); 


	delete spallnet;
	spallnet = NULL;
	while (xsecmodel.size()) {
		delete xsecmodel.back();
		xsecmodel.pop_back();
	}

	return ;
}

//*****************************************************************************************************************************************

void DRAGON::Print() {

	double Rmax = in->Rmax;

	unsigned int cr_irsun;
	if (in->gridtype == "2D")
		cr_irsun = (unsigned int) (in->robs/Rmax*(double)(in->numr-1)); /**< Radial Sun position. */
	else
		cr_irsun = (unsigned int) ((in->robs+Rmax)/(2.0*Rmax)*(double)(in->numx-1)); /**< Radial Sun position. */
	unsigned int cr_ixsun;
	unsigned int cr_iysun;
	if (in->gridtype == "3D") { //DG28.11.2013 added ixsun and iyzun in the 3D case
		cr_ixsun = (unsigned int) ((in->xobs+Rmax)/(2.0*Rmax)*(double)(in->numx-1));
		cr_iysun = (unsigned int) ((in->yobs+Rmax)/(2.0*Rmax)*(double)(in->numy-1));
	}   
	unsigned int cr_izsun = (unsigned int) ((in->zobs+in->zmax)/(2.0*in->zmax)*(double)(in->numz-1)), /**< Vertical Sun position. */

			cr_gas_model = (unsigned int)in->gas_model,
			cr_SNR_model = (unsigned int)in->SNR_model;

	double cr_rmin = in->Rmin,
			cr_rmax = Rmax,
			cr_zmin = -in->zmax,
			cr_zmax =  in->zmax,
			cr_rB = rB,
			cr_zt = in->zt,
			cr_zr = zr,
			cr_Bh = Bh,
			cr_D0 = in->D0*kpc*kpc/Myr,
			cr_D_ref_rig = in->D_ref_rig,
			cr_ind_diff = in->delta,
			cr_He_abundance = He_abundance,
			cr_sp_ref_rig = in->sp_ref_rig,
			cr_sp_ref_rig_norm = in->sp_ref_rig_norm,
			cr_spect_norm = in->spect_norm,
			cr_Ekmax = in->Ekmax,
			cr_Ekmin = in->Ekmin,
			cr_Ekfac = in->Ekfact,
			cr_u = u,
			//cr_tol = tolerance,
			//cr_alpha = alpha,
			cr_p = p,
			cr_vA = in->vAlfven*kpc/Myr/km,
			cr_v0 = in->v0*kpc/Myr/km,
			cr_dvdz = in->dvdz*kpc/Myr/km,
			cr_robs = in->robs,
			cr_zobs = in->zobs;

	double cr_xmin, cr_xmax, cr_ymin, cr_ymax, cr_xobs, cr_yobs; //DG28.11.2013
	if (in->gridtype == "3D") {
		cr_xmin = -in->Rmax;
		cr_xmax =  in->Rmax;
		cr_ymin = -in->Rmax;
		cr_ymax =  in->Rmax;
		cr_xobs =  in->xobs;
		cr_yobs =  in->yobs;
	}   

	double crmx = in->mx;
	double crtaudec = in->taudec;
	double crsigmav = in->sigmav;
	int crdmmode = in->dmmode;
	int crdmprof = in->dmprof;

	const long naxis = (in->gridtype == "2D") ? 3 : 4;
	int dimE = int(log(in->Ekmax/in->Ekmin)/log(in->Ekfact) + 1.9);
	long size_axes[naxis];
	dimE = gal->GetCoordinates()->GetDimE();


	size_axes[0] =dimE;
	if (in->gridtype == "2D") {
		size_axes[1] = in->numr;
		size_axes[2] = in->numz;
	}
	else {
		size_axes[1] = in->numx;
		size_axes[2] = in->numy;
		size_axes[3] = in->numz;
	}
	long nelements = 1;
	for (int i = 0; i < naxis; ++i) nelements *= size_axes[i];

	long fpixel = 1;
	int bitpix = FLOAT_IMG;

	/**
	 * Create the first HDU, which only contains header data.
	 */

	if (output_ptr) {
		if (fits_create_img(output_ptr, bitpix, naxis, size_axes, &status))
			fits_report_error(stderr, status);
		//
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "rmin",      &cr_rmin,            NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*)  "rmax",      &cr_rmax,            NULL, &status))
			fits_report_error(stderr, status);
		if (in->gridtype == "3D") {//DG28.11.2013
			if (fits_write_key(output_ptr, TDOUBLE, (char*) "xmin",      &cr_xmin,            NULL, &status))
				fits_report_error(stderr, status);
			if (fits_write_key(output_ptr, TDOUBLE, (char*) "xmax",      &cr_xmax,            NULL, &status))
				fits_report_error(stderr, status);
			if (fits_write_key(output_ptr, TDOUBLE, (char*) "ymin",      &cr_ymin,            NULL, &status))
				fits_report_error(stderr, status);
			if (fits_write_key(output_ptr, TDOUBLE, (char*) "ymax",      &cr_ymax,            NULL, &status))
				fits_report_error(stderr, status);
		}
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "zmin",      &cr_zmin,            NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "zmax",      &cr_zmax,            NULL, &status))
			fits_report_error(stderr, status);
		// 
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "robs",      &cr_robs,            NULL, &status))
			fits_report_error(stderr, status);
		if (in->gridtype == "3D") {//DG28.11.2013
			if (fits_write_key(output_ptr, TDOUBLE, (char*) "xobs",      &cr_xobs,            NULL, &status))
				fits_report_error(stderr, status);
			if (fits_write_key(output_ptr, TDOUBLE, (char*) "yobs",      &cr_yobs,            NULL, &status))
				fits_report_error(stderr, status);
		}
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "zobs",      &cr_zobs,            NULL, &status))
			fits_report_error(stderr, status);
		//
		if (fits_write_key(output_ptr, TUINT,   (char*) "irsun",     &cr_irsun,           NULL, &status))
			fits_report_error(stderr, status);
		if (in->gridtype == "3D") {//DG28.11.2013
			if (fits_write_key(output_ptr, TINT, (char*) "ixsun",      &cr_ixsun,            NULL, &status))
				fits_report_error(stderr, status);
			if (fits_write_key(output_ptr, TINT, (char*) "iysun",      &cr_iysun,            NULL, &status))
				fits_report_error(stderr, status);
		}
		if (fits_write_key(output_ptr, TUINT,   (char*) "izsun",     &cr_izsun,           NULL, &status))
			fits_report_error(stderr, status);
		//
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "u",         &cr_u,               NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "rB",        &cr_rB,              NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "zt",        &cr_zt,              NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "zr",        &cr_zr,              NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "Bh",        &cr_Bh,              NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "D0",        &cr_D0,              NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "Drefrig",   &cr_D_ref_rig,       NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "ind_diff",  &cr_ind_diff,        NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "He_ab",     &cr_He_abundance,    NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "sprefrig",  &cr_sp_ref_rig,      NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "sprignor",  &cr_sp_ref_rig_norm, NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "specnorm",  &cr_spect_norm,      NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "Ekmax",     &cr_Ekmax,           NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "Ekmin",     &cr_Ekmin,           NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "Ekin_fac",  &cr_Ekfac,           NULL, &status))
			fits_report_error(stderr, status);
		//    if (fits_write_key(output_ptr, TDOUBLE, (char*) "tol",       &cr_tol,             NULL, &status)) fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TUINT,   (char*) "gasmod",    &cr_gas_model,       NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TUINT,   (char*) "SNRmod",    &cr_SNR_model,       NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TINT,    (char*) "dimE",      &dimE,               NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TINT,    (char*) "dimr",      &(in->numr),         NULL, &status))
			fits_report_error(stderr, status);
		if (in->gridtype == "3D") {//DG28.11.2013
			if (fits_write_key(output_ptr, TINT, (char*) "dimx",      &(in->numx),            NULL, &status))
				fits_report_error(stderr, status);
			if (fits_write_key(output_ptr, TINT, (char*) "dimy",      &(in->numy),            NULL, &status))
				fits_report_error(stderr, status);
		}
		if (fits_write_key(output_ptr, TINT,    (char*) "dimz",      &(in->numz),         NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "vAlfven",   &cr_vA,              NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "v0_conv",   &cr_v0,              NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "dvdz_conv", &cr_dvdz,            NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "mx",        &crmx,               NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "taudec",    &crtaudec,           NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TDOUBLE, (char*) "sigmav",    &crsigmav,           NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TINT,    (char*) "dmmode",    &crdmmode,           NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr, TINT,    (char*) "dmprof",    &crdmprof,           NULL, &status))
			fits_report_error(stderr, status);
	}
	if (output_ptr_sp) {
		if (fits_create_img(output_ptr_sp, bitpix, 1, size_axes, &status))
			fits_report_error(stderr, status);

		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "rmin",      &cr_rmin,            NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "rmax",      &cr_rmax,            NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "zmin",      &cr_zmin,            NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "zmax",      &cr_zmax,            NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "robs",      &cr_robs,            NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "zobs",      &cr_zobs,            NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TUINT,   (char*) "irsun",     &cr_irsun,           NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TUINT,   (char*) "izsun",     &cr_izsun,           NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "u",         &cr_u,               NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "rB",        &cr_rB,              NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "zt",        &cr_zt,              NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "zr",        &cr_zr,              NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "Bh",        &cr_Bh,              NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "D0",        &cr_D0,              NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "Drefrig",   &cr_D_ref_rig,       NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "ind_diff",  &cr_ind_diff,        NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "He_ab",     &cr_He_abundance,    NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "sprefrig",  &cr_sp_ref_rig,      NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "sprignor",  &cr_sp_ref_rig_norm, NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "specnorm",  &cr_spect_norm,      NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "Ekmax",     &cr_Ekmax,           NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "Ekmin",     &cr_Ekmin,           NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "Ekin_fac",  &cr_Ekfac,           NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TUINT,   (char*) "gasmod",    &cr_gas_model,       NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TUINT,   (char*) "SNRmod",    &cr_SNR_model,       NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TINT,    (char*) "dimE",      &dimE,               NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TINT,    (char*) "dimr",      &(in->numr),         NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TINT,    (char*) "dimz",      &(in->numz),         NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "vAlfven",   &cr_vA,              NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "v0_conv",   &cr_v0,              NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "dvdz_conv", &cr_dvdz,            NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "mx",        &crmx,               NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "taudec",    &crtaudec,           NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "sigmav",    &crsigmav,           NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TINT,    (char*) "dmmode",    &crdmmode,           NULL, &status))
			fits_report_error(stderr, status);
		if (fits_write_key(output_ptr_sp, TINT,    (char*) "dmprof",    &crdmprof,           NULL, &status))
			fits_report_error(stderr, status);
	}


	if (output_ptr) {
		float* array = new float[nelements]();
		if (fits_write_img(output_ptr, TFLOAT, fpixel, nelements, array, &status))
			fits_report_error(stderr, status);
		delete [] array;
	}
	if (output_ptr_sp) {
		float* array = new float[size_axes[0]]();
		if (fits_write_img(output_ptr_sp, TFLOAT, fpixel, size_axes[0], array, &status))
			fits_report_error(stderr, status);
		delete [] array;
	}
}


