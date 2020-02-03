/**
 * @file spectrum.cc
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief In this file all the classes related to the CR injection spectrum are implemented.
 */

#include "spectrum.h"
#include "grid.h"
#include "input.h"
#include "utilities.h"

#ifdef HAVE_ROOT
#include "TFile.h"
#include "TTree.h"
#endif

#include "errorcode.h"
#include <cstdlib>

/**
 * @fn extern "C" double dmspec_(double&, int&, double&, double[], int&, double[], int&)
 * @brief Call to fortran interface to DarkSUSY
 */
extern "C" double dmspec_(double&, int&, double&, double[], int&, double[], int&);

using namespace std;

TSpectrum::TSpectrum(TGrid* coord) {
  spectrum = vector<double>(coord->GetDimE(), 0.0);
}

/*
//The old TSpectrum constructor without multiple breaks (before Simon's changes)
TSpectrum::TSpectrum(TGrid* Coord, Input* in, double inj, double lowinj, bool El) {
    
vector<double> pp;
if (El) pp = Coord->GetMomentumEl();
else pp = Coord->GetMomentum();
if (El) {
for (int i = 0; i < pp.size(); ++i) {
            
if (in->Spectrum == PowerLawBreak) {
                
if (pp[i] < in->sp_ref_rig_break_el) spectrum.push_back(pow(pp[i]/in->sp_ref_rig_break_el,-in->spect_ind_el_low));
else spectrum.push_back(pow(pp[i]/in->sp_ref_rig_break_el,-in->spect_ind_el)*exp(-pp[i]/in->cutoff_rig));
}
else spectrum.push_back(exp(-pow(in->sp_ref_rig_exp_el/pp[i],in->exp_cut_index_el))*pow(pp[i]/in->sp_ref_rig_break_el,-in->spect_ind_el)*exp(-pp[i]/in->cutoff_rig));
}
}
else {
for (int i = 0; i < pp.size(); ++i) {
if (in->Spectrum == PowerLawBreak) {
if (pp[i] < in->sp_ref_rig) spectrum.push_back(pow(pp[i]/in->sp_ref_rig,-lowinj));
else spectrum.push_back(pow(pp[i]/in->sp_ref_rig,-inj));
}
else {
//if (pp[i] < in->sp_ref_rig) spectrum.push_back(pow(pp[i]/in->sp_ref_rig,-lowinj));
//else spectrum.push_back(pow(pp[i]/in->sp_ref_rig,-inj));
spectrum.push_back(exp(-pow(in->sp_ref_rig_exp/pp[i],1.))*pow(pp[i]/in->sp_ref_rig,-inj));
}
}
}
}


//MW130621: see above, this HAS to be made consistent.
//DG130902: this is the new TSpectrum constructor, in my opinion this has to be called only when DMCMC block is used
//
//
*/

TSpectrum::TSpectrum(TGrid* Coord, Input* in, const vector<double>& break_positions, const vector<double>& slopes, const double CutoffRig, bool El, bool ExtraComponent) {
  
  if (in->feedback >0) cout << "Spectrum constructor " << endl;
  //DG: add possibility of exp cutoff again as in the old version!!
  
  vector<double> pp = (El) ? Coord->GetMomentumEl() : Coord->GetMomentum();
  vector<double> Ek = Coord->GetEk();

  if (El) {
    // DG30.09.2013 NEW implementation of injection slopes for nuclei and electrons with arbitrary number of breaks 
    if (El) cout << "Electrons" << endl;
    if (ExtraComponent) cout << "Extra component" << endl;
	
    vector<double> el_breaks;  
    vector<double> el_injections; 
    double el_cutoff ;
    if (ExtraComponent) {
      el_breaks     = in->inp_break_extra_positions;
      el_injections = in->inp_inj_extra_indexes;
      el_cutoff     = in->cutoff_rig_extra; cout << "Extra component cutoff: " << el_cutoff << endl; 
    }
    else {
      el_breaks     = in->inp_break_el_positions;
      el_injections = in->inp_inj_el_indexes;
      el_cutoff     = in->cutoff_rig_el; cout << "Electron cutoff: " << el_cutoff << endl;
    }

    int n_breaks = el_breaks.size();
    if (in->feedback >0) 
      cout << "Building ELECTRON spectrum with " << n_breaks << " breaks " << endl;

    for (int i = 0; i < pp.size(); ++i) {

      int slope_number = 0;
      if (n_breaks > 0) {
	for (int j = 0; j < n_breaks; j++) {
	  if (el_breaks[j] < pp[i])
	    slope_number++;
	}
      }
		
      if (n_breaks == 0)
	spectrum.push_back( pow(pp[i]/pp[0], -el_injections[0] ) * exp(-Ek[i]/el_cutoff) )	;    

      else {

	if (slope_number == 0)
	  spectrum.push_back( pow(pp[i]/el_breaks[0], -el_injections[0] )*exp(-Ek[i]/el_cutoff)  );
	else if (slope_number == 1)
	  spectrum.push_back( pow(pp[i]/el_breaks[0], -el_injections[1] )*exp(-Ek[i]/el_cutoff)  );			
	else {

	  double number = pow(pp[i]/el_breaks[slope_number-1], -el_injections[slope_number] );
	  for (int j=1; j<slope_number; j++)	
	    number *= pow(el_breaks[j]/el_breaks[j-1], -el_injections[j])  ;
		   
	  spectrum.push_back(number*exp(-Ek[i]/el_cutoff));	
	}
      }

      if (in->feedback >0) 
	cout << "energy = " << pp[i] << "; alpha = " <<  -el_injections[slope_number] << "; spectrum = " << spectrum.back() << endl;

    }
    cout << endl; 
  }

  else {
    // DG30.09.2013 NEW implementation of injection slopes for nuclei and electrons with arbitrary number of breaks 
    int n_breaks = break_positions.size();
    if (in->feedback >0)  {
      cout << "Building NUCLEAR spectrum with " << n_breaks << " breaks " << endl;
      cout << "Number of slopes: " << slopes.size() << endl;
    }

    for (int i = 0; i < pp.size(); ++i) {

      int slope_number = 0;
      if (n_breaks > 0) {
	for (int j = 0; j < n_breaks; j++) {
	  if (break_positions[j] < pp[i])
	    slope_number++;
	}
      }
		
      if (n_breaks == 0)
	spectrum.push_back( pow(pp[i]/pp[0], -slopes[0] ) )	;    
      
      else {
	
	if (slope_number == 0)
	  spectrum.push_back( pow(pp[i]/break_positions[0], -slopes[0] ) )	;
	else if (slope_number == 1)
	  spectrum.push_back( pow(pp[i]/break_positions[0], -slopes[1] ) )	;			
	else {
	  //cout<<i<<endl;
	  double number = pow(pp[i]/break_positions[slope_number-1], -slopes[slope_number] );
	  for (int j=1; j<slope_number; j++)	
	    number *= pow( break_positions[j]/break_positions[j-1], -slopes[j]  );
	  
	  spectrum.push_back(number);	
	}	
      }

      //if (CutoffRig > 0) { 
      //	for (int i = 0; i < pp.size(); ++i) spectrum[i] *= exp(-pp[i]/CutoffRig);
      //}

      if (in->feedback >0) 
	cout << "energy = " << pp[i] << "; alpha = " <<  -slopes[slope_number] << "; spectrum = " << spectrum.back() << endl;
    }
  }

  // Spectrum integral over energy: Ek*dN/dEk
  double integral = 0.;
  for ( size_t i = 0; i < Ek.size(); ++i )
    integral += Ek[i]*Ek[i]*spectrum[i];
  integral *= Coord->GetDeltaE();	

  if (in->feedback >0) 
    cout << "Spectrum has been integrated over energy. I = " << integral << endl;
        
 
  // DG30.09.2013 NEW implementation of injection slopes for nuclei and electrons with arbitrary number of breaks

  /*
    int n_breaks = break_positions.size();
    if (in->feedback >0){
    cout << "Building ";
    if (El==true) cout << "ELECTRON ";
    else cout << "NUCLEI ";
    cout << "spectrum with " << n_breaks << " breaks " << endl;
    }
    for (int i = 0; i < pp.size(); ++i) {
    
    int slope_number = 0;
    for (int j = 0; j < n_breaks; j++) {
    if (break_positions[j] < pp[i])
    slope_number++;
    }
    if (in->feedback >0) cout << "energy = " << pp[i] << " alpha = " <<  -slopes[slope_number] << endl;
      
    if (slope_number == 0)
    spectrum.push_back( pow(pp[i]/break_positions[0], -slopes[0] ) )	;
    else if (slope_number == 1)
    spectrum.push_back( pow(pp[i]/break_positions[0], -slopes[1] ) )	;
    else {
         
    double number = pow(pp[i]/break_positions[slope_number-1], -slopes[slope_number] );
    for (int j=1; j<slope_number; j++)
    number *= pow( break_positions[j]/break_positions[j-1], -slopes[j]  );
         
    spectrum.push_back(number);
         
    }
    }
   
    if (CutoffRig > 0) for (int i = 0; i < pp.size(); ++i) spectrum[i] *= exp(-pp[i]/CutoffRig);

    // Spectrum integral over energy: Ek*dN/dEk
    double integral = 0.;
    for ( size_t i = 0; i < Ek.size(); ++i )
    integral += Ek[i]*Ek[i]*spectrum[i];
    integral *= Coord->GetDeltaE();	

    if (in->feedback >0) 
    cout << "Spectrum has been integrated over energy. I = " << integral << endl;
  */
  return ;
}


TSpectrum::TSpectrum(TGrid* Coord, Input* in, int yieldk) {
  vector<double> Ek = Coord->GetEk();

  if (in->MOVING_CLUMP == true) {

    for (int i = 0; i < Ek.size(); ++i) {
	
      double clump_sp = in->clump_norm * pow(Ek[i],in->clump_inj) * exp(-Ek[i]/in->clump_cutoff); 

      if (in->DMs == Delta) {
	if (Ek[i] <= in->EkDelta && Ek[i+1] > in->EkDelta) spectrum.push_back(1.0/(Ek[i+1]-Ek[i]));
	else spectrum.push_back(0.0);
      }
      else
	spectrum.push_back(clump_sp);

    }

  } 
  else {
    if (in->DMs == Delta) {
      for (int i = 0; i < Ek.size(); ++i) {
	if (Ek[i] <= in->EkDelta && Ek[i+1] > in->EkDelta) spectrum.push_back(1.0/(Ek[i+1]-Ek[i]));
	else spectrum.push_back(0.0);
      }
    }
    else {
   
      int nE = Ek.size();
      double DSsp[nE];
      double Etot[nE];
      for (int i = 0; i < nE; ++i) Etot[i] = Ek[i];   
      dmspec_(in->mx, in->dmmode, in->sigmav, Etot, yieldk, DSsp, nE);
         
      //cout << "Injected spectrum" << endl;
      for (int i = 0; i < nE; ++i) {
	spectrum.push_back(DSsp[i]);
	//cout << Ek[i] << " " << DSsp[i] << endl;
      }
      if (in->feedback >0) cout << "Used DS spectrum" << endl;
    }
  }
}

TSpectrum::TSpectrum(TGrid* coord, string filename, int yieldk) {

  //DG17.10.2013
  //revised implementation	
  /*if (in->feedback >0)*/ cout << "TSpectrum constructor " << endl;

  spectrum = vector<double>(coord->GetDimE(), 0.0);
  vector<double> Ek = coord->GetEk();

  double logEk[Ek.size()];
  for (size_t i = 0; i < Ek.size(); ++i) {
    //cout << i << " " << Ek[i] << endl;
    logEk[i] = log10(Ek[i]); 
  }
    
  ifstream infile(filename.c_str(), ios::in);

  cout << "Reading data file " << filename.c_str() << endl;
  vector<double> logE_from_file, logsp_from_file;
  double x,y;
  while (infile >> x >> y) {
    logE_from_file.push_back(log10(x));
    logsp_from_file.push_back(log10(y));
    cout << x << " " << y << " " << log10(x) << " " << log10(y) <<endl;
  }
  cout << "Read "<< logE_from_file.size() << " points from data file" << endl;
  cout << "Interpolation to DRAGON energy grid " << endl;
  for (size_t i = 0; i < Ek.size(); ++i) { 
    if ( (logEk[i] > logE_from_file.back()+0.0001) || (logEk[i] < logE_from_file.front()) ) {
      spectrum[i] = 0.0;
      cout << Ek[i] << " - this energy is outside range provided by the user: spectrum here will set to zero. " << endl;	
      //cout << logEk[i] << " " << logE_from_file.back() << " " << logE_from_file.front() << endl;
    }
    else {
      int j = 0;
      while (logE_from_file[j] < logEk[i]) j++;
      j = 0;
      for (int jj=0; jj<logE_from_file.size(); jj++) {
	if(logEk[i] >=  logE_from_file[jj] ) {
	  j ++;
	}
      }
      cout << i << " " << Ek[i] << " " << j << " " << logsp_from_file[j-1] << " " << logsp_from_file[j] << endl;
      if(j>=logsp_from_file.size()) {
         cout << "Warning: j above the size of logsp_from_file vector " << endl;
         j = logsp_from_file.size()-1; 
      }
      double r1 = (logE_from_file[j]-logEk[i])/(logE_from_file[j] - logE_from_file[j-1]);
      spectrum[i] = pow(10, logsp_from_file[j-1]*r1 + (1.-r1)*logsp_from_file[j]);
      cout << Ek[i] << " " << " interpolated spectrum: "<< spectrum[i] << endl; 
    } 
  }
  infile.close();
  cout << "Spectrum was correctly interpolated" << endl;
  
  return  ;
}

TSpectrum::TSpectrum(TGrid* Coord, int mass, int dmmode, int yieldk) {
#ifdef HAVE_ROOT
  cout << "Mass = " << mass << endl;
    
  vector<double> Ek = Coord->GetEk();
    
  const char* possible_channels[29] = {"eL", "eR", "e", "MuL", "MuR", "Mu", "TauL", "TauR", "Tau", "q", "c", "b", "t", "WL", "WT", "W", "ZL", "ZT", "Z", "g", "Gamma", "h", "Nue", "NuMu", "NuTau", "Ve", "VMu", "VTau"};
    
  string channel = find_channel(dmmode);
    
  bool found = false;
  for (int i = 0; i < 29 && !found; ++i) found = !strcmp(channel.c_str(), possible_channels[i]);
    
  if (!found) {
    cerr << "Channel " << channel << " is not present in our database." << endl;
    exit(CHANNEL_NOT_PRESENT);
  }
    
  if (mass < 5 || mass > int(1e5)) {
    cerr << "DM mass out of range." << endl;
    exit(MASS_OUT_OF_RANGE);
  }
    
  TFile* rootfileDM = new TFile(EWdatafile.c_str(), "READONLY");
  TTree* tr;
  if (yieldk == 154) rootfileDM->GetObject("antiprotons", tr);
  else if (yieldk == 151) rootfileDM->GetObject("positrons", tr);
  tr->SetBranchStatus("*", 0);
  tr->SetBranchStatus("mDM", 1);
  tr->SetBranchStatus("log10E_over_mDM", 1);
  tr->SetBranchStatus(channel.c_str(), 1);
    
  ostringstream drawstr;
  drawstr << channel << ":log10E_over_mDM";
  ostringstream constrstr;
  constrstr << "mDM == " << mass;
  int Nap = tr->Draw(drawstr.str().c_str(), constrstr.str().c_str(), "goff");
  if (Nap == 0) {
    cerr << "No particles are available for this mass. Choose an integer value." << endl;
    exit(NO_MASS);
  }
  vector<double> energy;
  vector<double> orig_spectrum;
  for (int i = 0; i < Nap; ++i) {
    energy.push_back(double(mass)*pow(10, tr->GetV2()[i]));
    orig_spectrum.push_back(tr->GetV1()[i]/log(10.)/energy[i]);
        
    //cout << energy.back() << " " << orig_spectrum.back() << endl;
  }
  rootfileDM->Close();
    
  //cout << "Injected spectrum" << endl;
  for (size_t i = 0; i < Ek.size(); ++i) {
    if (Ek[i] > mass) spectrum.push_back(0);
    else spectrum.push_back(InterpSpectrum(Ek[i], energy, orig_spectrum));
    //cout << Ek[i] << " " << spectrum.back() << endl;
  }
  cout << "Used Ciafaloni spectrum" << endl;
#endif
  return ;
}

std::string TSpectrum::find_channel(int dmmode) {
  std::string channel;
  switch(dmmode) {
  case 3:
    channel = "e";
    break;
  case 13:
    channel = "W";
    break;
  case 12:
    channel="Z";
    break;
  case 17:
    channel="Mu";
    std::cerr << "Warning: did you mean VMu ??" << std::endl;
    break;
  case 19:
    channel="Tau";
    std::cerr << "Warning: did you mean VTau ??" << std::endl;
    break;
  case 22:
    channel = "c";
    break;
  case 24:
    channel = "t";
    break;
  case 25:
    channel = "b";
    break;
  case 26:
    channel = "g";
    break;
  default:
    std::cerr << "The mode you requested was not implemented. You requested mode " << dmmode << std::endl;
    exit(INVALIDDMMODE);
    break;
  }
    
  return channel;
}

double TSpectrum::InterpSpectrum(const double& en, const std::vector<double>& energy, const std::vector<double>& orig_spectrum) const { 
  int i = 0;
  while (energy[i] < en) ++i;
  if (i > energy.size()) std::cout << "Index out of range: Ek = " << en << std::endl;
  double r1 = log10(energy[i]/en)/log10(energy[i]/energy[i-1]);
  return pow(10., r1*log10(orig_spectrum[i-1]) + (1.0-r1)*log10(orig_spectrum[i]));
}

/*
  TSpectrumExtra::TSpectrumExtra(TGrid* grid, Input* in) : TSpectrum(grid) {
  
  vector<double> pp = grid->GetMomentumEl();
  
  for (int i = 0; i < pp.size(); ++i) {
  if (pp[i] < in->sp_ref_rig_break_el_extra) 
  spectrum[i] = (pow(pp[i]/in->sp_ref_rig_break_el_extra,-in->spect_ind_el_low_extra));
  else spectrum[i] = (pow(pp[i]/in->sp_ref_rig_break_el_extra,-in->spect_ind_el_extra)*exp(-pp[i]/in->cutoff_rig_extra));
  }}
*/
