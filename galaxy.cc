/**
 * @file galaxy.cc
 * @author Luca Maccione, Daniele Gaggero
 * @email luca.maccione@desy.de
 * @email daniele.gaggero@sissa.it
 * @brief In this file all the classes related to the model of the galaxy are implemented.
 */

#include "geometry.h"
#include "galaxy.h"
#include "grid.h"
#include "gas.h"
#include "input.h"
#include "nucleilist.h"
#include "sources.h"
#include "eloss.h"
#include "fitsio.h"
#include "errorcode.h"

#include "bfield.h"

#include "diffusion.h"

#include <fstream>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#define DIFFTHRESHOLD 0.2

using namespace std;

TConvectionVelocity::TConvectionVelocity(TGrid* coord, TGeometry* geom, Input* in, TSource* SourceTerm) {
    
  nrn_sn = SourceTerm->GetSource(in->xobs,in->yobs,in->zobs);
    
  if (in->feedback > 0) cout << "Called ConvectionVelocity"  << endl;
  dvdz = in->dvdz;
  conv_index_radial = in->conv_index_radial;
  set_profile_conv = in->set_profile_conv;
  conv_threshold = in->conv_threshold;
  vector<double> zgrid = coord->GetZ();
  dimz = zgrid.size();
  double rem_vel[dimz];
  double velocity=0;
  
  char buff[1000];
  sprintf(buff,"ASCII_spectra/%s/convection.dat",in->run_id.c_str());
  ofstream datafile;
  if(in->write_flag) datafile.open(buff);

  if (coord->GetType() == "2D") {
    vector<double> rgrid = coord->GetR();
    dimr = rgrid.size();
        
    for (int ir=0; ir<dimr; ir++) {
      double radius = rgrid[ir];

      for (int iz = 0; iz<dimz; iz++) {
	double zeta = zgrid[iz];
                
	if(fabs(zeta)<=in->z_k) velocity = ((((in->v0)-(in->vb))*pow(zeta,2.0)/pow(in->z_k,2.0))+(in->vb))* GetProfile(radius,0,zeta,SourceTerm);
	if(fabs(zeta)>in->z_k) velocity = ((in->v0) + dvdz*(fabs(zeta)-in->z_k)) * GetProfile(radius,0,zeta,SourceTerm);

	if((in->vb)>(in->v0)) cerr << "WARNING: vb > v0!" << endl;
         	
	//smoothen drop to avoid non-numerical values in fluxes
	if(set_profile_conv == Radial && conv_index_radial>0.){
	  if(ir==dimr-1) velocity=0.; //last bin  = zero
	  if(ir==dimr-3) rem_vel[iz]=velocity;
	  if(ir==dimr-2)velocity=0.5*rem_vel[iz]; //second to last bin is set to 0.5 * third to last bin
	}
      
	if(in->write_flag) datafile << radius << " " << zeta << " " << velocity/km/Myr*kpc << endl;

	vc.push_back(velocity);

	/*                cout << "[MW-DEBUG-VC] " << radius << " " << zeta << " | " << velocity << " " << GetProfile(radius,0,zeta,SourceTerm) << endl;*/
      }
    }

    //MW130705: CN coefficients now here for 2D, too
    for (unsigned int i = 0; i < dimr; ++i) {
      for (unsigned int k = 0; k < dimz; ++k) {

	double vCk = 0.0; // vC(i)
	double vCk1 = 0.0; // vC(i+1)
	double vC1k = 0.0; // vC(i-1)

	if ( coord->GetZ().at(k) > 0 ) {
	  vCk1 = 0.0;
	  vCk  = vc[conv_index(i,k)]/coord->GetDeltaZ_down(k);
	  vC1k = vc[conv_index(i,k-1)]/coord->GetDeltaZ_down(k);
	}
	else if ( coord->GetZ().at(k) < 0) {
	  vCk  = vc[conv_index(i,k)]/coord->GetDeltaZ_up(k);
	  vC1k = 0.0;
	  vCk1 = vc[conv_index(i,k+1)]/coord->GetDeltaZ_up(k);
	}
	else {
	  vCk  = vc[conv_index(i,k)]*2/(coord->GetDeltaZ_up(k) + coord->GetDeltaZ_down(k));
	  vC1k = -0.5*vc[conv_index(i,k-1)]/coord->GetDeltaZ_down(k);
	  vCk1 = -0.5*vc[conv_index(i,k+1)]/coord->GetDeltaZ_up(k);
	}

	CNconv_alpha1_z.push_back( vC1k );
	CNconv_alpha2_z.push_back( vCk );
	CNconv_alpha3_z.push_back( vCk1 );

	/*    cout << " [MW-DEBUG-CONV] " << i << " " << k << " " << conv_index(i,k) << " | " << vc[conv_index(i,k)] << " | " << vC1k << " " << vCk << " " << vCk1 << " | " << coord->GetDeltaZ_up(k) << " " << coord->GetDeltaZ_down(k) << " " << coord->GetDeltaZ(k) << endl;*/
      }
    }

  } // end 2D
  else {
    vector<double> xgrid = coord->GetX();
    dimx = xgrid.size();
    vector<double> ygrid = coord->GetY();
    dimy = ygrid.size();
        
    for (int ix=0; ix<dimx; ix++) {
      double x = xgrid[ix];
      for (int iy=0; iy<dimy; iy++) {
	double y = ygrid[iy];
	for (int iz = 0; iz<dimz; iz++) {
	  double z = zgrid[iz];
                    
	  if(fabs(z)<=in->z_k) velocity = ((((in->v0)-(in->vb))*pow(z,2.0)/pow(in->z_k,2.0))+(in->vb))* GetProfile(x,y,z,SourceTerm);
	  if(fabs(z)>in->z_k) velocity = ((in->v0) + dvdz*(fabs(z)-in->z_k)) * GetProfile(x,y,z,SourceTerm);
        
	  if((in->vb)>(in->v0)) cerr << "WARNING: vb > v0!" << endl;

	  velocity*= max( min( pow(geom->GetPattern(ix,iy,iz), in->SA_convec), in->SA_cut_convec), 1./in->SA_cut_convec );
	  velocity *= pow( in->LB_convec, coord->IsInLocalBubble(xgrid[ix],ygrid[iy],zgrid[iz]) );

	  vc.push_back(velocity);
	}
      }
    }

    if(set_profile_conv == Radial && conv_index_radial>0.){
      for (int ix=0; ix<dimx; ix++) for (int iy=0; iy<dimy; iy++) for (int iz = 0; iz<dimz; iz++){
	    //smoothen drop to avoid non-numerical values in fluxes, SK 06/13
	    if(ix==dimx-2 && iy==dimy-2)  vc[conv_index(ix,iy,iz)]=0.5*vc[conv_index(ix-1,iy,iz)];//second to last bin is set to 0.5 * third to last bin
	    if(iy==dimy-2 && ix==dimx-2)  vc[conv_index(ix,iy,iz)]=0.5*vc[conv_index(ix,iy-1,iz)];//second to last bin is set to 0.5 * third to last bin
	    if(ix==1 && iy==1)  vc[conv_index(ix,iy,iz)]=0.5*vc[conv_index(ix+1,iy,iz)];//second to last bin is set to 0.5 * third to last bin
	    if(iy==1 && ix==1)  vc[conv_index(ix,iy,iz)]=0.5*vc[conv_index(ix,iy+1,iz)];//second to last bin is set to 0.5 * third to last bin
	    if(ix==dimx-1 && iy==dimy-1) vc[conv_index(ix,iy,iz)]=0.; //border bins==0
	    if(ix==0 && iy==0) vc[conv_index(ix,iy,iz)]=0.; //border bins==0
	    if(ix==dimx-1 && iy==0) vc[conv_index(ix,iy,iz)]=0.; //border bins==0
	    if(ix==0 && iy==dimy-1) vc[conv_index(ix,iy,iz)]=0.; //border bins==0
	  }
    }

    //MW130624: CN coefficients are now here, the Evolutor just calls these vectors
    for (unsigned int i = 0; i < dimx; ++i) {
      for (unsigned int j = 0; j < dimy; ++j) {
	for (unsigned int k = 0; k < dimz; ++k) {

	  double vCk = 0.0; // vC(i)
	  double vCk1 = 0.0; // vC(i+1)
	  double vC1k = 0.0; // vC(i-1)

	  if ( coord->GetZ().at(k) > 0 ) {
	    vCk1 = 0.0;
	    vCk  = vc[conv_index(i,j,k)]/coord->GetDeltaZ_down(k);
	    vC1k = vc[conv_index(i,j,k-1)]/coord->GetDeltaZ_down(k);
	  }
	  else if ( coord->GetZ().at(k) < 0) {
	    vCk  = vc[conv_index(i,j,k)]/coord->GetDeltaZ_up(k);
	    vC1k = 0.0;
	    vCk1 = vc[conv_index(i,j,k+1)]/coord->GetDeltaZ_up(k);
	  }
	  else {
	    vCk  = vc[conv_index(i,j,k)]*2/(coord->GetDeltaZ_up(k) + coord->GetDeltaZ_down(k));
	    vC1k = -0.5*vc[conv_index(i,j,k-1)]/coord->GetDeltaZ_down(k);
	    vCk1 = -0.5*vc[conv_index(i,j,k+1)]/coord->GetDeltaZ_up(k);
	  }

	  CNconv_alpha1_z.push_back( vC1k );
	  CNconv_alpha2_z.push_back( vCk );
	  CNconv_alpha3_z.push_back( vCk1 );

	}
      }            
    }
  }

  if(in->write_flag) datafile.close();

}

double TConvectionVelocity::GetProfile(double x, double y, double zeta, TSource* SourceTerm) {
	
  //double nrn_sn  = 50.0*(0.79*exp(-pow(zobs/0.212,2.))+0.21*exp(-pow(zobs/0.636,2.)))*exp(-(robs*robs-robs*robs)/(6.8*6.8)) + 7.3*exp(-(robs-robs)/4.5-fabs(zobs)/0.325);
  double radial = 1;
  double result = 1.;
  //    cout << " nrn_sn in GetProfile ->" << nrn_sn << endl;
    
  switch(set_profile_conv) {
  case Constant :
    return 1.0;
    break;
            
    //case Exp :
    //    return exp(fabs(zeta)/zt);
    //    break;
    /*
      case Blasi :
      if (fabs(zeta) < 1.) return exp(fabs(zeta)/0.5);
      else return exp(1.0/0.5)*exp(fabs(zeta)/zt)/exp(1.0/zt);
      break;
             
      case Expr :
      return exp(fabs(zeta/zt))*(0.5*(tanh((radius-3.0)/0.25)+1.001));///cosh((radius-r0)/rd);
      break;
    */
  case Radial : //MW130621: Qtau is doing the same as Radial from Convection_new, just not evaluating it again, but instead taking the value of the SourceTerm.

    //MW130711: z-dependency should not be accounted for!
    zeta = 0;

    //MW130711: for comparison with DRAGON-KIT, don't normalize.
    radial = SourceTerm->GetSource(x,y,zeta);
    if (radial < conv_threshold) radial=conv_threshold;
    return pow(radial, conv_index_radial);

  case Qtau :
    /*
      if (radius > 3.7) radial = 50.0*(0.79*exp(-pow(zobs/0.212,2.))+0.21*exp(-pow(zobs/0.636,2.)))*exp(-((radius)*(radius)-robs*robs)/(6.8*6.8)) + 7.3*exp(-(radius-robs)/4.5-fabs(zobs)/0.325);
      else radial = 177.5*(0.79*exp(-pow(zobs/0.212,2.))+0.21*exp(-pow(zobs/0.636,2.)))*exp(-pow((radius-3.7)/2.1, 2.)) + 7.3*exp(-(radius-robs)/4.5-fabs(zobs)/0.325);
    */
    radial = SourceTerm->GetSource(x,y,zeta);
    radial /= nrn_sn;
            
    if (radial < conv_threshold) radial=conv_threshold;

    return pow(radial, conv_index_radial);
    break;
    /*
      case ExpRadial :
      {
      double rbulge = 3.;
      radial = (1. * pow ( radius/robs , 1.25 ) * exp ( -3.56*(radius-robs)/robs  ));
      if (radius <= rbulge)
      radial = (1. * pow ( rbulge/robs , 1.25 ) * exp ( -3.56*(rbulge-robs)/robs  ));
      result = exp(fabs(zeta/zt)) * (pow(radial, index_radial) + 0.01*pow(radial, -index_radial));
      //if (result < .2) result = .2;
      return result;
      break;
      }
    */
  default :
    return -1;
  }
    
}


TReaccelerationCoefficient::TReaccelerationCoefficient(vector<double> pp, TDiffusionCoefficient* dperp, TGeometry* geom, Input* in) {
    
  //sk: consider break in diffusion 
  double a;
  double Dpp_constant[pp.size()];

  for (unsigned int i = 0; i < pp.size(); ++i){
    if(pp[i] < in->rho_b){
      if (in->VariableDelta == true)
	a = (in->DiffT == Anisotropic) ? in->DeltaPar : in->delta_B + in->delta_A*in->robs;
      else
	a = (in->DiffT == Anisotropic) ? in->DeltaPar : dperp->GetDelta(); //MW130711: integrate Anisotropic Diffusion
    }
    else
      a = (in->DiffT == Anisotropic) ? in->DeltaPar : dperp->GetDelta_h();

    if (in->VariableDelta == true)
      Dpp_constant[i] = 1.0;
    else
      Dpp_constant[i]= 1.0/(a*(4.-a)*(4.-a*a));       // Ptuskin-2003
    if(in->diff_reacc == 1) Dpp_constant[i] *= 4.0/3.0; // Seo & Ptuskin
  }
    
  vector<double> DiffSpectrum = dperp->GetSpectrum();
  for (unsigned int i = 0; i < pp.size(); ++i) sp.push_back(Dpp_constant[i]*in->vAlfven*in->vAlfven*pp[i]*pp[i]/DiffSpectrum[i]);

  //     for (unsigned int i = 0; i < pp.size(); ++i) cout << "[MW-DEBUG REACC G] " << " " << i << " " << Dpp_constant[i] << " " << in->vAlfven << " " << pp[i] << " " << DiffSpectrum[i] << " | " << in->DiffT << " " << in->DeltaPar << " " << " " << dperp->GetDelta() << " " << dperp->GetDelta_h() << " " << endl;


  dimr = dperp->GetDimR();
  dimx = dperp->GetDimX();
  dimy = dperp->GetDimY();
  dimz = dperp->GetDimZ();
  vector<double> DiffProfile = (in->DiffT == Anisotropic) ? dperp->GetDPar() : dperp->GetDiffusionCoefficient();
    
  if(dperp->GetCoord()->GetType() == "3D")
    {
      int index = 0;
      for (vector<double>::iterator i = DiffProfile.begin(); i != DiffProfile.end(); ++i)
        {
	  int ix = dperp->GetCoord()->GetXFromIndexD_3D(index);
	  int iy = dperp->GetCoord()->GetYFromIndexD_3D(index);
	  int iz = dperp->GetCoord()->GetZFromIndexD_3D(index);

	  double xx = dperp->GetCoord()->GetX()[ix];
	  double yy = dperp->GetCoord()->GetY()[iy];
	  double zz = dperp->GetCoord()->GetZ()[iz];
	  
	  double reacc_spatial = 1.0/(*i);
	  
	  double spiral_factor_dperp = max( min( pow(geom->GetPattern(ix,iy,iz), in->SA_diff), in->SA_cut_diff), 1./in->SA_cut_diff );
	  double spiral_factor_dpp = max( min( spiral_factor_dperp * pow(geom->GetPattern(ix,iy,iz), 2*in->SA_vA), in->SA_cut_vA), 1./in->SA_cut_vA );
	  
	  reacc_spatial *= spiral_factor_dpp; //mw 130422
	  if (dperp->GetCoord()->IsInLocalBubble(xx,yy,zz)) reacc_spatial *= pow(in->LB_vA, 2*dperp->GetCoord()->IsInLocalBubble(xx,yy,zz)) * pow(in->LB_diff, dperp->GetCoord()->IsInLocalBubble(xx,yy,zz));
	  dpp.push_back(reacc_spatial);
	  
	  index++;
        }
    }
  else
    {
      for (vector<double>::iterator i = DiffProfile.begin(); i != DiffProfile.end(); ++i)
        {
	  double reacc_spatial = 1.0/(*i);
	  dpp.push_back(reacc_spatial);
        }
    }


vector<double> r_vec = dperp->GetCoord()->GetR();
if (in->VariableDelta == true){
  for (int ir=0; ir<dimr; ir++) {
    for (unsigned int k = 0; k < dimz; ++k) {
      double izr = index(ir,k);
      double aVar = in->delta_A*r_vec[ir]+in->delta_B;
      dpp[izr] *= 1.0/(aVar*(4.-aVar)*(4.-aVar*aVar));
    }
  }
  
 }

}

//********************************************************************************************************************************************************************************
  //************************************************************* THE GALAXY CONSTRUCTOR *******************************************************************************************
    //********************************************************************************************************************************************************************************


Galaxy::Galaxy(Input* in, TNucleiList* l) {

    if (in->feedback > 0) cout << "Welcome to the Galaxy constructor " << endl;
         
    if (in == NULL) {
      cerr << "No Input specified!" << endl;
      return ;
    }
    inp = in;
    
    ifstream infile(in->sourcedata.c_str(),ios::in);

    cout << " %%% in->sourcedata.c_str() " << in->sourcedata.c_str() << endl;

    if (!infile.is_open()) {
      cerr << "File " << in->sourcedata << " does not exist. Using config_files/template.source.param!" << endl;
      infile.open("config_files/template.source.param",ios::in);
      if (!infile.is_open()) {
	cerr << "config_files/template.source.param does not exist, either! Exiting." << endl;
	exit(NOSOURCEDATA);
      }

    }
    
    //new implementation DG29.09.2013
    int particle_ID;
    map<int,double> abundances_map;
    map<int, vector<double> > inj_indexes;
    map<int, vector<double> > break_positions;
    //the code reads the .source.param. First column: nucleus ID; second column: abundance; other columns: inj_slope - break rigidity - inj slope - break rigidity - (...) - highest energy inj_slope; the number of breaks is arbitrary
    //reads the .source.param to a table	
    typedef vector<double> Row;
    vector<Row> table;	 
    while (infile) {
      cout << "reading line" << endl;
      string line; getline(infile, line);
      istringstream temp_string(line);
      Row row;
      while (temp_string) {
	double data;
	temp_string >> data;
	row.push_back(data);
	cout << data << endl;
      }
      table.push_back(row);
    }
    //fills abundances_map inj_indexes and break_positions for each nucleus i using the table
    if (in->feedback > 0) cout << "Reading .source.param table with " << table.size() << " rows" << endl;
    for (unsigned i=0; i<table.size(); i++) {
      Row row; row = table[i];
      if (in->feedback > 0) cout << "Reading line in .source.param of size: " << row.size() << endl;
      int nid = 0;	
      for (unsigned j=0; j<row.size()-1; j++) {
	if (in->feedback > 0) cout << table[i][j] << ", ";
	if (j==0) {
	  nid = table[i][0]; cout << "$$$ " << nid << endl;
	}
	if (j==1) abundances_map[nid] = table[i][1];
	if (j>0 && j%2==0)
	  inj_indexes[nid].push_back(table[i][j]);	
	if (j>1 && j%2!=0)
	  break_positions[nid].push_back(table[i][j]);	
      }
      if (in->feedback >0) cout << endl;
      cout << nid << " <- nucleus ID | abundance -> " << abundances_map[nid] << endl;
    }    

    vector<int> list = l->GetList();

    if (in->feedback >0) cout << "Setting abundances, inj slopes and break positions for each nucleus in the list" << endl;
    for (vector<int>::iterator it_current_nucleus = list.begin(); it_current_nucleus != list.end(); ++it_current_nucleus) {
      // loop over NucleiList from nucleilist.cc

      map<int,double>::iterator it_current_nucleus_abundance = abundances_map.find(*it_current_nucleus);

      _fSourceAbundances[*it_current_nucleus] = (*it_current_nucleus_abundance).second;
      if (in->feedback >0) cout << " ** Nucleus id -> " << (*it_current_nucleus_abundance).first << ". Abundance found in .source.param -> " << _fSourceAbundances[*it_current_nucleus] << endl;


      if (it_current_nucleus_abundance != abundances_map.end()) {

	_fSourceAbundances[*it_current_nucleus] = (*it_current_nucleus_abundance).second;
	if (in->feedback >0) cout << "Nucleus id -> " << (*it_current_nucleus_abundance).first << ". Abundance found in .source.param -> " << _fSourceAbundances[*it_current_nucleus] << endl;

	if (inp->UseInjectionIndexAllNuclei == false) {	
	  if (in->feedback > 0) cout << "Break positions and slopes are NOT specified in the XML and are taken from .source.param file!" << endl;
	  _fInjSpectrum_rho[*it_current_nucleus]   = break_positions[*it_current_nucleus];	
	  _fInjSpectrum_alpha[*it_current_nucleus] = inj_indexes[*it_current_nucleus];
  		

	}
	else {

	  if (in->feedback >0) cout << "Break positions and slopes are taken from xml file!" << endl;

	  if (in->feedback >0) cout << "Number of slopes: " << inp->inp_inj_indexes.size() << endl;
	  if (in->feedback >0) cout << "Number of breaks: " << inp->inp_break_positions.size() << endl;

	  for (int j=0; j<inp->inp_break_positions.size(); j++)
	    _fInjSpectrum_rho[*it_current_nucleus] = inp->inp_break_positions;
	  for (int j=0; j<inp->inp_inj_indexes.size(); j++)
	    _fInjSpectrum_alpha[*it_current_nucleus] = inp->inp_inj_indexes;

	}	
      }
      else {

	if (in->feedback >0) cout << "Nucleus id -> " << (*it_current_nucleus_abundance).first << ". Abundance NOT found in .source.param! " << endl;
 
	_fSourceAbundances[*it_current_nucleus] = 0.0;
	_fInjSpectrum_rho[*it_current_nucleus].push_back(1.);
	_fInjSpectrum_alpha[*it_current_nucleus].push_back(0.);
	_fInjSpectrum_alpha[*it_current_nucleus].push_back(0.);
      }
    }
    
    _fSourceAbundances[-1000] = 1.0;
    _fSourceAbundances[1000] = 0.0;
    _fSourceAbundances[-999] = 0.0;
    _fSourceAbundances[-1998] = 0.0;

    //if (inp->TESTMODE == false) TESTMODE = false;
    //else TESTMODE = true;	  

    TESTMODE = inp->TESTMODE;

    if (inp->MOVING == false) MOVING = false;
    else MOVING = true;		
    //if (feedback >0) cout << "Source is moving? " << MOVING << endl; 		

    if (MOVING)  {

      source_x0 = inp->source_x0;
      source_y0 = inp->source_y0;
      source_z0 = inp->source_z0;
      source_vx = inp->source_vx;
      source_vy = inp->source_vy;
      source_vz = inp->source_vz;
    } 

    if (inp->MOVING_CLUMP == false) MOVING_CLUMP = false;
    else MOVING_CLUMP = true;		

    if (MOVING_CLUMP) { 

      clump_x0 = inp->clump_x0;
      clump_y0 = inp->clump_y0;
      clump_z0 = inp->clump_z0;
      clump_vx = inp->clump_vx;
      clump_vy = inp->clump_vy;
      clump_vz = inp->clump_vz;
      clump_deltat = inp->clump_deltat;
    } 

    //#ifdef DEBUG
   
    //for (vector<int>::iterator it = list.begin(); it != list.end(); ++it) cout <<"injection " <<  *it << " " << _fSourceAbundances[*it] << " " << _fInjSpectrum_rho_0[*it] << " " << _fInjSpectrum_rho_1[*it] << " " << _fInjSpectrum_rho_2[*it] << " " << _fInjSpectrum_alpha_0[*it] << " " << _fInjSpectrum_alpha_1[*it] << " " << _fInjSpectrum_alpha_2[*it] << " " << _fInjSpectrum_alpha_3[*it] <<endl;
    
    //#endif

    if (in->feedback > 0) cout << "Preparing the grid... " << endl; 

    if (inp->gridtype == "2D") _fCoordinates =  new TGrid2D(in);
    else _fCoordinates = new TGrid3D(in);
    if (in->feedback > 0) cout << "Grid done" << endl;
    
    if (in->feedback > 0) cout << "Preparing the geometry... " << endl;
    if (inp->SA_type == "None") _fGeometry = new TUniformGeometry(_fCoordinates, in);
    else _fGeometry = new TSpiralGeometry(_fCoordinates, in);
    if (in->feedback > 0) cout << "Geometry done" << endl;

    if (in->feedback > 0) cout << "Preparing the gas... " << endl;
    _fGas.push_back(new TH2Gas(_fCoordinates, in, _fGeometry));
    _fGas.push_back(new THIGas(_fCoordinates, in, _fGeometry));
    _fGas.push_back(new THIIGas(_fCoordinates, in, _fGeometry));

    _fTotalGas = new TGas(_fCoordinates, in);
    if (in->feedback > 0) cout << "[MW] Sum of TotalGas vector before filling(should be zero): " << _fTotalGas->GetTotalContent() << endl;
    //MW 130429: Construct TotalGas from other components
    *_fTotalGas += *_fGas[0];
    *_fTotalGas += *_fGas[1];
    *_fTotalGas += *_fGas[2];
    if (in->feedback > 0) cout << "[MW] Sum of TotalGas vector after filling (should be " << _fGas[0]->GetTotalContent()+_fGas[1]->GetTotalContent()+_fGas[2]->GetTotalContent() << "): " << _fTotalGas->GetTotalContent() << endl;

    // SETTING GAS ABUNDANCES

    _fGasAbundances[1001] = 1.;
    _fGasAbundances[2004] = 0.11;
    _fGasAbundances[6012] = 0.05;

    if (in->feedback > 0) cout << "Gas done" << endl;
    
    _fDMSource = new TDMSource(_fCoordinates, in);
    if (in->feedback > 0) cout << "Creating astrophysical source... " << endl;
    _fSource   = new TAstrophysicalSource(_fCoordinates, in, _fGeometry, in->SNR_model);
    _fSourceExtra = (in->prop_extracomp) ? new TAstrophysicalSource(_fCoordinates, in, _fGeometry, in->SNR_model_Extra) : NULL; //CAREFUL! hard coded model = model_extra
    if (in->feedback > 0) cout << "Sources done " << endl;
    
    switch (in->BM) {
    case Pshirkov:
      _fB = new TPshirkovField(in->B0disk, 5., 1., 10., in->B0halo, 8., 1.3, in->B0turb, 8.5, in->zt, in->robs, _fCoordinates, _fGeometry);
      break;
    case Farrar:
      _fB = new TFarrarField(in->betaFarrar, _fCoordinates, _fGeometry);
      break;
    case Uniform:
      _fB = new TUniformField(in->B0turb, _fCoordinates, _fGeometry);
      break;
    case Simple:
      _fB = new TSimpleField(in->b0, in->robs, _fCoordinates, _fGeometry);
      break;
    case ToyModel:
      cout << "ToyModel mag field was specified!" <<endl;
      _fB = new ToyModelField(in->bx, in->by, in->bz, in->bturb, _fCoordinates, _fGeometry); 
      break;
    default :
      _fB = NULL;
    }
    
    _fDperp = NULL;
    _fDpp = NULL;
    _fDperpEl = NULL;
    _fDppEl = NULL;
    
    if (in->DiffT != Anisotropic) {
        
      if (in->gridtype == "2D") _fDperp = new TDiffusionCoefficient2D(_fCoordinates, in, _fSource, _fB);
      else _fDperp = new TDiffusionCoefficient3D(_fCoordinates, in, _fSource, _fB, _fGeometry, 0, 0);
        
      if (in->feedback > 0) cout << "Diffusion coefficient done " << endl;
        
      _fDpp = (in->REACC) ? new TReaccelerationCoefficient(_fCoordinates->GetMomentum(), _fDperp, _fGeometry, in) : NULL;
    }
    
    _fVC = (in->CONVECTION) ? new TConvectionVelocity(_fCoordinates, _fGeometry, in, _fSource) : NULL;
    
    if (in->prop_lep || in->prop_extracomp || in->prop_DMel) {
        
      if (in->DiffT != Anisotropic) {
            
	if (in->gridtype == "2D") _fDperpEl = new TDiffusionCoefficient2D(_fCoordinates, in, _fSource, _fB, 1);
	else _fDperpEl = new TDiffusionCoefficient3D(_fCoordinates, in, _fSource, _fB, _fGeometry, 0, 0, 1);
            
	_fDppEl = (in->REACC) ? new TReaccelerationCoefficient(_fCoordinates->GetMomentumEl(), _fDperpEl, _fGeometry, in) : NULL;
      }
      _fISRF = new TISRF(_fCoordinates, ISRFfile, _fGeometry, in);
    }
    else {
      _fDperpEl = NULL;
      _fDppEl = NULL;
      //_fB = NULL;
      _fISRF = NULL;
    }
    }

void Galaxy::Delete() {
  if (_fTotalGas) delete _fTotalGas;
  for (vector<TGas*>::iterator i = _fGas.begin(); i != _fGas.end(); ++i) {
    if (*i) delete *i;
  }
  _fGas.clear();
  if (_fISRF) delete _fISRF;
  if (_fB) delete _fB;
  if (_fGeometry) delete _fGeometry;
}

Galaxy::~Galaxy() {
  if (_fCoordinates) delete _fCoordinates;
  if (_fSource) delete _fSource;
  if (_fSourceExtra) delete _fSourceExtra;
  if (_fDMSource) delete _fDMSource;
  if (_fDperp) delete _fDperp;
  if (_fDpp) delete _fDpp;
  if (_fDperpEl) delete _fDperpEl;
  if (_fDppEl) delete _fDppEl;
  if (_fVC) delete _fVC;
  if (_fISRF) delete _fISRF;
  if (_fB) delete _fB;
  if (_fTotalGas) delete _fTotalGas;
  for (vector<TGas*>::iterator i = _fGas.begin(); i != _fGas.end(); ++i) {
    if (*i) delete *i;
  }
  _fGas.clear();

  /*_fInjSpectrum_rho_0.clear();
    _fInjSpectrum_rho_1.clear();
    _fInjSpectrum_rho_2.clear();

    _fInjSpectrum_alpha_0.clear();
    _fInjSpectrum_alpha_1.clear();
    _fInjSpectrum_alpha_2.clear();
    _fInjSpectrum_alpha_3.clear();*/

  _fSourceAbundances.clear();
  _fInjSpectrum_rho.clear();
  _fInjSpectrum_alpha.clear();
}
