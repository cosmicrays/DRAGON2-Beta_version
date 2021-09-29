/**
 * @file    dragon.cc
 * @author  Daniele Gaggero, Carmelo Evoli, Pedro De La Torre Luque
 * @email   daniele.gaggero@gmail.com
 * @brief   Implementation of the DRAGON class. See the .h file
 */


#include "diffusion.h"

#include "geometry.h"
#include "grid.h"
#include "input.h"
#include "sources.h"
#include "errorcode.h"

#include "bfield.h"

using namespace std;

TDiffusionCoefficient::TDiffusionCoefficient(TGrid* Coord, Input* in, TSource* SourceTerm, TBField* bf, bool El) {
    
  _fCoordinates = Coord; 
    
  zt = in->zt;
  nrn_sn = SourceTerm->GetSource(8.5,0.0,0.0);
  index_radial = in->index_radial;
  set_profile = in->set_profile;

  delta = in->delta;
  delta_h = in->delta_h;
  rho_b = in->rho_b;
}

double TDiffusionCoefficient::GetProfile(double x, double y, double zeta, TSource* SourceTerm) {
	
  //double nrn_sn  = 50.0*(0.79*exp(-pow(zobs/0.212,2.))+0.21*exp(-pow(zobs/0.636,2.)))*exp(-(robs*robs-robs*robs)/(6.8*6.8)) + 7.3*exp(-(robs-robs)/4.5-fabs(zobs)/0.325);
  double radial = 1;
  double result = 1.;
  //    cout << " nrn_sn in GetProfile ->" << nrn_sn << endl;
    
  switch(set_profile) {
  case Constant :
    return 1.0;
    break;
            
  case Exp :
    return exp(fabs(zeta)/zt);
    break;
    /*
      case Blasi :
      if (fabs(zeta) < 1.) return exp(fabs(zeta)/0.5);
      else return exp(1.0/0.5)*exp(fabs(zeta)/zt)/exp(1.0/zt);
      break;
             
      case Expr :
      return exp(fabs(zeta/zt))*(0.5*(tanh((radius-3.0)/0.25)+1.001));///cosh((radius-r0)/rd);
      break;
    */
  case Qtau :
    /*
      if (radius > 3.7) radial = 50.0*(0.79*exp(-pow(zobs/0.212,2.))+0.21*exp(-pow(zobs/0.636,2.)))*exp(-((radius)*(radius)-robs*robs)/(6.8*6.8)) + 7.3*exp(-(radius-robs)/4.5-fabs(zobs)/0.325);
      else radial = 177.5*(0.79*exp(-pow(zobs/0.212,2.))+0.21*exp(-pow(zobs/0.636,2.)))*exp(-pow((radius-3.7)/2.1, 2.)) + 7.3*exp(-(radius-robs)/4.5-fabs(zobs)/0.325);
    */
    radial = SourceTerm->GetSource(x,y,zeta);
    radial /= nrn_sn;
    return pow(radial, index_radial);
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
  case Test :
    //if (y < 0.) return 1.;
    //else 
    return 1.;
    //	return exp(-y/2.0);
    break;		
  default :
    return -1;
  }
    
}


double TDiffusionCoefficient::GetXDerivative(double x, double y, double z, TSource* SourceTerm) {
    
  double delta_x = 0.00002; //Coord->GetDeltaX();
  double xup = x + 0.5*delta_x;
  double xdo = x - 0.5*delta_x;
  return ( GetProfile(xup,y,z, SourceTerm) - GetProfile(xdo,y,z, SourceTerm) )/delta_x ;
}

double TDiffusionCoefficient::GetYDerivative(double x, double y, double z, TSource* SourceTerm) {
    
  double delta_y = 0.00002; //Coord->GetDeltaY();
  double yup = y + 0.5*delta_y;
  double ydo = y - 0.5*delta_y;
  return ( GetProfile(x,yup,z, SourceTerm) - GetProfile(x,ydo,z, SourceTerm) )/delta_y;
}

double TDiffusionCoefficient::GetZDerivative(double x, double y, double z, TSource* SourceTerm) {
    
  double delta_z = 0.00002; //Coord->GetDeltaZ();
  double zup = z + 0.5*delta_z;
  double zdo = z - 0.5*delta_z;
  return ( GetProfile(x,y,zup, SourceTerm) - GetProfile(x,y,zdo, SourceTerm) )/delta_z;
}

// 2D

TDiffusionCoefficient2D::TDiffusionCoefficient2D(TGrid* Coord, Input* in, TSource* SourceTerm, TBField* bf, bool El) : TDiffusionCoefficient(Coord, in, SourceTerm, bf, El){
    
  delta = in->delta;

  vector<double> pp;
  vector<double> beta;
  vector<double> Ek;
  if (El) {
    pp = Coord->GetMomentumEl();
    beta = Coord->GetBetaEl();
    Ek = Coord->GetEk();
  }
  else {
    pp = Coord->GetMomentum();
    beta = Coord->GetBeta();
    Ek = Coord->GetEk();
  }
  //cout << in->D0 << " " << in->etaT << endl;
    
  if (in->VariableDelta == false) {

    delta = in->delta;

    cout << "building diff. coeff. with constant delta" << endl;
    for (int i = 0; i < pp.size(); ++i) {

      if ( delta > 0.0 ){
	
	if(in->D_ref_rig<=rho_b){
	  if(pp[i]<=rho_b) sp.push_back(in->D0*pow(beta[i], in->etaT)*pow(pp[i]/in->D_ref_rig, delta));
	  if(pp[i]>rho_b)  sp.push_back(in->D0*pow(beta[i], in->etaT)*pow(pp[i]/rho_b, delta_h)*pow(rho_b/in->D_ref_rig, delta));
	}
	if(in->D_ref_rig>rho_b){
	  if(pp[i]<=rho_b) sp.push_back(in->D0*pow(beta[i], in->etaT)*pow(rho_b/in->D_ref_rig, delta_h)*pow(pp[i]/rho_b, delta));
	  if(pp[i]>rho_b)  sp.push_back(in->D0*pow(beta[i], in->etaT)*pow(pp[i]/in->D_ref_rig, delta_h));
	}
      } else {
	const double lgEk = log10(Ek[i]);
	sp.push_back(in->D0*pow(10.0,-0.8990+1.1781*lgEk-0.1800*lgEk*lgEk+0.0141*lgEk*lgEk*lgEk));
      }
      
    }
  } else {

    double delta_A = in->delta_A;
    double delta_B = in->delta_B;

    cout << "building diff. coeff. with variable delta(R) = Ar + B " << endl;
    cout << delta_A << "r + " << delta_B << endl;

    vector<double> r_vec = Coord->GetR();
    for (int i = 0; i < r_vec.size(); i++)
      cout << r_vec[i] << " " ;
    cout << endl;
    vector<double> z_vec = Coord->GetZ();
    for (int i = 0; i < z_vec.size(); i++)
      cout << z_vec[i] << " " ;
    cout << endl;
    dimr = r_vec.size();
    cout << dimr << endl;
    dimz = z_vec.size();
    cout << dimz << endl;
    unsigned int irsun = (unsigned int) ((in->robs-r_vec.front())/(r_vec.back()-r_vec.front())*(double)(dimr-1));
    unsigned int izsun = (unsigned int) ((in->zobs-z_vec.front())/(z_vec.back()-z_vec.front())*(double)(dimz-1));

    for (int ip = 0; ip < pp.size(); ++ip) {
      //cout << ip << endl;
      for (int iz = 0; iz < z_vec.size(); iz++) {
	for (int ir = 0; ir < r_vec.size(); ir++) {

	  double current_delta = delta_A*r_vec[ir]+delta_B;

	  if (r_vec[ir] > in->diff_threshold)
	    current_delta = delta_A*in->diff_threshold + delta_B;

	  double current_delta_hi = current_delta;

	  if(in->D_ref_rig<=rho_b){
	    if(pp[ip]<=rho_b) spectrum_extended.push_back(in->D0*pow(beta[ip], in->etaT)*pow(pp[ip]/in->D_ref_rig, current_delta));
	    if(pp[ip]>rho_b)  spectrum_extended.push_back(in->D0*pow(beta[ip], in->etaT)*pow(pp[ip]/rho_b, current_delta_hi)*pow(rho_b/in->D_ref_rig, current_delta));
	  }
	  if(in->D_ref_rig>rho_b){
	    if(pp[ip]<=rho_b) spectrum_extended.push_back(in->D0*pow(beta[ip], in->etaT)*pow(rho_b/in->D_ref_rig, current_delta_hi)*pow(pp[ip]/rho_b, current_delta));
	    if(pp[ip]>rho_b)  spectrum_extended.push_back(in->D0*pow(beta[ip], in->etaT)*pow(pp[ip]/in->D_ref_rig, current_delta_hi));
	  }	

	  if (ir == irsun && iz == izsun) {
	    if(in->D_ref_rig<=rho_b){
	      if(pp[ip]<=rho_b) sp.push_back(in->D0*pow(beta[ip], in->etaT)*pow(pp[ip]/in->D_ref_rig, current_delta));
	      if(pp[ip]>rho_b)  sp.push_back(in->D0*pow(beta[ip], in->etaT)*pow(pp[ip]/rho_b, current_delta_hi)*pow(rho_b/in->D_ref_rig, current_delta));
	    }
	    if(in->D_ref_rig>rho_b){
	      if(pp[ip]<=rho_b) sp.push_back(in->D0*pow(beta[ip], in->etaT)*pow(rho_b/in->D_ref_rig, current_delta_hi)*pow(pp[ip]/rho_b, current_delta));
	      if(pp[ip]>rho_b)  sp.push_back(in->D0*pow(beta[ip], in->etaT)*pow(pp[ip]/in->D_ref_rig, current_delta_hi));
	    }
	  }		
	}
      }
    }
  }

  // sp.push_back(in->D0*pow(beta[i], in->etaT)*pow(pp[i]/D_ref_rig, delta)); //MW130621: replace
    
  vector<double> r = Coord->GetR();
  vector<double> z = Coord->GetZ();
    
  dimr = r.size();
  dimz = z.size();
    
  //MW130705
  vector<double> dr_up = Coord->GetDeltaR_up();
  vector<double> dz_up = Coord->GetDeltaZ_up();
  vector<double> dr_down = Coord->GetDeltaR_down();
  vector<double> dz_down = Coord->GetDeltaZ_down();

  for (unsigned int i = 0; i < dimr; ++i) {
    double radius = r[i];
    double Deltar = Coord->GetDeltaR(i);
        
    for (unsigned int j = 0; j < dimz; ++j) {
      double zeta = z[j];
      double Deltaz = Coord->GetDeltaZ(j);

      dperp.push_back(GetProfile(radius, 0, zeta, SourceTerm));
      //phi.push_back(0.5/Deltar*(dperp.back()/max(u,radius)));  // +GetRDerivative(radius, 0, zeta, SourceTerm)  // Variable Phi not used in this implementation
      //             psi.push_back(0.5/Deltaz*GetZDerivative(radius, 0, zeta, SourceTerm));
    }
  }


  //MW130625: order correct?
  double sp_ = 0;
  double sp_rup_ = 0;
  double sp_rdown_ = 0;
  
  for (unsigned int i = 0; i < dimr; ++i) {
    double dr_central = 0.5 * (dr_up[i]+dr_down[i]);
            
    for (unsigned int k = 0; k < dimz; ++k) {
      double dz_central = 0.5 * (dz_up[k]+dz_down[k]);
                    
      //MW130624: Trying to gain speed at the cost of memory...
      //Tests gave me a speed up of 45% with no considerable memory loss!
      double indspat       = index(i,k);
      double indspat_rup   = index(i+1,k);
      double indspat_rdown = index(i-1,k);
      double indspat_zup   = index(i,k+1);
      double indspat_zdown = index(i,k-1);

      for  (int ip = 0; ip < pp.size(); ip++) {

	if (in->VariableDelta == false){
	  sp_ = sp[ip];
	  sp_rup_ = sp_;
	  sp_rdown_ = sp_;
	  sp_zup_ = sp_;
	  sp_zdown_ = sp_;
	}
	
	else{
	  sp_ = spectrum_extended[index_rzp(ip,k,i)];
	    
	    if (i == dimr-1)
	      sp_rup_ = sp_;
	    else
	      sp_rup_ = spectrum_extended[index_rzp(ip,k,i+1)];
	    
	    if (i == 0)
	      sp_rdown_ = sp_;
	    else
	      sp_rdown_ = spectrum_extended[index_rzp(ip,k,i-1)];
	    
	    if (k == dimz-1)
	      sp_zup_ = sp_;
	    else
	      sp_zup_ = spectrum_extended[index_rzp(ip,k+1,i)];
	    
	    if (k == 0)
	      sp_zdown_ = sp_;
	    else
	      sp_zdown_ = spectrum_extended[index_rzp(ip,k-1,i)];
	}
	
	
	double D       = dperp[indspat]*sp_;
	double D_rup   = dperp[indspat_rup]*sp_rup_;
	double D_rdown = dperp[indspat_rdown]*sp_rdown_;
	double D_zup   = dperp[indspat_zup]*sp_zup_;
	double D_zdown = dperp[indspat_zdown]*sp_zdown_;

	CNdiff_alpha1_r.push_back(D/(dr_central*dr_down[i]) - (D_rup-D_rdown)/(4*dr_central*dr_central) - 0.5*D/dr_central/max(u,r[i]) ); 
	CNdiff_alpha2_r.push_back(D/(dr_central*dr_up[i]) + D/(dr_down[i]*dr_central));
	double temp2 = 0.;
	if (i > 1) temp2 += (0.5/dr_central * D/max(u,r[i]));
	CNdiff_alpha3_r.push_back((D_rup-D_rdown)/(4*dr_central*dr_central) + D/(dr_up[i]*dr_central) + temp2); 
				  
	CNdiff_alpha1_z.push_back(D/(dz_central*dz_down[k]) - (D_zup-D_zdown)/(4*dz_central*dz_central));
	CNdiff_alpha2_z.push_back(D/(dz_central*dz_up[k]) + D/(dz_down[k]*dz_central));
	CNdiff_alpha3_z.push_back((D_zup-D_zdown)/(4*dz_central*dz_central) + D/(dz_up[k]*dz_central));

      }
    }
  }
 
}



// 3D

TDiffusionCoefficient3D::TDiffusionCoefficient3D(TGrid* Coord, Input* in, TSource* SourceTerm, TBField* bf, TGeometry* geom, double A, double Z, bool El) : TDiffusionCoefficient(Coord, in, SourceTerm, bf, El) {
    
  diffT = in->DiffT;
    
  vector<double> pp;
  vector<double> beta;
  if (El) {
    pp = Coord->GetMomentumEl();
    beta = Coord->GetBetaEl();
  }
  else {
    pp = Coord->GetMomentum();
    beta = Coord->GetBeta();
  }
  //cout << in->D0 << " " << in->etaT << endl;
    
  vector<double> x = Coord->GetX();
  vector<double> y = Coord->GetY();
  vector<double> z = Coord->GetZ();
    
  dimx = x.size();
  dimy = y.size();
  dimz = z.size();
  dimE = pp.size();

  //MW130624
  vector<double> dx_up = Coord->GetDeltaX_up();
  vector<double> dy_up = Coord->GetDeltaY_up();
  vector<double> dz_up = Coord->GetDeltaZ_up();
  vector<double> dx_down = Coord->GetDeltaX_down();
  vector<double> dy_down = Coord->GetDeltaY_down();
  vector<double> dz_down = Coord->GetDeltaZ_down();

  if (diffT == Isotropic) {
    
    for ( int i = 0; i < pp.size(); ++i ){
      double rho_n = in->D_ref_rig;
      
      //sp.push_back(in->D0*pow(beta[i], in->etaT)*pow(pp[i]/D_ref_rig, delta)); //MW130621
      if(rho_n<=rho_b){
	if(pp[i]<=rho_b) sp.push_back(in->D0*pow(beta[i], in->etaT)*pow(pp[i]/rho_n, delta));
	if(pp[i]>rho_b)  sp.push_back(in->D0*pow(beta[i], in->etaT)*pow(pp[i]/rho_b, delta_h)*pow(rho_b/rho_n, delta));
      }
      if(rho_n>rho_b){
	if(pp[i]<=rho_b) sp.push_back(in->D0*pow(beta[i], in->etaT)*pow(rho_b/rho_n, delta_h)*pow(pp[i]/rho_b, delta));
	if(pp[i]>rho_b)  sp.push_back(in->D0*pow(beta[i], in->etaT)*pow(pp[i]/rho_n, delta_h));
      }
    }
    
    for (unsigned int i = 0; i < dimx; ++i) {
      for (unsigned int j = 0; j < dimy; ++j) {
	//double radius = sqrt(x[i]*x[i]+y[j]*y[j]);
	for (unsigned int k = 0; k < dimz; ++k) {
	  //double zeta = z[j];
	  double dperp_profile = GetProfile(x[i], y[j], z[k], SourceTerm);
	  
	  //mw, 130326, 130415
	  double spiral_factor_dperp = max( min( pow(geom->GetPattern(i,j,k), in->SA_diff), in->SA_cut_diff), 1./in->SA_cut_diff );
	  dperp.push_back( dperp_profile * pow( in->LB_diff, Coord->IsInLocalBubble(x[i],y[j],z[k]) ) * spiral_factor_dperp);
	  //                     phix.push_back(GetXDerivative(x[i], y[j], z[k], SourceTerm));
	  //                     phiy.push_back(GetYDerivative(x[i], y[j], z[k], SourceTerm));
	  //                     phiz.push_back(GetZDerivative(x[i], y[j], z[k], SourceTerm));
	}
      }
    }
        
    //MW130625: order correct?
    for (unsigned int i = 0; i < dimx; ++i) {
      double dx_central = 0.5 * (dx_up[i]+dx_down[i]);
            
      for (unsigned int j = 0; j < dimy; ++j) {
	double dy_central = 0.5 * (dy_up[j]+dy_down[j]);
                
	for (unsigned int k = 0; k < dimz; ++k) {
	  double dz_central = 0.5 * (dz_up[k]+dz_down[k]);
                    
	  //MW130624: Trying to gain speed at the cost of memory...
	  //Tests gave me a speed up of 45% with no considerable memory loss!
	  double indspat       = index(i,j,k);
	  double indspat_xup   = index(i+1,j,k);
	  double indspat_xdown = index(i-1,j,k);
	  double indspat_yup   = index(i,j+1,k);
	  double indspat_ydown = index(i,j-1,k);
	  double indspat_zup   = index(i,j,k+1);
	  double indspat_zdown = index(i,j,k-1);

	  for  (int ip = 0; ip < dimE; ip++) {
       
	    double sp_ = sp[ip];

	    //MW130624: This would be the place to account for a spatial dependency of Delta
	    // Question: Which places are there that use sp[] and dperp[]?
	    //                         double rho_n = D_ref_rig;
	    //                         double delta_loc, delta_h_loc;
	    // 
	    //                         delta_loc = delta * pow( in->LB_delta, Coord->IsInLocalBubble(x[i],y[j],z[k]) );
	    //                         delta_h_loc = delta_h * pow( in->LB_delta, Coord->IsInLocalBubble(x[i],y[j],z[k]) );
	    //                         
	    //                         if(rho_n<=rho_b){
	    //                             if(pp[ip]<=rho_b) sp = in->D0*pow(beta[ip], in->etaT)*pow(pp[ip]/rho_n, delta_loc);
	    //                             if(pp[ip]>rho_b)  sp = in->D0*pow(beta[ip], in->etaT)*pow(pp[ip]/rho_b, delta_h_loc)*pow(rho_b/rho_n, delta_loc);
	    //                         }
	    //                         if(rho_n>rho_b){
	    //                             if(pp[ip]<=rho_b) sp = in->D0*pow(beta[ip], in->etaT)*pow(rho_b/rho_n, delta_h_loc)*pow(pp[i]/rho_b, delta_loc);
	    //                             if(pp[ip]>rho_b)  sp = in->D0*pow(beta[ip], in->etaT)*pow(pp[ip]/rho_n, delta_h_loc);
	    //                         }

	    double D       = dperp[indspat]*sp_;
	    double D_xup   = dperp[indspat_xup]*sp_;
	    double D_xdown = dperp[indspat_xdown]*sp_;
	    double D_yup   = dperp[indspat_yup]*sp_;
	    double D_ydown = dperp[indspat_ydown]*sp_;
	    double D_zup   = dperp[indspat_zup]*sp_;
	    double D_zdown = dperp[indspat_zdown]*sp_;
        
	    CNdiff_alpha1_x.push_back(D/(dx_central*dx_down[i]) - (D_xup-D_xdown)/(4*dx_central*dx_central));
	    CNdiff_alpha2_x.push_back(D/(dx_central*dx_up[i]) + D/(dx_down[i]*dx_central));
	    CNdiff_alpha3_x.push_back((D_xup-D_xdown)/(4*dx_central*dx_central) + D/(dx_up[i]*dx_central));
                        
	    CNdiff_alpha1_y.push_back(D/(dy_central*dy_down[j]) - (D_yup-D_ydown)/(4*dy_central*dy_central));
	    CNdiff_alpha2_y.push_back(D/(dy_central*dy_up[j]) + D/(dy_down[j]*dy_central));
	    CNdiff_alpha3_y.push_back((D_yup-D_ydown)/(4*dy_central*dy_central) + D/(dy_up[j]*dy_central));
                        
	    CNdiff_alpha1_z.push_back(D/(dz_central*dz_down[k]) - (D_zup-D_zdown)/(4*dz_central*dz_central));
	    CNdiff_alpha2_z.push_back(D/(dz_central*dz_up[k]) + D/(dz_down[k]*dz_central));
	    CNdiff_alpha3_z.push_back((D_zup-D_zdown)/(4*dz_central*dz_central) + D/(dz_up[k]*dz_central));
                        
	    //                         cout << "lets see " << i << " " << j << " " << k << " " << ip << " " << sp_ << " " << index(i,j,k) << " " << dz_up[k] << " " << dz_down[k] << " " << dz_central << " " << D << " " << D_zup << " " << D_zdown << " " << CNdiff_alpha1_z.back() << "  " << CNdiff_alpha2_z.back() << "  " << CNdiff_alpha3_z.back() << "  " << endl;
	    //                         CNdiff_alpha1_x.back() << "  " << CNdiff_alpha2_x.back() << "  " << CNdiff_alpha3_x.back() << "  " << CNdiff_alpha1_y.back() << "  " << CNdiff_alpha2_y.back() << "  " << CNdiff_alpha3_y.back() << "  "
	  }
	}
      }
    }
  }
  else if (diffT == Anisotropic) {
        
    cout << "In TDiffusionCoefficient3D, Anisotropic" << endl;
        
    for (unsigned int ip = 0; ip < dimE; ip++)
      //MW: this still contains no break. 130725
      sp.push_back( (El) ? in->Dpar*pow(beta[ip], in->etaTpar)*pow(pp[ip]/fabs(Z)/in->D_ref_rig, in->DeltaPar) : in->Dpar*pow(beta[ip], in->etaTpar)*pow(A*pp[ip]/fabs(Z)/in->D_ref_rig, in->DeltaPar) );
        
    for (int i = 0; i < dimx; ++i) {
      for (int k = 0; k < dimy; ++k) {
	//double radius = sqrt(x[i]*x[i]+y[k]*y[k]);
                
	for (int j = 0; j < dimz; ++j) {
	  //  double zeta = z[j];
                    
	  double profDPerp = GetProfileDPerp(x[i], y[k], z[j], SourceTerm);
	  double profDPar = GetProfileDPar(x[i], y[k], z[j], SourceTerm);
                    
	  Dpar_profile.push_back(profDPar);
                    
	  //cout << "calling GetBregVersors(int,int,int) " << endl;	
	  std::vector<double> versors = bf->GetBregVersors(i,k,j);
	  //cout << versors[0] << " " << versors[1] << " " << versors[2] << endl;	
                    
	  std::vector<double> GradProfDPerp = GetGradProfileDPerp(x[i], y[k], z[j], SourceTerm);
                    
	  std::vector<double> GradProfDPar = GetGradProfileDPar(x[i], y[k], z[j], SourceTerm);
                    
	  std::map<char, std::vector<double> > GradVersors = GetGradVersors(x[i], y[k], z[j], bf);
                    
	  for (unsigned int ip = 0; ip < dimE; ip++) {
                        
	    double Dpar = sp[ip];
                        
	    double Dparprof = Dpar*profDPar;
	    dpar_total.push_back(Dparprof);
                        
	    double DPerp = (El) ? in->Dperp*pow(beta[ip], in->etaTperp)*pow(pp[ip]/fabs(Z)/in->D_ref_rig, in->DeltaPerp) : in->Dperp*pow(beta[ip], in->etaTperp)*pow(A*pp[ip]/fabs(Z)/in->D_ref_rig, in->DeltaPerp);
                        
	    double Dperpprof = DPerp*profDPerp;
	    dperp_total.push_back(Dperpprof);
                        
	    double axx = (Dparprof - Dperpprof)*pow(versors[0],2) + Dperpprof;
	    double ayy = (Dparprof - Dperpprof)*pow(versors[1],2) + Dperpprof;
	    double azz = (Dparprof - Dperpprof)*pow(versors[2],2) + Dperpprof;

	    /*if (j == dimz/2 && ip == dimE/2 && i%5 == 0 && k%5 == 0) {
	      cout << i << " <--- i --- j ---> " << j << endl;
	      cout << "Dpar  = " << Dparprof  << endl;
	      cout << "Dperp = " << Dperpprof << endl;
	      cout << "versor x = " << versors[0] <<  " axx = " << axx << endl;
	      cout << "versor y = " << versors[1] <<  " ayy = " << ayy << endl;
	      cout << "versor z = " << versors[2] <<  " azz = " << azz << endl;
	      }*/
			
                        
	    double axy = 2.0*(Dparprof - Dperpprof)*versors[0]*versors[1];
	    double ayz = 2.0*(Dparprof - Dperpprof)*versors[2]*versors[1];
	    double axz = 2.0*(Dparprof - Dperpprof)*versors[0]*versors[2];
                        
	    alpha_xx.push_back(axx);
	    alpha_yy.push_back(ayy);
	    alpha_zz.push_back(azz);
                        
	    alpha_xy.push_back(axy);
	    alpha_yz.push_back(ayz);
	    alpha_xz.push_back(axz);
                        
	    double ux =
	      (Dpar*GradProfDPar[0] - DPerp*GradProfDPerp[0])*pow(versors[0],2)
	      + DPerp*GradProfDPerp[0]
	      + (Dparprof - Dperpprof)*2.0*versors[0]*GradVersors['x'][0]    // dx alpha_xx
	      + (Dpar*GradProfDPar[1] - DPerp*GradProfDPerp[1])*versors[0] * versors[1]
	      + (Dparprof - Dperpprof)*versors[0]*GradVersors['y'][1]
	      + (Dparprof - Dperpprof)*versors[1]*GradVersors['y'][0]   // dy alpha_xy
	      + (Dpar*GradProfDPar[2] - DPerp*GradProfDPerp[2])*versors[2] * versors[0]
	      + (Dparprof - Dperpprof)*versors[2]*GradVersors['z'][0]
	      + (Dparprof - Dperpprof)*versors[0]*GradVersors['z'][2]   // dz alpha_xz
	      ;
                        
	    double uy =
	      (Dpar*GradProfDPar[1] - DPerp*GradProfDPerp[1])*pow(versors[1],2)
	      + DPerp*GradProfDPerp[1]
	      + (Dparprof - Dperpprof)*2.0*versors[1]*GradVersors['y'][1]    // dy alpha_yy
	      + (Dpar*GradProfDPar[0] - DPerp*GradProfDPerp[0])*versors[0] * versors[1]
	      + (Dparprof - Dperpprof)*versors[0]*GradVersors['x'][1]
	      + (Dparprof - Dperpprof)*versors[1]*GradVersors['x'][0]   // dx alpha_xy
	      + (Dpar*GradProfDPar[2] - DPerp*GradProfDPerp[2])*versors[2] * versors[1]
	      + (Dparprof - Dperpprof)*versors[2]*GradVersors['z'][1]
	      + (Dparprof - Dperpprof)*versors[1]*GradVersors['z'][2]   // dz alpha_yz
	      ;
                        
	    double uz =
	      (Dpar*GradProfDPar[2] - DPerp*GradProfDPerp[2])*pow(versors[2],2)
	      + DPerp*GradProfDPerp[2]
	      + (Dparprof - Dperpprof)*2.0*versors[2]*GradVersors['z'][2]    // dz alpha_zz
	      + (Dpar*GradProfDPar[1] - DPerp*GradProfDPerp[1])*versors[2] * versors[1]
	      + (Dparprof - Dperpprof)*versors[2]*GradVersors['y'][1]
	      + (Dparprof - Dperpprof)*versors[1]*GradVersors['y'][2]   // dy alpha_yz
	      + (Dpar*GradProfDPar[0] - DPerp*GradProfDPerp[0])*versors[2] * versors[0]
	      + (Dparprof - Dperpprof)*versors[2]*GradVersors['x'][0]
	      + (Dparprof - Dperpprof)*versors[0]*GradVersors['x'][2]   // dx alpha_xz
	      ;
                        
	    phix.push_back(ux);
	    phiy.push_back(uy);
	    phiz.push_back(uz);

	    //     cout << "[MW-DEBUG DIFF] DiffT " << diffT << " and so at " << i << " " << k << " " << j << " | " << ux << " " << uy << " " << uz << " | " << axx << axy << axz << ayy << ayz << azz << " | " << endl;

	    /*if (j == dimz/2 && ip == dimE/2 && i%5 == 0 && k == dimy/2) {
	      cout << i << " <--- i --- j ---> " << j << endl;
	      cout << "Dpar  = " << Dparprof  << endl;
	      cout << "Dperp = " << Dperpprof << endl;
	      cout << "versor x = " << versors[0] <<  " axx = " << axx << endl;
	      cout << "versor y = " << versors[1] <<  " ayy = " << ayy << endl;
	      cout << "versor z = " << versors[2] <<  " azz = " << azz << endl;
	      cout << "axx = " << axx << endl;
	      cout << "ayy = " << ayy << endl;
	      cout << "azz = " << azz << endl;
	      cout << "ux =  " << ux << endl;
	      cout << "uy =  " << uy << endl; 
	      }*/
	  }
	}
      }
    }
    
  }
  
  return ;
}

double TDiffusionCoefficient3D::GetProfileDPar(double x, double y, double z, TSource* SourceTerm) {
  double radial = 1;
  double result = 1.;
  //    cout << " nrn_sn in GetProfile ->" << nrn_sn << endl;
    
  switch(set_profile) {
  case Constant :
    return 1.0;
    break;
            
  case Exp :
    return exp(fabs(z)/zt);
    break;
    /*
      case Blasi :
      if (fabs(zeta) < 1.) return exp(fabs(zeta)/0.5);
      else return exp(1.0/0.5)*exp(fabs(zeta)/zt)/exp(1.0/zt);
      break;
             
      case Expr :
      return exp(fabs(zeta/zt))*(0.5*(tanh((radius-3.0)/0.25)+1.001));///cosh((radius-r0)/rd);
      break;
    */
  case Qtau :
    /*
      if (radius > 3.7) radial = 50.0*(0.79*exp(-pow(zobs/0.212,2.))+0.21*exp(-pow(zobs/0.636,2.)))*exp(-((radius)*(radius)-robs*robs)/(6.8*6.8)) + 7.3*exp(-(radius-robs)/4.5-fabs(zobs)/0.325);
      else radial = 177.5*(0.79*exp(-pow(zobs/0.212,2.))+0.21*exp(-pow(zobs/0.636,2.)))*exp(-pow((radius-3.7)/2.1, 2.)) + 7.3*exp(-(radius-robs)/4.5-fabs(zobs)/0.325);
    */
    radial = SourceTerm->GetSource(x,y,z);
    radial /= nrn_sn;
    return pow(radial, -index_radial);
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
  case Test:
    //if (x < 0.) return 1.;
    //else 
    //return 1.;
    return exp(x/5.0);
    break; 
  default :
    return -1;
  }
}

double TDiffusionCoefficient3D::GetXDerivativeDPar(double x, double y, double z, TSource* SourceTerm) {
    
  double delta_x = 0.00002; //Coord->GetDeltaX();
  double xup = x + 0.5*delta_x;
  double xdo = x - 0.5*delta_x;
  return ( GetProfileDPar(xup,y,z, SourceTerm) - GetProfileDPar(xdo,y,z, SourceTerm) )/delta_x ;
}

double TDiffusionCoefficient3D::GetYDerivativeDPar(double x, double y, double z, TSource* SourceTerm) {
    
  double delta_y = 0.00002; //Coord->GetDeltaY();
  double yup = y + 0.5*delta_y;
  double ydo = y - 0.5*delta_y;
  return ( GetProfileDPar(x,yup,z, SourceTerm) - GetProfileDPar(x,ydo,z, SourceTerm) )/delta_y;
}

double TDiffusionCoefficient3D::GetZDerivativeDPar(double x, double y, double z, TSource* SourceTerm) {
    
  double delta_z = 0.00002; //Coord->GetDeltaZ();
  double zup = z + 0.5*delta_z;
  double zdo = z - 0.5*delta_z;
  return ( GetProfileDPar(x,y,zup, SourceTerm) - GetProfileDPar(x,y,zdo, SourceTerm) )/delta_z;
}

vector<double> TDiffusionCoefficient3D::GetXDerivativeBF(double x, double y, double z, TBField* SourceTerm) {
    
  double delta_x = 0.00002; //Coord->GetDeltaX();
  double xup = x + 0.5*delta_x;
  double xdo = x - 0.5*delta_x;
  vector<double> vers1 = SourceTerm->GetBregVersors(xup,y,z);
  vector<double> vers2 = SourceTerm->GetBregVersors(xdo,y,z);
  vector<double> final(3,0);
  final[0] = (vers1[0]-vers2[0])/delta_x;
  final[1] = (vers1[1]-vers2[1])/delta_x;
  final[2] = (vers1[2]-vers2[2])/delta_x;
  return final ;
}

vector<double> TDiffusionCoefficient3D::GetYDerivativeBF(double x, double y, double z, TBField* SourceTerm) {
    
  double delta_y = 0.00002; //Coord->GetDeltaY();
  double yup = y + 0.5*delta_y;
  double ydo = y - 0.5*delta_y;
  vector<double> vers1 = SourceTerm->GetBregVersors(x,yup,z);
  vector<double> vers2 = SourceTerm->GetBregVersors(x,ydo,z);
  vector<double> final(3,0);
  final[0] = (vers1[0]-vers2[0])/delta_y;
  final[1] = (vers1[1]-vers2[1])/delta_y;
  final[2] = (vers1[2]-vers2[2])/delta_y;
  return final ;
}

vector<double> TDiffusionCoefficient3D::GetZDerivativeBF(double x, double y, double z, TBField* SourceTerm) {
    
  double delta_z = 0.00002; //Coord->GetDeltaZ();
  double zup = z + 0.5*delta_z;
  double zdo = z - 0.5*delta_z;
  vector<double> vers1 = SourceTerm->GetBregVersors(x,y,zup);
  vector<double> vers2 = SourceTerm->GetBregVersors(x,y,zdo);
  vector<double> final(3,0);
  final[0] = (vers1[0]-vers2[0])/delta_z;
  final[1] = (vers1[1]-vers2[1])/delta_z;
  final[2] = (vers1[2]-vers2[2])/delta_z;
  return final ;
}


std::vector<double> TDiffusionCoefficient3D::GetGradProfileDPerp(double x, double y, double z, TSource* ST) {
  vector<double> prof(3,0);
  prof[0] = GetXDerivative(x,y,z,ST);
  prof[1] = GetYDerivative(x,y,z,ST);
  prof[2] = GetZDerivative(x,y,z,ST);
  return prof;
}

std::vector<double> TDiffusionCoefficient3D::GetGradProfileDPar(double x, double y, double z, TSource* ST) {
  vector<double> prof(3,0);
  prof[0] = GetXDerivativeDPar(x,y,z,ST);
  prof[1] = GetYDerivativeDPar(x,y,z,ST);
  prof[2] = GetZDerivativeDPar(x,y,z,ST);
  return prof;
    
}

std::map<char, std::vector<double> > TDiffusionCoefficient3D::GetGradVersors(double x, double y, double z, TBField* bf) {
  map<char, std::vector<double> > gradver;

  gradver['x'] = GetXDerivativeBF(x,y,z,bf);
  gradver['y'] = GetYDerivativeBF(x,y,z,bf);
  gradver['z']  = GetZDerivativeBF(x,y,z,bf);
    
  return gradver;
}
