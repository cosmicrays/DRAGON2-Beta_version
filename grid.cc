/**
 * @file grid.cc
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief See the .h file.
 */

#include "grid.h"
#include "input.h"
#include "constants.h"
#include <cmath>

#include <fstream>
#include <limits>

using namespace std;
TGrid::TGrid(Input* in) {
   
  if (in->feedback >0) cout << "Welcome to the grid constructor" << endl;
   
   dimx = dimy = dimr = 1; // 2D<->3D compatibility
   dimz = in->numz;
   zmax = in->zmax;
   if (in->feedback >0) cout << "DimZ = "<<dimz<<endl;
   
   string divz;
   stringstream temp;
   double next_point;
   
   if (in->dimz_equidistant)
   {
      if (in->feedback >0) cout << "// [MW] z bins equidistant" << endl;
      zmin = -zmax;
      for (int i=0; i<(dimz-1); ++i)
      {
         double val = zmin + (double)i/(dimz-1)*(zmax-zmin);
         //mw: round to 1e-8, solves precision issues
         if(val!=0) val = (long)((1e8*val)+(val/fabs(val))*.5)/1e8;
         z.push_back( ( fabs(val) < 1e-8 ) ? 0. : val);
      }
      z.push_back(zmax);
   }
   else
   {
      divz = in->divz;
      int i=0,j;
      while(i<dimz)
      {
         j = divz.find(';');
         temp << divz.substr(0,j);
         temp >> next_point;
         z.push_back(in->zmax * next_point);
         divz = divz.substr(j+1);
         i++;
         temp.clear();
      }
      //rescale, so that in->zmax is really the maximum absolute value
      double z_rescale = in->zmax/max(abs(z.front()), abs(z.back()));
      for( vector<double>::iterator it = z.begin(); it!=z.end(); ++it) *it *= z_rescale;
   }
   if (in->feedback >0){
     cout << "// [MW] z[] grid; setting DimZ = " << dimz << endl;
     for (int i=0;i<dimz;i++)
       {
	 cout<<"z["<<i<<"] = " << z[i];
	 if(in->gridtype=="3D" && in->LB_shape == "Cuboid" && (abs(z[i]-in->zobs)<=in->LB_az))cout << " /* Local Cuboid */ ";
	 cout << "; ";
       }
     cout << endl;
   }
   //MW130618: fill dz_up and dz_down vectors to speed up the evolutor
   dz_up.push_back(z[1]-z[0]);
   dz_down.push_back(z[1]-z[0]);
   for (int i=1;i<dimz-1;i++)
   {
      dz_up.push_back(z[i+1]-z[i]);
      dz_down.push_back(z[i]-z[i-1]);
   }
   dz_up.push_back(z[dimz-1]-z[dimz-2]);
   dz_down.push_back(z[dimz-1]-z[dimz-2]);
   if (in->feedback >0) cout << "dimz is " << dimz << " and z holds " << z.size() << " elements and dz_up " << dz_up.size() << " and dz_down " << dz_down.size() << endl;
   
   Ekin_factor = in->Ekfact;
   Ekmin = in->Ekmin;
   Ekmax = in->Ekmax;
   DeltalogE = log(Ekin_factor);
   dimE = int(log(Ekmax/Ekmin)/DeltalogE + 1.9);
   if( in->spallationxsec == Fluka) {
     dimE = 160;
     dimE = int(log(Ekmax/Ekmin)/DeltalogE)+2;
     cout << "grid.cc: Fluka Xsec , dimE = " << dimE << " " << Ekmin << " " << Ekmax << " " << Ekin_factor << endl; 
   }
   if( in->propDM) {
     cout << "Propagate DM? " << in->propDM << endl;
     cout << "Mass = " << in->mx << endl;
     dimE += 4;
   }
   double Ei;
   for (int i = 0; i < dimE; ++i) {
     if (dimE == 1) {
       Ek.push_back(Ekmin);
     } else {
       Ei = exp(log(Ekmin)+(double)i*log(Ekin_factor));
       Ek.push_back(Ei);
     }
      cout << "grid.cc: " << i << " " << Ek[i] << endl; 
      
      gamma.push_back(1.0+Ek.back()/mp);
      beta.push_back(sqrt(1.0-1.0/pow(gamma.back(),2)));
      momentum.push_back(gamma.back()*mp*beta.back());
      
      gammael.push_back(1.0+Ek.back()/MeleGeV);
      betael.push_back(sqrt(1.0-1.0/pow(gammael.back(),2)));
      momentumel.push_back(gammael.back()*MeleGeV*betael.back());
      //         cout << "[MW-DEBUG GRID] p: " << i << "/" << dimE << " ... " << gamma.back() << " " <<momentum.back() << " " << beta.back() << " " << mp << endl;
      //         cout << "[MW-DEBUG GRID] e: " << i << "/" << dimE << " ... " << gammael.back() << " " <<momentumel.back() << " " << betael.back() << " " << MeleGeV << endl;
   }
}

TGrid2D::TGrid2D(Input* in) : TGrid(in) {
   
   type = "2D";
   
   dimr = in->numr;
   
   string divr;
   stringstream temp;
   double next_point;
   
   if (in->dimr_equidistant) {
      //cout << "// [MW] r bins equidistant" << endl;
      if (in->feedback >0) cout << "Building the 2D equidistant grid " << endl;
      for (int i=0; i<(dimr-1); ++i)
      {
         double val = in->Rmin + (double)i/(dimr-1)*(in->Rmax - in->Rmin);
         //             cout << i << " val " << val << " rmin " << in->Rmin << " " << i/(dimr-1) << " " << (double)i/(dimr-1)*(in->Rmax - in->Rmin) << endl;
         //mw: round to 1e-8, solves precision issues
         if(val!=0) val = (long)((1e8*val)+(val/abs(val))*.5)/1e8;
         r.push_back( ( abs(val) < 1e-8 ) ? 0. : val);
      }
      r.push_back(in->Rmax);
   }
   else
   {
      divr = in->divr;
      int i=0,j;
      while(i<dimr)
      {
         j = divr.find(';');
         temp << divr.substr(0,j);
         temp >> next_point;
         r.push_back(in->Rmax * next_point);
         divr = divr.substr(j+1);
         i++;
         temp.clear();
      }
   }
   //rescale
   double r_rescale = in->Rmax/max(abs(r.front()), abs(r.back()));
   for( vector<double>::iterator it = r.begin(); it!=r.end(); ++it) *it *= r_rescale;
   
   
   x = r;
   dimx = dimr;
   dimy = 1;
   
   //MW130709: now also 2D
   dr_up.push_back(r[1]-r[0]);
   dr_down.push_back(r[1]-r[0]);
   for (int i=1;i<dimr-1;i++)
   {
      dr_up.push_back(r[i+1]-r[i]);
      dr_down.push_back(r[i]-r[i-1]);
   }
   dr_up.push_back(r[dimr-1]-r[dimr-2]);
   dr_down.push_back(r[dimr-1]-r[dimr-2]);
   
   //     cout.precision(8);
   if (in->feedback >0){
     cout << "// [MW] r[] grid; setting DimR = " << dimr << endl;
     for (int i=0;i<dimr;i++) cout<<"r["<<i<<"] = " << r[i] << "; ";
     cout << endl;
   }
}

TGrid3D::TGrid3D(Input* in) : TGrid(in) {
   
   type = "3D";
   
   //local bubble shape
   if( (LB_ax = in->LB_ax) == -1 ) LB_ax = 0.3;
   if( (LB_ay = in->LB_ay) == -1 ) LB_ay = 0.3;
   if( (LB_az = in->LB_az) == -1 ) LB_az = 0.3;
   
   //precision issues
   LB_ax += 1e-8;
   LB_ay += 1e-8;
   LB_az += 1e-8;
   
   LB_shape       = in->LB_shape;
   LB_smearing    = in->LB_smearing;
   
   in_SA_type         = in->SA_type;
   in_LB_MagField     = in->LB_MagField;
   in_SA_MagField     = in->SA_MagField;
   in_SA_cut_MagField = in->SA_cut_MagField;
   in_LB_ISRF         = in->LB_ISRF;
   in_SA_ISRFStar     = in->SA_ISRFStar;
   in_SA_ISRFDust     = in->SA_ISRFDust;
   in_SA_cut_ISRF     = in->SA_cut_ISRF;
   
   dimx = in->numx;
   dimy = in->numy;
   
   string divx, divy;
   stringstream temp;
   double next_point;
   
   //matze
   if (in->dimx_equidistant) {
     if (in->feedback >0) cout << "// [MW] x bins equidistant" << endl;
      for (int i=0; i<(dimx-1); ++i)
      {
         double val = -in->Rmax + (double)i/(dimx-1)*(2.0*in->Rmax);
         //mw: round to 1e-8, solves precision issues
         if(val!=0) val = (long)((1e8*val)+(val/abs(val))*.5)/1e8;
         x.push_back( ( abs(val) < 1e-8 ) ? 0. : val);
      }
      x.push_back(in->Rmax);
   }
   else
   {
      divx = in->divx;
      int i=0,j;
      while(i<dimx)
      {
         j = divx.find(';');
         temp << divx.substr(0,j);
         temp >> next_point;
         x.push_back(in->Rmax * next_point);
         divx = divx.substr(j+1);
         i++;
         temp.clear();
      }
   }
   
   if (in->dimy_equidistant) {
       if (in->feedback >0) cout << "// [MW] y bins equidistant" << endl;
      for (int i=0; i<(dimy-1); ++i)
      {
         double val = -in->Rmax + (double)i/(dimy-1)*(2.0*in->Rmax);
         //mw: round to 1e-8, solves precision issues
         if(val!=0) val = (long)((1e8*val)+(val/abs(val))*.5)/1e8;
         y.push_back( ( abs(val) < 1e-8 ) ? 0. : val);
      }
      y.push_back(in->Rmax);
   }
   else
   {
      divy = in->divy;
      int i=0,j;
      while(i<dimy)
      {
         j = divy.find(';');
         temp << divy.substr(0,j);
         temp >> next_point;
         y.push_back(in->Rmax * next_point);
         divy = divy.substr(j+1);
         i++;
         temp.clear();
      }
   }
   
   //rescale, so that in->Rmax is really the maximum absolute value of both directions
   double xy_rescale = in->Rmax/max( max(abs(x.front()), abs(x.back())) , max(abs(y.front()), abs(y.back())) );
   for( vector<double>::iterator it = x.begin(); it!=x.end(); ++it) *it *= xy_rescale;
   for( vector<double>::iterator it = y.begin(); it!=y.end(); ++it) *it *= xy_rescale;
   
   if (in->feedback >0){
     cout << "// [MW] x[] grid; setting DimX = " << dimx << endl;
     for (int i=0;i<dimx;i++)
       {
	 cout<<"x["<<i<<"] = " << x[i];
	 if(LB_shape == "Cuboid" && (abs(x[i]-in->xobs)<=LB_ax))cout << " /* Local Cuboid */";
	 cout << "; ";
       }
     cout << endl << "// [MW] y[] grid; setting DimY = " << dimy << endl;
     for (int i=0;i<dimy;i++)
       {
	 cout<<"y["<<i<<"] = " << y[i];
	 if(LB_shape == "Cuboid" && (abs(y[i]-in->yobs)<=LB_ay))cout << " /* Local Cuboid */";
	 cout << "; ";
       }
     cout << endl;
   }
   //MW130618: fill d..._up and d..._down vectors to speed up the evolutor
   dx_up.push_back(x[1]-x[0]);
   dx_down.push_back(x[1]-x[0]);
   for (int i=1;i<dimx-1;i++)
   {
      dx_up.push_back(x[i+1]-x[i]);
      dx_down.push_back(x[i]-x[i-1]);
   }
   dx_up.push_back(x[dimx-1]-x[dimx-2]);
   dx_down.push_back(x[dimx-1]-x[dimx-2]);
   
   dy_up.push_back(y[1]-y[0]);
   dy_down.push_back(y[1]-y[0]);
   for (int i=1;i<dimy-1;i++)
   {
      dy_up.push_back(y[i+1]-y[i]);
      dy_down.push_back(y[i]-y[i-1]);
   }
   dy_up.push_back(y[dimy-1]-y[dimy-2]);
   dy_down.push_back(y[dimy-1]-y[dimy-2]);
   
   
   //write Bubble shape to file
   std::ofstream bub;
   bub.open("data/BUBBLE.DAT");
   for(int ix=0; ix<dimx; ix++)
      for(int iy=0; iy<dimy; iy++)
         for(int iz=0; iz<dimz; iz++)
            if(IsInLocalBubble_Indexed(ix, iy, iz)!=0) bub << x[ix] << " " << y[iy] << " " << z[iz] << " " << IsInLocalBubble_Indexed(ix, iy, iz) << " " << endl;
   bub.close();
    if (in->feedback >0) cout << "data/BUBBLE.DAT written, smearing is " << LB_smearing << endl;
   
}

double TGrid::IsInLocalBubble(double xx, double yy, double zz)
{
  double xobs = 8.3; //in->xobs;
  double yobs = 0.;  //in->yobs;
  double zobs = 0.;  //in->zobs;

   if(type=="3D")
   {
      double IsInLocal = 0;
      
      if(LB_shape == "Cuboid")
	IsInLocal = ( (fabs(xx-xobs)<=LB_ax) && (fabs(yy-yobs)<=LB_ay) && (fabs(zz-zobs)<=LB_az) );
      else if(LB_shape == "Ellipsoid")
         IsInLocal = ( pow((xx-xobs)/LB_ax,2) + pow((yy-yobs)/LB_ay,2) + pow((zz-zobs)/LB_az,2) <= 1 );
      
      
      if (LB_smearing == "step" || LB_smearing == "Step")
      {
         
         const double lim_out = 4;
         const double lim_in = 2;
         const double step_out = 0.25;
         const double step_in = 0.5;
         
         IsInLocal = ( (fabs(xx-xobs)<lim_out*LB_ax) ? ( (fabs(xx-xobs)<lim_in*LB_ax) ? ( (fabs(xx-xobs) < LB_ax) ? 1 : step_in ) : step_out ) : 0 )
         * ( (fabs(yy-yobs)<lim_out*LB_ay) ? ( (fabs(yy-yobs)<lim_in*LB_ay) ? ( (fabs(yy-yobs) < LB_ay) ? 1 : step_in ) : step_out ) : 0 )
         * ( (fabs(zz-zobs)<lim_out*LB_az) ? ( (fabs(zz-zobs)<lim_in*LB_az) ? ( (fabs(zz-zobs) < LB_az) ? 1 : step_in ) : step_out ) : 0 );
      }
      else if (LB_smearing == "linear" || LB_smearing == "Linear")
      {
         const double lim = 4;
         
         IsInLocal = ( (fabs(xx-xobs)<lim*LB_ax) ? ( (fabs(xx-xobs) < LB_ax) ? 1 : (fabs(xx/LB_ax) - lim)/(1-lim) ) : 0 )
         * ( (fabs(yy-yobs)<lim*LB_ay) ? ( (fabs(yy-yobs) < LB_ay) ? 1 : (fabs(yy/LB_ay) - lim)/(1-lim) ) : 0 )
         * ( (fabs(zz-zobs)<lim*LB_az) ? ( (fabs(zz-zobs) < LB_az) ? 1 : (fabs(zz/LB_az) - lim)/(1-lim) ) : 0 );
      }
      else if (LB_smearing == "gauss" || LB_smearing == "Gauss" || LB_smearing == "gaussian" || LB_smearing == "Gaussian")
      {
         const double lim = 3;
         
         IsInLocal = ( (fabs(xx-xobs)<1.5*lim*LB_ax) ? ( (fabs(xx-xobs) < LB_ax) ? 1 : (exp(-pow((2./lim)*(fabs(xx-xobs)/LB_ax-1),2))) ) : 0 )
         * ( (fabs(yy-yobs)<1.5*lim*LB_ay) ? ( (fabs(yy-yobs) < LB_ay) ? 1 : (exp(-pow((2./lim)*(fabs(yy-yobs)/LB_ay-1),2))) ) : 0 )
         * ( (fabs(zz-zobs)<1.5*lim*LB_az) ? ( (fabs(zz-zobs) < LB_az) ? 1 : (exp(-pow((2./lim)*(fabs(zz-zobs)/LB_az-1),2))) ) : 0 );
         
         //cout << "[MW-DEBUG-BUBBLE] LB is Gaussian and values are " << xx << " " << yy << " " << zz << " | " << IsInLocal << endl;
      }
      
      return IsInLocal;
   }
   else
   {
     cout << "[MW] WARNING: tried to call IsInLocalBubble for a Non-3D grid. There is no bubble." << endl;
     return 0;
   }
}

double TGrid::IsInLocalBubble_Indexed(int ix, int iy, int iz) //MW: just in case it's called from some scope where the X,Y,Z are not known.
{
   if( ix<0 or ix>=x.size() or iy<0 or iy>=y.size() or iz<0 or iz>=z.size() ) return 0;
   return IsInLocalBubble(x[ix], y[iy], z[iz]);
}
