/**
 * @file grid.h
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief File in which the TGrid class is defined.
 */

#ifndef _GRID_H
#define _GRID_H

#include <iostream>
#include <vector>
#include <math.h>
#include <sstream>

using namespace std;
class Input;

/**
 * @class TGrid
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class to handle kinematics variables and the geometry of the galaxy.
 */

class TGrid {
 public:
  TGrid() { }
  TGrid(Input*);
  
  virtual ~TGrid() {
    z.clear();
    Ek.clear();
    beta.clear();
    gamma.clear();
    momentum.clear();

    betael.clear();
    gammael.clear();
    momentumel.clear();

    r.clear();
    x.clear();
    y.clear();

    dx_up.clear(); /**< X grid spacing, up. */
    dx_down.clear(); /**< X grid spacing, down. */
    dy_up.clear(); /**< Y grid spacing, up. */
    dy_down.clear(); /**< Y grid spacing, down. */
    dr_up.clear(); /**< R grid spacing, up. */
    dr_down.clear(); /**< R grid spacing, down. */
    dz_up.clear(); /**< Z grid spacing, up. */
    dz_down.clear(); /**< Z grid spacing, down. */
  }

  inline vector<double>& GetZ() { return z; } /**< Returns the Z grid. */

  inline vector<double>& GetEk() { return Ek; } /**< Returns the kinetic energy grid. */
  inline vector<double>& GetBeta() { return beta; } /**< Returns the beta grid. */
  inline vector<double>& GetGamma() { return gamma; } /**< Returns the boost factor grid. */
  inline vector<double>& GetMomentum() { return momentum; } /**< Returns the momentum grid. */

  inline vector<double>& GetBetaEl() { return betael; } /**< Returns the beta grid for electrons. */
  inline vector<double>& GetGammaEl() { return gammael; } /**< Returns the boost factor grid for electrons. */
  inline vector<double>& GetMomentumEl() { return momentumel; } /**< Returns the momentum grid for electrons. */

  //MW130619: up- and downward steps, in vectors
  inline vector<double>& GetDeltaX_up() { return dx_up; }
  inline vector<double>& GetDeltaY_up() { return dy_up; }
  inline vector<double>& GetDeltaZ_up() { return dz_up; }
  inline vector<double>& GetDeltaR_up() { return dr_up; }
  inline vector<double>& GetDeltaX_down() { return dx_down; }
  inline vector<double>& GetDeltaY_down() { return dy_down; }
  inline vector<double>& GetDeltaZ_down() { return dz_down; }
  inline vector<double>& GetDeltaR_down() { return dr_down; }
  inline double GetDeltaX_up(int k) const   { return dx_up.at(k);   } /**< Returns Delta x, upwards. */
  inline double GetDeltaY_up(int k) const   { return dy_up.at(k);   } /**< Returns Delta y, upwards. */
  inline double GetDeltaR_up(int k) const   { return dr_up.at(k);   } /**< Returns Delta x, upwards. */
  inline double GetDeltaZ_up(int k) const   { return dz_up.at(k);   } /**< Returns Delta z, upwards. */
  inline double GetDeltaX_down(int k) const { return dx_down.at(k); } /**< Returns Delta x, downwards. */
  inline double GetDeltaY_down(int k) const { return dy_down.at(k); } /**< Returns Delta y, downwards. */
  inline double GetDeltaR_down(int k) const { return dr_down.at(k); } /**< Returns Delta z, downwards. */
  inline double GetDeltaZ_down(int k) const { return dz_down.at(k); } /**< Returns Delta z, downwards. */
  inline double GetDeltaX_central(int k) const { return 0.5 * (dx_up.at(k) + dx_down.at(k)); }
  inline double GetDeltaY_central(int k) const { return 0.5 * (dy_up.at(k) + dy_down.at(k)); }
  inline double GetDeltaR_central(int k) const { return 0.5 * (dr_up.at(k) + dr_down.at(k)); }
  inline double GetDeltaZ_central(int k) const { return 0.5 * (dz_up.at(k) + dz_down.at(k)); }
  
  //MW130619: symmetric scheme
  inline double GetDeltaX(int i) {
                                   if(i==0)           return (x[1]-x[0]);
                                   else if(i==dimx-1) return (x[dimx-1] - x[dimx-2]);
                                   else               return 0.5*(x[i+1] - x[i-1]);
                                 }
  inline double GetDeltaY(int i) {
                                   if(i==0)           return (y[1]-y[0]);
                                   else if(i==dimy-1) return (y[dimy-1] - y[dimy-2]);
                                   else               return 0.5*(y[i+1] - y[i-1]);
                                 }
  inline double GetDeltaZ(int i) {
                                   if(i==0)           return (z[1]-z[0]);
                                   else if(i==dimz-1) return (z[dimz-1] - z[dimz-2]);
                                   else               return 0.5*(z[i+1] - z[i-1]);
                                 }

  inline double GetDeltaR(int i) {
                                   if(i==0)           return (r[1]-r[0]);
                                   else if(i==dimr-1) return (r[dimr-1] - r[dimr-2]);
                                   else               return 0.5*(r[i+1] - r[i-1]);
                                 }
  
  inline double GetEkMin() const { return Ekmin; } /**< Returns Delta z. */
  inline double GetEkMax() const { return Ekmax; } /**< Returns Delta z. */
  inline double GetDeltaE() const { return DeltalogE; } /**< Returns log_{10}(Delta E/E). */
  inline int GetDimZ() const { return dimz; } /**< Returns dimz. */
  inline int GetDimE() const { return dimE; } /**< Returns dimE. */

  virtual vector<double>& GetR() { return r; }  /**< Returns the R grid. */
  virtual vector<double>& GetX() { return x; } /**< Returns the X grid. */
  virtual vector<double>& GetY() { return y; } /**< Returns the Y grid. */
  
  //MW
  virtual int GetXFromIndexD_3D(int indexD) { return (((indexD)/dimz)/dimy); }
  virtual int GetYFromIndexD_3D(int indexD) { return (((indexD)/dimz)%dimy); }
  virtual int GetZFromIndexD_3D(int indexD) { return  ((indexD)%dimz); }
  virtual int GetXFromIndexD_4D(int indexD) { return (((indexD/dimE)/dimz)/dimy); }
  virtual int GetYFromIndexD_4D(int indexD) { return (((indexD/dimE)/dimz)%dimy); }
  virtual int GetZFromIndexD_4D(int indexD) { return  ((indexD/dimE)%dimz); }
  
  virtual int GetDimR() { return dimr; } /**< Returns dimr. */
  virtual int GetDimX() { return dimx; } /**< Returns dimr. */
  virtual int GetDimY() { return dimy; } /**< Returns dimr. */

  inline int index(int ir, int iz) { return ir*dimz+iz; }
  /**
   * @fn inline int index(int ir, int iz)
   * @brief Compute linearized index for 2-D spatial matrices. 
   * @param ir radial index
   * @param iz vertical index
   * @return ir*dimz+iz
   */
  inline int index(int ir, int iz, int ip) { return (ir*dimz+iz)*dimE+ip; }
  /**
   * @fn inline int index(int ir, int iz, int ip)
   * @brief Compute linearized index for 3-D momentum-spatial matrices. 
   * @param ir radial index
   * @param iz vertical index
   * @param ip momentum index
   * @return (ir*dimz+iz)*dimE+ip
   */

  virtual int indexD(int ix, int iy, int iz) {
        if(ix<0) ix = 0;
        if(iy<0) iy = 0;
        if(iz<0) iz = 0;
        if(ix>=dimx) ix = dimx - 1;
        if(iy>=dimy) iy = dimy - 1;
        if(iz>=dimz) iz = dimz - 1;
        
        return (ix*dimy + iy)*dimz+iz;
  }
    /**
   * @fn inline int indexD(int ix, int iy, int iz)
   * @brief Compute linearized index for 3-D spatial matrices. 
     */
  virtual long int indexD(int ix, int iy, int iz, int ip) { return (indexD(ix,iy,iz))*dimE+ip; }
  /**
   * @fn inline int indexD(int ix, int iy, int iz, int ip)
   * @brief Compute linearized index for 4-D momentum-spatial matrices. 
   */
   inline string GetType() const { return type; }
   
   //MW: Local Bubble accessoires
   double IsInLocalBubble(double xx, double yy, double zz); //MW
   double IsInLocalBubble_Indexed(int ix, int iy, int iz); //MW
   string in_SA_type; //MW
   double in_SA_MagField, in_SA_cut_MagField, in_LB_MagField, in_SA_ISRFStar, in_SA_ISRFDust, in_SA_cut_ISRF, in_LB_ISRF; //MW: because the Input object is not available anywhere.

 protected:
  int dimE; /**< Dimension of energy grid. */
  int dimz; /**< Dimension of vertical grid. */
  double zmax; /**< Vertical upper edge of simulation box. */
  double zmin; /**< Vertical lower edge of simulation box. */
  double Ekmin;
  double Ekmax;
  double Deltaz; /**< Vertical grid spacing. */
  double Ekin_factor; /**< Energy grid spacing. */
  double DeltalogE; /**< Radial grid spacing (logarithmic). */

  vector<double> z; /**< Array of points in the vertical grid. */ 
  vector<double> Ek; /**< Array of points in the kinetic energy grid. */ 
  vector<double> beta; /**< Array of points in the beta grid. */ 
  vector<double> gamma; /**< Array of points in the boost grid. */ 
  vector<double> momentum; /**< Array of points in the momentum grid. */ 

  vector<double> betael; /**< Array of points in the beta grid (electrons). */
  vector<double> gammael; /**< Array of points in the boost grid (electrons). */ 
  vector<double> momentumel; /**< Array of points in the momentum grid (electrons). */ 
  string type;
  
  int dimr; /**< Dimension of radial grid. */
//  double Deltar; /**< Radial grid spacing. */
  vector<double> r; /**< Array of points in the radial grid. */ 
  
  int dimx; /**< Dimension of x grid. */
  int dimy; /**< Dimension of y grid. */
  vector<double> x; /**< Array of points in the radial grid. */
  vector<double> y; /**< Array of points in the radial grid. */

  vector<double> dx_up; /**< X grid spacing, up. */
  vector<double> dx_down; /**< X grid spacing, down. */
  vector<double> dy_up; /**< Y grid spacing, up. */
  vector<double> dy_down; /**< Y grid spacing, down. */
  vector<double> dr_up; /**< R grid spacing, up. */
  vector<double> dr_down; /**< R grid spacing, down. */
  vector<double> dz_up; /**< Z grid spacing, up. */
  vector<double> dz_down; /**< Z grid spacing, down. */
  
  double LB_ax, LB_ay, LB_az; //dimensions of Local Bubble
  string LB_shape, LB_smearing;
};

class TGrid2D : public TGrid {

 public:
  TGrid2D() { } /**< Default constructor. */
  TGrid2D(Input*); 
  /**< Constructor given some input.
   * @param Ekfact Delta E / E.
   * @param numr Dimension of the radial grid.
   * @param numz Dimension of the vertical grid.
   * @param zmax_ Height of the simulation box.
   */

  ~TGrid2D() { } /**< Destructor. */

 protected:

};

class TGrid3D : public TGrid {

 public:
  TGrid3D() { } /**< Default constructor. */
  TGrid3D(Input*);
  /**< Constructor given some input.
   * @param Ekfact Delta E / E.
   * @param numr Dimension of the radial grid.
   * @param numz Dimension of the vertical grid.
   * @param zmax_ Height of the simulation box.
   */

  ~TGrid3D() { } /**< Destructor. */


};
#endif
