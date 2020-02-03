/**
 * @file gas.h
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief File where gas classes are defined.
 */

#ifndef _GAS_H
#define _GAS_H

#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

class TGrid;
class TGeometry;
class Input;

/**
 * @class TGas
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief General class describing a gas distribution in the galaxy. Can also describe the total distribution of the gas, sumed over all the components.
 */
class TGas {

 public:
  TGas() { } /**< Default constructor. */
  TGas(TGrid*, Input*); /**< Constructor given the geometry of the galaxy. Used to reproduce exactly the Galprop gas distribution (for comparison). @sa TGrid */
  TGas(const TGas& gas) :
    dimr(gas.dimr),
    dimx(gas.dimx),
    dimy(gas.dimy),
    dimz(gas.dimz),
    density(gas.density)
      { }
    /**< Copy constructor. */

  inline int GetDimr() { return dimr; }
  /**< Returns dimension of radial grid. */
  inline int GetDimx() { return dimx; }
  /**< Returns dimension of X grid. */
  inline int GetDimy() { return dimy; }
  /**< Returns dimension of Y grid. */
  inline int GetDimz() { return dimz; }
  /**< Returns dimension of vertical grid. */
  inline vector<double>& GetGas()  { return density; }
  /**< Returns gas density (as a vector). */
  inline double GetGas(int i)  { return density[i]; }
  /**< Returns gas density at a given point (linearized). */
  inline double GetGas(int i, int j)  { return density[this->index(i,j)]; }
  inline double GetGas(int i, int j, int k)  { return density[this->index(i,j,k)]; }
  /**
   * @fn inline double GetGas(int i, int j)
   * @param i radial position of the given point.
   * @param j vertical position of the given point.
   * @return Gas density at a given point (i,j). 
   */
  TGas& operator += (const TGas& g) {
    for (int i = 0; i < density.size(); ++i) density[i] += g.density[i];
    return *this;
  }
  /**< Sum up several gas components. */

  inline double GetTotalContent() {
    double GasSum=0;
    for (int i = 0; i < dimx * dimy * dimz ; i++) GasSum+=density[i];
    return GasSum;
  }
  /**< MW: Sum over the whole density vector. */

  void Check() {
    for (int i = 0; i < density.size(); ++i) density[i] = max(3e-6,density[i]);
  }
  /**< Make sure that gas density is larger than 3e-6 (for comparison with Galprop). */
  virtual ~TGas() { density.clear(); }
  /**< Destructor. */

    /**
     * @fn double nH2_Gal(double,double)
     * @brief H2 Gas density Galprop models.
     * @param r Radial position
     * @param z Vertical position
     */
    double nH2_Gal(double,double, Input*);
    
    /**
     * @fn double nHI_Gal(double,double)
     * @brief HI Gas density Galprop models.
     * @param r Radial position
     * @param z Vertical position
     */  
    double nHI_Gal(double,double);
    
    /**
     * @fn double nHII_Gal(double,double)
     * @brief HII Gas density Galprop models.
     * @param Rkpc Radial position
     * @param Zkpc Vertical position
     */
    double nHII_Gal(double,double);
    
    /**
     * @fn double nHII_av(double,double,double,double)
     * @brief HII Gas density averaged Galprop models.
     * @param R Radial position
     * @param z Vertical position
     * @param dz Vertical size of region where to average 
     * @param dzz Substep size
     */
    double nHII_av(double, double,double,double);
    
    /**
     * @fn double nH2_av(double,double,double,double)
     * @brief H2 Gas density averaged Galprop models.
     * @param R Radial position
     * @param z Vertical position
     * @param dz Vertical size of region where to average 
     * @param dzz Substep size
     */
    double nH2_av(double, double,double,double, Input*);
    

    /**
     * @fn double nHI_av(double,double,double,double)
     * @brief HI Gas density averaged Galprop models.
     * @param R Radial position
     * @param z Vertical position
     * @param dz Vertical size of region where to average 
     * @param dzz Substep size
     */
    double nHI_av(double, double,double,double);

    double X_CO(double, Input*);

 protected:
  int dimr; /**< Radial dimension of the spatial grid. */
  int dimx; /**< Radial dimension of the spatial grid. */
  int dimy; /**< Radial dimension of the spatial grid. */
  int dimz; /**< Vertical dimension of the spatial grid. */
  vector<double> density; /**< Gas density (linearized). */
  inline int index(int ir, int iz) { return ir*dimz+iz; }
  inline int index(int ix, int iy, int iz) { return (ix*dimy+iy)*dimz+iz; }
  /**
   * @fn inline virtual int index(int ir, int iz)
   * @brief Convert between matrix and linear representation of the gas density.
   * @param ir Radial position
   * @param iz Vertical position
   * @return ir*dimz+iz
   */
};

/**
 * @class TH2Gas
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class describing the H2 distribution in the galaxy. Derived from TGas.
 */
class TH2Gas : public TGas {

 public:
  TH2Gas() {} /**< Default constructor. */
  TH2Gas(TGrid* Coord, Input* in, TGeometry*); /**< Constructor given the galactic geometry. */

  virtual ~TH2Gas() { } /**< Destructor. */

 protected:
};

/**
 * @class THIGas
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class describing the HI distribution in the galaxy. Derived from TGas.
 */
class THIGas : public TGas {

 public:
  THIGas() {} /**< Default constructor. */
  THIGas(TGrid* Coord, Input* in, TGeometry*); /**< Constructor given the galactic geometry. */

  virtual ~THIGas() { } /**< Destructor. */

 protected:

};

/**
 * @class THIIGas
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class describing the HII distribution in the galaxy. Derived from TGas.
 */
class THIIGas : public TGas {

 public:
  THIIGas() {} /**< Default constructor. */
  THIIGas(TGrid* Coord, Input* in, TGeometry*); /**< Constructor given the galactic geometry. */
 
  virtual ~THIIGas() { } /**< Destructor. */

 protected:

};
#endif

