/**
 * @file spectrum.h
 * @author Luca Maccione
 * @email luca.maccione@lmu.de
 * @brief In this file all the classes related to the model of the sources are defined.
 */

#ifndef _SOURCE_H
#define _SOURCE_H

#include <iostream>
#include <vector>
#include <string>
#include <math.h>

#include "constants.h"

class Input;
class TGrid;
class TGeometry;

using namespace std;

/**
 * @class TSource
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief General class describing the spatial distribution of CR sources.
 */
class TSource {
    
public:
    TSource() {} /**< Default constructor. */
    
    virtual ~TSource() { source.clear(); }
    /**< Destructor. */
    
    inline vector<double>& GetSource() { return source; }
    /**< Get source vector. */
    inline double GetSource(int i) { return source[i]; }
    /**< Get source density at given position in linearized representation. */
    inline double GetSource(int ir /**< radial position */, int iz /**< vertical position */) { return source[index(ir,iz)]; }
    inline double GetSource(int ix /**< radial position */, int iy, int iz /**< vertical position */) { return source[index(ix,iy,iz)]; }
    /**< Get source density at given position in matrix representation. */
   virtual double GetSource(double x, double y, double z)=0;
 
protected:
    SNRType SNR_model;
    bool isDmsource;

    double Alpha;
    double Beta;
    double Gamma;
    double Rs;
    double rc;
    double rhoc;
    
    Input* in;
    TGrid* Coord;
    double ringmin;
    double ringmax;
    double rings_period;
    double rings_phase;
    double pointsrc_x;
    double pointsrc_y;
    double pointsrc_z;
   	
    vector<double> source; /**< source distribution */
    int dimx; /**< radial dimension of simulation box. */
    int dimy; /**< radial dimension of simulation box. */
    int dimz; /**< vertical dimension of simulation box. */
    inline int index(int ir /**< radial position */, int iz /**< vertical position */) { return ir*dimz+iz; }
    inline int index(int ix /**< radial position */, int iy, int iz /**< vertical position */) { return (ix*dimy+iy)*dimz+iz; }
    /**
     * @fn inline virtual int index(int ir, int iz)
     * @brief convert from matrix to linear representation.
     * @return ir*dimz+iz
     */

};

/**
 * @class TAstrophysicalSource
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class describing the distribution of astrophysical CR sources
 */
class TAstrophysicalSource : public TSource {
    
public:
    TAstrophysicalSource(TGrid* /**< Geometry. */, Input*, TGeometry*, SNRType);
    /**< Constructor given a geometry. */
    
    ~TAstrophysicalSource() { }
    /**< Destructor. */
   virtual double GetSource(double x, double y, double z) { return SourceDistribution(x,y,z); }

protected:
   double SourceDistribution(double x, double y, double z);

};

/**
 * @class TDMSource
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class describing the distribution of DM sources of CRs.
 */
class TDMSource : public TSource {
    
public:
    TDMSource(TGrid* /**< Geometry. */, Input*);
    /**< Constructor given a geometry and properties of DM particle. */
    
    ~TDMSource() { }
    /**< Destructor. */
    //double DM_profile(double, double);
    /**< DM distribution. */
   virtual double GetSource(double x, double y, double z) {double radius = sqrt(x*x+y*y); return DM_profile(radius,z); }
protected:
    double DM_profile_av(double, double, double, double);
    /**< DM distribution averaged over some region, to avoid problems with steep profiles. */
   double DM_profile(double, double);
   /**< DM distribution. */


};

#endif
