/**
 * @file galaxy.h
 * @author Luca Maccione, Daniele Gaggero
 * @email luca.maccione@desy.de
 * @email daniele.gaggero@sissa.it
 * @brief In this file all the classes related to the model of the galaxy are defined.
 */

#ifndef _GALAXY_H
#define _GALAXY_H

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>

#include "constants.h"

using namespace std;

class TGrid;
class TGeometry;
class TGas;
class Input;
class TNucleiList;
class TISRF;
class TBField;
class TSource;
class TDiffusionCoefficient;

/**
 * @class TConvectionVelocity
 * @author Luca Maccione & Daniele Gaggero
 * @email luca.maccione@lmu.de
 * @brief Initialize Convection velocity
 */

class TConvectionVelocity {
    
public:
    TConvectionVelocity() {}
    TConvectionVelocity(TGrid* coord, TGeometry*, Input*, TSource*);
    ~TConvectionVelocity() { vc.clear(); }
    
    
    inline vector<double>& GetVC() { return vc; }
    /**< Get the velocity array. */
    inline double GetVC(int ind) { return vc[ind]; }
    /**< Get velocity at fixed vectorized position */
    
    inline double GetVC(int ir, int iz) { return vc[conv_index(ir,iz)];  }
    inline double GetVC(int ix, int iy, int iz) { return vc[conv_index(ix,iy,iz)];  } //MW130620
    inline double GetDvdz() { return dvdz; }
    double GetProfile(double, double, double, TSource*);
    /**< Get increment */
    
    //MW130624
    inline double GetCNconv_alpha1_z(int indspat) { return CNconv_alpha1_z[indspat]; }
    inline double GetCNconv_alpha2_z(int indspat) { return CNconv_alpha2_z[indspat]; }
    inline double GetCNconv_alpha3_z(int indspat) { return CNconv_alpha3_z[indspat]; }

protected:
    int dimr;
    int dimx;
    int dimy;
    int dimz;
    vector<double> vc;
    DPerpType set_profile_conv;
    double dvdz;
    double nrn_sn;
    double conv_index_radial;
    double conv_threshold;
    inline int conv_index(int ir /**< Radial index. */, int iz /**< Vertical index. */)
    {
        if(ir<0) ir = 0;
        if(iz<0) iz = 0;
        if(ir>=dimr) ir = dimr - 1;
        if(iz>=dimz) iz = dimz - 1;
        
        return ir*dimz + iz;
    } 
    inline int conv_index(int ix /**< radial index 1. */, int iy /**< radial index 1. */, int iz /**< Vertical index. */)
    {
        if(ix<0) ix = 0;
        if(iy<0) iy = 0;
        if(iz<0) iz = 0;
        if(ix>=dimx) ix = dimx - 1;
        if(iy>=dimy) iy = dimy - 1;
        if(iz>=dimz) iz = dimz - 1;
        
        return (ix*dimy + iy)*dimz+iz; //MW130711
    }

    //MW130624
    vector<double> CNconv_alpha1_z;
    vector<double> CNconv_alpha2_z;
    vector<double> CNconv_alpha3_z;
};

/**
 * @class TReaccelerationCoefficient
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Position dependent diffusion coefficient in momentum space (reacceleration). Part of it independent of the nucleus.
 */
class TReaccelerationCoefficient {
    
public:
    TReaccelerationCoefficient() {} /**< Default constructor. */
    TReaccelerationCoefficient(vector<double> /**< Momentum. */, TDiffusionCoefficient* /**< Diffusion coefficient. */, TGeometry*, Input*);
    /**< Constructor given vA and the diffusion coefficient. */
    ~TReaccelerationCoefficient() { dpp.clear(); sp.clear(); }
    /**< Destructor. */
    
    inline const vector<double>& GetReaccelerationCoefficient() const { return dpp; }
    /**< Obtain dpp. */
    inline double GetReaccelerationCoefficient(int i /**< Linearized position index. */) { return dpp[i]; }
    /**< Obtain dpp at given linearized position. */
    inline double GetReaccelerationCoefficient(int isp /**< Linearized position index. */, int ip /**< Energy index. */) { return dpp[isp]*sp[ip]; }
    /**< Obtain the diffusion coefficient in momentum space at given linearized position and energy. */
    inline double GetReaccelerationCoefficient(int ir /**< Radial index. */, int iz /**< Vertical index. */, int ip /**< Energy index. */) { return dpp[index(ir,iz)]*sp[ip]; }
    inline double GetReaccelerationCoefficient(int ix /**< Radial index. */, int iy, int iz /**< Vertical index. */, int ip /**< Energy index. */) { return dpp[index(ix,iy,iz)]*sp[ip]; }
    /**< Obtain the diffusion coefficient in momentum space at given position and energy. */
    inline const vector<double>& GetSpectrum() const { return sp; }
    /**< Obtain the energy spectrum. */
    inline double GetSpectrum(int i /**< Energy index. */) { return sp[i]; }
    /**< Obtain the energy spectrum at given energy. */
    
protected:
    vector<double> dpp; /**< Spatial profile. */
    vector<double> sp; /**< Energy spectrum. */
    int dimr; /**< Radial dimension of simulation box. */
    int dimx; /**< Radial dimension of simulation box. */
    int dimy; /**< Radial dimension of simulation box. */
    int dimz; /**< Vertical dimension of simulation box. */
    
    inline int index(int ir /**< Radial index. */, int iz /**< Vertical index. */)
    {
        if(ir<0) ir = 0;
        if(iz<0) iz = 0;
        if(ir>=dimr) ir = dimr - 1;
        if(iz>=dimz) iz = dimz - 1;
        
        return ir*dimz + iz;
    } 
    inline int index(int ix, int iy, int iz) //MW130711: fixing border issues
    {
        if(ix<0) ix = 0;
        if(iy<0) iy = 0;
        if(iz<0) iz = 0;
        if(ix>=dimx) ix = dimx - 1;
        if(iy>=dimy) iy = dimy - 1;
        if(iz>=dimz) iz = dimz - 1;
        
        return (ix*dimy + iy)*dimz+iz;
    } 
    
    /**
     * @fn inline int index(int ir, int iz)
     * @brief convert from matrix to linear representation.
     * @return ir*dimz+iz
     */
};

/**
 * @class Galaxy
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class containing a general description of the galaxy.
 */
class Galaxy {
    
public:
    
    Galaxy() {} /**< Default constructor. */
    Galaxy(Input* /**< User input. */, TNucleiList*); /**< Constructor given user input. */

    bool GetTestMode() { return TESTMODE; };
    bool IsSourceMoving() { return MOVING; };	
    bool IsMovingClump() { return MOVING_CLUMP; };

    double GetSourceX0() { return source_x0; };
    double GetSourceY0() { return source_y0; };
    double GetSourceZ0() { return source_z0; };

    double GetSourceVX() { return source_vx; };
    double GetSourceVY() { return source_vy; };
    double GetSourceVZ() { return source_vz; }; 

    double GetClumpX0() { return clump_x0; };
    double GetClumpY0() { return clump_y0; };
    double GetClumpZ0() { return clump_z0; };

    double GetClumpVX() { return clump_vx; };
    double GetClumpVY() { return clump_vy; };
    double GetClumpVZ() { return clump_vz; }; 
	
    double GetClumpDeltaT() { return clump_deltat; }; 

    TGrid* GetCoordinates() { return _fCoordinates; }
    TGeometry* GetGeometry() { return _fGeometry; }
    /**< Obtain the Geometry description. */
    TSource* GetSource() { return _fSource; }
    /**< Obtain the CR source distribution. */
    TSource* GetSourceExtra() { return _fSourceExtra; }
    /**< Obtain the source distribution for the extra component. */
    TSource* GetDMSource() { return _fDMSource; }
    /**< Obtain the DM source distribution. */
    vector<TGas*>& GetGas() { return _fGas; }
    /**< Obtain the ISM gas components. */
    TGas* GetTotalGas() { return _fTotalGas; }
    /**< Obtain the distribution of total gas in the galaxy. */

    map<int, double> GetGasAbundances() { return _fGasAbundances; }    
    /**< Obtain the relative abundances of the various gas components. */

    TDiffusionCoefficient* GetDiffCoeff() { return _fDperp; }
    /**< Obtain the diffusion coefficient. */
    TReaccelerationCoefficient* GetReaccCoeff() { return _fDpp; }
    /**< Obtain the diffusion coefficient in momentum space. */
    
    TDiffusionCoefficient* GetDiffCoeffEl() { return _fDperpEl; }
    /**< Obtain the diffusion coefficient for electrons/positrons. */
    TReaccelerationCoefficient* GetReaccCoeffEl() { return _fDppEl; }
    /**< Obtain the diffusion coefficient for electrons/positrons in momentum space. */
    TConvectionVelocity* GetVC() { return _fVC; }
    /**< Obtain the convection velocity. */
    
    TISRF* GetISRF() { return _fISRF; }
    /**< Obtain the description of InterStellar Radiation Fields. */
    TBField* GetBField() { return _fB; }
    /**< Obtain the Magnetic Field. */
    
    Input* GetInput() { return inp; }
    void Delete();
    double GetSourceAbundance(int uid /**< Unique ID of nucleus. */) {
        map<int,double>::iterator found = _fSourceAbundances.find(uid);
        if (found != _fSourceAbundances.end()) return (*found).second;
        else return -1;
    }
    /**
     * @fn double GetSourceAbundance(int uid)
     * @brief Obtain source abundance of specified nucleus
     */
    //  const double GetInjLowSpectrum() const { return inj_low_nucl; }
    /**< Obtain the injection spectral index at low energy (below some ref rigidity. */


    //DG.29.09.2013
    vector<double> GetInjSpectrum_rho(int uid /**< Unique ID of nucleus. */) {
	vector<double> result;	
        map<int, vector<double> >::iterator found = _fInjSpectrum_rho.find(uid);
        if (found != _fInjSpectrum_rho.end()) result = (*found).second;
        return result;
    }
    vector<double> GetInjSpectrum_alpha(int uid /**< Unique ID of nucleus. */) {
	vector<double> result;	
        map<int, vector<double> >::iterator found = _fInjSpectrum_alpha.find(uid);
        if (found != _fInjSpectrum_alpha.end()) result = (*found).second;
        return result;
    }
    

    //double GetInjSpectrum_rho_0(int uid /**< Unique ID of nucleus. */) {
    //    map<int,double>::iterator found = _fInjSpectrum_rho_0.find(uid);
    //    if (found != _fInjSpectrum_rho_0.end()) return (*found).second;
    //    else return 0;
    //}
    //double GetInjSpectrum_rho_1(int uid /**< Unique ID of nucleus. */) {
    //    map<int,double>::iterator found = _fInjSpectrum_rho_1.find(uid);
    //    if (found != _fInjSpectrum_rho_1.end()) return (*found).second;
    //    else return 0;
    //}
    //double GetInjSpectrum_rho_2(int uid /**< Unique ID of nucleus. */) {
    //    map<int,double>::iterator found = _fInjSpectrum_rho_2.find(uid);
    //    if (found != _fInjSpectrum_rho_2.end()) return (*found).second;
    //    else return 0;
    //}
    //double GetInjSpectrum_alpha_0(int uid /**< Unique ID of nucleus. */) {
    //    map<int,double>::iterator found = _fInjSpectrum_alpha_0.find(uid);
    //    if (found != _fInjSpectrum_alpha_0.end()) return (*found).second;
    //    else return 0;
    //}
    //double GetInjSpectrum_alpha_1(int uid /**< Unique ID of nucleus. */) {
    //    map<int,double>::iterator found = _fInjSpectrum_alpha_1.find(uid);
    //    if (found != _fInjSpectrum_alpha_1.end()) return (*found).second;
    //    else return 0;
    //}
    //double GetInjSpectrum_alpha_2(int uid /**< Unique ID of nucleus. */) {
    //    map<int,double>::iterator found = _fInjSpectrum_alpha_2.find(uid);
    //    if (found != _fInjSpectrum_alpha_2.end()) return (*found).second;
    //    else return 0;
    //}
    //double GetInjSpectrum_alpha_3(int uid /**< Unique ID of nucleus. */) {
    //    map<int,double>::iterator found = _fInjSpectrum_alpha_3.find(uid);
    //    if (found != _fInjSpectrum_alpha_3.end()) return (*found).second;
    //    else return 0;
    //}



    /**
     * @fn double GetInjSpectrum(int uid)
     * @brief Obtain high energy injection index of specified nucleus
     */

/*    
    double GetInjLowSpectrum(int uid) {
        map<int,double>::iterator found = _fInjSpectrum_low.find(uid);
        if (found != _fInjSpectrum_low.end()) return (*found).second;
        else return 0;
    }
*/
    /**
     * @fn double GetInjLowSpectrum(int uid)
     * @brief Obtain low energy injection index of specified nucleus
     */
    
    ~Galaxy();
    /**< Destructor. */
    
protected:

    //OLD structure with 2 breaks ----------------------------------------------------------------
    //map<int /**< Unique ID of nucleus. */, double /**< Injection index. */> _fInjSpectrum_rho_0;
    //map<int /**< Unique ID of nucleus. */, double /**< Injection index. */> _fInjSpectrum_rho_1;
    //map<int /**< Unique ID of nucleus. */, double /**< Injection index. */> _fInjSpectrum_rho_2;
    //map<int /**< Unique ID of nucleus. */, double /**< Injection index. */> _fInjSpectrum_alpha_0;
    //map<int /**< Unique ID of nucleus. */, double /**< Injection index. */> _fInjSpectrum_alpha_1;
    //map<int /**< Unique ID of nucleus. */, double /**< Injection index. */> _fInjSpectrum_alpha_2;
    //map<int /**< Unique ID of nucleus. */, double /**< Injection index. */> _fInjSpectrum_alpha_3;
    //end of OLD structure with 2 breaks -----------------------------------------------------------

    //DG.29.09.2013
    //Each nucleus is associated to a vector of injection slopes and another vector containing the break rigidities
    //The number of breaks is arbitrary!!
    map<int, vector<double> > _fInjSpectrum_rho;  
    map<int, vector<double> > _fInjSpectrum_alpha;
    
    map<int /**< Unique ID of nucleus. */, double /**< Source abundance. */> _fSourceAbundances;
    /**
     * @fn map<int, double> _fSourceAbundances;
     * @brief Map associating to each nucleus its corresponding source abundance.
     */
    
    TGrid* _fCoordinates; /**< Geometry and kinematics. */
    TGeometry* _fGeometry;
    TSource* _fSource; /**< CR sources. */
    TSource* _fSourceExtra; /**< CR sources. */
    TSource* _fDMSource; /**< DM sources. */

    vector<TGas*> _fGas; /**< Array of gas components. */
    TGas* _fTotalGas; /**< Total gas. */

    map<int, double> _fGasAbundances; /*relative abundances with respect to H of various gas components */

    TDiffusionCoefficient* _fDperp; /**< Diffusion coefficient. */
    TReaccelerationCoefficient* _fDpp; /**< Reacceleration coefficient. */
    
    // Electron quantities
    TDiffusionCoefficient* _fDperpEl; /**< Diffusion coefficient for electrons/positrons. */
    TReaccelerationCoefficient* _fDppEl; /**< Reacceleration coefficient for electrons/positrons. */
    
    TConvectionVelocity* _fVC; /**< Convection velocity. */
    TISRF* _fISRF; /**< InterStellar Radiation Field. */
    TBField* _fB; /**< Magnetic Field. */
    Input* inp;
    bool TESTMODE;
    bool MOVING;	
    bool MOVING_CLUMP;
    double source_x0;
    double source_y0;
    double source_z0;
    double source_vx;
    double source_vy;
    double source_vz;
    double clump_x0;	
    double clump_y0;	
    double clump_z0;	
    double clump_vx;
    double clump_vy;
    double clump_vz;	
    double clump_deltat;

};

#endif
