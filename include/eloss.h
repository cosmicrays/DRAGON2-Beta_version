/**
 * @file eloss.h
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief In this file all the classes related to energy losses are defined.
 */

#ifndef _ELOSS_H
#define _ELOSS_H

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

using namespace std;

class TGrid;
class TGas;
class TInput;
class TNucleiList;
class TBField;

/**
 * @class TISRF
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Description of Interstellar Radiation Field, according to Galprop model. CMB, IR and starlight are included. See ApJ...
 */
class TISRF {
    
public:
    TISRF() {} /**< Default constructor. */
    
    TISRF(TGrid* /**< Geometry description */, string /**< File containing the ISRF*/, TGeometry*, Input*);
    /**
     * @fn TISRF(TGrid*, string)
     * @brief Read the ISRF contained in file and interpolate to fit simulation box. Taken from Galprop package.
     */
    
    virtual ~TISRF() { ISRField.clear(); e_loss_Compton.clear(); nu_array.clear(); }
    /**< Destructor. */
    
    inline vector<double>& GetISRF() { return ISRField; } /**< Return the ISRF field distribution. */
    inline double  GetISRF(int ir /**< radial position */, int iz /**< vertical position */, int inu /**< frequency position */) { return ISRField[index(inu,ir,iz)]; }
    inline double  GetISRF(int ix /**< radial position */, int iy, int iz /**< vertical position */, int inu /**< frequency position */) { return ISRField[index(inu,ix,iy,iz)]; }
    /**
     * @fn inline double  GetISRF(int ir, int iz, int inu)
     * @brief Get ISRF at given position and frequency, with wrapping between matrix and linear representation.
     * @return ISRF density.
     */
    inline vector<double>& GetNuArray() { return nu_array; } /**< Get array of frequencies. */
    inline double GetDnu() const { return Dnu; } /**< Get (logarithmic) spacing in frequency space. */
    inline double GetElossCompton(int inu /**< Frequency position */, int iE /**< Energy position */) { return e_loss_Compton[iE*dimnu+inu]; }
    /**
     * @fn inline double GetElossCompton(int inu, int iE)
     * @brief Get the energy loss rate of electrons/positrons of energy index iE for Compton scattering onto ISRF at frequency index inu.
     */
    
protected:
    vector<double> ISRField; /**< ISRF density. */
    vector<double> e_loss_Compton; /**< Energy loss for Compton scattering. */
    vector<double> nu_array; /**< Array of frequency of ISRF. */
    double Dnu; /**< Logarithmic spacing of frequency. */
    int dimx; /**< radial dimension of simulation box. */
    int dimy; /**< radial dimension of simulation box. */
    int dimz; /**< vertical dimension of simulation box. */
    int dimr_isrf; /**< radial dimension of ISRF representation in file. */
    int dimz_isrf; /**< vertical dimension of ISRF representation in file. */
    int dimnu;  /**< frequency dimension of ISRF representation in file. */
    int ncomp; /**< Components accounted for in ISRF. */
    inline int isrf_index(int ir, int iz, int inu, int icomp) { return ((icomp*dimnu +inu)*dimz_isrf + iz)*dimr_isrf + ir; }
    /**
     * @fn inline int index(int ir, int iz, int inu, int icomp)
     * @brief Convert matrix to linear representation.
     * @return ((icomp*dimnu +inu)*dimz_isrf + iz)*dimr_isrf + ir
     */
    
    inline int index(int inu, int ir, int iz) { return (inu*dimx+ir)*dimz+iz; }
    /**
     * @fn inline int index(int inu, int ir, int iz)
     * @brief Convert matrix to linear representation.
     * @return (inu*dimr+ir)*dimz+iz
     */
    inline int index(int inu, int ix, int iy, int iz) { return ((inu*dimx+ix)*dimy+iy)*dimz+iz; }
    /**
     * @fn inline int index(int inu, int ir, int iz)
     * @brief Convert matrix to linear representation.
     * @return (inu*dimr+ir)*dimz+iz
     */
    
    inline long double f1(double z /**< */) {
        
        return 
        ((z+6.0+3.0/z)*gsl_sf_log_1plusx(2.0*z)      // !!!!! Using log(1+2*z) gives largely wrong results !!!!!
         - (22.0*gsl_pow_3(z)/3.0 + 24.0*gsl_pow_2(z) + 18.0*z + 4)/gsl_pow_2(1.0+2.0*z)
         -2.0
         + 2.0*gsl_sf_dilog(-2.0*z));
    }
    /**
     * @fn  inline long double f1(double z)
     * @brief Function appearing in computing IC energy losses. See Jones...
     * @warning Must use GSL high precision functions to get the correct result.
     */
    inline long double f2(double z /**< */) {
        return 
        (z+31.0/6.0+5.0/z+1.5/z/z)*gsl_sf_log_1plusx(2.0*z)
        -(22.0/3.0*gsl_pow_3(z) + 28.0*gsl_pow_2(z) + 103.0/3.0*z + 17.0 + 3.0/z)/gsl_pow_2(1.0+2.0*z)
        -2.0
        +gsl_sf_dilog(-2.0*z);
    }
    /**
     * @fn  inline long double f2(double z)
     * @brief Function appearing in computing IC energy losses. See Jones...
     * @warning Must use GSL high precision functions to get the correct result.
     */
};

/**
 * @class TEnergyLoss
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief General class describing energy losses in DRAGON.
 */
class TEnergyLoss {
public:
    virtual ~TEnergyLoss() { dpdt.clear(); }
    /**< Destructor. */
    TEnergyLoss(TGrid* coord, Input* in);
    /**< Constructor given a geometry. */
    
    TEnergyLoss(const TEnergyLoss& el) :
    dimE(el.dimE),
    dimx(el.dimx),
    dimy(el.dimy),
    dimz(el.dimz),
    dpdt(el.dpdt)
    { }
    /**< Copy constructor. */
    
    inline int GetDimr() { return dimx; } /**< Get radial dimension of simulation box. */
    inline int GetDimx() { return dimx; } /**< Get radial dimension of simulation box. */
    inline int GetDimy() { return dimy; } /**< Get radial dimension of simulation box. */
    inline int GetDimz() { return dimz; } /**< Get vertical dimension of simulation box. */
    inline int GetDimE() { return dimE; } /**< Get energy dimension of simulation box. */
    inline vector<double>& GetDpdt() { return dpdt; }  /**< Get energy loss distribution. */
    inline double GetDpdt(int ir /**< Radial position */, int iz /**< Vertical position */, int ip /**< Energy position */) { return dpdt[index(ir,iz,ip)]; }
    inline double GetDpdt(int ix /**< Radial position */, int iy, int iz /**< Vertical position */, int ip /**< Energy position */) { return dpdt[index(ix,iy,iz,ip)]; }
    /**
     * @fn inline virtual double GetDpdt(int ix, int iy, int iz, int ip)
     * @brief Get energy loss at given position and energy with matrix to linear representation conversion.
     */
    inline double GetDpdt(int i) { return dpdt[i]; } /**< Get energy loss at given position and energy in linear representation. */
    
    inline void Set_dpdt_Zero() {for(int i=0;i<dpdt.size(); ++i) dpdt[i]=0.;} //fk 130701

    TEnergyLoss& operator += (const TEnergyLoss& el) {
        for (int i = 0; i < dpdt.size(); ++i) dpdt[i] += el.dpdt[i];
        return *this;
    }
    /**< Sum up several energy loss terms. */
    
protected:
    int dimE; /**< Energy dimension of simulation box. */ 
    int dimx; /**< Radial dimension of simulation box. */
    int dimy; /**< Radial dimension of simulation box. */
    int dimz; /**< Vertical dimension of simulation box. */ 
    vector<double> dpdt; /**< Energy loss distribution. */
    inline int index(int ir, int iz, int ip) { return (ir*dimz+iz)*dimE+ip; }
    /**
     * @fn inline int index(int ir, int iz, int ip)
     * @brief Convert matrix to linear representation.
     * @return (ir*dimz+iz)*dimE+ip
     */
    inline int index(int ix, int iy, int iz, int ip) { return ((ix*dimy + iy)*dimz+iz)*dimE+ip; }
    /**
     * @fn inline int index(int ix, int iy, int iz, int ip)
     * @brief Convert matrix to linear representation.
     */
    Input* inp;
};

/**
 * @class TIonizationLoss
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Ionization energy losses.
 */
class TIonizationLoss : public TEnergyLoss {
public:
    TIonizationLoss(TGrid* /**< geometry */, vector<TGas*> /**< ISM gas components */, TGas* /**< Total ISM gas */, Input*, double /**< Nucleus mass */, double /**< Nucleus charge */);
    /**
     * @fn TIonizationLoss(TGrid*, vector<TGas*>, TGas*, double, double)
     * @brief Initialize ionization losses.
     */
    
    virtual ~TIonizationLoss() { }
    /**< Destructor. */
    
protected: 
    
};

/**
 * @class TCoulombLoss
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Coulomb energy losses
 */
class TCoulombLoss : public TEnergyLoss {
public:
    TCoulombLoss(TGrid* /**< geometry */, vector<TGas*> /**< ISM gas components */, Input*, double /**< Nucleus mass */, double /**< Nucleus charge */);
    /**
     * @fn TCoulombLoss(TGrid*, vector<TGas*>, double, double)
     * @brief Initialize Coulomb losses.
     */ 
    
    virtual ~TCoulombLoss() { }
    /**< Destructor. */
    
protected:
    
};

/**
 * @class TBremsstrahlungLoss
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Bremsstrahlung energy losses
 */
class TBremsstrahlungLoss : public TEnergyLoss {
public:
    TBremsstrahlungLoss(TGrid* /**< geometry */, vector<TGas*> /**< ISM gas components */, TGas* /**< Total ISM gas */, Input*);
    /**
     * @fn TBremsstrahlungLoss(TGrid*, vector<TGas*>, double, double)
     * @brief Initialize Bremsstrahlung losses (only leptons).
     */ 
    
    virtual ~TBremsstrahlungLoss() { }
    /**< Destructor. */
    
protected:
    
};

/**
 * @class TSynchrotronLoss
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Synchrotron energy losses
 */
class TSynchrotronLoss : public TEnergyLoss {
public:
    TSynchrotronLoss(TGrid* /**< geometry */, TBField* /**< Magnetic Field*/, Input*);
    /**
     * @fn TSynchrotronLoss(TGrid*, TBField*)
     * @brief Initialize Synchrotron losses (only leptons).
     */ 
    
    virtual ~TSynchrotronLoss() { }
    /**< Destructor. */
    
protected: 
};

/**
 * @class TICSLoss
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Inverse Compton Scattering energy loss
 */
class TICSLoss : public TEnergyLoss {
public:
    TICSLoss(TGrid* /**< geometry */, TISRF* /**< ISRF */, Input*);
    /**
     * @fn TICSLoss(TGrid*, TISRF*)
     * @brief Initialize ICS losses (only leptons).
     */ 
    
    virtual ~TICSLoss() { }
    /**< Destructor. */
    
protected: 
};

#endif
