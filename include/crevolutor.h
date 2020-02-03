/**
 * @file crevolutor.h
 * @author Luca Maccione
 * @author Daniele Gaggero
 * @email luca.maccione@desy.de
 * @email daniele.gaggero@sissa.it
 * @brief Classes for the solution of the transport equation are defined.
 */

#ifndef _CREVOLUTOR_H
#define _CREVOLUTOR_H

#include <iostream>
#include <vector>

#include "fitsio.h"
#include "galaxy.h"

using namespace std;

class TGas;
class TGrid;
class TSource;
class TConvectionVelocity;
class TDiffusionCoefficient;
class TReaccelerationCoefficient;
class Input;
class TInelasticCrossSection;
class TXSecBase;

/**
 * @class TCREvolutorBasis
 * @author Luca Maccione
 * @email luca.maccione@lmu.de
 * @brief Abstract class for the solution of the transport equation.
 *
 * This class defines all the methods needed to solve the transport equation, and that a concrete implementation of a solution algorithm should implement. This structure allows to use several algorithm in cascade. For example, one can use the fast Cranck-Nicholson method a la Galprop implemented in TCREvolutor to quickly obtain a rough approximation to the solution, followed by a more reliable, but slower, ADI version of this method, implemented in TCREvolutorADI. In such a way one is able to obtain reliable solutions in a reasonable amount of time.
 */

class TCREvolutorBasis {
    
public:
    
    //modified
    virtual void Run(vector<double>&, vector<double>&, TInelasticCrossSection*, const vector<double>&, const vector<double>&, const vector<double>&, double, double, double, int, bool SecEl=false, int K_electron=0, bool isDM=false, bool isextra=false)=0;

    virtual void Run(vector<double>&, vector<double>&, TInelasticCrossSection*, const vector<double>&, const vector<double>&, const vector<double>&, TDiffusionCoefficient*, TReaccelerationCoefficient*, double, double, double, int, bool SecEl=false, int K_electron=0, bool isDM=false, bool isextra=false)=0;

    virtual void Run(vector<double>&, vector<double>&, TInelasticCrossSection*, const vector<double>&, const vector<double>&, const vector<double>&, TDiffusionCoefficient*, TReaccelerationCoefficient*, double, double, double, int, bool SecEl=false, int K_electron=0, bool isDM=false, bool isextra=false, Galaxy* gal=NULL)=0;
    
    /**< Actually solve the transport equation and return the propagated density. */
    virtual ~TCREvolutorBasis() {
        riac1.clear();
        riac2.clear();
        riac3.clear();
        Pdotup.clear();
        Pdotdown.clear();
}
    /**< Destructor. */
    
    double CalcInjFactor(int K_electron, bool SecEl, bool isDM, double A, double Z) {
        if (K_electron > 0) return 0.; // the nucleus that already attached an electron only has a secondary source term!!
        if (SecEl) return 0.0;
        else if (isDM) return 1.0;
        else if (A > 0) return gal->GetSourceAbundance(int(A+1000*Z))/A;
        else return 1.0;
    }

    TSource* GetSourceTerm(bool isDM, bool isextra) {
        if (isDM) return gal->GetDMSource();
        else if (isextra) return gal->GetSourceExtra();
        else return gal->GetSource();
    }

protected:
    int dimr; /**< Radial dimension of spatial grid. */
    int dimx; /**< Radial dimension of spatial grid. */
    int dimy; /**< Radial dimension of spatial grid. */
    int dimz; /**< Vertical dimension of spatial grid. */
    int dimE; /**< Dimension of energy grid. */
    
    vector<double> riac1; /**< Matrix related to reacceleration. */
    vector<double> riac2; /**< Matrix related to reacceleration. */
    vector<double> riac3; /**< Matrix related to reacceleration. */
    vector<double> Pdotup; /**< Matrix related to energy loss. */
    vector<double> Pdotdown; /**< Matrix related to energy loss. */
    
    
    TGrid* coord; /**< Pointer to kinematics and geometry. */
    Galaxy* gal; /**< Pointer to galaxy structure. */
    Input* in;
    inline int index(int ir, int iz, int ip) { return (ir*dimz+iz)*dimE+ip; }
    inline long index(int ix, int iy, int iz, int ip) { return ((ix*dimy+iy)*dimz+iz)*dimE+ip; }
    /**
     * @fn inline virtual int index(int ir, int iz, int ip)
     * @brief Convert from matrix to linearized indexes.
     * @param ir Radial position
     * @param iz Vertical position
     * @param ip Energy position
     * @return (ir*dimz+iz)*dimE+ip
     */
};

/**
 * @class TCREvolutor
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class for the solution of the transport equation, derived from TCREvolutorBasis.
 *
 * This class implements a second order operator splitting implicit Cranck-Nicholson method. In order to speed up the solution, and also to compare with Galprop, we start with a very large time step, and after some iterations we reduce the time step. Minimum time step has to be of order of 100 yr to have reasonable solutions for nuclei, and of order of 1 yr for leptons. \sa Nrept \sa dtfactor \sa dtmin \sa dtmax
 */
class TCREvolutor : public TCREvolutorBasis {
    
public:
    TCREvolutor(Galaxy*);  /**< Constructor, given a model for the galaxy. */
    virtual ~TCREvolutor() {    } /**< Destructor. */
    
    //modified
    virtual void Run(vector<double>&, vector<double>&, TInelasticCrossSection* /**< Inelastic cross section */, const vector<double>& /**< Energy losses */, const vector<double>& /**< Secondary source contribution */, const vector<double>& /**< Injection spectrum */, double /**< Mass number */, double /**< Charge */, double /**< Lifetime */, int /**< Unique ID of daughter nucleus */, bool SecEl=false /**< Whether it is a secondary electron, or proton, in which case the primary source term has to be 0. Default is false. */, int K_electron=0/*K_electron*/, bool isDM=false, bool isextra=false);
    virtual void Run(vector<double>&, vector<double>&, TInelasticCrossSection* /**< Inelastic cross section */, const vector<double>&, const vector<double>&, const vector<double>&, TDiffusionCoefficient*, TReaccelerationCoefficient*, double, double, double, int, bool SecEl=false, int K_electron=0, bool isDM=false, bool isextra=false) { /*...*/ }
    virtual void Run(vector<double>&, vector<double>&, TInelasticCrossSection* /**< Inelastic cross section */, const vector<double>&, const vector<double>&, const vector<double>&, TDiffusionCoefficient*, TReaccelerationCoefficient*, double, double, double, int, bool SecEl=false, int K_electron=0, bool isDM=false, bool isextra=false, Galaxy* gal = NULL) { /*...*/ }
};

/**
 * @class TCREvolutorADI
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class for the solution of the transport equation, derived from TCREvolutorBasis.
 *
 * This class implements a second order operator splitting alternate-direction implicit (ADI) Cranck-Nicholson method. It essentially reimplements the method TCREvolutor::Run(...).
 */
class TCREvolutorADI : public TCREvolutorBasis {
    
public:
    TCREvolutorADI(Galaxy*); /**< Constructor, given a model for the galaxy. */
    
    virtual ~TCREvolutorADI() {
        
        N1.clear();
        N2.clear();
        riac1.clear();
        riac2.clear();
        riac3.clear();
        Pdotup.clear();
        Pdotdown.clear();
        /*
         Rzz.clear();
         dzz.clear();
         lodzz.clear();
         uodzz.clear();
         yy.clear();
         
         Re.clear();
         de.clear();
         odeu.clear();
         oded.clear();
         ee.clear();
         
         Rrr.clear();
         drr.clear();
         lodrr.clear();
         uodrr.clear();
         xx.clear();
         */
    } /**< Default destructor. */
    
    virtual void Run(vector<double>&, vector<double>&, TInelasticCrossSection* /**< Inelastic cross section */, const vector<double>& /**< Energy losses */, const vector<double>& /**< Secondary source contribution */, const vector<double>& /**< Injection spectrum */, double /**< Mass number */, double /**< Charge */, double /**< Lifetime */, int /**< Unique ID of daughter nucleus */, bool SecEl=false /**< Whether it is a secondary electron, or proton, in which case the primary source term has to be 0. Default is false. */, int K_electron=0 /** K electron */, bool isDM=false, bool isextra=false);
    
    virtual void Run(vector<double>&, vector<double>&, TInelasticCrossSection* /**< Inelastic cross section */, const vector<double>&, const vector<double>&, const vector<double>&, TDiffusionCoefficient*, TReaccelerationCoefficient*, double, double, double, int, bool SecEl=false, int K_electron=0, bool isDM=false, bool isextra=false) { /*...*/ }

    
    virtual void Run(vector<double>&, vector<double>&, TInelasticCrossSection* /**< Inelastic cross section */, const vector<double>&, const vector<double>&, const vector<double>&, TDiffusionCoefficient*, TReaccelerationCoefficient*, double, double, double, int, bool SecEl=false, int K_electron=0, bool isDM=false, bool isextra=false, Galaxy* gal = NULL) { /*...*/ }

protected:
    double tolerance;
    double dtADI;
    vector<double> N1;
    vector<double> N2;

    //TGas* totalgas; /**< Pointer to total gas density. */
    //TSource* source; /**< Pointer to primary source spatial distribution. */
    //TConvectionVelocity* vC; /**< Convection velocity. */
    
//    vector<double> riac1; /**< Matrix related to reacceleration. */
//    vector<double> riac2; /**< Matrix related to reacceleration. */
//    vector<double> riac3; /**< Matrix related to reacceleration. */
//    vector<double> Pdotup; /**< Matrix related to energy loss. */
//    vector<double> Pdotdown; /**< Matrix related to energy loss. */
    
    //  vector<double> Rzz; /**< Propagation in Z direction. R term. */
    //  vector<double> dzz; /**< Propagation in Z direction. Diagonal term. */
    //  vector<double> lodzz; /**< Propagation in Z direction. Lower diagonal term. */
    //  vector<double> uodzz; /**< Propagation in Z direction. Upper diagonal term. */
    //  vector<double> yy; /**< Propagation in Z direction. Solution term. */
    
    //  vector<double> Re; /**< Propagation in Momentum direction. R term. */
    //  vector<double> de; /**< Propagation in Momentum direction. Diagonal term. */
    //  vector<double> odeu; /**< Propagation in Momentum direction. Upper diagonal term. */
    //  vector<double> oded; /**< Propagation in Momentum direction. Lower diagonal term. */
    //  vector<double> ee; /**< Propagation in Momentum direction. Solution term. */
    
    //  vector<double> Rrr; /**< Propagation in R direction. R term. */
    //  vector<double> drr; /**< Propagation in R direction. Diagonal term. */
    //  vector<double> lodrr; /**< Propagation in R direction. Lower diagonal term. */
    //  vector<double> uodrr; /**< Propagation in R direction. Upper diagonal term. */
    //  vector<double> xx; /**< Propagation in R direction. Solution term. */
    
    double FindMax(TDiffusionCoefficient* dperp /**< Diffusion coefficient. */);
    /**
     * @fn double FindMax(TDiffusionCoefficient* dperp)
     * @brief Find the maximum of the diffusion coefficient, to set a reasonable dt.
     */
    double FindTLoss(const vector<double>& dpdt, vector<double>& totmomentum);
    double FindTRiacc(TReaccelerationCoefficient* dpp, vector<double>& totmomentum);
    double FindTInt(TGas* gas, const vector<double>& xsec);
};


// MODIFIED 17-09-2012

// *********************************************************************************************************************************************************
// *********************************************************************************************************************************************************
// *********************************************************************************************************************************************************

/**
 * @class TCREvolutor3D
 * @author Daniele Gaggero
 * @author Luca Maccione
 * @email daniele.gaggero@sissa.it
 * @brief 3D solution of the transport equation.
 *
 * Description to be added
 */

class TCREvolutor3D : public TCREvolutorBasis {
    
public:
    TCREvolutor3D(Galaxy*); /**< Constructor, given a model for the galaxy. */
    
    virtual void Run(vector<double>&, vector<double>&, TInelasticCrossSection* /**< Inelastic cross section */, const vector<double>&, const vector<double>&, const vector<double>&, double, double, double, int, bool SecEl=false, int K_electron=0 /** K electron */, bool isDM=false, bool isextra=false);  /**< Actually solve the transport equation and return the propagated density. */

    virtual void Run(vector<double>&, vector<double>&, TInelasticCrossSection* /**< Inelastic cross section */, const vector<double>&, const vector<double>&, const vector<double>&, TDiffusionCoefficient*, TReaccelerationCoefficient*, double, double, double, int, bool SecEl=false, int K_electron=0, bool isDM=false, bool isextra=false) { /*****/ }

    virtual void Run(vector<double>&, vector<double>&, TInelasticCrossSection* /**< Inelastic cross section */, const vector<double>&, const vector<double>&, const vector<double>&, TDiffusionCoefficient*, TReaccelerationCoefficient*, double, double, double, int, bool SecEl=false, int K_electron=0, bool isDM=false, bool isextra=false, Galaxy* gal = NULL);// { /**< Solve the 3D Anisotropic equation, return the propagated density. */ }
    
    virtual ~TCREvolutor3D() {
        
      }
};


// *********************************************************************************************************************************************************
// *********************************************************************************************************************************************************
// *********************************************************************************************************************************************************


class TCREvolutor3DADI : public TCREvolutorBasis {
    
public:
    TCREvolutor3DADI(Galaxy*) { } /**< Constructor, given a model for the galaxy. */
    
    virtual void Run(vector<double>&, vector<double>&, TInelasticCrossSection* /**< Inelastic cross section */, const vector<double>&, const vector<double>&, const vector<double>&, double, double, double, int, bool SecEl=false, int K_electron=0 /** K electron */, bool isDM=false, bool isextra=false) { } /**< Actually solve the transport equation and return the propagated density. */
    virtual void Run(vector<double>&, vector<double>&, TInelasticCrossSection* /**< Inelastic cross section */, const vector<double>&, const vector<double>&, const vector<double>&, TDiffusionCoefficient*, TReaccelerationCoefficient*, double, double, double, int, bool SecEl=false, int K_electron=0, bool isDM=false, bool isextra=false) { /*...*/ }
    virtual void Run(vector<double>&, vector<double>&, TInelasticCrossSection* /**< Inelastic cross section */, const vector<double>&, const vector<double>&, const vector<double>&, TDiffusionCoefficient*, TReaccelerationCoefficient*, double, double, double, int, bool SecEl=false, int K_electron=0, bool isDM=false, bool isextra=false, Galaxy* gal=NULL) { /*...*/ }
    virtual ~TCREvolutor3DADI() { }
};

#endif
