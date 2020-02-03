/**
 * @file particle.h
 * @author Luca Maccione
 * @author Daniele Gaggero
 * @email luca.maccione@desy.de
 * @email daniele.gaggero@sissa.it
 * @brief Definition of TParticle class is given.
 */


#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>

#include "fitsio.h"
#include "nucleilist.h"

using namespace std;

class TXSecBase;
class TNucleiList;
class Galaxy;
class TGeometry;
class TCREvolutorBasis;
class Input;
class TEnergyLoss;
class TSpectrum;
class TSpallationNetwork;
class TInelasticCrossSection;
class TDiffusionCoefficient;
class TReaccelerationCoefficient;

/**
 * @class TParticle
 * @author Luca Maccione
 * @email luca.maccione@lmu.de
 * @brief This class is how the properties of a nucleus, its spatial distribution and energy spectrum are described in DRAGON.
 
 * A nucleus is defined to have mass and charge, a unique ID computed from them, possibly some decay channel and daughter nuclei. Moreover, it has a distribution in physical space and in momentum which has to be found after propagation, and its propagation occurs in a given galaxy, where the nucleus suffers energy losses and inelastic scattering. Of course, nuclei are injected by sources in the ISM with some spectrum.
 */

//modified
class TParticle {
   
protected:
   int A; /**< Mass number */
   int Z; /**< Charge */
   int uid; /**< Unique ID: = 1000*Z+A */
   int dimz; /**< vertical dimension of the box */
   int dimE; /**< energy dimension of the box */
   double lifetime; /**< life time for unstable nuclei (=-1 if stable) */
   DECMODE decmode;
   int daughter; /**< UID of daughter nucleus */
   int issec; /**< Whether is a secondary species (already exist a primary one). */
   int isDM;
   int isExtra;
   int isTPP;
   int K_electron; //modified
   Galaxy* _fGalaxy; /**< galaxy model */
   TInelasticCrossSection* _fInXSec; /**< Inelastic cross section */
   TDiffusionCoefficient* _fDiff;
   TReaccelerationCoefficient* _fDpp;
   TSpectrum* sp; /**< Injection spectrum */
   vector<TEnergyLoss*> eloss; /**< Energy losses (several components) */
   vector<double> density; /**< Spatial-energy distribution of propagated nucleus. */
   vector<double> density_previous;
   //modified
   virtual vector<double> ComputeSecondarySource(vector<TParticle*> /**< Array of already propagated heavier nuclei, that can possibly generate present nucleus during spallation*/, TSpallationNetwork* spnet, vector<TXSecBase*>)=0;
   /**
    * @fn vector<double> ComputeSecondarySource(vector<TParticle*>)
    * @return energy-spatial distribution of secondary source.
    * @brief Compute the contribution from spallation of heavier nuclei to the source distribution of present nucleus.
    * @todo Include also inverse beta decay.
    */
   
   inline int index(int ir /**< Radial position */, int iz /**< Vertical position */, int ip /**< Energy position */) {
      return (ir*dimz+iz)*dimE+ip;
   }
   virtual int index(int ix /**< Radial position */, int iy, int iz /**< Vertical position */, int ip /**< Energy position */)=0;// { /*...*/   }
   /**
    * @fn inline int index(int ir, int iz, int ip)
    * @brief Convert matrix to linear representation.
    * @return (ir*dimz+iz)*dimE+ip
    */
   
   Input* in;
   
public:
   TParticle() {} /**< Default constructor. */
   TParticle(TParticle& part);
   /**< Copy constructor. */
   
   TParticle(int A_ /**< Mass number */, int Z_ /**< Charge */, Galaxy* gal /**< A model for the galaxy */, Input* in_ /**< User input */, bool issec, vector<TXSecBase*> xsecmodel, TNucleiList* l, int K_electron, bool isDM, bool isextra, bool isTPP_=false); //modified
   /**< The constructor normally used in DRAGON */
   
   virtual ~TParticle();
   /**< Destructor */
   
   TParticle& operator += (const TParticle& addendum) {
      if (uid != addendum.uid) {
         cerr << "You are summing apples with bananas. Check your routines!" << endl;
         return *this;
      }
      for (unsigned int i = 0; i < density.size(); ++i)
         density[i] += addendum.density[i];
      return *this;
   }
   /**< To add nuclei of the same species but with different origin (e.g. primary and secondary protons, to obtain normalization). */
   
   //modified
   void Evolve(vector<TParticle*> /**< Array of already propagated nuclei, to obtain the secondary source contribution. */, vector<TCREvolutorBasis*> /**< Array of solution algorithms, to use them in cascade. */, TSpallationNetwork* spnet, vector<TXSecBase*>, bool isSecondIteration=false);
   /**< Propagate nuclear density from sources. */
   
   virtual void Print(fitsfile* /**< Output file */, double /**< Normalization factor*/)=0; /**< Print propagated density and spectrum to file. */
   virtual void PrintSpectrum(fitsfile* /**< Output file */, double /**< Normalization factor*/)=0; /**< Print propagated spectrum to file. */
   //virtual vector <double> GetSpectra(double /**< Normalization factor*/,string,int,double* ) {//SK to grab the spectra/}

   virtual vector<double> GetSpectrumAtSunPosition()=0;
   virtual double GetFluxAtSunPosition(int)=0;
   
   virtual double FindNormalization(const double&, const double&)=0;
   /**
    * @fn double FindNormalization()
    * @brief Find normalization factor, by comparing with observed proton flux.
    * @return the normalization factor.
    * @warning Use GSL cspline interpolator.
    */
   
   //virtual double FindNormalization(TParticle*, const double&)=0;
   /**
    * @fn double FindNormalization()
    * @brief Find normalization factor, by comparing with observed proton flux.
    * @return the normalization factor.
    * @warning Use GSL cspline interpolator.
    */
   
   void PrintELosses();
   
   inline TSpectrum* GetSpectrum() {
      return sp;
   } /**< Get injection spectrum. */
   inline double GetDensity(int ir /**< Radial position */, int iz /**< Vertical position */, int ip /**< Energy position */) {
      return density[index(ir,iz,ip)];
   }
   inline double GetDensity(int ix /**< Radial position */, int iy, int iz /**< Vertical position */, int ip /**< Energy position */) {
      return density[index(ix,iy,iz,ip)];
   }
   /**< Return density at given position and energy. Wrap from matrix to linear representation. */
   inline vector<double>& GetDensity() {
      return density;
   } /**< Get density/spectrum array. */
   inline double GetDensity(int i /**< Linearized index */) {
      return density[i];
   }  /**< Return density at given position and energy in linear representation. */
   
   inline int GetUid() {
      return uid;
   } /**< Get unique nucleus ID = 1000*Z+A */
   inline double GetLifetime() const {
      return lifetime;
   } /**< Get Life time */
   inline int GetDaughter() const {
      return daughter;
   } /**< Get UID of daughter nucleus. */
   inline int GetA() {
      return A;
   } /**< Get Mass number. */
   inline int GetZ() {
      return Z;
   } /**< Get charge. */
   inline int GetIsSec() {
      return issec;
   }
   inline int IsDM() {
      return isDM;
   }
   inline int IsExtra() {
      return isExtra;
   }
   inline int IsTPP() {
      return isTPP;
   }
   vector<double> GetEnergyGrid() const;
   
};

class TParticle2D : public TParticle {
   
public:
   TParticle2D() {} /**< Default constructor. */
   TParticle2D(int A_ /**< Mass number */, int Z_ /**< Charge */, Galaxy* gal /**< A model for the galaxy */, Input* in_ /**< User input */, bool issec, vector<TXSecBase*> xsecmodel, TNucleiList* l, int K_electron, bool isDM, bool isextra, bool isTPP_=false); /**< The constructor normally used in DRAGON */
   
   ~TParticle2D() { }
   /**< Destructor */
   
   virtual void Print(fitsfile* /**< Output file */, double /**< Normalization factor*/); /**< Print propagated density and spectrum to file. */
   virtual void PrintSpectrum(fitsfile* /**< Output file */, double /**< Normalization factor*/); /**< Print propagated spectrum to file. */

   virtual vector<double> GetSpectrumAtSunPosition();
   virtual double GetFluxAtSunPosition(int);
   
   //virtual vector <double> GetSpectra(double /**< Normalization factor*/,string,int,double*);
   
   virtual double FindNormalization(const double&, const double&);
   /**
    * @fn double FindNormalization()
    * @brief Find normalization factor, by comparing with observed proton flux.
    * @return the normalization factor.
    * @warning Use GSL cspline interpolator.
    */
   //virtual double FindNormalization(TParticle*, const double&);
   
protected:
   int dimr;
   virtual vector<double> ComputeSecondarySource(vector<TParticle*> /**< Array of already propagated heavier nuclei, that can possibly generate present nucleus during spallation*/, TSpallationNetwork* spnet, vector<TXSecBase*>);
   virtual int index(int ix /**< Radial position */, int iy, int iz /**< Vertical position */, int ip /**< Energy position */) { throw -1; }
   /**
    * @fn vector<double> ComputeSecondarySource(vector<TParticle*>)
    * @return energy-spatial distribution of secondary source.
    * @brief Compute the contribution from spallation of heavier nuclei to the source distribution of present nucleus.
    * @todo Include also inverse beta decay.
    */
};


class TParticle3D : public TParticle {
   
public:
   TParticle3D() {} /**< Default constructor. */
   TParticle3D(int A_ /**< Mass number */, int Z_ /**< Charge */, Galaxy* gal /**< A model for the galaxy */, Input* in_ /**< User input */, bool issec, vector<TXSecBase*> xsecmodel, TNucleiList* l, int K_electron, bool isDM, bool isextra, bool isTPP_=false); //modified
   /**< The constructor normally used in DRAGON */
   
   ~TParticle3D() {}
   /**< Destructor */
   
   virtual void Print(fitsfile* /**< Output file */, double /**< Normalization factor*/); /**< Print propagated density and spectrum to file. */
   virtual void PrintSpectrum(fitsfile* /**< Output file */, double /**< Normalization factor*/); /**< Print propagated spectrum to file. */

   virtual vector<double> GetSpectrumAtSunPosition();
   virtual double GetFluxAtSunPosition(int);

   virtual double FindNormalization(const double&, const double&);
   /**
    * @fn double FindNormalization()
    * @brief Find normalization factor, by comparing with observed proton flux.
    * @return the normalization factor.
    * @warning Use GSL cspline interpolator.
    */
   
   //virtual vector <double> GetSpectra(double /**< Normalization factor*/,string,int,double*);
   
   //virtual double FindNormalization(TParticle*, const double&);
   
protected:
   int dimx;
   int dimy;
   
   //modified
   virtual vector<double> ComputeSecondarySource(vector<TParticle*> /**< Array of already propagated heavier nuclei, that can possibly generate present nucleus during spallation*/, TSpallationNetwork* spnet, vector<TXSecBase*>);
   /**
    * @fn vector<double> ComputeSecondarySource(vector<TParticle*>)
    * @return energy-spatial distribution of secondary source.
    * @brief Compute the contribution from spallation of heavier nuclei to the source distribution of present nucleus.
    * @todo Include also inverse beta decay.
    */
   
   virtual int index(int ix /**< Radial position */, int iy, int iz /**< Vertical position */, int ip /**< Energy position */) {
      return ((ix*dimy+iy)*dimz+iz)*dimE+ip;
   }
   /**
    * @fn inline int index(int ir, int iz, int ip)
    * @brief Convert matrix to linear representation.
    * @return (ir*dimz+iz)*dimE+ip
    */
};

#endif
