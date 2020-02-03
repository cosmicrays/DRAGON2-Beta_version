/**
 * @file constants.h
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 *
 * @brief Definition of many constants and settings that are used in the code.
 */

#ifndef _CONSTANTS_H
#define _CONSTANTS_H

//#include "../config.h"
//#define DEBUG

#include <string>
#include <math.h>

// DM

/**
 * @enum DMreaction
 * @brief Chose whether to propagate products of DM decay or annihilation.
 */
typedef enum {
  Decay, /**< Set DM decay. */
  Annihilation /**< SetDM annihilation. */
} DMreaction;

/**
 * @enum DMspectrum
 * @brief Chose whether to use DarkSUSY spectra or a Delta function.
 */
typedef enum {
  Delta, /**< Set Delta function at kinetic energy \sa EkDelta. */
  DarkSUSY, /**< Set DarkSUSY spectrum. */
  SelfTable,
  EWCorrections
} DMspectrum;

/**
 * @enum DMprofile
 * @brief Chose the spatial DM profile: Isothermal, NFW, Kra, Moore, Einasto. \sa include/input.h . For the DM density at solar position \sa rhos.
 */
typedef enum {
  ISO=0, /**< Set Isothermal profile. */
  NFW=1, /**< Set NFW profile. */
  Kra=2, /**< Set Kra profile. */
  Moore=3, /**< Set Moore profile. */
  Einasto=4 /**< Einasto profile. */
} DMprofile;

typedef enum {SM96, galprop_2004, galprop_2010, constant, dragon} xco_modes;
typedef enum {Azimuthal, Spiral, Prouda, Wmap, Tkachev, Pshirkov, Farrar, Uniform, Simple, ToyModel} MagFieldModel;
//ToyModel: magnetic field directed along x axis. Useful to test anisotropic diffusion

/**
 * @enum GasType
 * @brief Chose the gas distribution: Bronfmann + Ferriere, Nakanishi \& Sophue or Galprop.
 */
typedef enum {
  BronfFerr, /**< Set Bronfmann + Ferriere distribution. */
  NS, /**< Set Nakanishi \& Sophue distribution. */
  Galprop, /**< Set Galprop distribution. */
  UniformGas
  //ToyModel
} GasType; 

typedef enum {
  GalpropISRF,
  UniformISRF
} ISRFType;		

/**
 * @enum SNRType
 * @brief Chose the SNR (CR source) distribution: Galprop, Ferriere, a Point Source at Galactic Center, a Ring (\sa ringmin, ringmax).
 */
typedef enum {
  Galprop_,
  Const,
  Lorimer, /**< based on PSR catalogue: Lorimer et al., Mon. Not. R. Astron. Soc. 372, 777–800 (2006)  */
  Ferriere,   /**<  based on PSR catalogue + disk stars: K. Ferriere, Rev.Mod.Phys. 73, 1031-1066 (2001) */
  CaseBhattacharya,   /**<  based on SNR catalogue: Case and  Bhattacharya, Astronomy and Astrophysics Supplement, v.120, p.437-440 (1996) */
  GiguereKaspi, /**< Set GiguereKaspi distribution */
  PointSource, /**< Set Point Source at GC. */
  Ring, /**< Set Ring of sources: \sa ringmin \sa ringmax . */
  Rings,
  BlasiSmooth,
  OnlyExtra,
  Custom
} SNRType;

/**
 * @enum DPerpType
 * @brief Spatial properties of the diffusion coefficient. 
 */
typedef enum {
  Constant, /**< constant diffusion coefficient. */
  Exp, /**< Exponentially growing with vertical scale. \sa zt . */
  Blasi, /**< Suggested by Pasquale. */
  Expr, /**< Expr . */
  Radial,
  Qtau, /**< radial behavior according to SNR distribution. */
  Test,
  ExpRadial /**< Radial behavior according to SNR distribution and exponential growth with vertical scale. */
} DPerpType;


/**
 * @enum DiffusionType
 * @brief Chose whether diffusion should be isotropic (default) or anisotropic, in 3D.
 */
typedef enum {
    Isotropic, /**< Set DM decay. */
    Anisotropic /**< SetDM annihilation. */
} DiffusionType;

/**
 * @enum SpallXSec
 * @brief Choose the spallation cross section model.
 */
typedef enum {
  GalpropXSec, /**< Use Galprop cross sections. */ 
  Webber03, /**< Use Webber 2003 cross sections where available. */
  Fluka,
  DRAGON2
} SpallXSec;

/**
 * @enum LeptonYields
 * @brief Choose model for lepton yields in secondary electron and positron production.
 */
typedef enum {
  GalpropTable, /**< Tabulated Galprop model. */ 
  Kamae, /**< Kamae 2007 yields. \sa include/cparam.h . */
  Pohl,
  FlukaLep
} LeptonYields;

typedef enum {
  GalpropFunction,
  QGSJET,
  FlukaAp,
  Winkler
} AntiprotonYields;

typedef enum {
  Ptuskin94,
  Ptuskin03
} ReaccType;

typedef enum {
  PowerLawBreak,
  ExpSuppressed
} SpectralShape;

// Data files
const std::string ISRFfile = "data/MilkyWay_DR0.5_DZ0.1_DPHI10_RMAX20_ZMAX5_galprop_format.fits.gz"; /**< ISRF model. */
//const std::string ElTablefile = "data/Electron_production_new.dat";  /**< Tabulated electron yields from Galprop. */   //NOT MANTAINED ANYMORE -- use Kamae option
//const std::string PosTablefile = "data/Positron_production_new.dat"; /**< Tabulated positron yields from Galprop. */   //NOT MANTAINED ANYMORE -- use Kamae option
//const std::string FlukaElTablefile = "data/Fluka_Electron_production_new.dat";  /** Fluka model D.G. & N.M. 2014 **/ 
//const std::string FlukaPosTablefile = "data/Fluka_Positron_production_new.dat";  /** Fluka model D.G. & N.M. 2014 **/
//const std::string FlukaNucleiTablefile = "data/Fluka_Nuclei_production.dat";  /** Fluka model D.G. & N.M. 2014 **/
//const std::string FlukaLightNucleiTablefile = "data/Fluka_LightNuclei_production.dat";  /** Fluka model D.G. & N.M. 2014 **/
//const std::string FlukaApTablefile = "data/Fluka_Antiproton_production_new.dat" ;
//const std::string FlukaTertiaryApTablefile = "data/Fluka_TertiaryAntiproton_production_new.dat" ;
//const std::string FlukaProtTablefile = "data/Fluka_Proton_production_new.dat" ;
//
//const std::string Fluka_inelastic_datafile = "data/FlukaIneXsec_new.dat";  /** Fluka Inelastic Xsec D.G. & N.M. 2014 **/
//
//const std::string Webber03Data = "data/webber_xsec.dat"; /**< Tabulated Webber 2003 spallation cross sections. */
//const std::string Webber03DataTotal = "data/webber_xsec_total.dat"; /**< Tabulated Webber 2003 inelastic cross sections. */
//const std::string BNLdata = "data/bnl_data_clean.dat"; /**< Nuclear properties database. */
const std::string BNLdata = "data/isotope_list.dat"; /**< Nuclear properties database. */
//const std::string sourcedata = "data/Source.param";
//const std::string EWdatafile = "data/DMspectrum_EWcontribution.root";
const std::string EWdatafile = "data/DMspectra_EWcontribution_2.0.root"; //NOT MANTAINED ANYMORE

// Fundamental constants
const double C         = 2.99792458e10 ; /**< [cm s^-1] Speed of light in vacuum. */
const double h_planck  = 6.6260755e-27 ; /**< [erg sec] Planck's constant. */

const double year      = 365.25*24.*3600. ;
const double Myr       = 1e6*year ;
const double kpc       = 3.08568e21 ;
const double km        = 1e5 ;

const double Clight    = C * Myr ; /**< [cm Myr^-1] */

//modified
const double mp        = 0.938 ; //272029          ; /**< [GeV] Proton mass. */
const double Pi        = 3.14159265358979312 ;
const double eV_to_erg = 1.60217733e-12 ;
const double erg_to_eV = 1./eV_to_erg ;

const double  M_Sun_in_GeV = 1.115e57; /**< Solar mass in GeV */ //EDIT HANI

//const double  robs      = 8.3 ; /**< [kpc] value of r at Sun's position. */
// Matze: hard coded for now, ToDo as input
//const double  xobs      = 8.3 ; /**< [kpc] value of x at Sun's position. */
//const double  yobs      = 0. ; /**< [kpc] value of y at Sun's position. */
//const double  robs      = sqrt(xobs*xobs + yobs*yobs); /**< calculated [kpc] value of r at Sun's position. */
//const double  zobs      = 0. ; /**< [kpc] value of z at Sun's position. */

// Utils
const double u                 = 1./200. ; /**< [kpc] just a little value to regularize diffusion behavior near GC. */
const double p                 = 3.0 ;     /**< Number of differential operators to be updated. */
const int cross_section_option = 12 ; /**< Option related to Galprop cross section model. */
const double t_half_limit      = 1.e4 ; /**< [yr] Lifetime below which nuclei are considered completely unstable. */
//modified
// const double minlifetime       = 1.e3; // [seconds!!] Lifetime below which a nucleus is not considered
//modified by Nicola March 24, 2015
const double minlifetime       = 1.e-17; // [seconds!!] Lifetime below which a nucleus is not considered

// Galactic Magnetic Field
const double rB        = 10. ; /**< [kpc] MF radial scale length. */
const double zr        = 2.0 ; /**< [kpc] regular MF vertical scale. */
const double Bh        = 6.1 ; /**< [muG] intensity of MF at Sun's position. */
//const double D_ref_rig = 3. ; /**< [GV] reference rigidity for D. */
const double rd        = 2. ; /**< [kpc] radial length scale of D. */
const double r0        = 2. ; /**< [kpc] radial position of maximum of D. */
const double dzzGal    = 0.01 ; /**< [kpc] value to average over the gas distribution, used in \sa GasType::Galprop. */

const double  He_abundance = 0.11 ;/**< Helium abundance in the ISM. */

// Energy losses
const double ALPHAf = 1.0/137.035989561; /**< Fine structure constant. */
const double Mele = 0.5109990615;         /**< [MeV/c^2] electron rest mass. */
const double MeleGeV = 1e-3*Mele;         /**< [GeV/c^2] electron rest mass. */
const double MpGeV = 0.938; 		  /**< [GeV/c^2] proton rest mass in GeV */
const double Rele = 2.8179409238e-13;      /**< [cm] =e^2/mc^2 classical electron radius. */
const double H2PiC   = 0.19732705359e-10; /**< [MeV*cm] =hc/(2Pi) conversion const. */
const double PIR0H2C2  = Pi*Rele*H2PiC*H2PiC;           /**< [MeV^2*cm^3] */
const double PIR02MC2C = Pi*Rele*Rele*C*(Mele*1.e6);    /**< [eV*cm^3/s] = Pi*e^4/mc */
const double PIR02C    = Pi*Rele*Rele*C;                /**< [cm^3/s] */
const double AFR02MC2C = ALPHAf*Rele*Rele*C*(Mele*1.e6); /**< [eV*cm^3/s] = e^6/mc/hc */
 
const double EiH = 13.6e-6;   /**< [MeV] Hydrogen (in electron EL) */
const double EiHe = 24.59e-6;  /**< [MeV] Helium (in electron EL) */
const double TH = 62.8;    /**< [g/cm^2] Hydrogen */
const double THe = 93.1;   /**< [g/cm^2] Helium */

const double ZHe = 2.; /**< He charge. */

// Ionization
const double EH = 19.e-6; /**< [MeV] H  eff. ioniz. potential (in nucl. loss) */
const double EHe= 44.e-6; /**< [MeV] He eff. ioniz. potential (in nucl. loss) */

// Coulomb
const double BK = 1.38066e-23/1.60218e-13;     /**< [MeV/K] Bolzmann constant */
const double Te = 1.0e4;                       /**< [K] temperature of warm ionized medium */
const double bet_e = sqrt(2.*BK*Te/Mele); /**< Some parameter */
const double xm = pow(pow(3.*sqrt(Pi)/4.,1./3.) *bet_e,3); /**< Some other parameter */

// Bremsstrahlung
const double MH = 1.67e-24;                            /**< [grams] the mass of H atom. */
const double MHe= 1.66e-24 *4.;                        /**< [grams] the mass of He atom. */
const double gam1 = 100.; /**< Boosts */
const double gam2 = 800.; /**< Boosts */

const double sigmaT  = 6.6524e-25/(kpc*kpc) ; /**< [kpc^2] Thomson cross section */

const double amu = 0.931494;  // GeV, mass 12C / 12 = atomic mass unit

#endif
