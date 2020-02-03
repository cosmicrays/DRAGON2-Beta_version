/**
 * @file input.h
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief File where class Input is defined.
 */

#ifndef _INPUT_H
#define _INPUT_H

#include <iostream>
#include <string>
#include <vector>

#include "constants.h"

class TiXmlElement;

using namespace std;

/**
 * @class Input
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class to read user's input.
 */
class Input {
   
private:
   int QueryIntAttribute(string, TiXmlElement*);
   int QueryIntAttributeWithDefault(string, TiXmlElement*, int);
   double QueryDoubleAttribute(string, TiXmlElement*);
   double QueryDoubleAttributeWithDefault(string, TiXmlElement*, double);
   string QueryStringAttribute(string, TiXmlElement*);
   vector<double> ReadInjectionStructure(string name, TiXmlElement* el1);
   
public:
   Input() { } /**< The default constructor. */
   int LoadFile(const string);
   void Print();
   
   //  Input(char* run, double& Ekfact_, int& Nr, int& Nz, double& D0_, double& zt_, double& zmax_, double& delta_, double& index_radial_, double& ab_C_, double& ab_N_, double& vA_, double& vC_, double& dvdz_conv_, double etaT_, char* profile);
   /**< Constructor given individual data. This is used when doing MCMC. */
   
   ~Input() { }
   /**< Destructor. Does nothing. */
   
   // Grid
   double Ekfact; /**< Delta E/E. */
   double Ekmax; /**< [Gev/n] Maximal energy. */
   double Ekmin; /**< [GeV/n] Minimal energy. */
   int numr; /**< Dimension of R grid. */
   int numx; /**< Dimension of X grid. */
   int numy; /**< Dimension of Y grid. */
   bool dimx_equidistant, dimy_equidistant, dimz_equidistant, dimr_equidistant; //mw
   string divx, divy, divz, divr;                                               //mw
   int numz; /**< Dimension of Z grid */
   double xobs; /**< [kpc] value of x at Sun's position. */
   double yobs; /**< [kpc] value of y at Sun's position. */
   double robs; /**< calculated [kpc] value of r at Sun's position. */
   double zobs; /**< [kpc] value of z at Sun's position. */
   double Rmin; /**< [kpc] */
   double Rmax; /**< [kpc] Radial end of the Galaxy. */
   double zmax; /**< [kpc] Half height of the propagation box. */
   string gridtype; /**< 2D or 3D */
   
   // Limits for the nuclear chain
   int Zmax; /**< Maximal propagated charge. */
   int Zmin ; /**< Minimal propagated charge. */
   bool prop_nuclei;
   bool prop_ap; /**< Whether to propagate antiprotons. */
   bool prop_deuteron; /**< Whether to propagate antideuterons. */
   bool prop_lep; /**< Whether to propagate electrons and positrons. */
   bool prop_secap; /**< Whether to propagate antiprotons. */
   bool prop_secdeuteron; /**< Whether to propagate antideuterons. */
   //bool prop_secel; /**< Whether to propagate electrons and positrons. */
   bool prop_TPP; /**< Whether to propagate positrons from triple pair production */
   bool prop_DMap; /**< Whether to propagate antiprotons. */
   bool prop_DMdeuteron; /**< Whether to propagate antideuterons. */
   bool prop_DMel; /**< Whether to propagate electrons and positrons. */
   bool prop_extracomp;  /**< @brief Whether to propagate the extra component (electron/positron).  */
   bool Kcapture;
   bool DoubleRun;
   
   // MW: local bubble and spiral arms for 3D nonequi usage
   double LB_ax, LB_ay, LB_az; //dimensions of Local Bubble
   string LB_shape, LB_smearing;
   double LB_source, LB_DMsource, LB_gastotal, LB_gasH2, LB_gasHI, LB_gasHII, LB_diff, LB_delta, LB_convec, LB_vA, LB_MagField, LB_ISRF;
   string SA_type;
   double SA_diff, SA_convec, SA_vA, SA_source, SA_DMsource, SA_gasH2, SA_gasHI, SA_gasHII, SA_MagField, SA_ISRFStar, SA_ISRFDust;
   double SA_cut_diff, SA_cut_convec, SA_cut_vA, SA_cut_source, SA_cut_DMsource, SA_cut_gasH2, SA_cut_gasHI, SA_cut_gasHII, SA_cut_MagField, SA_cut_ISRF;
   
   int num_arms;
   double spiral_width;
   vector<double> arms_Kvec;
   vector<double> arms_r0vec;
   vector<double> arms_theta0vec;
   
   // I/O
   bool fullstore; /**< Set whether to store full spatial density of propagated CRs. */
   bool partialstore; /**< Set whether to store only CR spectra at Solar system position. */
   string filename; /**< Run name. */
   bool timestepstore; //MW130801 write each Evolutor dt loop to file
   int feedback;

   // Algorithm
   int Nrept; /**< Number of iterations for a single timestep value \sa dtfactor \sa dtmin \sa dtmat. Only used in class \sa TCREvolutor .*/
   double dtfactor; /**< Factor of reduction of timestep after \sa Nrept iterations. */
   double dtmin; /**< [Myr] Minimum value of timestep. */
   double dtmax; /**< [Myr] Maximum (starting) value of timestep. */
   
   bool ADI; /**< Controls whether to refine the solution with the ADI scheme. */
   bool OpSplit; /**< Controls whether to use the OS (Galprop) scheme. */
   double tol; /**< Tolerance for convergence of the solution. Used only when ADI == true */
   double alpha; /**< Courant parameter for the ADI scheme */
   
   bool TESTMODE;
   bool MOVING;
   
   double source_x0;
   double source_y0;
   double source_z0;
   double source_vx;
   double source_vy;
   double source_vz;
   
   // Galaxy
   GasType gas_model; /**< Set Gas model. */
   ISRFType ISRF_model;
   SNRType SNR_model; /**< Set Source distribution model. */
   SNRType SNR_model_Extra; /**< Set Source distribution model. */
   DiffusionType DiffT;
   xco_modes xco_mode;
   double xco_constant;
   double xco_inner;
   double xco_outer;
   double ringmin; /**< [kpc] Inner radius of the ring if option \sa SNRType::Rings is selected. */
   double ringmax; /**< [kpc] Outer radius of the ring if option \sa SNRType::Rings is selected. */
   double ringmin_extra; /**< [kpc] Inner radius of the ring if option \sa SNRType::Rings is selected. */
   double ringmax_extra; /**< [kpc] Outer radius of the ring if option \sa SNRType::Rings is selected. */
   double rings_period	;
   double rings_phase  ;
   double pointsrc_x;
   double pointsrc_y;
   double pointsrc_z;
   double pointsrc_Extra_x;
   double pointsrc_Extra_y;
   double pointsrc_Extra_z;
   
   // Reacceleration
   double vAlfven; /**< [km/s] Alfven velocity. */
   ReaccType diff_reacc; /**< Use diffusion reacceleration model: 1 = Seo & Ptuskin 94; 2 = Ptuskin 2003. */
   bool REACC;
   
   // Convection
   double v0; /**< [km/s] Convection velocity at z_k (base velocity). */
   double vb; /**< [km/s] Convection velocity at z=0. */
   double f_b; /**< [0,1] vb=v0*fb. */
   double dvdz; /**< [km/s/kpc] Vertical derivative of convection velocity. */
   double z_k;  /**< [kpc] height above which we have linear increase in convection. */
   bool CONVECTION;
   double conv_index_radial;
   double conv_threshold;
   DPerpType set_profile_conv;
      
   // Diffusion
   DPerpType set_profile; /**< Type of diffusion. */
   double D0; /**< [1e28 cm^2/s] Normalization of diffusion coefficient. */
   double D_ref_rig;
   double zt; /**< [kpc] Vertical scale of diffusion coefficient. */
   double etaT; /**< index of dependence of the diffusion coefficient with \beta. */
   double delta; /**< Power law index of energy dependence of diffusion coefficient. */
   double delta_h; /**< Power law index of energy dependence of diffusion coefficient. */
   int VariableDelta; /**< if this flag is activated, delta is a function of R */
   double diff_threshold;
   double delta_A;
   double delta_B;
   double rho_b; /**< Power law index of energy dependence of diffusion coefficient. */
   double index_radial; /**< Power law index of the scaling of diffusion coefficient with sources in the radial direction. */
   bool DIFFUSION; //fk 130701
   
   //Daniele's block for anisotropic diffusion
   double Dpar;
   double etaTpar;
   double DeltaPar;
   double etaTperp;
   double DeltaPerp;
   double Dperp;
   
   //Energy loss & radioactive decay fk 130701
   bool ELOSS;
   bool RDECAY;
   
   // Cross sections
   SpallXSec spallationxsec; /**< Set spallation cross section model. */
   LeptonYields ly; /**< Set lepton yields. */
   AntiprotonYields apy;
   bool SPALL; //fk 130701
   
   // Magnetic Field
   MagFieldModel BM;
   double B0disk;
   double B0halo;
   double B0turb;
   double betaFarrar;
   double b0;
   double bx,by,bz,bturb;
   
   // CR normalization
   double sp_ref_rig_norm; /**< [GeV] energy for normalization of the final hadron flux. */
   double sp_ref_rig_el; /**< [GeV] energy for normalization of the final electron flux. */
   double sp_ref_rig_el_extra; /**< [GeV] energy for normalization of the final electron flux. */
   double spect_norm; /**< [m^-2 sr^-1 s^-1 GeV^-1] normalization at \sa sp_ref_rig_norm  from Biermann A&A 330, 389 (1998). */
   double spect_norm_el; /**< [m^-2 sr^-1 s^-1 GeV^-1] normalization at \sa sp_ref_rig_el  from Biermann A&A 330, 389 (1998). */
   double spect_norm_el_extra;
   string sourcedata;
   SpectralShape Spectrum;
   double sp_ref_rig_exp_el;
   double sp_ref_rig_exp;
   double exp_cut_index_el;
   
   // not used now!
   double sp_ref_rig; /**< [GV] energy for abundance normalization or breaks in the sources. */
   
   //DMCMC
   string run_id;
   int write_flag;
   int fit_mod;
   
   bool   UseInjectionIndexAllNuclei;
   //OLD implementation of injection slopes for the nuclei -----------
   /*
    double rho_0; //position of first  break in GV
    double rho_1; //position of second break in GV
    double rho_2; //position of third  break in GV
    double alpha_0; // lowest energy slope
    double alpha_1;
    double alpha_2;
    double alpha_3; // highest energy slope*/
   //-------------------------------------------
   
   //NEW implementation of injection slopes for the nuclei with arbitrary number of breaks ---
   vector<double> inp_inj_indexes;
   vector<double> inp_break_positions;
   double cutoff_rig;
   //-----------------------------------------------------------------------------------------
   
   //NEW implementation of injection slopes for the electrons with arbitrary number of breaks ---
   vector<double> inp_inj_el_indexes;
   vector<double> inp_break_el_positions;
   double cutoff_rig_el; /**< [GV] Rigidity of cut off of primary electron injection spectra. */
   //-----------------------------------------------
   //NEW implementation of injection slopes for the extra component with arbitrary number of breaks ---
   vector<double> inp_inj_extra_indexes;
   vector<double> inp_break_extra_positions;
   double cutoff_rig_extra; /**< [GV] Rigidity of cut off of primary electron injection spectra. */
   
   // not used now!
   //double sp_ref_rig; /**< [GV] energy for abundance normalization or breaks in the sources. */
   // for Leptons
   // not used now! -----------------------------------------------------------------------------------------------------
   //double cutoff_rig; /**< [GV] Rigidity of cut off of primary electron injection spectra. */
   //double spect_ind_el_low; /**< Source injection spectral index for primary electrons below \sa sp_ref_rig_break_el. */
   //double spect_ind_el; /**< Source injection spectral index for primary electrons above \sa sp_ref_rig_break_el. */
   //double sp_ref_rig_break_el;  /**< [GV] energy for break in electron injection spectrum. */
   // -------------------------------------------------------------------------------------------------------------------
   
   //double cutoff_rig_extra; /**< [GV] Rigidity of cut off of primary electron injection spectra. */
   //double spect_ind_el_low_extra; /**< Source injection spectral index for primary electrons below \sa sp_ref_rig_break_el. */
   //double spect_ind_el_extra; /**< Source injection spectral index for primary electrons above \sa sp_ref_rig_break_el. */
   //double sp_ref_rig_break_el_extra;  /**< [GV] energy for break in electron injection spectrum. */
   
   // Antiprotons
   int antiproton_cs; /**< antiproton cross section option. */
   bool scaling; /**< use the propagated He spectrum to compute He contribution to antiproton spectrum, or use a scaling relation from the proton spectrum. */
   
   // DM
   bool propDM;
   bool MOVING_CLUMP;
   bool analytical_refinement;
   
   double clump_x0;
   double clump_y0;
   double clump_z0;
   double clump_vx;
   double clump_vy;
   double clump_vz;
   double clump_norm;
   double clump_size;
   double clump_inj;
   double clump_cutoff;
   double clump_deltat;
   int stop_after_timestep;
   
   double mx; /**< DM mass. */
   double taudec; /**< Lifetime of DM particle. */
   double sigmav; /**< Thermally averaged DM annihilation cross section. */
   DMreaction DMr; /**< Set DM reaction channel. */
   DMspectrum DMs; /**< Set DM reaction model. */
   double rhos;    /**< DM solar system density   GeV/cm^3 */
   double EkDelta; /**< Ek in case of Delta injection for DM */
   string MySelfTableDMap;
   string MySelfTableDMdbar;
   string MySelfTableDMel;
   
   int dmmode; /**< DM channel.
                Ch No  Particles       Old Ch No
                5     h10 h30         7
                6     h20 h30         11
                8     z0 h10          8
                9     z0 h20          9
                11    w+ h- / w- h+   10
                12    z0 z0           6
                13    w+ w-           5
                17    mu+ mu-         13
                19    tau+ tau-       4
                22    cc-bar          1
                24    tt-bar          3
                25    bb-bar          2
                26    gluon gluon     12
                29    z gamma         14
                */
   
   DMprofile dmprof; /**
                      @var DMprofile dmprof
                      @brief DM spatial profile.
                      Iso      NFW     Kra    Moore    Einasto
                      Alpha                   2        1       2      1.5
                      Beta                    2        3       3      3.0
                      Gamma                   0        1      0.4     1.5
                      Rs                     3.5      20.     10.     28.
                      DMprof_                 0        1       2       3         4
                      */
};

#endif
