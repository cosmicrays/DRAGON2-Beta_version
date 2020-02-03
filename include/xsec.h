#ifndef _XSEC_H
#define _XSEC_H


#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <vector>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>
#include "constants.h"
#include "input.h"
#include "grid.h"
#include "utilities.h"
#include "nucleilist.h"

class TGrid;
class TNucleiList;
class Input;


extern "C" void yieldx_(int*,int*,int*,int*,float*,float*);
/**
 * @fn extern "C" void yieldx_(int*,int*,int*,int*,float*,float*)
 * @brief Silberberg & Tsao isotopic production cross section
 */
extern "C" double wsigma_(int*,int*,int*,int*,double*);
/**
 * @fn extern "C" double wsigma_(int*,int*,int*,int*,double*)
 * @brief Wrapper to Webber's (1993) code, from Galprop
 */
extern "C" void set_sigma_(int*); /**< Initialization of Webber's code. */
extern "C" void sigtap2_(int*); /**< initialization of the Barashenkov & Polanski cross section code. */
extern "C" double sighad_(int*,double*,double*,double*,double*,double*); /**< Barashenkov & Polanski pA total cross section. */


//namespace DRAGON {
  
  class TXSecBase {
    
  public:
    TXSecBase() { }
    virtual ~TXSecBase() { }
    
    /* virtual std::vector<double> GetXSec(int /\**< Charge of primary. *\/, int /\**< Mass of primary. *\/, */
    /* 				      int /\**< Charge of secondary. *\/, int /\**< Mass of secondary. *\/)=0; */
    
    //virtual std::vector<double> GetXSec(int /**< uid of nucleus. */)=0;
    
    /* virtual void InitTableInelastic() { */
    /*   throw -1; */
    /* } */
    
    virtual std::vector<double> GetXSec(int, int)=0;
    
    virtual void nucleon_cs(int, double /**< kinetic energy per nucleon of the beam momentum particles (GeV) */,
    			  int /**< =+/-1 is the charge of beam and At is the atomic number of target nuclei */,
    			  int /**< =1 for proton, = -1 for antiproton collisions. */, int /**< =1 for proton collisions. */,
    			  double* /**< [mbarn] is the proton-proton inelastic cross sect. */,
    			  double* /**< [mbarn] is the proton-nucleus inelastic cross sect. (p+A2) */,
    			  double* /**< [mbarn] is the antiproton-proton inelastic cross sect.(nonannihil.) */,
    			  double* /**< [mbarn] put =PA_inel, is the antiproton-nucleus inelastic cross sect.(nonan.) */,
    			  double* /**< [mbarn] is the antiproton-proton annihilation cross sect. */,
    			  double* /**< [mbarn] is the antiproton-nucleus annihilation cross sect. */)=0; // { /*...*/ }
    
    virtual void Kcapture_cs(double, int, int, double*, double*)=0; // { /* */ } /* Cross section of electron attachment/stripping process */
    
    virtual double GetXSec(std::pair<int, int> /**< Pair of parent-daughter nuclei uid. */, double /**< Energy. */)=0; // { /*...*/ }
    /* virtual double GetTotalXSec(int /\**< Nucleus uid. *\/, double /\**< Energy. *\/)=0; // { /\*...*\/ } */
    /* /\* virtual double GetTotalApHXSec(double)=0; *\/ */
    /* /\* virtual double GetTotalApHeXSec(double)=0; *\/ */
    /* /\* virtual std::vector<double> GetTotalXSec_vec(int /\\**< Nucleus uid. *\\/)=0; *\/ */
    /* /\* virtual std::map<std::pair<int, int>, std::vector<double> > GetTotalXSec_for_each_gas_type()=0; *\/ */
    /* /\* virtual bool IsPresent(int uid /\\**< Nucleus uid. *\\/)=0; *\/ */
    /* /\* virtual bool IsPresent(std::pair<int, int> uid /\\**< Nucleus parent-daughter uid. *\\/)=0; *\/ */
    
    
    
    inline double GetHefactor() { 
      return 1.26 / He_abundance; 
    } 
    
    
    virtual double antiproton_cc1(gsl_integration_workspace*, size_t, int, double, double, int, int, int, int)=0; // { /\*...*\/ } */
    
  };
  
  class TSpallationNetwork {
    
  public:
    //TSpallationNetwork(double*, std::string*);
    TSpallationNetwork(TGrid*  /**< Kinematics. */, Input*, std::vector<TXSecBase*>  /**< Model of the cross section */,
    		      std::vector<int>&  /**< List of propagating nuclei. */);
    
    ~TSpallationNetwork() {
            
      spall.clear();
      spall_apel.clear();
      DRAGONEnergyVector.clear();
      spallationXsections.clear();
      spallationXsectionsInterpolated.clear();
      
    }
    
    std::vector<std::pair<double, double> > GetXSecPair(int /**< Parent nucleus uid */, int /**< Daughter nucleus uid. */);
    
    virtual std::vector<double> GetXSec(int /**< Parent nucleus uid */, int /**< Daughter nucleus uid. */);
    
    std::pair<double, double> GetXSecPair(int /**< Parent nucleus uid */, int /**< Daughter nucleus uid. */, double /** Energy. */);
    
    std::pair<double, double> GetXSecPair(int /**< Parent nucleus uid */, int /**< Daughter nucleus uid. */, int /** Energy index. */);

    virtual std::vector<double> GetXSecTPP(std::vector<double> /** ISRF frequency vector*/);

    //void spec_ini();
    
    //double spec_int(double ep, double es, int id, int reac);
    
    //inline std::vector<double> GetXSec(int /**< Parent nucleus uid */, int /**< Daughter nucleus uid. */) {return std::vector<double>(0);}
    //inline double GetXSec(int /**< Parent nucleus uid */, int /**< Daughter nucleus uid. */, double /** Energy. */) {return 0.;}
    //inline double GetXSec(int /**< Parent nucleus uid */, int /**< Daughter nucleus uid. */, int /** Energy index. */) {return 0.;}
    
    virtual std::vector<double> GetXSecApEl(int /**< Parent nucleus uid */, int /**< Daughter particle uid. */, int /** Energy index of parent particle. */);
    void InitXSecWinkler(double);

    
  private:
    
   std::string spallXsecFile, spallXsecFile_pos, spallXsecFile_el, spallXsecFile_p, spallXsecFile_ap, spallXsecFile_Tap;
   
   std::ifstream datafileSpall; 
   
   std::map< std::pair<int, int>, std::vector< std::pair<double, double> > > spallationXsections;	
   std::map< std::pair<int, int>, std::vector< std::pair<double, double> > > spallationXsectionsInterpolated;
   std::map< std::pair<int, int>, std::vector< double> > spall;	
   //std::map< std::pair<int, int>, std::vector< double> > spall_ApEl;	

   std::map< std::pair<int, int>, std::vector< std::vector<double> > > spall_apel;
   //std::map< std::pair<int, int>, std::vector< std::vector<double> > > spallationXsectionsInterpolated_apel;
   
   //std::map<std::pair<int, int>, std::vector<std::pair<double, double> > > spallationXsections_fullIntegral;	
   //std::map<std::pair<int, int>, std::vector<std::pair<double, double> > > spallationXsectionsInterpolated_fullIntegral;
   
   std::vector<double> spallationKineticEnergyVector; /**< Energy array. */
   
   std::vector<double> DRAGONEnergyVector; /**< Energy array. */
   //int index_matrix(int i /**< Energy index of daughter particle. */, int j /**< Energy index of parent particle. */) { return i*Nprotons+j; }
   /**
    * @fn int index_matrix(int i, int j)
    * @brief Convert matrix to linearized index, for data tables. \sa InitDataTables.
    * @return i*801+j
    */

   int dimEn;
   int Nelectrons, Nprotons;
   void InitXSecKamae(double);
   //void InitXSecGalprop(double);
   void InitXSecPohl(double);
   void InitDataTablesPohl();   

   void Retrieve_XSecs(Input*, std::vector<int>&, TGrid* );

   std::vector<double> ProdXsec;
   std::vector< std::vector<double> > ElppPohl;
   std::vector< std::vector<double> > ElHepPohl;
   std::vector< std::vector<double> > PosppPohl;
   std::vector< std::vector<double> > PosHepPohl;
   
   std::vector<double> Matrix_El_pp;    /**< Table for the electron spectrum from pp interactions. */
   std::vector<double> Matrix_El_pHe;   /**< Table for the electron spectrum from pHe interactions. */
   std::vector<double> Matrix_El_Hep;   /**< Table for the electron spectrum from Hep interactions. */
   std::vector<double> Matrix_El_HeHe;  /**< Table for the electron spectrum from HeHe interactions. */
   std::vector<double> Matrix_Pos_pp;   /**< Table for the positron spectrum from pp interactions. */
   std::vector<double> Matrix_Pos_pHe;  /**< Table for the positron spectrum from pHe interactions. */
   std::vector<double> Matrix_Pos_Hep;  /**< Table for the positron spectrum from Hep interactions. */
   std::vector<double> Matrix_Pos_HeHe; /**< Table for the positron spectrum from HeHe interactions. */
   std::vector<double> Matrix_Ap_pp;    /**< Table for the ap spectrum from pp interactions. */
   std::vector<double> Matrix_Ap_pHe;   /**< Table for the ap spectrum from pHe interactions. */
   std::vector<double> Matrix_Ap_Hep;   /**< Table for the ap spectrum from Hep interactions. */
   std::vector<double> Matrix_Ap_HeHe;  /**< Table for the ap spectrum from HeHe interactions. */
   std::vector<double> Matrix_Ap_apHe;
   std::vector<double> Matrix_Ap_app;
   std::vector<double> Matrix_p_pp;    /**< Table for the p spectrum from pp interactions. */
   std::vector<double> Matrix_p_pHe;   /**< Table for the p spectrum from pHe interactions. */
   std::vector<double> Matrix_p_Hep;   /**< Table for the p spectrum from Hep interactions. */
   std::vector<double> Matrix_p_HeHe;  /**< Table for the p spectrum from HeHe interactions. */
   std::vector<double> Matrix_3Ap_app;    /**< Table for the p spectrum from app interactions. */
   std::vector<double> Matrix_3Ap_apHe;   /**< Table for the p spectrum from apHe interactions. */

   double cs_pp;
   double cs_pHe;
   double cs_Hep;
   double cs_HeHe;
   
   double PP_inel;
   double PA_inel;
   double aPP_non;
   double aPA_non;
   double aPP_ann;
   double aPA_ann;
   unsigned int ind;
   int i_ap, j_pr, index, index1;
   double u, valueup, valuefix, Eap_i, Epr, Ep_i;
 };
 
 class TInelasticCrossSection {
   
 public:
   TInelasticCrossSection() {} /**< Default constructor. */
   
   TInelasticCrossSection(TGrid*, Input*, int /**< Nucleus uid. */, int /**< K_electron */, std::vector<TXSecBase*>);
   /*TInelasticCrossSection(double* , std::string* , int  , int , std::vector<TXSecBase*>){}*/
   
   virtual ~TInelasticCrossSection() {
     totalxsec.clear();
     //xsec.clear();
     //std::cout << "Destroyed entity" << std::endl;
     //xsec_for_each_gas_type.clear();
   } /**< Destructor. */
   
   inline std::vector<double>& GetXSec() {
     return totalxsec; } /**< Obtain the cross section.*/ 
   
   inline std::map<std::pair<int, int>, std::vector<double> >& GetXSec_extended() {
   std::cout << "returning extended inelastic cross section" << std::endl;
   return xsec_for_each_gas_type;
   } /**only for Fluka model */
   
   inline double GetXSec(int i /**< Energy index. */) {
     return totalxsec[i];
   } /**< Obtain cross section at given energy. */
   
   
   double GetTotalIne (std::vector<double> , std::vector<double> , double);
   double GetIneXSec (double );
   
   
 protected:

   std::vector<double> e_vec; /*< Energy array. */ 
   
   //std::vector<double> xsec; /**< Array of Ine cross section. */
   
   std::map<std::pair<int, int>, std::vector<double> > xsec_for_each_gas_type; /**< Array of cross section for each nucleus against each gas type, only for Fluka model. */
   
   std::map<int, std::vector<double> > mapInexsec;
   std::vector<double> totalxsec;
   std::ifstream dataFileIne; 
   
   
 };
 
 
 
 /* class TWebber03: public TXSecBase {  */
 
 /* public:  */
 
 /*   TWebber03();  */
 
 /*   ~TWebber03();  */
   
 /*   double GetXSec(std::pair<int, int> /\*< Pair of parent-daughter nuclei uid. *\/, double /\*< Energy. *\/);  */
   
 /*   double GetTotalXSec(int /\*< Nucleus uid. *\/, double /\*< Energy. *\/);  */
   
 /*   /\* 	virtual std::vector<double> GetTotalXSec_vec(int /\\**< Nucleus uid. *\\/) { *\/ */
 /*   /\* 		throw -1; *\/ */
 /*   /\* 	} *\/ */
   
 /*   /\* 	virtual void InitTableInelastic() { *\/ */
 /*   /\* 		throw -1; *\/ */
 /*   /\* 	} *\/ */
   
 /*   /\* 	virtual double GetTotalApHXSec(double) { *\/ */
 /*   /\* 		throw -1; *\/ */
 /*   /\* 	} *\/ */
   
 /*   /\* 	virtual double GetTotalApHeXSec(double) { *\/ */
 /*   /\* 		throw -1; *\/ */
 /*   /\* 	} *\/ */
   
 /*   /\* 	virtual std::map<std::pair<int, int>, std::vector<double> > GetTotalXSec_for_each_gas_type() { *\/ */
 /*   /\* 		throw -1; *\/ */
 /*   /\* 	} *\/ */
   
 /*   //virtual std::vector<double> GetXSec(int /\*< uid of nucleus. *\/) {  */
 /*   //throw -2;  */
 /*   //}  */
   
 /*   /\* 	virtual std::vector<double> GetXSec(int, int) { *\/ */
 /*   /\* 		throw -2; *\/ */
 /*   /\* 	} *\/ */
   
 /*   /\* 	virtual std::vector<double> GetXSec(int /\\**< Charge of primary. *\\/, int /\\**< Mass of primary. *\\/, *\/ */
 /*   /\* 			int /\\**< Charge of secondary. *\\/, int /\\**< Mass of secondary. *\\/) { *\/ */
 /*   /\* 		throw -2; *\/ */
 /*   /\* 	} *\/ */
   
	
 /*   /\* 	virtual void Kcapture_cs(double, int, int, double*, double*) { *\/ */
 /*   /\* 		throw -2; *\/ */
 /*   /\* 	} *\/ */
   
 /*   //virtual bool IsPresent(int uid /\*< Nucleus uid. *\/) {  */
 /*   //std::map<int, std::vector<double> >::iterator it = totalxsec.find(uid);  */
 /*   //		if (it != totalxsec.end())  */
 /*   //		return true;  */
 /*   //		else  */
 /*   //		return false;  */
 /*   //}  */
   
 /*   /\* 	virtual bool IsPresent(std::pair<int, int> uid /\\**< Nucleus parent-daughter uid. *\\/) { *\/ */
 /*   /\* 		std::map<std::pair<int, int>, std::vector<double> >::iterator it = xsec.find(uid); *\/ */
 /*   /\* 		if (it != xsec.end()) *\/ */
 /*   /\* 			return true; *\/ */
 /*   /\* 		else *\/ */
 /*   /\* 			return false; *\/ */
 /*   /\* 	} *\/ */
   
 /*   /\* 	//void Print(TGrid*); *\/ */
   
 /*   /\* 	virtual double antiproton_cc1(gsl_integration_workspace*, size_t, int, double, double, int, int, int, int) { *\/ */
 /*   /\* 		return -1; *\/ */
 /*   /\* 	} *\/ */
   
 /* protected: */
 /*   /\* 	std::map<std::pair<int, int>, std::vector<double> > xsec; *\/ */
   
 /*   std::vector<double> totalxsec;  */
   
 /*   /\* 	std::vector<double> energy; /\\**< Energy array. *\\/ *\/ */
   
 /*   /\* 	void convert(const int channel /\\**< Database ID. *\\/, int& uid1 /\\**< uid of parent nucleus (it will be modified). *\\/, *\/ */
 /*   /\* 			int& uid2 /\\**< uid of daughter nucleus (it will be modified). *\\/) { *\/ */
 /*   /\* 		int A1 = channel / 1000000; *\/ */
 /*   /\* 		int Z2 = channel % 100; *\/ */
 /*   /\* 		int Z1 = (channel - A1 * 1000000) / 10000; *\/ */
 /*   /\* 		uid1 = A1 + Z1 * 1000; *\/ */
 /*   /\* 		int A2 = (channel - A1 * 1000000 - Z1 * 10000) / 100; *\/ */
 /*   /\* 		uid2 = A2 + Z2 * 1000; *\/ */
 /*   /\* 		return; *\/ */
 /*   /\* 	} *\/ */
   
 /*   /\* 	void convert_single(const int channel /\\**< Database ID. *\\/, *\/ */
 /*   /\* 			int& uid1 /\\**< uid of nucleus (it will be modified). *\\/) { *\/ */
 /*   /\* 		int A1 = channel / 100; *\/ */
 /*   /\* 		int Z1 = channel % 100; *\/ */
 /*   /\* 		uid1 = A1 + Z1 * 1000; *\/ */
 /*   /\* 		return; *\/ */
 /*   /\* 	} *\/ */
 /* }; */
 
 class TGalpropXSec : public TXSecBase { 
   
 public: 
   TGalpropXSec(TGrid* co, Input* in){
     coord = co;
     inn = in;
     energy = co->GetEk();
     beta = co->GetBeta();
     dimEn = energy.size();
}
   
   
   ~TGalpropXSec(){
     energy.clear();
     beta.clear();
   }
   
   virtual std::vector<double> GetXSec(int /*< Charge of primary. */, int /*< Mass of primary. */, 
				       int /*< Charge of secondary. */, int /*< Mass of secondary. */); 
   
   //virtual std::vector<double> GetXSec(int , int , int , int , std::vector<std::pair<double, double> >, std::vector<double>);
   
   /*virtual std::vector<double> GetTotalXSec_vec(int /*< Nucleus uid. *\/) { */
   /*  throw -1; */
   /*}*/ 
   
   /* 	virtual void InitTableInelastic() { */
   /* 		throw -1; */
   /* 	} */
   
   /* 	virtual bool IsPresent(int uid /\**< Nucleus uid. *\/) { */
   /* 		return true; */
   /* 	} */
   
   /* 	virtual bool IsPresent(std::pair<int, int> uid /\**< Nucleus parent-daughter uid. *\/) { */
   /* 		return true; */
   /* 	} */
   
   /* 	//virtual void Print(TGrid*) { */
   /* 	//} */
   
   /* 	virtual std::vector<double> GetXSec(int /\**< uid of nucleus. *\/) { */
   /* 		throw -1; */
   /* 	} */
   
   /* 	virtual double GetTotalApHXSec(double) { */
   /* 		throw -1; */
   /* 	} */
   
   /* 	virtual double GetTotalApHeXSec(double) { */
   /* 		throw -1; */
   /* 	} */
   
   virtual double GetXSec(std::pair<int, int> /*< Pair of parent-daughter nuclei uid. */, double /*< Energy. */){ throw -1; }
   
   virtual void InitXSec();
   
   /* 	virtual double GetTotalXSec(int /\**< Nucleus uid. *\/, double /\**< Energy. *\/) { */
   /* 		throw -1; */
   /* 	} */
   
   virtual std::vector<double> GetXSec(int /*< uid of nucleus. */, int /*< K electron */); 
   
   /* 	virtual std::map<std::pair<int, int>, std::vector<double> > GetTotalXSec_for_each_gas_type() { */
   /* 		throw -1; */
   /* 	} */
   
   /* 	~TGalpropXSec() { */
   /* 	} */
   
   virtual void nucleon_cs(int, double /**< kinetic energy per nucleon of the beam momentum particles (GeV) */, int /**< =+/-1 is the charge of beam and At is the atomic number of target nuclei */, int /**< =1 for proton, = -1 for antiproton collisions. */, int /**< =1 for proton collisions. */, double* /**< [mbarn] is the proton-proton inelastic cross sect. */, double* /**< [mbarn] is the proton-nucleus inelastic cross sect. (p+A2) */, double* /**< [mbarn] is the antiproton-proton inelastic cross sect.(nonannihil.) */, double* /**< [mbarn] put =PA_inel, is the antiproton-nucleus inelastic cross sect.(nonan.) */, double* /**< [mbarn] is the antiproton-proton annihilation cross sect. */, double* /**< [mbarn] is the antiproton-nucleus annihilation cross sect. */); 
   
   virtual double antiproton_cc1(gsl_integration_workspace*, size_t, int, double, double, int, int, int, int); 

   virtual void Kcapture_cs(double Ek, int Zp, int Zt, double *attach, double *strip);
   
   void set_sigma_cc(); /**< initialization of Webber's code */
   double wsigma_cc(int,int,int,int,double); /**< Webber's isotopic production cross section */
   double yieldx_cc(int,int,int,int,float); /**< Silberberg & Tsao isotopic production cross section */
   void sigtap_cc(int); /**< initialization of the Barashenkov & Polanski cross section code */
   double sighad_cc(int,double,double,double,double,double); /**< Barashenkov & Polanski pA total cross section */

   double He_to_H_CS_ratio(double /*< energy of the primary (GeV/nucleon). */, int /*< primary charge. */, 
   			   int /*< primary mass. */, int /*< secondary charge. */, int /*< secondary mass. */); 
 protected:

   int dimEn;
   std::vector<double> energy;
   std::vector<double> beta;
   TGrid* coord;
   Input* inn;

   std::vector<std::string> data_filename; /**< data files to read. */
   int n_data[3][3]; /**< their dimensions n1,n2,n3. times number of files */
   float *nucdat, *csdat, *protdat;
   
   
   double nucdata(int, int, int, int, int, int, int*, int*, double*); 
   double isotope_cs(double, int, int, int, int, int, int*); 
   double eval_cs(double, int, int, int*);
   
   inline double FI(double X, double X1, double X2, double F1, double F2) { 
    return ((F1 - F2) * X + X1 * F2 - X2 * F1) / (X1 - X2); 
   }
   
   inline double He_to_H_CS_totratio(int IAI /* Mass of primary */) { 
    return 2.10 / pow((double) IAI, 0.055); 
   }
   //void set_sigma_cc(); /**< initialization of Webber's code */
   //double wsigma_cc(int,int,int,int,double); /**< Webber's isotopic production cross section */
   //double yieldx_cc(int,int,int,int,float); /**< Silberberg & Tsao isotopic production cross section */
   //void sigtap_cc(int); /**< initialization of the Barashenkov & Polanski cross section code */
   //double sighad_cc(int,double,double,double,double,double); /**< Barashenkov & Polanski pA total cross section */
   inline int fnuc(int z,int a) { return (100 * (z) + (a)); }
   inline int inuc(float b) { return (int)(100 * (b) + 0.1); }
 }; 



double tan_ng(double, void*);
 
   //} /* namespace DRAGON */

#endif
