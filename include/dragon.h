/**
 * @file dragon.h
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 */

#include <iostream>
#include <vector>
#include <string>

#include "fitsio.h"

class Input;
class Galaxy;
class TNucleiList;
class TParticle;
class TCREvolutorBasis;
class TGrid;
class TXSecBase;
class TSpectrum;
class TSource;

using namespace std;

/**
 * @class DRAGON
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief %DRAGON main class.
 *
 * This class takes care of reading user input and settings, 
 * creating the necessary galactic structure, the algorithm to solve the diffusion equation, 
 * the structure to contain the propagated nuclear densities and the output.
 *
 * \sa Input \sa constants.h \sa TNucleiList \sa Galaxy \sa TCREvolutorBasis \sa TParticle 
 */

class DRAGON {
    
public:
    DRAGON() { } /**< The default constructor */
    
    DRAGON(Input* in_ /**< Wrapper to user input. */);
    /**< Constructor given some input. */ 
    
    ~DRAGON();
    /**< Destructor, that takes care also of closing the output files. */
    
    void Run(); /**< Main routine of DRAGON, that actually propagates the nuclear chain. */
    void Print(); /**< Print the first HDU, with general information on the run. */
    // void CalChi2();

    //  void PrintChi2(); /**< Compute and print chi2, for MCMC runs */
    /**
     * @fn void PrintRatios(double modpotential = 0.0)
     * @brief Print B/C, N/O, C/O and 10Be/9Be, for MCMC runs.
     * @param modpotential [GV] Modulation potential
     */
    //void PrintRatios(double modpotential = 0.0);
    //void PrintRatiosSubFe(double modpotential = 0.0);
    //void PrintAnisotropyElectrons(double modpotential=0.0);
    //void PrintAnisotropyProtons(double modpotential=0.0);
    
    /**
     * @fn vector<double> particle_modulated_spectrum(TParticle* particle, double potential)
     * @brief Modulate particle spectra at Sun position.
     * @param particle Pointer to a nuclear species.
     * @param potential Modulation potential.
     * @return Modulated spectrum at Sun position.
     * @warning Make use of GSL cspline.
     */
    //vector<double> particle_modulated_spectrum(TParticle* particle, double potential);
    
    /**
     * @fn TParticle* FindParticle(int uid)
     * @brief Find particle in particle array.
     * @param uid Unique ID of nuclear species.
     * @return Particle density of specified type
     */
    TParticle* FindParticle(int uid, bool sec = false, int isDM = 0, int isextra=0);
    vector<TParticle*> FindSpecies(int Z);
   
   inline double GetNorm() const { return norm; }
   inline double GetNormEl() const { return normel; }
   inline double GetNormExtra() const { return normextra; }
   inline double GetNormDM() const { return normDM; }
   
protected:
    Galaxy* gal; /**< Pointer to a Galaxy object, where all the ambient information is stored. */
    Input* in; /**< Pointer to a Input object, an interface to user input. */
    vector<TCREvolutorBasis*> alg; /**< Container of possibly many resolution algorithms to be used in cascade. */
    TNucleiList* list; /**< List of nuclei to be propagated. */
    vector<TParticle*> particles; /**< Container of nuclei densities and spectra. */
    fitsfile* output_ptr; /**< Output structure for spectra and galactic densities of CRs. Used if fullstore == true. \sa fullstore */
    fitsfile* output_ptr_sp; /**< Output structure only for solar CR spectra. Used if partialstore == true. \sa partialstore */
    int status; /**< FITS error status. */
    double normDM;
    double norm;
    double normel;
    double normextra;
   /*
  std::vector <double> e_prim;
  std::vector <double> e_sec;
  std::vector <double> p_prim;
  std::vector <double> p_sec;
  std::vector <double> pb_sec;
  std::vector <double> pb_ter;
  std::vector <double> ratio_pb_p;
  std::vector <double> eb;
  std::vector <double> C_12;
  std::vector <double> C_13;
  std::vector <double> C_14;
  std::vector <double> B_10;
  std::vector <double> B_11;
  std::vector <double> ratio_b_c;
  std::vector <double> pos_frac;
  std::vector <double> Be_9;
  std::vector <double> Be_10;
  std::vector <double> Be_10_9;
  std::vector <double> Al_26;
  std::vector <double> Al_27;
  std::vector <double> Al_26_27;
  std::vector <double> e_unmod;
  std::vector <double> eb_unmod;
  std::vector <double> p_unmod;
*/

     
    /**
     * @fn string make_filename(const char* argv, const char* mode) 
     * @brief Produce file name from run name.
     * @param argv Run name
     * @param mode extension
     */
    string make_filename(const char* argv, const char* mode) {
        string filename = "!output/";
        filename += argv;   
        return filename + mode;
    }
   
   TParticle* CreateParticle(const string&, int A_ /**< Mass number */, int Z_ /**< Charge */, Galaxy* gal /**< A model for the galaxy */, Input* in /**< User input */, bool issec, vector<TXSecBase*> xsecmodel, TNucleiList* l, int K_electron, bool isDM, bool isextra, bool isTPP_=false);
};
