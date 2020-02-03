/**
 * @file utilities.h
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Definition of some useful functions.
 */

#ifndef _UTILITIES_H
#define _UTILITIES_H

#include <vector>
#include <string>
#include "grid.h"

class Input;

/**
 * @namespace Utility
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Namespace grouping useful functions.
 */
namespace Utility {
    
    const int dimDATA = 38; /**< Dimension of B/C, N/O and C/O databasis. */
    const int dimDATABe = 11; /**< Dimension of 10Be/9Be databasis. */
    
    /**
     * @struct data
     * @author Carmelo Evoli
     * @email evoli@sissa.it
     * @brief Structure of a data point for B/C, N/O and C/O joint analysis.
     */
    struct data {
        double E; /**< Energy of experimental point. */
        double BC; /**< Observed B/C. */
        double sigma_BC; /**< Experimental error on B/C. */
        double CO; /**< Observed C/O. */
        double sigma_CO; /**< Experimental error on C/O. */
        double NO; /**< Observed N/O. */
        double sigma_NO; /**< Experimental error on N/O. */
        double phi;  /**< Modulation potential to be adopted for the specific experiment. */      
        int experiment; /**< Experiment ID. */
    };
    
    /**
     * @struct data_be
     * @author Carmelo Evoli
     * @email evoli@sissa.it
     * @brief Structure of a data point for 10Be/9Be analysis.
     */
    struct data_be {
        double E; /**< Energy of experimental point. */
        double Be; /**< Observed 10Be/9Be. */
        double sigma_Be; /**< Experimental error on 10Be/9Be. */
        double phi; /**< Modulation potential to be adopted for the specific experiment. */
        int experiment ;/**< Experiment ID. */
    };
    
    /**
     * @struct chi2_struct
     * @author Carmelo Evoli
     * @email evoli@sissa.it
     * @brief Structure for Chi2 computation.
     */
    struct chi2_struct{
        double all[3]; /**< Computed for all experiments. all[0] = B/C, all[1] = C/O, all[2] = N/O.*/
        int N_all; /**< Number of degrees of freedom*/
        double HEAO3[3];  /**< Computed for HEAO3 experiment. all[0] = B/C, all[1] = C/O, all[2] = N/O.*/
        int N_HEAO3; /**< Number of degrees of freedom*/
        double CREAM[3];  /**< Computed for CREAM experiment. all[0] = B/C, all[1] = C/O, all[2] = N/O.*/
        int N_CREAM; /**< Number of degrees of freedom*/
        double ATIC[3];  /**< Computed for ATIC experiment. all[0] = B/C, all[1] = C/O, all[2] = N/O.*/
        int N_ATIC; /**< Number of degrees of freedom*/
    };
    
    /**
     * @fn void solve_tridag(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&)
     * @author Luca Maccione
     * @email luca.maccione@desy.de
     * @brief Solve a tridiagonal system, from Numerical Recipes.
     * @param a Lower diagonal
     * @param b Diagonal
     * @param c Upper diagonal
     * @param r Termini noti
     * @param u Solution
     * @return Nothing, but updates u
     */
    //  void solve_tridag(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, std::vector<double>&);
    //void solve_tridag(double* a, double* b, double* c, double* r, double* u, int n);
    void solve_tridag(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& r, std::vector<double>& u, int n);
    
    /**
     * @fn  void insert_data(int, double, double, double, double, double, double, double, double, int, Utility::data*);
     * @author Carmelo Evoli
     * @email evoli@sissa.it
     * @brief Fill B/C, N/O, C/O data structure. 
     * @param i
     * @param E
     * @param BC
     * @param sigma_BC
     * @param CO
     * @param sigma_CO
     * @param NO
     * @param sigma_NO
     * @param phi
     * @param experiment
     * @param data_set
     * @return Nothing, but updates data_set
     */
    void insert_data(int, double, double, double, double, double, double, double, double, int, Utility::data*);
    
    /**
     * @fn  void insert_data_be(int, double, double, double, double, int, Utility::data_be*);
     * @author Carmelo Evoli
     * @email evoli@sissa.it
     * @brief Fill 10Be/9Be data structure. 
     * @param i
     * @param E
     * @param Be
     * @param sigma_Be
     * @param phi
     * @param experiment
     * @param data_set
     * @return Nothing, but updates data_set
     */
    void insert_data_be(int, double, double, double, double, int, Utility::data_be*);
    
    /**
     * @fn void id_nuc(int, int&, int&)
     * @author Luca Maccione
     * @email luca.maccione@desy.de
     * @brief Compute A,Z from Unique ID.
     * @param uid Unique ID of the nucleus
     * @param A Mass number of the nucleus
     * @param Z Charge of the nucleus
     * @return Nothing, but updates A and Z
     */

    void id_nuc(int, int&, int&);
    
    //#ifdef HAVE_DS
    /**
     * @fn double dmspec_cc(double, int, double, double, int)
     * @author Luca Maccione
     * @email luca.maccione@desy.de
     * @brief Wrapper to DarkSUSY
     * @param mx Mass of DM particle
     * @param mode Annihilation channel
     * @param sigmav Thermally averaged cross section
     * @param egev Energy [GeV]
     * @param yieldk Particles to be output (151 = positrons, 154 = antiprotons)
     * @return Particle spectrum at energy egev.
     */
    //  double dmspec_cc(double, int, double, double[], int, double[]);
    //#endif
    
};

#endif
