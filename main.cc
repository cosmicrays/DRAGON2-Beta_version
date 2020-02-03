/**
 * @mainpage %DRAGON Documentation: Index Page
 *
 * @section intro_sec Introduction
 *
 * The CR propagation equation from a continuos distribution of sources can be written in the general form
 * \f[ \frac{\partial N^i}{\partial t} - {\nabla}\cdot \left( D\,{\nabla}
 -v_{c}\right)N^{i} + \frac{\partial}{\partial p} \left(\dot{p}-\frac{p}{3}{\nabla}\cdot{v}_{c}\right) N^i -\frac{\partial}{\partial p} p^2 D_{pp} \frac{\partial}{\partial p} \frac{N^i}{p^2} =  Q^{i}(p,r,z) + \sum_{j>i}c\beta n_{\rm gas}(r,z)
 \sigma_{ji}N^{j} -  c\beta n_{\rm gas}\sigma_{\rm in}(E_{k})N^{i} \f]
 *
 * Here \f$N^i(p,r,z)\f$ is the number density of the \f$i\f$-th atomic species; \f$p\f$ is its momentum; 
 * \f$\beta\f$ its velocity in units of the speed of light \f$c\f$; 
 * \f$\sigma_{in}\f$ is the total inelastic cross section onto the ISM gas, whose density is \f$n_{\rm gas}\f$; 
 * \f$\sigma_{ij}\f$ is the production cross-section of a nuclear species \f$j\f$ by the fragmentation of the \f$i\f$-th one; 
 * \f$D\f$ is the spatial diffusion coefficient; \f$v_{c}\f$ is the convection velocity. 
 * The last term on the l.h.s. describes diffusive reacceleration of CRs in the turbulent galactic magnetic field. 
 * 
 * %DRAGON adopts a second-order Cranck-Nicholson scheme with Operator Splitting and time overrelaxation to solve the diffusion 
 * equation. This provides fast a solution that is enough accurate for the average user. Occasionally, users may want to have
 * very accurate solutions to their problem. To enable this feature, users may get close to the accurate solution by using the fast 
 * method, and then switch to a more accurate solution scheme, featuring the Alternating-Direction-Implicit (ADI) Cranck-Nicholson 
 * scheme.
 *
 * Some parts of %DRAGON are built following <a href="http://galprop.stanford.edu/">GALPROP</a>, v50p. The first reason is that it is a waste of time to reimplement standard parts, 
 * like energy losses, in which nothing new has to be found. The second reason is that it is essential to be able to compare our predictions with that of the Galprop code, and this 
 * can be done only by following the details of its implementation. Therefore, we kept in the code some features and models used in Galprop, like nuclear cross-sections, the gas 
 * distribution, the convergence \"technique\". However, each of these models is accompanied by other models, which can be selected by setting the appropriate switch. This is done 
 * very easily using the well known C++ structure of abstract/derived classes. The code is then very flexible and easy to manage and to modify or update. 
 *
 * The code was built having in mind a few motivations:
 * - in order to find good propagation models one needs to run the code thousand times. Therefore we wrote the code aiming at performances and with an efficient memory management.
 * - we wanted to propagate DM originated cosmic rays. Therefore we wrote %DRAGON as a library that can be coupled to, e.g., DarkSUSY.
 * - from the physics point of view, we wanted to have a position dependent diffusion model, which requires a substantial modification of the discretization scheme.
 *
 * @section install_sec Installation
 *
 * %DRAGON comes with one library and one executable. The library contains the whole structure that is used to solve 
 * the CR propagation equation, and can be linked against other programs exploiting %DRAGON classes. 
 * The executable is the result of coupling the %DRAGON library with a driver routine, which reads user's input 
 * and solves the transport equation.
 *
 * @subsection step1 Step 1:  Download and unpack the code.
 * %DRAGON can be found <a href="OFFICIAL_DRAGON.tar.gz">here</a>.
 *
 * Then just tar xfz %DRAGON.tar.gz
 * @subsection step2 Step 2: Run the configure macro. 
 * configure needs some input from the user. User must tell configure where some external libraries are located.
 * %DRAGON requires <a href="http://www.gnu.org/software/gsl/">GSL libraries</a> 
 * and the <a href="http://heasarc.gsfc.nasa.gov/fitsio/">CFITSIO library</a>.
 * It also optionally requires <a href="http://www.physto.se/~edsjo/darksusy/">DarkSUSY</a>, if user wants to propagate DM originated cosmic rays.
 *
 * A typical command line is:
 * 
 * ./configure --with-gsl-path=<PATH_TO_GSL_EXEC> --with-cfitsio-include=<PATH_TO_CFITSIO_INCLUDE_FILES> --with-cfitsio-library=<PATH_TO_CFITSIO_LIBRARY> FFLAGS='-fPIC'
 *
 * This will preset %DRAGON to propagate astrophysically generated cosmic rays. 
 * If user wants instead to propagate DM originated cosmic rays, then he should configure with the following command line:
 *
 * ./configure --with-gsl-path=<PATH_TO_GSL_EXEC> --with-cfitsio-include=<PATH_TO_CFITSIO_INCLUDE_FILES> --with-cfitsio-library=<PATH_TO_CFITSIO_LIBRARY> FFLAGS='-fPIC' 
 * --with-ds-include=<PATH_TO_DARKSUSY_INCLUDE_FILES> --with-ds-library=<PATH_TO_DARKSUSY_LIBRARY_FILES>
 *
 * In this case %DRAGON will only propagate DM originated cosmic rays. 
 * The default installation path is in the same folder as the 
 * source code is (the program automatically creates the bin/ and lib/ subfolders). 
 * It can be set via --prefix=<NEW_INSTALLATION_PATH>
 * @subsection step3 Step 3: Compile the code.
 * Just run make and make install. The executable will be in $PREFIX/bin.
 * @subsection step4 Step 4: Enjoy DRAGON.
 * Run $PREFIX/bin/%DRAGON <argument_list>
 *
 * @section Example_sec Examples
 *  
 */

/**
 * @file main.cc
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief This is the driver routine.
 *
 * It takes in input the name for the run and the dimensions of the grid.
 * Then it instantiates an object of type DRAGON, makes the initializations
 * and takes care of the printout of the final results. \sa Input
 */

#include <iostream>
#include<ctime>
#include <math.h>

#include "input.h"
#include "dragon.h"
#include "constants.h"
#include "errorcode.h"

using namespace std;

int main(int argc, char** argv) {
   
    if (argc != 2) {
        cerr << "Usage: DRAGON <xml-file>" << endl;
        exit(NOFILE);
    }
    
    time_t time_s;
    time_t time_e;
    
    time(&time_s); // fix initial time

    Input* inp = new Input();
    inp->LoadFile( argv[1] );
    if (inp->feedback >0) inp->Print();
        
    if (inp->feedback >0) cout << "**************************" << endl;
    if (inp->feedback >0) cout << "*** Welcome to DRAGON! ***" << endl;
    if (inp->feedback >0) cout << "**************************" << endl << endl;
    
    //cout << "*** Reading xml input file... " << argv[1] << " ***" << endl;
    
    if (inp->feedback >0) cout << "*** Creating DRAGON  ***" << endl << endl;
    DRAGON* dr = new DRAGON(inp);
    
    if (inp->feedback >0) cout << "*** End of DRAGON constructor ***" << endl << endl;

    dr->Print();
    
    dr->Run();
    
    //if(inp->write_flag) dr->CalChi2();
 
    delete dr;
    delete inp;
    
    time(&time_e); // fix final time
    cout << "Solution found in " << (double)(time_e-time_s) << " s." << endl;
    
    return 0;
}
