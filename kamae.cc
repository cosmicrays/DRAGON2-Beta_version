/**
 * @file kamae.cc
 * @brief Implementation file for cparamlib
 *
 * It implements Kamae et al model for secondary electron/positron production in pp interactions. Obtained from ... It will not be documented.
 */

#include "kamae.h"
#include "errorcode.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

PARAMSET parameters;

double KamaeYields::GetSigma(double Eg, double Ekinbullet, PARTICLE_IDS par) {
  if (Ekinbullet >= 4.e+5) Ekinbullet = 4.e+5;
    switch (par) {
        case ID_GAMMA :
            gamma_param(Ekinbullet, &parameters);
            break;
        case ID_ELECTRON :
            elec_param(Ekinbullet, &parameters);
            break;
        case ID_POSITRON :
            posi_param(Ekinbullet, &parameters);
            break;
        default:
            cerr << "Wrong particle ID in KamaeYields" << endl;
            exit(WRONGID);
    }
    
    return 1.e-3*sigma_incl_tot(par,Eg, Ekinbullet, &parameters)/Eg;
}

