/**
 * @file utilities.cc
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Implementation of useful functions
 */
#include "utilities.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>

#include "constants.h"
#include "input.h"
#include "grid.h"

using namespace std;

void Utility::solve_tridag(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& r, vector<double>& u, int n) {
 int j = 0;
 double bet = 0.0; 
 vector<double> gam(n,0.);  //double gam[n];
 //	One vector of workspace, gam, is needed. 
 if (b[0] == 0.0) cerr << "Error 1 in tridag: the first diagonal term is 0!! " << endl; 
 //If this happens, then you should rewrite your equations as a set of order N-1, with u1 trivially eliminated.
 bet = b[0];  
 u[0] = r[0] / bet; 
 for (j = 1; j < n; j++) {	//Decomposition and forward substitution.
   //double* gm = gam+j;
   //(*gm) = c[j-1]/bet;  
   gam[j] = c[j-1]/bet;	
   //bet = b[j] - a[j]*(*gm);
   bet = b[j] - a[j]*gam[j];
   if (bet == 0.0){ 
        cout << "j = 0 " << " --> diagonal term b[0] = " << b[0] << " off diagonal term a[0] = " << a[0] << " c[0] = " << c[0] << " u[0] = " << u[0] << " bet = b[0] " << endl;
        cout << "j = " << j << " --> diagonal term b[j] = " << b[j] << " off diagonal term a[j] = " << a[j] << " gam[j] = " << gam[j] << " bet = b[j] - a[j]*c[j-1]/bet " << bet << endl;
	cerr << "Error 2 in tridag: bet = 0!" << endl; 
    }   
    u[j] = (r[j] - a[j]*u[j-1])/bet;
 } 
 for (j = (n-2); j >= 0; j--)
 u[j] -= gam[j+1]*u[j+1];	//Backsubstitution.
 return ;
}

/*
void Utility::solve_tridag(double* a, double* b, double* c, double* r, double* u, int n) {
    
    int j = 0;
    double bet = 0.0; 
    
    double gam[n];
    //	One vector of workspace, gam, is needed. 
    if (b[0] == 0.0) cerr << "Error 1 in tridag" << endl; 
    //If this happens, then you should rewrite your equations as a set of order N-1, with u1 trivially eliminated. 
    u[0] = r[0] / (bet = b[0]); 
    double* gm;
    for (j = 1; j < n; j++) {	//Decomposition and forward substitution.
        a++;
        r++;
        u++;
        b++;
        gm = gam+j;
        (*gm) = (*c)/bet; 
        c++;
        bet = *b - (*a)*(*gm);
        if (bet == 0.0) cerr << "Error 2 in tridag" << endl; 
        *(u) = (*r - (*a)*(*(u-1)))/bet;
    } 
    u--;
    for (j = (n-2); j >= 0; j--) { 
        *u -= (*gm)*(*(u+1));	//Backsubstitution.
        u--;
        gm--;
    }
    return ;
}
*/


void Utility::insert_data(
                          int i, 
                          double E, 
                          double BC, double sigma_BC, 
                          double CO, double sigma_CO, 
                          double NO, double sigma_NO, 
                          double phi,
                          int experiment, 
                          Utility::data* data_set
                          )
{
    data_set[i].E = E;
    data_set[i].BC = BC;
    data_set[i].sigma_BC = sigma_BC;
    data_set[i].CO = CO;
    data_set[i].sigma_CO = sigma_CO;
    data_set[i].NO = NO;
    data_set[i].sigma_NO = sigma_NO;
    data_set[i].phi = phi;
    data_set[i].experiment = experiment;
    return;
}

void Utility::insert_data_be(
                             int i, 
                             double E, 
                             double Be, double sigma_Be, 
                             double phi,
                             int experiment, 
                             Utility::data_be* data_set
                             )
{
    data_set[i].E = E;
    data_set[i].Be = Be;
    data_set[i].sigma_Be = sigma_Be;
    data_set[i].phi = phi;
    data_set[i].experiment = experiment;
    return;
}

void Utility::id_nuc(int uid, int& A, int& Z) {
    
    if (uid == -999) {
        A = 1;
        Z = -1;
        return ;
    }
    if (uid == -998) {
        A = 2;
        Z = -1;
        return ;
    }
    
    A = int(uid%1000);
    Z = int(uid/1000);
    
    return ;
}

