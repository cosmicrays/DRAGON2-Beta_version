/**
 * @file nucleilist.cc
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class TNucleiList is implemented. @sa nucleilist.h
 */

#include "nucleilist.h"
#include "constants.h"
#include "utilities.h"
#include "input.h"

#include <algorithm>
#include <stdlib.h>
#include <fstream>
#include <string>

using namespace std;

/**
 * @fn bool comparator(int i, int j)
 * @brief Function defining a weak ordering of nuclei, so that they are propagated from the heaviest to the lightest, with the appropriate inversions to take into account beta decays.
 * @param i uid of the first nucleus.
 * @param j uid of the second nucleus.
 */
bool comparator(int i, int j) {
   int Ai = -1000;
   int Zi = -1000;
   Utility::id_nuc(i, Ai, Zi);
   
   int Aj = -1000;
   int Zj = -1000;
   Utility::id_nuc(j, Aj, Zj);
   
   if ( (Ai == Aj) && ( (Zi == 4 && Zj == 5) || (Zi == 6 && Zj == 7) || (Zi == 1 && Zj == 2) ||(Zi == 14 && Zj == 16)) ) return true;
   if ( ( (Ai==14) && (Aj==15) && (Zi == 6) && (Zj == 6) ) || ( (Ai==14) && (Aj==16) && (Zi == 6) && (Zj == 6) ) ) return true; //To make 6014 --> 7014 in the correct order (first 6014, then 7014). The rest of short lived beta decays must be ghost nuclei!    modified by Pedro 12/02/2020.
   if (Zi == Zj && (Ai>Aj)) return true;
   if (Zi>Zj) return true;
   return false;
}

DECMODE strtoDecMode(const string& str) {
   
   if (str == "EC") return EC;
   if (str == "B-") return BM;
   //if (str == "B+") return BP; //only two channels have b+ decay, and they are very short lived. They are just included as ghost -- 12/02/2020
   if (str == "BB") return BB;
   if (str == "IT") return IT;
   if (str == "ECB-") return ECBM;
   if (str == "ECB+") return ECBP;
   return STABLE;
}

TNucleiList::TNucleiList(Input* in) {
   ifstream infile(BNLdata.c_str(), ios::in); // Read databasis

   cout << "TNucleiList --> " << BNLdata.c_str() << endl;
   char s[3000];
   infile.getline(s,3000);
   
   int A = 0;
   int Z = 0;
   double tau = 0.;
   double tau_naked = 0.;
   string Mode;
   string name;
   
   int id;
   while (infile >> A >> Z >> tau >> tau_naked >> Mode >> name) {
      
      if (Z > in->Zmax || Z < in->Zmin) continue; /**< Select only nuclei in the wished charge range. \sa constants.h */
      id = 1000*Z+A; // Compute uid
      cout << "^^ id " << A << " " << Z << " " << tau << " " << tau_naked << " " << Mode << " " << name << endl;
      if (tau > minlifetime || Mode == "STABLE") {
         uid.push_back(id);
         lifetime[id]       = tau/Myr;
         lifetime_naked[id] = tau_naked/Myr;
         cout << "...adding to the list this nucleus: id = " << id << "; lifetime in Myr = " << tau/Myr << endl;
         betadec[id] = strtoDecMode(Mode);
      }

      else id_short.push_back(id); // Short lived nucleus; THIS LIST IS NOT USED. Short-lived nuclei are not propagated but are treated properly in the spallation network
   }
   
   infile.close();
   
   sort(uid.begin(), uid.end(), comparator); // Sort nuclei according to ordering defined by function comparator.
   sort(id_short.begin(), id_short.end(), comparator);
   
   if (in->prop_ap || in->prop_DMap) {
      uid.push_back(1000*(-1)+1); // Secondary antiprotons
      lifetime[1000*(-1)+1] = -1;
      betadec[1000*(-1)+1] = STABLE;
      
   }
   if (in->prop_lep || in->prop_extracomp || in->prop_DMel) {
      
      uid.push_back(1000*(-1)); // Primary electrons
      lifetime[1000*(-1)] = -1;
      betadec[1000*(-1)] = STABLE;
      
      uid.push_back(1000*(1)); // Secondary positrons
      lifetime[1000*(1)] = -1;
      betadec[1000*(1)] = STABLE;
      
   }
   
   if (in->prop_deuteron) {
      uid.push_back(1000*(-1)+2); // secondary antideuterons
      lifetime[1000*(-1)+2] = -1;
      betadec[1000*(-1)+2] = STABLE;
   }
   
   return ;
}
