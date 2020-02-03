/**
 * @file nucleilist.h
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class TNucleiList is defined.
 */

#ifndef _NUCLEILIST_H
#define _NUCLEILIST_H

#include <iostream>
#include <vector>
#include <map>

/**
 * @class TNucleiList
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class to read the nuclei database and compute the list of nuclei to be propagated. Also lifetimes and decay channels are read in the database. 
 */

typedef enum {STABLE, BM, EC, BB, IT, BP, ECBM, ECBP} DECMODE;

class Input;
class TNucleiList {

 public:

  TNucleiList(Input*);
  inline std::vector<int>& GetList() { return uid; } 
  /**< Get list of nuclei (through Unique ID). */
  inline std::vector<int>& GetListShort() { return id_short; }
  /**< Get list of short lived nuclei (not used in the program). */
  inline std::map<int, double>& GetLifeTime() { return lifetime; }
  /**< Get lifetimes. */
  inline double GetLifeTime(int i) { return lifetime[i]; }
  inline double GetLifeTime_naked(int i) { return lifetime_naked[i]; }
  /**
   * @fn inline double GetLifeTime(int i)
   * @brief Obtain life time of specified nucleus.
   * @param i Unique ID of the nucleus.
   */
  inline DECMODE GetDecayMode(int i) { return betadec[i]; }
  /**
   * @fn inline int GetDecayMode(int i)
   * @brief Obtain decay mode of specified nucleus.
   * @param i Unique ID of the nucleus.
   */

  ~TNucleiList() {
    uid.clear();
    id_short.clear();
    lifetime.clear();
  }

 private: 
  std::vector<int> uid; /**< Array of unique IDs = 1000*Z+A .*/
  std::vector<int> id_short; /**< Array of short lived nuclei. */
  std::map<int, double> lifetime; /**< Lifetimes of nuclei (stored according to their uid). */
  std::map<int, double> lifetime_naked; /**< Lifetimes of nuclei (stored according to their uid). */
  std::map<int, DECMODE> betadec; /**< Whether nuclei decay beta+ or beta- . */

};
#endif
