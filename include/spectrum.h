/**
 * @class TSpectrum
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Injection spectrum of CRs.
 */

#include <vector>
#include <iostream>
#include <string>
#include <math.h>
#include <fstream>
#include <sstream>

class TGrid;
class Input;

class TSpectrum {
    
 public:
  TSpectrum() {} /**< Default constructor. */
  TSpectrum(TGrid* coord); /**< Default constructor with some geometry. */
  TSpectrum(TGrid*, Input*, const std::vector<double>& break_positions, const std::vector<double>& slopes, double CutoffRig = -1, bool El = false, bool ExtraComponent = false);
  TSpectrum(TGrid*, Input*, int /**< Outgoing particles: 151 = antiprotons, 154 = positrons. */); 
  TSpectrum(TGrid*, std::string /**< File with the tabulated model. */, int yieldk);
  TSpectrum(TGrid* Coord, int mass, int dmmode, int yieldk);
    
  virtual ~TSpectrum() { }
  /**< Destructor. */
    
  inline const std::vector<double>& GetSpectrum() const { return spectrum; }
  /**< Get the total spectrum. */
  inline double GetSpectrum(int i /**< Energy index. */) { return spectrum[i]; }
  /**
   * @fn inline double GetSpectrum(int i)
   * @brief Get spectrum at some energy.
   */
    
 protected:
  std::vector<double> spectrum; /**< A CR spectrum. */
    
  std::string find_channel(int dmmode);
  /*{
    std::string channel;
    switch(dmmode) {
    case 13: 
    channel = "W";
    break;
    case 12:
    channel="Z";
    break;
    case 17:
    channel="Mu";
    std::cerr << "Warning: did you mean VMu ??" << std::endl;
    break;
    case 19:
    channel="Tau";
    std::cerr << "Warning: did you mean VTau ??" << std::endl;
    break;
    case 22:
    channel = "c";
    break;
    case 24:
    channel = "t";
    break;
    case 25:
    channel = "b";
    break;
    case 26:
    channel = "g";
    break;
    default:
    std::cerr << "The mode you requested was not implemented. You requested mode " << dmmode << std::endl;
    break;
    }
     
    return channel;
    }
  */
  double InterpSpectrum(const double& en, const std::vector<double>& energy, const std::vector<double>& orig_spectrum) const ;
  /*{ 
    int i = 0;
    while (energy[i] < en) ++i;
    if (i > energy.size()) std::cout << "index out of range" << std::endl;
    double r1 = log10(energy[i]/en)/log10(energy[i]/energy[i-1]);
    return pow(10., r1*log10(orig_spectrum[i-1]) + (1.0-r1)*log10(orig_spectrum[i]));
    }
  */
};

class TSpectrumExtra : public TSpectrum {
  
 public : 
  TSpectrumExtra(TGrid*, Input*);
  virtual ~TSpectrumExtra() { }
    
};
