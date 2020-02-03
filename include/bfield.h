//
//  bfield.h
//  DRAGON
//
//  Created by Luca Maccione on 04/06/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef DRAGON_bfield_h
#define DRAGON_bfield_h

/**
 * @class TBField
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Description of a magnetic field. Only a 2D, smooth magnetic field is modeled here. \sa rB \sa zr \sa Bh
 */

#include <vector>
#include <math.h>
#include <map>

class TGrid;
class TGeometry;

class TBField {
    
public:
    TBField() {} /**< Default constructor */
    TBField(TGrid*); /**< Constructor in a given geometry. */
    virtual ~TBField() { MagField.clear(); }
    /**< Destructor. */
    
    inline std::vector<double>& GetBField() { return MagField; }
    /**< Get the magnetic field distribution. */
    //inline double GetTurb() { return eta; }
    /**< Get turbulence factor. */
    double GetEnDensity(int i /**< linearized position index */);
    /**< Get energy density (linearized representation). */
    double GetEnDensity(int ir /**< radial position */, int iz /**< vertical position */);
    double GetEnDensity(int ix /**< radial position */, int iy, int iz /**< vertical position */);
    /**
     * @fn double GetEnDensity(int ir, int iz)
     * @brief Get energy density with matrix to linear representation wrapping.
     * @return [eV/cm^3] Magnetic energy density.
     */
    inline std::map<int, std::vector<double> >& GetBreg() {return Breg; }
    
    inline std::map<int, std::vector<double> >& GetBregVersors() {return Breg_versors; }

    std::vector<double> GetBreg(int counter) {
        std::vector<double> Breturn(3, 0);
        Breturn[0] = Breg[0][counter];
        Breturn[1] = Breg[1][counter];
        Breturn[2] = Breg[2][counter];
        return Breturn;
    }

    inline std::vector<double> GetBreg(int ix, int iy, int iz) { return GetBreg(index(ix,iy,iz)); }

    std::vector<double> GetBregVersors(int counter) {
        std::vector<double> Breturn(3, 0);
        Breturn[0] = Breg_versors[0][counter];
        Breturn[1] = Breg_versors[1][counter];
        Breturn[2] = Breg_versors[2][counter];
        return Breturn;
    }
    
    inline std::vector<double> GetBregVersors(int ix, int iy, int iz) { return GetBregVersors(index(ix,iy,iz)); }
    virtual std::vector<double> GetField(double x, double y, double z)=0;
    std::vector<double> GetBregVersors(double x, double y, double z);
    
protected:
    std::vector<double> MagField; /**< Magnetic field distribution, in muG. */
    std::map<int, std::vector<double> > Breg;
    std::map<int, std::vector<double> > Breg_versors;

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    int dimx; /**< Dimension of radial box. */
    int dimy; /**< Dimension of radial box. */
    int dimz; /**< Dimension of vertical box. */
    inline int index(int ir /**< radial position */, int iz /**< vertical position */) { return ir*dimz+iz; }
    /**
     * @fn inline int index(int ir, int iz)
     * @brief Convert matrix to linear representation.
     * @return ir*dimz+iz
     */
    inline int index(int ix /**< radial position */, int iy, int iz /**< vertical position */) { return (ix*dimy+iy)*dimz+iz; }
    /**
     * @fn inline int index(int ir, int iz)
     * @brief Convert matrix to linear representation.
     * @return ir*dimz+iz
     */
};

class TUniformField : public TBField {
  
public:
    TUniformField() { }
    TUniformField(double B, TGrid* Coord, TGeometry*);
    virtual std::vector<double> GetField(double x, double y, double z) { return std::vector<double>(3,0); }
};


class TPshirkovField : public TBField {

public :
    TPshirkovField() { }
    TPshirkovField(double B0, double RC, double Z0, double R0, double B0H, double R0H, double Z0H, double B0turb, double rscale_turb, double zscale_turb, double robs, TGrid*, TGeometry*);
    virtual std::vector<double> GetField(double x, double y, double z);
    
    protected :
    double B0;
    double RC;
    double Z0;
    double R0;
    double B0H;
    double R0H;
    double Z0H;
    double B0turb;
    double rscale_turb;
    double zscale_turb;
    double robs;
};


class TFarrarField : public TBField {
    
    public :
    TFarrarField() :
    bring(0.1),
    hdisk(0.40),
    wdisk(0.27),
    Bn(1.4),
    Bs(-1.1),
    rn(9.22),
    rs(16.7),
    wh(0.2),
    z0(5.3),
    BX(4.6),
    Theta0X(49.0*M_PI/180.0),
    rcX(4.8),
    rX(2.9),
    p(11.5*M_PI/180.0)
    {
        bj = std::vector<double>(8,0);
        bj[0] = 0.1;
        bj[1] = 3.0;
        bj[2] = -0.9;
        bj[3] = -0.8;
        bj[4] = -2.0;
        bj[5] = -4.2;
        bj[6] = 0.0;
        bj[7] = 0.0;

       fj = std::vector<double>(8,0);
       fj[0] = 0.130;
       fj[1] = 0.165;
       fj[2] = 0.094;
       fj[3] = 0.122;
       fj[4] = 0.13;
       fj[5] = 0.118;
       fj[6] = 0.084;
       fj[7] = 0.156;
       
       
       for (int i = 0; i < 7; i++) bj[7] -= (fj[i]*bj[i]/fj[7]);

       
        rj = std::vector<double>(8,0);
       rj[0] = 5.1;
       rj[1] = 6.3;
       rj[2] = 7.1;
       rj[3] = 8.3;
       rj[4] = 9.8;
       rj[5] = 11.4;
       rj[6] = 12.7;
       rj[7] = 15.5;

    }
    TFarrarField(double,TGrid*, TGeometry*);
    virtual std::vector<double> GetField(double x, double y, double z);
    
    protected :
    std::vector<double> GetFieldCyl(double r, double phi, double z);
    std::vector<double> GetAvField(double r, double z);
    
    std::vector<double> bj;
    std::vector<double> fj;
    std::vector<double> rj;
    double bring;
    double hdisk;
    double wdisk;
    
    double Bn;
    double Bs;
    double rn;
    double rs;
    double wh;
    double z0;
    
    double BX;
    double Theta0X;
    double rcX;
    double rX;
    
    double p;

};

class ToyModelField : public TBField {

public :
    ToyModelField() { }
    ToyModelField(double Bx_, double By_, double Bz_, double B0turb_, TGrid* Coord, TGeometry*);
    virtual std::vector<double> GetField(double x, double y, double z);
    
    protected :
    double Bx, By, Bz;
    double B0turb;
};

class TSimpleField : public TBField {
  
public:
    TSimpleField() { }
    TSimpleField(double B, double rsun, TGrid* Coord, TGeometry*);
//     virtual std::vector<double> GetField(double x, double y, double z) { return std::vector<double>(3,0); } //MW130710: Compatibility with Daniele's and Luca's structures, when merging Bfields
    virtual std::vector<double> GetField(double x, double y, double z); //MW130710: Compatibility with Daniele's and Luca's structures, when merging Bfields

protected:
    double B0;
    double robs;
    
    double r_scale;
    double z_scale;
    double r_temp;

};

#endif
