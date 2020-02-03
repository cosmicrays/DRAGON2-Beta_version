//
//  geometry.h
//  DRAGON
//
//  Created by Luca Maccione on 17/07/13.
//
//

#ifndef __DRAGON__geometry__
#define __DRAGON__geometry__

#include <string>
#include <vector>

class Input;
class TGrid;

using namespace std;

class TGeometry {
    
public:
    TGeometry(TGrid* coord_, Input* in);
    virtual ~TGeometry() { }
    double GetPattern(int x, int y, int z);
    void ApplySpiral(std::vector<double>& dens, double index, double cut);
    
protected:
    TGrid* Coord;
    Input* inp;
    vector<double> density;
};

class TUniformGeometry : public TGeometry {
    
public:
    TUniformGeometry(TGrid* coord, Input* in) : TGeometry(coord, in) { }
    virtual ~TUniformGeometry() { }
    
};

class TSpiralGeometry : public TGeometry {
    
public:
    TSpiralGeometry(TGrid*, Input*);
    virtual ~TSpiralGeometry() { }
};


#endif /* defined(__DRAGON__geometry__) */
