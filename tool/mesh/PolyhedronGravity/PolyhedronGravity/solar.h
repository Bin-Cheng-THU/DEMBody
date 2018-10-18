//
//  Created by Yang Yu on 22-11-11.
//
//

#ifndef _SOLAR_H_
#define _SOLAR_H_

#include "constant.h"
#include "mat3d.h"
#include "odeset.h"
//using namespace std;

typedef struct keporb {
    
    double SemiMajorAxis; // unit: m
    double Eccentricity;
    double LongAscendNode; // unit: rad
    double Inclination; // unit: rad
    double ArgPeriapsis; // unit: rad
    double MeanAnomaly; // unit: rad
    
    Vector Position, Velocity; // SXYZ
    
} KEPORB;

typedef struct solar {
    
    int SolarTide;
    int SolarPressure;
    double Reflection;
    double SOI;         // radius of the sphere of influence of the binary system with respect to the sun
    // defined as: a*(mA/MS)^(2/5)
    
} SOLAR;

void KepCopy(KEPORB &kep, KEPORB &kep1);
void KepStepping(ODE_OPTION &odeopt, KEPORB &kep);
void KEP2RV(KEPORB &kep);

//void TrackOrbit(double t, ELLIPELEM &orbelem, Vector R);

#endif
