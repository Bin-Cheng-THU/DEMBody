//
//  solar.cpp
//
//
//  Created by Yang Yu on 10-11-16.
//
//

#include "solar.h"

void KepCopy(KEPORB &kep, KEPORB &kep1)
{
    kep1.SemiMajorAxis = kep.SemiMajorAxis;
    kep1.Eccentricity = kep.Eccentricity;
    kep1.LongAscendNode = kep.LongAscendNode;
    kep1.Inclination = kep.Inclination;
    kep1.ArgPeriapsis = kep.ArgPeriapsis;
    kep1.MeanAnomaly = kep.MeanAnomaly;
    //
    vectorCopy(kep.Position, kep1.Position);
    vectorCopy(kep.Velocity, kep1.Velocity);
}

void KepStepping(ODE_OPTION &odeopt, KEPORB &kep)
{
    double n;
    n = sqrt(G*SM/(kep.SemiMajorAxis*kep.SemiMajorAxis*kep.SemiMajorAxis));
    kep.MeanAnomaly = odeopt.StepSize*n + kep.MeanAnomaly;
    
    KEP2RV(kep); 
}

void KEP2RV(KEPORB &kep)
{
    double E0, M0, err, E, theta, r, v, mu, w, gamma;
    double tmpd, cw, sw, cwg, swg, cl, sl, ci, si;
    
    mu = G*SM;
    
    E0 = kep.MeanAnomaly;
    while (1)
    {
        M0 = E0 - kep.Eccentricity*sin(E0);
        err = kep.MeanAnomaly - M0;
        if (fabs(err)<PrecM2E)
        {
            E = E0;
            break;
        }
        E0 = E0 + err/(1.0 - kep.Eccentricity*cos(E0));
    }
    
    tmpd = tan(E/2.0)/sqrt((1.0-kep.Eccentricity)/(1.0+kep.Eccentricity));
    theta = atan(tmpd)*2.0;
    
    r = kep.SemiMajorAxis*(1.0-kep.Eccentricity*kep.Eccentricity)/(1.0+kep.Eccentricity*cos(theta));
    v = sqrt(2.0*mu/r - mu/kep.SemiMajorAxis);
    
    gamma = atan(kep.Eccentricity*sin(theta)/(1.0+kep.Eccentricity*cos(theta)));
    
    w = theta + kep.ArgPeriapsis;
    
    cw = cos(w);
    sw = sin(w);
    cl = cos(kep.LongAscendNode);
    sl = sin(kep.LongAscendNode);
    ci = cos(kep.Inclination);
    si = sin(kep.Inclination);
    cwg = cos(w - gamma);
    swg = sin(w - gamma);
    
    kep.Position[0] = r*(cw*cl - sw*ci*sl);
    kep.Position[1] = r*(cw*sl + sw*ci*cl);
    kep.Position[2] = r*sw*si;
    
    kep.Velocity[0] = v*(-swg*cl - cwg*ci*sl);
    kep.Velocity[1] = v*(-swg*sl + cwg*ci*cl);
    kep.Velocity[2] = v*cwg*si;
}











