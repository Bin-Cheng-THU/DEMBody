//
//  
//
//  Created by Yang Yu on 14-11-11.
//
//

#ifndef _CONSTANT_H_
#define _CONSTANT_H_

// input & output const settings
const int WidthOutput = 7; // The output width of the step number in the file name
const int NumParameter = 11; // The number of input parameters from PAR file
const int WidthInt = 22; // The output width of integer data
const int WidthDouble = 22; // The output width of double data
const int PrecDouble = 12; // The output precession of double data

// solver constants
const int EqnDimAsteroid = 7; // The dimension of binary equations
const int EqnDimDebris = 6; // The dimension of debris equations
const double PrecM2E = 1.0e-15; // The iterational precession of Mean Anomaly to Eccentric Anomaly
//const double ZeroVec[3] = {0.0, 0.0, 0.0};

// solar system physical constants
const double PI = 3.14159265359; // The value of PI
const double G = 6.674184e-11; // The gravitational constant, unit: m^3/kg/s^2
const double SM = 1.9891e30; // The solar mass, unit: kg
const double AU = 1.49597871e11; // The astronomical unit, unit: m
const double SRF = 1.361e3; // The solar radiation flux at 1 AU, unit: kg/s^3
const double LS = 2.988e8; // The light speed value, unit: m/s

// IDs of asteroid components
const int ID_Asteroid = -1; // The ID of the asteroid

// global constants for R-K-F solver
const double RKFa[13]=
{ 0.0, 2.0/27.0,  1.0/9.0,  1.0/6.0,  5.0/12.0, 1.0/2.0, 5.0/6.0, 1.0/6.0, 2.0/3.0, 1.0/3.0, 1.0,  0.0, 1.0};
const double RKFb[13][12]=
{
    {  0.0,            0.0,       0.0,         0.0,           0.0,           0.0,         0.0,           0.0,       0.0,        0.0,       0.0, 0.0, },
    {  2.0/27.0,       0.0,       0.0,         0.0,           0.0,           0.0,         0.0,           0.0,       0.0,        0.0,       0.0, 0.0, },
    {  1.0/36.0,       1.0/12.0,  0.0,         0.0,           0.0,           0.0,         0.0,           0.0,       0.0,        0.0,       0.0, 0.0, },
    {  1.0/24.0,       0.0,       1.0/8.0,     0.0,           0.0,           0.0,         0.0,           0.0,       0.0,        0.0,       0.0, 0.0, },
    {  5.0/12.0,       0.0,       -25.0/16.0,  25.0/16.0,     0.0,           0.0,         0.0,           0.0,       0.0,        0.0,       0.0, 0.0, },
    {  1.0/20.0,       0.0,       0.0,         1.0/4.0,       1.0/5.0,       0.0,         0.0,           0.0,       0.0,        0.0,       0.0, 0.0, },
    {  -25.0/108.0,    0.0,       0.0,         125.0/108.0,   -65.0/27.0,    125.0/54.0,  0.0,           0.0,       0.0,        0.0,       0.0, 0.0, },
    {  31.0/300.0,     0.0,       0.0,         0.0,           61.0/225.0,    -2.0/9.0,    13.0/900.0,    0.0,       0.0,        0.0,       0.0, 0.0, },
    {  2.0,            0.0,       0.0,         -53.0/6.0,     704.0/45.0,    -107.0/9.0,  67.0/90.0,     3.0,       0.0,        0.0,       0.0, 0.0, },
    {  -91.0/108.0,    0.0,       0.0,         23.0/108.0,    -976.0/135.0,  311.0/54.0,  -19.0/60.0,    17.0/6.0,  -1.0/12.0,  0.0,       0.0, 0.0, },
    {  2383.0/4100.0,  0.0,       0.0,         -341.0/164.0,  4496.0/1025.0, -301.0/82.0, 2133.0/4100.0, 45.0/82.0, 45.0/164.0, 18.0/41.0, 0.0, 0.0, },
    {  3.0/205.0,      0.0,       0.0,         0.0,           0.0,           -6.0/41.0,   -3.0/205.0,    -3.0/41.0, 3.0/41.0,   6.0/41.0,  0.0, 0.0, },
    {  -1777.0/4100.0, 0.0,       0.0,         -341.0/164.0,  4496.0/1025.0, -289.0/82.0, 2193.0/4100.0, 51.0/82.0, 33.0/164.0, 12.0/41.0, 0.0, 1.0, }
};
const double RKFc[13]=
{ 0.0, 0.0, 0.0, 0.0, 0.0, 34.0/105.0, 9.0/35.0, 9.0/35.0, 9.0/280.0, 9.0/280.0, 0.0, 41.0/840.0, 41.0/840.0};

#endif
