//
//  Created by Yang Yu on 14-11-11.
//
//

#ifndef _VOLINT_H_
#define _VOLINT_H_

#include <iostream>
#include "mat3d.h"
#include "polyhedron.h"

// volint.cpp

void VolumeIntegrals(POLYHEDRON &p);
void CompWfNf(POLYHEDRON &p, double *&Wf, double **&Nf);
void FaceIntegrals(POLYHEDRON &p, double wf, double nf[3], int id, int a, int b, int c, double F[12]);
void ProjectionIntegral(POLYHEDRON &p, int a, int b, int id, double P[10]);
void VolumeInfoCalcu(double VolIntZeroth, double VolIntFirst[3], double VolIntSecond[6], POLYHEDRON &p);

#endif
