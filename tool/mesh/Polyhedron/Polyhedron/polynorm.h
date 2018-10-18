//
//  Created by Yang Yu on 14-11-11.
//
//

#ifndef _POLYNORM_H_
#define _POLYNORM_H_

#include <iostream>
#include "mat3d.h"
#include "polyhedron.h"
#include "jacobi_eigen.h"
#include "constant.h"

// polynorm.cpp
void PolyNormalize(POLYHEDRON &p, double &lenscale);
void InertSort(double d[3], Matrix m);
void normVertices(POLYHEDRON &p, double a, Matrix m);

#endif
