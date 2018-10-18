//
//  Created by Yang Yu on 10-11-17.
//
//

#ifndef _PARTICLES_H_
#define _PARTICLES_H_

#include "stdlib.h"
#include <sstream>
#include "body.h"
#include "mat3d.h"
#include "gravity.h"

//#include <omp.h>

using namespace std;

typedef struct pointsData {
    //
	int PointsNum;
    double **PointsPos, **PointsForce;
} POINTSDATA;

void LoadPoints(string pointsfile, POINTSDATA &m);
void AllocPoints(POINTSDATA &m);
void WritePoints(string pointsfile, POINTSDATA &m);
void ConductPoints(POINTSDATA &m, POLYHEDRON &p);

#endif
