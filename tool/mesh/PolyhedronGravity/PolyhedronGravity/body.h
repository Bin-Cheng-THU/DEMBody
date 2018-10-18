//
//  Created by Yang Yu on 22-11-11.
//
//

#ifndef _BODY_H_
#define _BODY_H_

#include <iostream>
#include <fstream>
#include <string>
#include "mat3d.h"
#include "constant.h"

using namespace std;

typedef struct polyhedron {
    
    Vector Inertia;
    Vector Ellipsoid;
    // the envelop tri-axial ellipsoid in the principal body-fixed frame for occultation checking
    int NumVerts, NumFaces, NumEdges;
    int **Faces, **Edges;
    double Volume, Density, Mass;
    double *EdgeLens;
    double **Vertices, **EdgeNormVecs, **FaceNormVecs;
    
} POLYHEDRON;

void LoadPolyhedron(string polyfile, double density, POLYHEDRON &p);
void PolyhedronAlloc(POLYHEDRON &p);

#endif
