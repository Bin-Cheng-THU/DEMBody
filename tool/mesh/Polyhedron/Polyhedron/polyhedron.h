//
//  Created by Yang Yu on 14-11-11.
//
//

#ifndef _POLYHEDRON_H_
#define _POLYHEDRON_H_

#include "mat3d.h"

typedef struct polyhedron {
    
    double Volume;
    Vector Centroid;
    Matrix InertiaTensor;
    int NumVerts, NumFaces, NumEdges;
    int **Faces, **Edges;
    double *EdgeLens;
    double **Vertices, **EdgeNormVecs, **FaceNormVecs;
    
} POLYHEDRON;

#endif
