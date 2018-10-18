//
//  Created by Yang Yu on 14-11-11.
//
//

#ifndef _MESH_H_
#define _MESH_H_

#include <iostream>
#include "mat3d.h"
#include "polyhedron.h"

// mesh.cpp

void VertFaceLoad(char* filename, POLYHEDRON &p);
void MeshParaAlloc(POLYHEDRON &p);
void MeshParaCalcu(POLYHEDRON &p);

#endif
