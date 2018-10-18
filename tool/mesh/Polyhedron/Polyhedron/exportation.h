//
//  exportation.h
//  
//
//  Created by Yang Yu on 14-11-15.
//
//

#ifndef _EXPORTATION_H_
#define _EXPORTATION_H_

#include <iostream>
#include <fstream>
#include <iomanip>
#include "mat3d.h"
#include "polyhedron.h"

void PolyhedronOutput(char* filename,POLYHEDRON p);
void VolumeParaOutput(char* filename, double vol0, Vector cntr0, Matrix tnsr0, double lenscale, POLYHEDRON p);

#endif 
