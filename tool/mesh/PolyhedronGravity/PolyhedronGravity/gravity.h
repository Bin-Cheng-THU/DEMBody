//
//  gravity.h
//  
//
//  Created by Yang Yu on 10/25/17.
//
//

#ifndef _GRAVITY_H_
#define _GRAVITY_H_

#include "body.h"

void GravAttraction(POLYHEDRON &p, Vector r, Vector F);
void GravPotential(POLYHEDRON &p, Vector r, double &U);

#endif /* defined(____gravity__) */
