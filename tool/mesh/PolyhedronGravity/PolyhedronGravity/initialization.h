//
//  initialization.h
//  
//
//  Created by Yang Yu on 10/25/17.
//
//

#ifndef _INITIALIZATION_H_
#define _INITIALIZATION_H_

#include <algorithm>
#include "asteroid.h"
#include "solar.h"
#include "particles.h"
#include "sysdyn.h"

using namespace std;

void LoadParas(char* parafile, string &polyfile, string &debrisfile, string &iniconfile, ODE_OPTION &odeopt, SOLAR &ss);
void LoadInicon(string iniconfile, DEBRIS &d, KEPORB &kep, PHASE &x, CLOUD &cl);

#endif /* */
