//
//  Created by Yang Yu on 14-11-24
//
//

#ifndef _KINEMATICS_H_
#define _KINEMATICS_H_

#include "constant.h"
#include "gravity.h"
#include <iomanip>

typedef double Quaternion[4];

void QuatSet(Quaternion q, double a0, double a1, double a2, double a3);
void QuatCopy(const Quaternion q1, Quaternion q2);
void QuatNorm(Quaternion q);
void QuatMultiply(const Quaternion q1, const Quaternion q2, Quaternion q);
void Quat2DCM(const Quaternion q, Matrix m);
void Quat2IvDCM(const Quaternion q, Matrix m);
//void DCM2Quat(const Matrix m, Quaternion q);
//void Euler2Quat(const double a[3], Quaternion q);

#endif
