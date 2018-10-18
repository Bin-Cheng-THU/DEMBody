#include "kinematics.h"

void QuatSet(Quaternion q, double a0, double a1, double a2, double a3)
{
    q[0] = a0;
    q[1] = a1;
    q[2] = a2;
    q[3] = a3;
}

void QuatCopy(const Quaternion q1, Quaternion q2)
{
    q2[0] = q1[0];
    q2[1] = q1[1];
    q2[2] = q1[2];
    q2[3] = q1[3];
}

void QuatNorm(Quaternion q)
{
    double mag = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    //assert(mag > 0.0);
    for (int i=0; i<4; i++) q[i] = q[i] / mag;
}

void QuatMultiply(const Quaternion q1, const Quaternion q2, Quaternion q)
{
    double tmpd;
    Vector tmpv, v1, v2;
    
    vectorSet(v1, q1[1], q1[2], q1[3]);
    vectorSet(v2, q2[1], q2[2], q2[3]);
    
    tmpd = vectorDot(v1,v2);
    tmpd = q1[0]*q2[0] - tmpd;
    
    vectorCross(v1,v2,tmpv);
    vectorScale(v1,q2[0],v1);
    vectorScale(v2,q1[0],v2);
    
    vectorAdd(tmpv,v1,tmpv);
    vectorAdd(tmpv,v2,tmpv);
    
    QuatSet(q,tmpd,tmpv[0],tmpv[1],tmpv[2]);  
}

/*
 Quat2DCM: transfer quaternion to direction cosine matrix
 Quat2IvDCM: transfer quaternion to the inverse of direction cosine matrix
           DCM - direction cosine matrix
           the input q must be unit quaternion
*/

void Quat2DCM(const Quaternion q, Matrix m)
{
    double qq[10];
    //
    qq[0] = q[0] * q[0];
    qq[1] = q[1] * q[1];
    qq[2] = q[2] * q[2];
    qq[3] = q[3] * q[3];
    qq[4] = 2.0 * q[0] * q[1];
    qq[5] = 2.0 * q[0] * q[2];
    qq[6] = 2.0 * q[0] * q[3];
    qq[7] = 2.0 * q[1] * q[2];
    qq[8] = 2.0 * q[1] * q[3];
    qq[9] = 2.0 * q[2] * q[3];
    //
    m[0][0] = qq[0] + qq[1] - qq[2] - qq[3];
    m[0][1] = qq[7] + qq[6];
    m[0][2] = qq[8] - qq[5];
    m[1][0] = qq[7] - qq[6];
    m[1][1] = qq[0] - qq[1] + qq[2] - qq[3];
    m[1][2] = qq[9] + qq[4];
    m[2][0] = qq[8] + qq[5];
    m[2][1] = qq[9] - qq[4];
    m[2][2] = qq[0] - qq[1] - qq[2] + qq[3];
}

void Quat2IvDCM(const Quaternion q, Matrix m)
{
    double qq[10];
    //
    qq[0] = q[0] * q[0];
    qq[1] = q[1] * q[1];
    qq[2] = q[2] * q[2];
    qq[3] = q[3] * q[3];
    qq[4] = 2.0 * q[0] * q[1];
    qq[5] = 2.0 * q[0] * q[2];
    qq[6] = 2.0 * q[0] * q[3];
    qq[7] = 2.0 * q[1] * q[2];
    qq[8] = 2.0 * q[1] * q[3];
    qq[9] = 2.0 * q[2] * q[3];
    //
    m[0][0] = qq[0] + qq[1] - qq[2] - qq[3];
    m[0][1] = qq[7] - qq[6];
    m[0][2] = qq[8] + qq[5];
    m[1][0] = qq[7] + qq[6];
    m[1][1] = qq[0] - qq[1] + qq[2] - qq[3];
    m[1][2] = qq[9] - qq[4];
    m[2][0] = qq[8] - qq[5];
    m[2][1] = qq[9] + qq[4];
    m[2][2] = qq[0] - qq[1] - qq[2] + qq[3];
}

/*
 Euler2Quat: transfer Euler angles to unit quaternion
 Euler angle uses the definition of 3-1-3 rotation:
 input a[0] = phi is the precession about z-axis
       a[1] = theta is the nutation about x'-axis
       a[2] = varphi is the spin about z"-axis


void Euler2Quat(const double a[3], Quaternion q)
{
    double Cp, Sp, Ct, St, Cv, Sv;
    
    Cp = cos(a[0]/2.0);
    Sp = sin(a[0]/2.0);
    Ct = cos(a[1]/2.0);
    St = sin(a[1]/2.0);
    Cv = cos(a[2]/2.0);
    Sv = sin(a[2]/2.0);
    //
    q[0] = Cp*Ct*Cv - Sp*Ct*Sv;
    q[1] = Cp*St*Cv + Sp*St*Sv;
    q[2] = Sp*St*Cv - Cp*St*Sv;
    q[3] = Sp*Ct*Cv + Cp*Ct*Sv;
    //
    QuatNorm(q);
}
 */
//

/*



 */

/*
void testConserv(double y[26], POLYHEDRON &p, CLUSTER &c)
{
    int i;
    double E;
    PHASE x;
    Vector tmpV1, tmpV2;
    
    for (i=0; i<3; i++)
    {
        x.AlphaPos[i] = y[i];
        x.AlphaVel[i] = y[i+3];
        x.AlphaOrien[i] = y[i+6];
        x.AlphaAngVel[i] = y[i+10];
        x.BetaPos[i] = y[i+13];
        x.BetaVel[i] = y[i+16];
        x.BetaOrien[i] = y[i+19];
        x.BetaAngVel[i] = y[i+23];
    }
    x.AlphaOrien[3] = y[9];
    x.BetaOrien[3] = y[22];
    
    getMomentumMoment(p, c, x, tmpV1);
    getMomentum(p, c, x, tmpV2);
    getMechanicEnergy(p, c, x, E);
    cout<<setiosflags(ios::scientific)<<setprecision(22);
    cout<<tmpV1[0]<<" "<<tmpV1[1]<<" "<<tmpV1[2]<<" | "<<tmpV2[0]<<" "<<tmpV2[1]<<" "<<tmpV2[2]<<" | "<<E<<endl;
}
*/









