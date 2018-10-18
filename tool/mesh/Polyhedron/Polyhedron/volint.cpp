//
//  volume.cpp
//  
//
//  Created by Yang Yu on 14-11-14.
//
//

#include "volint.h"

void VolumeIntegrals(POLYHEDRON &p)
{
    int i, j, a, b, c;
    double nx, ny, nz, tmp;
    double *Wf, **Nf;
    double F[12], wf, nf[3];
    double VolIntZeroth, VolIntFirst[3], VolIntSecond[6];
    
    VolIntZeroth = 0.0;
    for(i=0; i<3; i++) VolIntFirst[i] = 0.0;
    for(i=0; i<6; i++) VolIntSecond[i] = 0.0;
    
    CompWfNf(p, Wf, Nf);
    
    for (i=0; i<p.NumFaces; i++)
    {
        wf = Wf[i];
        for(j=0; j<3; j++) nf[j] = Nf[i][j];

        nx = fabs(nf[0]);
        ny = fabs(nf[1]);
        nz = fabs(nf[2]);
        
        if ((nx>ny) && (nx>nz)) c = 0;
        else if (ny>nz) c = 1;
        else c = 2;
        
        a = (c + 1) % 3;
        b = (a + 1) % 3;
        
        FaceIntegrals(p, wf, nf, i, a, b, c, F);
        
        if (a==0) tmp = F[0];
        else if (b==0) tmp = F[1];
        else tmp = F[2];
            
        VolIntZeroth = VolIntZeroth + Nf[i][0] * tmp;
        
        VolIntFirst[a] = VolIntFirst[a] + Nf[i][a] * F[3];
        VolIntFirst[b] = VolIntFirst[b] + Nf[i][b] * F[4];
        VolIntFirst[c] = VolIntFirst[c] + Nf[i][c] * F[5];
        
        VolIntSecond[a] = VolIntSecond[a] + Nf[i][a] * F[6];
        VolIntSecond[b] = VolIntSecond[b] + Nf[i][b] * F[7];
        VolIntSecond[c] = VolIntSecond[c] + Nf[i][c] * F[8];
        VolIntSecond[a+3] = VolIntSecond[a+3] + Nf[i][a] * F[9];
        VolIntSecond[b+3] = VolIntSecond[b+3] + Nf[i][b] * F[10];
        VolIntSecond[c+3] = VolIntSecond[c+3] + Nf[i][c] * F[11];
    }
    
    for(j=0; j<3; j++)
    {
        VolIntFirst[j] = VolIntFirst[j] / 2.0;
        VolIntSecond[j] = VolIntSecond[j] / 3.0;
        VolIntSecond[j+3] = VolIntSecond[j+3] / 2.0;
    }
    
    VolumeInfoCalcu(VolIntZeroth, VolIntFirst, VolIntSecond, p);
    
}

//
//  Compute the Wf and Nf
//

void CompWfNf(POLYHEDRON &p, double *&Wf, double **&Nf)
{
    int i, j;
    Vector ra, rb, rc, rab, rbc, tmpv;
    
    Wf = (double *) malloc(p.NumFaces*sizeof(double));
    Nf = (double **) malloc(p.NumFaces*sizeof(double *));
    for (i=0; i<p.NumFaces; i++) Nf[i] = (double *) malloc(3*sizeof(double));
    
    for(i=0; i<p.NumFaces; i++)
    {
        vectorSet(ra, p.Vertices[p.Faces[i][0]][0], p.Vertices[p.Faces[i][0]][1], p.Vertices[p.Faces[i][0]][2]);
        vectorSet(rb, p.Vertices[p.Faces[i][1]][0], p.Vertices[p.Faces[i][1]][1], p.Vertices[p.Faces[i][1]][2]);
        vectorSet(rc, p.Vertices[p.Faces[i][2]][0], p.Vertices[p.Faces[i][2]][1], p.Vertices[p.Faces[i][2]][2]);
        
        vectorSub(ra,rb,rab);
        vectorSub(rb,rc,rbc);
        
        vectorCross(rab, rbc, tmpv);
        vectorNorm(tmpv);
        
        Wf[i] = - vectorDot(ra, tmpv);
        
        for(j=0; j<3; j++) Nf[i][j] = tmpv[j];
    }
    
}

//
//  Compute the face integrals
//

void FaceIntegrals(POLYHEDRON &p, double wf, double nf[3], int id, int a, int b, int c, double F[12])
{
    double P[10], K[4];
    
    ProjectionIntegral(p, a, b, id, P);
    
    K[0] = 1.0 / nf[c];
    K[1] = K[0] * K[0];
    K[2] = K[1] * K[0];
    K[3] = K[2] * K[0];
    
    F[0] = K[0] * P[1];
    F[1] = K[0] * P[2];
    F[2] = - K[1] * (nf[a] * P[1] + nf[b] * P[2] + wf * P[0]);
    
    F[3] = K[0] * P[3];
    F[4] = K[0] * P[5];
    F[5] = K[2] * (nf[a] * nf[a] * P[3] + 2.0 * nf[a] * nf[b] * P[4] + \
                   nf[b] * nf[b] * P[5] + wf * (2.0 * (nf[a] * P[1] + \
                   nf[b] * P[2]) + wf * P[0]));
    
    F[6] = K[0] * P[6];
    F[7] = K[0] * P[9];
    F[8] = - K[3] * (nf[a] * nf[a] * nf[a] * P[6] + 3.0 * nf[a] * nf[a] * nf[b] * P[7] + \
                     3.0 * nf[a] * nf[b] * nf[b] * P[8] + nf[b] * nf[b] * nf[b] * P[9] + \
                     3.0 * wf * (nf[a] * nf[a] * P[3] + 2.0 * nf[a] * nf[b] * P[4] + \
                     nf[b] * nf[b] * P[5]) + wf * wf * (3.0 * (nf[a] * P[1] + nf[b] * P[2]) + \
                     wf * P[0]));
                     
    F[9] = K[0] * P[7];
    F[10] = - K[1] * (nf[a] * P[8] + nf[b] * P[9] + wf * P[5]);
    F[11] = K[2] * (nf[a] * nf[a] * P[6] + 2.0 * nf[a] * nf[b] * P[7] + nf[b] * nf[b] * P[8] + \
                    wf * (2.0 * (nf[a] * P[3] + nf[b] * P[4]) + wf * P[1]));
                     
}

//
//  Compute the projection integral
//

void ProjectionIntegral(POLYHEDRON &p, int a, int b, int id, double P[10])
{
    int i;
    double A0, A1, DA;
    double B0, B1, DB;
    double A02, A03, A04, B02, B03, B04;
    double A12, A13, B12, B13;
    double C1, CA, CAA, CAAA, CB, CBB, CBBB;
    double CAB, KAB, CAAB, KAAB, CABB, KABB;
    
    for(i=0; i<10; i++) P[i] = 0.0;
    
    for (i=0; i<3; i++)
    {
        A0 = p.Vertices[p.Faces[id][i]][a];
        B0 = p.Vertices[p.Faces[id][i]][b];
        A1 = p.Vertices[p.Faces[id][(i+1) % 3]][a];
        B1 = p.Vertices[p.Faces[id][(i+1) % 3]][b];
        
        DA = A1 - A0;
        DB = B1 - B0;
        
        A02 = A0 * A0;
        A03 = A02 * A0;
        A04 = A03 * A0;
        
        B02 = B0 * B0;
        B03 = B02 * B0;
        B04 = B03 * B0;
        
        A12 = A1 * A1;
        A13 = A12 * A1;
        B12 = B1 * B1;
        B13 = B12 * B1;
        
        C1 = A1 + A0;
        CA = A1 * C1 + A02;
        CAA = A1 * CA + A03;
        CAAA = A1 * CAA + A04;
        CB = B1 * (B1 + B0) + B02;
        CBB = B1 * CB + B03;
        CBBB = B1 * CBB + B04;
        CAB = 3.0 * A12 + 2.0 * A1 * A0 + A02;
        KAB = A12 + 2.0 * A1 * A0 + 3.0 * A02;
        CAAB = A0 * CAB + 4.0 * A13;
        KAAB = A1 * KAB + 4.0 * A03;
        CABB = 4.0 * B13 + 3.0 * B12 * B0 + 2.0 * B1 * B02 + B03;
        KABB = B13 + 2.0 * B12 * B0 + 3.0 * B1 * B02 + 4.0 * B03;
        
        P[0] = P[0] + DB * C1;
        P[1] = P[1] + DB * CA;
        P[3] = P[3] + DB * CAA;
        P[6] = P[6] + DB * CAAA;
        P[2] = P[2]+ DA * CB;
        P[5] = P[5]+ DA * CBB;
        P[9] = P[9]+ DA * CBBB;
        P[4] = P[4]+ DB * (B1 * CAB + B0 * KAB);
        P[7] = P[7]+ DB * (B1 * CAAB + B0 * KAAB);
        P[8] = P[8]+ DA * (A1 * CABB + A0 * KABB);
    }
    
    P[0] = P[0] / 2.0;
    P[1] = P[1] / 6.0;
    P[3] = P[3] / 12.0;
    P[6] = P[6] / 20.0;
    P[2] = - P[2] / 6.0;
    P[5] = - P[5] / 12.0;
    P[9] = - P[9] / 20.0;
    P[4] = P[4] / 24.0;
    P[7] = P[7] / 60.0;
    P[8] = - P[8] / 60.0;
    
}

//
// get the volume, centroid and normalized inertia tensor (density = 1.0)
//

void VolumeInfoCalcu(double VolIntZeroth, double VolIntFirst[3], double VolIntSecond[6], POLYHEDRON &p)
{
    p.Volume = VolIntZeroth;

    vectorSet(p.Centroid, VolIntFirst[0]/VolIntZeroth, \
                          VolIntFirst[1]/VolIntZeroth, \
                          VolIntFirst[2]/VolIntZeroth);
    
    p.InertiaTensor[0][0] = VolIntSecond[1] + VolIntSecond[2];
    p.InertiaTensor[1][1] = VolIntSecond[2] + VolIntSecond[0];
    p.InertiaTensor[2][2] = VolIntSecond[0] + VolIntSecond[1];
    p.InertiaTensor[0][1] = - VolIntSecond[3];
    p.InertiaTensor[1][0] = - VolIntSecond[3];
    p.InertiaTensor[1][2] = - VolIntSecond[4];
    p.InertiaTensor[2][1] = - VolIntSecond[4];
    p.InertiaTensor[0][2] = - VolIntSecond[5];
    p.InertiaTensor[2][0] = - VolIntSecond[5];
    
}
                     


