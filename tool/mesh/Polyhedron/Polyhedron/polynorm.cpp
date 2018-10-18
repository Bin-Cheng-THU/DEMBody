//
//  frames.cpp
//  
//
//  Created by Yang Yu on 14-11-15.
//
//

#include "polynorm.h"

void PolyNormalize(POLYHEDRON &p, double &lenscale)
{
    int i, dim, it_max, it_num, rot_num;
    double DVC[3], tmpL1[9], tmpL2[9];
    Matrix VMX, tmpM1, tmpM2, tmpM3, JMX, tranMX;
    double tmp, error_frobenius;
    
    //lenscale = 0.001;//pow(0.75 * p.Volume / PI, 1.0/3.0);
    lenscale = 1.0;
    tmp = vectorMagSq(p.Centroid);
    
    matrixZero(tmpM1);
    matrixZero(tmpM2);
    
    for (i=0; i<3; i++)
    {
        tmpM1[i][i] = tmp;
        tmpM2[i][i] = p.Centroid[i] * p.Centroid[i];
    }
    
    tmpM2[0][1] = p.Centroid[0] * p.Centroid[1];
    tmpM2[1][0] = tmpM2[0][1];
    tmpM2[0][2] = p.Centroid[0] * p.Centroid[2];
    tmpM2[2][0] = tmpM2[0][2];
    tmpM2[1][2] = p.Centroid[1] * p.Centroid[2];
    tmpM2[2][1] = tmpM2[1][2];
    
    matrixSub(tmpM1,tmpM2,tmpM3);
    matrixScale(tmpM3,p.Volume,tmpM1);
    matrixSub(p.InertiaTensor,tmpM1,JMX);
    
    dim = 3;
    it_max = 200;
    
    matrixExtend(JMX,tmpL1);
    
    jacobi_eigenvalue(dim, tmpL1, it_max, tmpL2, DVC, &it_num, &rot_num );
    
    matrixContract(tmpL2,VMX);

    printf( "\n" );
    printf( "  Eigen decomposing: " );
    printf( "  Number of iterations = %d\n", it_num );
    printf( "  Number of rotations  = %d\n", rot_num );
    
    printf( "  Eigenvalues D:" );
    r8vec_print(dim, DVC);
    
    printf( "  Eigenvector matrix V:" );
    r8mat_print(dim, dim, tmpL2);

    error_frobenius = r8mat_is_eigen_right(dim, dim, tmpL1, tmpL2, DVC);
    printf( "\n" );
    printf( "  Frobenius norm error in eigensystem A*V-D*V = %g\n", error_frobenius);
    
    matrixCopy(VMX, tranMX);
    InertSort(DVC, tranMX);
    
    matrixTranspose(tranMX);
    
    normVertices(p, 1.0/lenscale, tranMX);

}

void InertSort(double d[3], Matrix m)
{
    int i;
    double tmp;
    Vector tmpV;
    
    
    if ((d[0]>d[1])&&(d[0]>d[2]))
    {
        for (i=0; i<3; i++)
        {
            tmpV[i] = m[i][2];
            m[i][2] = m[i][0];
            m[i][0] = tmpV[i];
        }
        tmp = d[2];
        d[2] = d[0];
        d[0] = tmp;
        
        if (d[1]<d[0])
        {
            for (i=0; i<3; i++)
            {
                tmpV[i] = m[i][1];
                m[i][1] = m[i][0];
                m[i][0] = tmpV[i];
            }
            tmp = d[1];
            d[1] = d[0];
            d[0] = tmp;
        }
    }
    else if (d[1]>d[2])
    {
        for (i=0; i<3; i++)
        {
            tmpV[i] = m[i][2];
            m[i][2] = m[i][1];
            m[i][1] = tmpV[i];
        }
        tmp = d[2];
        d[2] = d[1];
        d[1] = tmp;
        
        if (d[1]<d[0])
        {
            for (i=0; i<3; i++)
            {
                tmpV[i] = m[i][1];
                m[i][1] = m[i][0];
                m[i][0] = tmpV[i];
            }
            tmp = d[1];
            d[1] = d[0];
            d[0] = tmp;
        }
    }
    else
    {
        if (d[1]<d[0])
        {
            for (i=0; i<3; i++)
            {
                tmpV[i] = m[i][1];
                m[i][1] = m[i][0];
                m[i][0] = tmpV[i];
            }
            tmp = d[1];
            d[1] = d[0];
            d[0] = tmp;
        }
    }
    
    tmp = matrixDet(m);
    
    if (tmp<0.0)
    {
        for (i=0; i<3; i++) m[i][0] = - m[i][0];
    }
}


void normVertices(POLYHEDRON &p, double a, Matrix m)
{
    int i, j;
    Vector tmpV, tmpU;
    
    for (i=0; i<p.NumVerts; i++)
    {
        vectorSet(tmpV, p.Vertices[i][0], p.Vertices[i][1], p.Vertices[i][2]);
        
        vectorSub(tmpV, p.Centroid, tmpU);
        
        vectorTransform(m, tmpU, tmpV);
        
        vectorScale(tmpV, a, tmpU);
        
        for (j=0; j<3; j++) p.Vertices[i][j] = tmpU[j];
    }
}

