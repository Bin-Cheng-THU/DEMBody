#include "mat3d.h"

void vectorCopy(const Vector u,Vector v)
{
	v[0] = u[0];
	v[1] = u[1];
	v[2] = u[2];
}

void vectorScale(const Vector u,const double s,Vector v)
{
	v[0] = u[0]*s;
	v[1] = u[1]*s;
	v[2] = u[2]*s;
}

void vectorAdd(const Vector v1,const Vector v2,Vector v)
{
	v[0] = v1[0] + v2[0];
	v[1] = v1[1] + v2[1];
	v[2] = v1[2] + v2[2];
}

void vectorSub(const Vector v1,const Vector v2,Vector v)
{
	v[0] = v1[0] - v2[0];
	v[1] = v1[1] - v2[1];
	v[2] = v1[2] - v2[2];
}

void vectorCross(const Vector v1,const Vector v2,Vector v)
{
    double a, b, c;
    
    a = v1[1]*v2[2] - v1[2]*v2[1];
    b = v1[2]*v2[0] - v1[0]*v2[2];
    c = v1[0]*v2[1] - v1[1]*v2[0];
    
    vectorSet(v, a, b, c);
}

void vectorMultiply(const Vector v1,const Vector v2,Vector v)
{
    v[0] = v1[0] * v2[0];
    v[1] = v1[1] * v2[1];
    v[2] = v1[2] * v2[2];
}

double vectorDot(const Vector v1,const Vector v2)
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

double vectorMagSq(const Vector v)
{
	return vectorDot(v,v);
}

double vectorMag(const Vector v)
{
	return sqrt(vectorMagSq(v));
	}

void vectorNorm(Vector v)
{
	double mag = vectorMag(v);
	assert(mag > 0.0);
	vectorScale(v,1.0/mag,v);
}

void vectorSet(Vector v,const double x,const double y,const double z)
{
	v[0] = x;
	v[1] = y;
	v[2] = z;
}

void vectorZero(Vector v)
{
	vectorSet(v,0.0,0.0,0.0);
}

void vectorTransform(const Matrix m,const Vector u,Vector v)
{
    double a, b, c;
    
    a = vectorDot(m[0],u);
    b = vectorDot(m[1],u);
    c = vectorDot(m[2],u);
    
    vectorSet(v, a, b, c);
}

void matrixCopy(const Matrix a,Matrix b)
{
	vectorCopy(a[0],b[0]);
	vectorCopy(a[1],b[1]);
	vectorCopy(a[2],b[2]);
}

void matrixZero(Matrix m)
{
    m[0][0] = 0.0;
    m[0][1] = 0.0;
    m[0][2] = 0.0;
    m[1][0] = 0.0;
    m[1][1] = 0.0;
    m[1][2] = 0.0;
    m[2][0] = 0.0;
    m[2][1] = 0.0;
    m[2][2] = 0.0;
}

void matrixIdentity(Matrix m)
{
	m[0][0] = 1.0;
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = 1.0;
	m[1][2] = 0.0;
	m[2][0] = 0.0;
	m[2][1] = 0.0;
	m[2][2] = 1.0;
}

void matrixDiagonal(const Vector v,Matrix m)
{
	m[0][0] = v[0];
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = v[1];
	m[1][2] = 0.0;
	m[2][0] = 0.0;
	m[2][1] = 0.0;
	m[2][2] = v[2];
}

void matrixScale(const Matrix a,const double s,Matrix b)
{
	vectorScale(a[0],s,b[0]);
	vectorScale(a[1],s,b[1]);
	vectorScale(a[2],s,b[2]);
}

void matrixSub(const Matrix a,const Matrix b,Matrix c)
{
    vectorSub(a[0],b[0],c[0]);
    vectorSub(a[1],b[1],c[1]);
    vectorSub(a[2],b[2],c[2]);
}

void matrixMultiply(const Matrix a,const Matrix b,Matrix c)
{
    Matrix d;
    
    d[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0];
	d[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1];
	d[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];
	d[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0];
	d[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1];
	d[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];
	d[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
	d[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
	d[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];
    
    matrixCopy(d, c);
}

void matrixTranspose(Matrix a)
{
    double tmp;
    
    tmp = a[0][1];
    a[0][1] = a[1][0];
    a[1][0] = tmp;
    
    tmp = a[0][2];
    a[0][2] = a[2][0];
    a[2][0] = tmp;
    
    tmp = a[1][2];
    a[1][2] = a[2][1];
    a[2][1] = tmp;
}

double matrixDet(const Matrix a)
{
    return a[0][0]*(a[1][1]*a[2][2]-a[2][1]*a[1][2]) \
          -a[0][1]*(a[1][0]*a[2][2]-a[2][0]*a[1][2]) \
          +a[0][2]*(a[1][0]*a[2][1]-a[2][0]*a[1][1]);
}

double matrixTrace(const Matrix a)
{
    return a[0][0]+a[1][1]+a[2][2];
}

/*
 matrixExtend and matrixContract implement the transformation between
 Matrix m (3x3 elements) and 1-dimensional array a (9 elements):
 
 m[0][0]  m[0][1]  m[0][2]     a[0]  a[3]  a[6]
 
 m[1][0]  m[1][1]  m[1][2]  =  a[1]  a[4]  a[7]
 
 m[2][0]  m[2][1]  m[2][2]     a[2]  a[5]  a[8]
 */

void matrixExtend(const Matrix m, double* a)
{
    a[0] = m[0][0];
    a[1] = m[1][0];
    a[2] = m[2][0];
    a[3] = m[0][1];
    a[4] = m[1][1];
    a[5] = m[2][1];
    a[6] = m[0][2];
    a[7] = m[1][2];
    a[8] = m[2][2];
}

void matrixContract(const double* a, Matrix m)
{
    m[0][0] = a[0];
    m[0][1] = a[3];
    m[0][2] = a[6];
    m[1][0] = a[1];
    m[1][1] = a[4];
    m[1][2] = a[7];
    m[2][0] = a[2];
    m[2][1] = a[5];
    m[2][2] = a[8];
}








