#include "points.h"
#include <fstream>
#include <iomanip>
using namespace std;

void LoadPoints(string pointsfile, POINTSDATA &m)
{
    int i;
	const char* file;
    FILE * infile;

    file = pointsfile.c_str();
    infile = fopen(file, "r");
	if(!infile) {printf("fopen error!\n");exit;}

	fscanf(infile, "%d", &m.PointsNum);

	AllocPoints(m);

	for(i=0; i<m.PointsNum; i++) fscanf (infile,"%lf%lf%lf",&m.PointsPos[i][0],&m.PointsPos[i][1],&m.PointsPos[i][2]);
}

void AllocPoints(POINTSDATA &m)
{
    int i;
    
    m.PointsPos = new double*[m.PointsNum];
    for (i=0; i<m.PointsNum; i++) m.PointsPos[i] = new double[3];
    
    m.PointsForce = new double*[m.PointsNum];
    for (i=0; i<m.PointsNum; i++) m.PointsForce[i] = new double[3];
}

void WritePoints(string pointsfile, POINTSDATA &m)
{
#define WidthInt 12
#define WidthDouble 22
#define PrecDouble 12
#define N 5

	ofstream forcefile;
    int i,n,tn;
	char* file;
    FILE * infile;
	const char *forcename;
	char force[N] = {'f','o','r','c','e'};

	forcename = pointsfile.c_str();

    n = strlen(forcename);
    for (i=0; i<n; i++)
    {
        if (forcename[n-i-1] == '.') break;
    }
    tn = i;
    
    file = (char *) malloc((n-tn+N)*sizeof(char));
    
    for (i=0; i<n-tn; i++) file[i] = *(forcename+i);
    for (i=0; i<N; i++) file[n-tn+i] = force[i];
	file[n-tn+N] = 0;
    
    forcefile.open(file);

	forcefile<<setiosflags(ios::scientific)<<setprecision(PrecDouble);

    for (i=0; i<m.PointsNum; i++)
    {
        forcefile<<setw(WidthDouble)<<m.PointsForce[i][0]<<setw(WidthDouble)<<m.PointsForce[i][1]\
                <<setw(WidthDouble)<<m.PointsForce[i][2]<<endl;
    }
}

void ConductPoints(POINTSDATA &m, POLYHEDRON &p)
{
	for (int i=0; i<m.PointsNum; i++)
	{
		Vector r, f;
		r[0] = m.PointsPos[i][0];
		r[1] = m.PointsPos[i][1];
		r[2] = m.PointsPos[i][2];

		GravAttraction(p,r,f);
		
		m.PointsForce[i][0]	= f[0];
		m.PointsForce[i][1]	= f[1];
		m.PointsForce[i][2]	= f[2];
	}
}