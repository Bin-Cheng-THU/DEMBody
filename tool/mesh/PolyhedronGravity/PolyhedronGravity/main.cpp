//
//  main.cpp
//  
//
//  Created by Yang Yu on 11/24/14.
//
//

#include <iostream>
#include "body.h"
#include "kinematics.h"
#include "points.h"
//

using namespace std;
//
int main(int argc,char *argv[])
{
    POLYHEDRON p;
	POINTSDATA m;

	Vector r, f;
	double density;

    r[0] = -2000.0;
	r[1] = -2000.0;
	r[2] = 0.0;

    setbuf(stdout,(char *)NULL);
    
    if (argc<3)
    {
        (void) fprintf(stderr,"Usage: %s file\n",argv[0]);
        exit(1);
    }
    
	density = atof(argv[3]);
	density = ((int)(density*1000))/1000.0;

	LoadPolyhedron(argv[1], density, p);

	cout << "Complete Polyhedron loading..." << endl;

	LoadPoints(argv[2], m);

	cout << "Complete Points loading..." << endl;

	ConductPoints(m, p);

	WritePoints(argv[2], m);

	GravAttraction(p,r,f);

	cout << "Test point: " << f[0] << " " << f[1] << " " << f[2];
            
    return 0;
    
}


