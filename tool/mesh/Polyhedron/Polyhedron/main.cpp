//
//  ESPI (Extracting Single Polyhedron Information)
//
//  Created by Yang Yu on 22 Nov 2014
//
//  improved by using 3d linear algebra library "mat3d" and "polyhedron" structure
//
//
#include <iostream>
#include "mat3d.h"
#include "polyhedron.h"
#include "mesh.h"
#include "volint.h"
#include "polynorm.h"
#include "exportation.h"
//

using namespace std;

//

int main(int argc,char *argv[])
{
    int i, j, k;
    double lenscale, vol0;
    Vector cntr0;
    Matrix tnsr0;
    POLYHEDRON p;
    
    setbuf(stdout,(char *)NULL);
    
    if (argc<=1)
    {
        (void) fprintf(stderr,"Usage: %s file [ file ... ]\n",argv[0]);
        exit(1);
    }
    
    for (i=1; i<argc; i++)
    {
        cout<<"~~~~~~~~~~~~~~~~~~~~~ Load Vertices & Faces ~~~~~~~~~~~~~~~~~~~~~"<<endl;
        
        VertFaceLoad(argv[i], p);
        
        cout<<"~~~~~~~~~~~~~~~ Calculate Volume Integrals ~~~~~~~~~~~~~~~"<<endl;
      
        VolumeIntegrals(p);
        
        cout<<"~~~~~~~~~~~~ Save Original Volume Information ~~~~~~~~~~~~"<<endl;

        vol0 = p.Volume;
        for (j=0; j<3; j++)
        {
            cntr0[j] = p.Centroid[j];
            for (k=0; k<3; k++) tnsr0[j][k] = p.InertiaTensor[j][k];
        }
        
        cout<<"~~~~~~~~~~~~~~~~~ Normalize the Polyhedron ~~~~~~~~~~~~~~~"<<endl;
        
        PolyNormalize(p, lenscale);
        
        cout<<"~~~~~~~~~~~~~~~ Calculate Volume Integrals ~~~~~~~~~~~~~~~"<<endl;
        
        VolumeIntegrals(p);
        
        cout<<"~~~~~~~~~~~~~~ Calculate the Mesh Parameters ~~~~~~~~~~~~~~"<<endl;
        
        MeshParaCalcu(p);
        
        cout<<"~~~~ Output the Mesh Parameters and Volume Information ~~~~"<<endl;

        PolyhedronOutput(argv[i], p);
        VolumeParaOutput(argv[i], vol0, cntr0, tnsr0, lenscale, p);
        
        cout<<"~~~~~~~~~~~~~~~~~~~~~~~ Completed ~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
 
    }

    return 0;
    
}
