//
//  exportation.cpp
//  
//
//  Created by Yang Yu on 14-11-15.
//
//

#include "exportation.h"

using namespace std;

void PolyhedronOutput(char* filename,POLYHEDRON p)
{
#define WidthInt 12
#define WidthDouble 22
#define PrecDouble 12
#define N 10
    
    ofstream polyfile;
    int i, n, tn;
    char *polyname;
    char polyn[N] = {'p','o','l','y','h','e','d','r','o','n'};
    
    n = strlen(filename);
    for (i=0; i<n; i++)
    {
        if (filename[n-i-1] == '.') break;
    }
    tn = i;
    
    polyname = (char *) malloc((n-tn+N)*sizeof(char));
    
    for (i=0; i<n-tn; i++) polyname[i] = *(filename+i);
    for (i=0; i<N; i++) polyname[n-tn+i] = polyn[i];
	polyname[n-tn+N] = 0;
    
    polyfile.open(polyname);
    
    polyfile<<setiosflags(ios::scientific)<<setprecision(PrecDouble);
    
    polyfile<<setw(WidthDouble)<<p.Volume<<endl;
    polyfile<<setw(WidthDouble)<<p.Centroid[0]<<setw(WidthDouble)<<p.Centroid[1]<<setw(WidthDouble)<<p.Centroid[2]<<endl;
    for (i=0; i<3; i++)
    {
        polyfile<<setw(WidthDouble)<<p.InertiaTensor[i][0]<<setw(WidthDouble)<<p.InertiaTensor[i][1]\
        <<setw(WidthDouble)<<p.InertiaTensor[i][2]<<endl;
    }
    
    polyfile<<setw(WidthInt)<<p.NumVerts<<setw(WidthInt)<<p.NumFaces<<setw(WidthInt)<<p.NumEdges<<endl;
    
 //   polyfile<<setiosflags(ios::scientific)<<setprecision(PrecDouble);
    
    for (i=0; i<p.NumVerts; i++)
    {
        polyfile<<setw(WidthDouble)<<p.Vertices[i][0]<<setw(WidthDouble)<<p.Vertices[i][1]\
                <<setw(WidthDouble)<<p.Vertices[i][2]<<endl;
    }
    
    for (i=0; i<p.NumFaces; i++)
    {
        polyfile<<setw(WidthInt)<<p.Faces[i][0]<<setw(WidthInt)<<p.Faces[i][1]<<setw(WidthInt)<<p.Faces[i][2]\
                <<setw(WidthDouble)<<p.FaceNormVecs[i][0]<<setw(WidthDouble)<<p.FaceNormVecs[i][1]\
                <<setw(WidthDouble)<<p.FaceNormVecs[i][2]<<endl;
    }

    for (i=0; i<p.NumEdges; i++)
    {
        polyfile<<setw(WidthInt)<<p.Edges[i][0]<<setw(WidthInt)<<p.Edges[i][1]<<setw(WidthInt)<<p.Edges[i][2]\
                <<setw(WidthInt)<<p.Edges[i][3]<<setw(WidthDouble)<<p.EdgeLens[i]<<setw(WidthDouble)\
                <<p.EdgeNormVecs[i][0]<<setw(WidthDouble)<<p.EdgeNormVecs[i][1]<<setw(WidthDouble)\
                <<p.EdgeNormVecs[i][2]<<setw(WidthDouble)<<p.EdgeNormVecs[i][3]<<setw(WidthDouble)\
                <<p.EdgeNormVecs[i][4]<<setw(WidthDouble)<<p.EdgeNormVecs[i][5]<<endl;
    }
    
    polyfile.close();

#undef WidthInt
#undef WidthDouble
#undef PrecDouble
#undef N
}

void VolumeParaOutput(char* filename, double vol0, Vector cntr0, Matrix tnsr0, double lenscale, POLYHEDRON p)
{
#define WidthDouble 21
#define PrecDouble 12
#define N 6
    
    ofstream volfile;
    int i,n,tn;
    char *volname;
    char voln[N] = {'v','o','l','u','m','e'};
 
    n = strlen(filename);
    for (i=0; i<n; i++)
    {
        if (filename[n-i-1] == '.') break;
    }
    tn = i;
    volname = (char *) malloc((n-tn+N)*sizeof(char));
    
    for (i=0; i<n-tn; i++) volname[i] = *(filename+i);
    for (i=0; i<N; i++) volname[n-tn+i] = voln[i];
	volname[n-tn+N] = 0;
    
    volfile.open(volname);
    
    volfile<<setiosflags(ios::scientific)<<setprecision(PrecDouble);
    
    volfile<<"The Volume Information in Original Frame (where the initial vertice coordinates are expressed): "<<endl;
    volfile<<endl;
    volfile<<"   1. The volume V ="<<setw(WidthDouble)<<vol0<<endl;
    volfile<<endl;
    volfile<<"   2. The centroid (xc, yc, zc) = ("<<setw(WidthDouble)<<cntr0[0]<<setw(WidthDouble)<<cntr0[1]<<setw(WidthDouble)<<cntr0[2]<<"   )"<<endl;
    volfile<<endl;
    volfile<<"   3. The inertia tensor matrix (density = 1)"<<endl;
    volfile<<endl;
    volfile<<"                [Ixx  Ixy  Ixz]   ["<<setw(WidthDouble)<<tnsr0[0][0]<<setw(WidthDouble)<<tnsr0[0][1]<<setw(WidthDouble)<<tnsr0[0][2]<<"   ]"<<endl;
    volfile<<endl;
    volfile<<"                [Ixy  Iyy  Iyz] = ["<<setw(WidthDouble)<<tnsr0[1][0]<<setw(WidthDouble)<<tnsr0[1][1]<<setw(WidthDouble)<<tnsr0[1][2]<<"   ]"<<endl;
    volfile<<endl;
    volfile<<"                [Ixz  Iyz  Izz]   ["<<setw(WidthDouble)<<tnsr0[2][0]<<setw(WidthDouble)<<tnsr0[2][1]<<setw(WidthDouble)<<tnsr0[2][2]<<"   ]"<<endl;
    volfile<<endl;
    volfile<<endl;
    volfile<<"The Volume Information in Normalized Frame (origin at the mass center, scaled by L, and x-, y-, z-axes are minimum, medium, maximum principle inertial axes): "<<endl;
    volfile<<endl;
    volfile<<"   1. The length scale used for normalization L ="<<lenscale<<endl;
    volfile<<endl;
    volfile<<"   2. The volume V ="<<setw(WidthDouble)<<p.Volume<<endl;
    volfile<<endl;
    volfile<<"   3. The centroid (xc, yc, zc) = ("<<setw(WidthDouble)<<p.Centroid[0]<<setw(WidthDouble)<<p.Centroid[1]<<setw(WidthDouble)<<p.Centroid[2]<<"   )"<<endl;
    volfile<<endl;
    volfile<<"   4. The inertia tensor matrix (density = 1)"<<endl;
    volfile<<endl;
    volfile<<"                [Ixx  Ixy  Ixz]   ["<<setw(WidthDouble)<<p.InertiaTensor[0][0]<<setw(WidthDouble)<<p.InertiaTensor[0][1]<<setw(WidthDouble)<<p.InertiaTensor[0][2]<<"   ]"<<endl;
    volfile<<endl;
    volfile<<"                [Ixy  Iyy  Iyz] = ["<<setw(WidthDouble)<<p.InertiaTensor[1][0]<<setw(WidthDouble)<<p.InertiaTensor[1][1]<<setw(WidthDouble)<<p.InertiaTensor[1][2]<<"   ]"<<endl;
    volfile<<endl;
    volfile<<"                [Ixz  Iyz  Izz]   ["<<setw(WidthDouble)<<p.InertiaTensor[2][0]<<setw(WidthDouble)<<p.InertiaTensor[2][1]<<setw(WidthDouble)<<p.InertiaTensor[2][2]<<"   ]"<<endl;
    volfile.close();

#undef PrecDouble
#undef WidthDouble
#undef N
}








