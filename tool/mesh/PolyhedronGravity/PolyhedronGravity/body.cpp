#include "body.h"

void LoadPolyhedron(string polyfile, double density, POLYHEDRON &p)
{
    int i, j;
    double tmpd, tmpp;
    const char* file;
    FILE * infile;
	double dump1,dump2,dump3;
    
    file = polyfile.c_str();
    
    infile = fopen(file, "r+");
	if(!infile) printf("fopen error!\n");

    fscanf(infile, "%lf", &p.Volume);
    fscanf(infile, "%lf%lf%lf", &dump1,&dump2,&dump3);

    fscanf(infile,"%lf%lf%lf",&p.Inertia[0],&dump1,&dump2);
	fscanf(infile,"%lf%lf%lf",&dump1,&p.Inertia[1],&dump2);
	fscanf(infile,"%lf%lf%lf",&dump1,&dump2,&p.Inertia[2]);
    
    fscanf(infile, "%d%d%d", &p.NumVerts, &p.NumFaces, &p.NumEdges);
    
    PolyhedronAlloc(p);
    
    for(i=0; i<p.NumVerts; i++) fscanf (infile,"%lf%lf%lf",&p.Vertices[i][0],&p.Vertices[i][1],&p.Vertices[i][2]);

    for(i=0; i<p.NumFaces; i++) fscanf (infile,"%d%d%d%lf%lf%lf",&p.Faces[i][0],&p.Faces[i][1],&p.Faces[i][2], \
                                        &p.FaceNormVecs[i][0],&p.FaceNormVecs[i][1],&p.FaceNormVecs[i][2]);
    
    for(i=0; i<p.NumEdges; i++) fscanf (infile,"%d%d%d%d%lf%lf%lf%lf%lf%lf%lf",&p.Edges[i][0],&p.Edges[i][1], \
                                        &p.Edges[i][2],&p.Edges[i][3],&p.EdgeLens[i],&p.EdgeNormVecs[i][0], \
                                        &p.EdgeNormVecs[i][1],&p.EdgeNormVecs[i][2],&p.EdgeNormVecs[i][3], \
                                        &p.EdgeNormVecs[i][4],&p.EdgeNormVecs[i][5]);
    
    fclose (infile);
	p.Density = density;
    p.Mass = p.Volume * p.Density;
    
    for (i=0; i<3; i++)
    {
        tmpd = 0.0;
        for (j=0; j<p.NumVerts; j++)
        {
            tmpp = fabs(p.Vertices[j][i]);
            if (tmpp>tmpd) tmpd = tmpp;
        }
        p.Ellipsoid[i] = tmpd;
    }
    
}

void PolyhedronAlloc(POLYHEDRON &p)
{
    int i;
    
    p.Vertices = new double*[p.NumVerts];
    for (i=0; i<p.NumVerts; i++) p.Vertices[i] = new double[3];
    
    p.Faces = new int*[p.NumFaces];
    p.FaceNormVecs = new double*[p.NumFaces];
    for (i=0; i<p.NumFaces; i++)
    {
        p.Faces[i] = new int[3];
        p.FaceNormVecs[i] = new double[3];
    }
    
    p.Edges = new int*[p.NumEdges];
    p.EdgeNormVecs = new double*[p.NumEdges];
    for (i=0; i<p.NumEdges; i++)
    {
        p.Edges[i] = new int[4];
        p.EdgeNormVecs[i] = new double[6];
    }
    
    p.EdgeLens = new double[p.NumEdges];
}




