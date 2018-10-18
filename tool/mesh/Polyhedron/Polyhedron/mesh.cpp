//
//
//  Created by Yang Yu on 11/14/14.
//
//

#include "mesh.h"

//
/*
 subroutine VertFaceLoad: Load the mesh data of the polyhedron (the size of
 the variant arranged here).
 */
//

using namespace std;

void VertFaceLoad(char* filename, POLYHEDRON &p)
{
    int i, n1, n2;
    FILE * infile;
    
    infile = fopen(filename, "r");
    fscanf(infile, "%d%d", &p.NumVerts, &p.NumFaces);
   
    MeshParaAlloc(p);
    
    for(i=0; i<p.NumVerts; i++)
    {
        fscanf (infile,"%lf%lf%lf",&p.Vertices[i][0],&p.Vertices[i][1],&p.Vertices[i][2]);
    }
    for(i=0; i<p.NumFaces; i++)
    {
        fscanf (infile,"%d%d%d",&p.Faces[i][0],&p.Faces[i][1],&p.Faces[i][2]);
    }
    fclose (infile);
}

//
/*
 subroutine MeshParaAlloc: alloc the size of mesh parameters
 */
//

void MeshParaAlloc(POLYHEDRON &p)
{
    int i;
    
    p.NumEdges = p.NumFaces * 3 / 2;
    
    p.Vertices = (double **) malloc(p.NumVerts*sizeof(double *));
    for (i=0; i<p.NumVerts; i++)
    {
        p.Vertices[i] = (double *) malloc(3*sizeof(double));
    }
    
    p.Faces = (int **) malloc(p.NumFaces*sizeof(int *));
    p.FaceNormVecs = (double **) malloc(p.NumFaces*sizeof(double *));
    for (i=0; i<p.NumFaces; i++)
    {
        p.Faces[i] = (int *) malloc(3*sizeof(int));
        p.FaceNormVecs[i] = (double *) malloc(3*sizeof(double));
    }
    
    p.Edges = (int **) malloc(p.NumEdges*sizeof(int *));
    p.EdgeNormVecs = (double **) malloc(p.NumEdges*sizeof(double *));
    for (i=0; i<p.NumEdges; i++)
    {
        p.Edges[i] = (int *) malloc(4*sizeof(int));
        p.EdgeNormVecs[i] = (double *) malloc(6*sizeof(double));
    }
    
    p.EdgeLens = new double[p.NumEdges];
    //(double *) malloc(p.NumEdges*sizeof(double));
}


//
/*
 subroutine MeshParaCalcu: to create the parameters of the polyhedral mesh
 */
//

void MeshParaCalcu(POLYHEDRON &p)
{
    bool flag;
    
    int i, j, k, count, fn;
    int **FaceEdgeMap;
	int m,n;
	m = 3*p.NumFaces;
	n = 2;
	FaceEdgeMap = (int**)malloc(sizeof(int*) *m);
	for(i=0;i<m;i++)
		FaceEdgeMap[i] = (int*)malloc(sizeof(int)*n);
    
    double lam,tmp,nap;
    Vector tmpv, tmpu, ra, rb, rc, rp, ap, rab, rac, rbc;
    
    //  Generate edge parameter Edges, format: VERTEX 1, VERTEX 2, FACE 1, FACE 2
    
    for(i=0; i<p.NumFaces; i++)
    {
        for(j=0; j<3; j++)
        {
            FaceEdgeMap[3*i+j][0] = p.Faces[i][j];
            FaceEdgeMap[3*i+j][1] = p.Faces[i][(j+1)%3];
        }
    }

    count = -1;
    
    for(i=0; i<3*p.NumFaces-1; i++)
    {
        for(j=i+1; j<3*p.NumFaces; j++)
        {
            flag = ((FaceEdgeMap[j][0] == FaceEdgeMap[i][0]) && (FaceEdgeMap[j][1] == FaceEdgeMap[i][1])) || \
                   ((FaceEdgeMap[j][0] == FaceEdgeMap[i][1]) && (FaceEdgeMap[j][1] == FaceEdgeMap[i][0]));
            
            if (flag)
            {
                count = count + 1;
                
                for(k=0; k<2; k++)
                {
                    p.Edges[count][k] = FaceEdgeMap[i][k];
                }
                
                p.Edges[count][2] = i / 3;
                p.Edges[count][3] = j / 3;
            }
        }
    }
    
    //  Generate edge parameters EdgeLength, EdgeNormVecs
    
    for(i=0; i<p.NumEdges; i++)
    {
        vectorSet(rb, p.Vertices[p.Edges[i][0]][0], p.Vertices[p.Edges[i][0]][1], p.Vertices[p.Edges[i][0]][2]);
        vectorSet(rc, p.Vertices[p.Edges[i][1]][0], p.Vertices[p.Edges[i][1]][1], p.Vertices[p.Edges[i][1]][2]);
        vectorSub(rb, rc, tmpv);
        
        p.EdgeLens[i] = vectorMag(tmpv);
        
        for(k=0; k<2; k++)
        {
            fn = p.Faces[p.Edges[i][k+2]][0] + p.Faces[p.Edges[i][k+2]][1] + p.Faces[p.Edges[i][k+2]][2] - \
                 p.Edges[i][0] - p.Edges[i][1];

            vectorSet(ra, p.Vertices[fn][0], p.Vertices[fn][1], p.Vertices[fn][2]);
            vectorSub(ra,rc,rac);
            vectorSub(rb,rc,rbc);
            
            lam = vectorDot(rac,rbc) / vectorMagSq(rbc);
            
            tmp = 1.0 - lam;
            vectorScale(rb,lam,tmpv);
            vectorScale(rc,tmp,tmpu);
            
            vectorAdd(tmpv,tmpu,rp);
            vectorSub(rp,ra,ap);
            vectorNorm(ap);
            
            for(j=0; j<3; j++) p.EdgeNormVecs[i][3*k+j] = ap[j];
        }
        
    }

    //  Generate face parameters FaceNormVecs
    
    for(i=0; i<p.NumFaces; i++)
    {
        vectorSet(ra, p.Vertices[p.Faces[i][0]][0], p.Vertices[p.Faces[i][0]][1], p.Vertices[p.Faces[i][0]][2]);
        vectorSet(rb, p.Vertices[p.Faces[i][1]][0], p.Vertices[p.Faces[i][1]][1], p.Vertices[p.Faces[i][1]][2]);
        vectorSet(rc, p.Vertices[p.Faces[i][2]][0], p.Vertices[p.Faces[i][2]][1], p.Vertices[p.Faces[i][2]][2]);
        
        vectorSub(ra,rb,rab);
        vectorSub(rb,rc,rbc);
        
        vectorCross(rab, rbc, tmpv);
        vectorNorm(tmpv);
        
        for(j=0; j<3; j++) p.FaceNormVecs[i][j] = tmpv[j];
    }
    
	for(i = 0; i < m; i ++)
		free(FaceEdgeMap[i]);
	free(FaceEdgeMap);
}

//















