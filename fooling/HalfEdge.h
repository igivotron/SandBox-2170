#ifndef HALFEDGE_H
#define HALFEDGE_H

#include <stdlib.h>
#include <stdio.h>

typedef struct HalfEdge HalfEdge;
typedef struct Face Face;
typedef struct Vertex Vertex;
typedef struct TriangularMesh TriangularMesh;


struct Vertex {
    double x, y, z;
    int index;
    HalfEdge* halfEdge;  
};

struct HalfEdge {
    HalfEdge* opposite;  
    HalfEdge* next;      
    HalfEdge* prev;   
    Vertex* vertex;      
    Face* face;    
    int index;
};

struct Face {
    HalfEdge* halfEdge;  
    int index;
};

struct TriangularMesh {
    Vertex* vertices;
    HalfEdge* halfEdges;
    Face* faces;
    int numVertices;
    int numHalfEdges;
    int numFaces;
};

TriangularMesh* createTriangularMesh(double** vertices, int numVertices, int** faces, int numFaces);
int freeTriangularMesh(TriangularMesh* mesh);
int saveMeshToOBJ(TriangularMesh* mesh, const char* filename);

#endif // HALFEDGE_H
