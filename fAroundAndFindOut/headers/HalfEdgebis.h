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
    int halfEdge;      // index of the HalfEdge in the mesh's halfEdges array
};

struct HalfEdge {
    int opposite;  // index of the opposite HalfEdge in the mesh's halfEdges array
    int next;      // index of the next HalfEdge in the mesh's halfEdges array  
    int prev;      // index of the previous HalfEdge in the mesh's halfEdges array
    int vertex;    // index of the Vertex in the mesh's vertices array
    int face;      // index of the Face in the mesh's faces array
    int index;
    int valid;      // 1 if valid, 0 if removed
};

struct Face {
    int halfEdge;  
    int index;
    int valid; // 1 if valid, 0 if removed
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

void get_opposite(TriangularMesh* mesh, int he_index);

#endif // HALFEDGE_H
