#ifndef HALFEDGE_H
#define HALFEDGE_H

#include <stddef.h>

typedef struct HalfEdge HalfEdge;
typedef struct Vertex Vertex;
typedef struct Face Face;

struct Vertex {
    double coord[2]; 
    int index;   // Index of the vertex
    HalfEdge *half_edge; // One of the half-edges originating from the vertex
};

struct Face {
    int index;
    HalfEdge *half_edge; // One of the half-edges bordering the face
};

struct HalfEdge {
    int index; // Index of the origin vertex
    HalfEdge *twin;   // the twin half-edge
    HalfEdge *next;   // the next half-edge in the face
    HalfEdge *prev;   // the previous half-edge in the face
    Face *face;   // face this half-edge belongs to
    Vertex *vertex; // Origin vertex of the half-edge
};

typedef struct {
    size_t vertex_count;
    Vertex *vertices;
    size_t face_count;
    Face *faces;
    size_t half_edge_count;
    HalfEdge *half_edges; // Array of half-edges
} TriangularMesh;

TriangularMesh *allocMesh(size_t vertex_count);
int initMesh(TriangularMesh *mesh, size_t vertex_count, size_t face_count, size_t half_edge_count);
void freeMesh(TriangularMesh *mesh);

void initVertex(Vertex *v, double x, double y, int index, HalfEdge *he);
void initFace(Face *f, int index, HalfEdge *he);
void initHalfEdge(HalfEdge *he, int index, HalfEdge *twin, HalfEdge *next, HalfEdge *prev, Face *face, Vertex *vertex);

#endif // HALFEDGE_H