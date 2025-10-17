#include <stdlib.h>
#include "HalfEdge.h"


TriangularMesh* allocMesh(size_t vertex_count) {
    TriangularMesh* mesh = (TriangularMesh*)malloc(sizeof(TriangularMesh));
    if (!mesh) return NULL;

    mesh->half_edges = (HalfEdge*)calloc(6*vertex_count+24, sizeof(HalfEdge));
    mesh->vertices = (Vertex*)calloc(vertex_count+4, sizeof(Vertex));
    mesh->faces = (Face*)calloc(2*vertex_count+8, sizeof(Face));

    // if (!mesh->half_edges || !mesh->vertices || !mesh->faces) {
    //     free(mesh->half_edges);
    //     free(mesh->vertices);
    //     free(mesh->faces);
    //     free(mesh);
    //     return NULL;
    // }

    mesh->vertex_count = 0;
    mesh->half_edge_count = 0;
    mesh->face_count = 0;

    return mesh;
}

void freeMesh(TriangularMesh* mesh) {
    if (!mesh) return;

    free(mesh->half_edges);
    free(mesh->vertices);
    free(mesh->faces);
    free(mesh);
}

void initVertex(Vertex* v, double x, double y, int index, HalfEdge* he) {
    if (!v) return;
    v->coord[0] = x;
    v->coord[1] = y;
    v->index = index;
    v->half_edge = he;
}

void initFace(Face* f, int index, HalfEdge* he) {
    if (!f) return;
    f->index = index;
    f->half_edge = he;
}

void initHalfEdge(HalfEdge* he, int index, HalfEdge* twin, HalfEdge* next, HalfEdge* prev, Face* face, Vertex* vertex) {
    if (!he) return;
    he->index = index;
    he->twin = twin;
    he->next = next;
    he->prev = prev;
    he->face = face;
    he->vertex = vertex;
}