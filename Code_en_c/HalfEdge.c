#include <stdlib.h>
#include <stdio.h>
#include "headers/HalfEdge.h"


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

int writeMeshToFile(const TriangularMesh *mesh, const char *filename) {
    if (!mesh || !filename) return -1;
    FILE *f = fopen(filename, "w");
    if (!f) return -2;

    // Header: counts
    fprintf(f, "%zu\n", mesh->face_count);
    // write for each face the coordinates of its vertices
    for (size_t i = 0; i < mesh->face_count; ++i) {
        const Face *face = &mesh->faces[i];
        HalfEdge *he = face->half_edge;
        for (int j = 0; j < 3; ++j) {
            if (he) {
                fprintf(f, "%g %g ", he->vertex->coord[0], he->vertex->coord[1]);
                he = he->next;
            }

        }
        fprintf(f, "\n");
    }

    fclose(f);
    return 0;
}