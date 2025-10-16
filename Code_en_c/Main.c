#include <stdio.h>
#include <stdlib.h>
#include "HalfEdge.h"
#include "robust_predicates_c.h"

double min(const double *arr, size_t n) {
    if (n == 0) return 0.0;
    double m = arr[0];
    for (size_t i = 1; i < n; ++i) {
        if (arr[i] < m) m = arr[i];
    }
    return m;
}

double max(const double *arr, size_t n) {
    if (n == 0) return 0.0;
    double m = arr[0];
    for (size_t i = 1; i < n; ++i) {
        if (arr[i] > m) m = arr[i];
    }
    return m;
}



TriangularMesh *del2d(double *x, double *y, size_t n) {
    TriangularMesh *mesh = allocMesh(n);
    if (!mesh) {
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }
    int nv = mesh->vertex_count;
    Vertex *vertices = mesh->vertices;
    int nh = mesh->half_edge_count;
    HalfEdge *half_edges = mesh->half_edges;
    int nf = mesh->face_count;
    Face *faces = mesh->faces;

    // vertices=[[x1,y1],[x2,y2],...]
    // triangles=[[4,36,24],[12,97,53],...]
    double xmin = min(x,n);
    double ymin = min(y,n);
    double xmax = max(x,n);
    double ymax = max(y,n);

    // add the four points
    double L = (xmax - xmin + ymax - ymin) / 20;
    nv=4;
    initVertex(&vertices[0], xmin - L, ymin - L, 0, &half_edges[0]);
    initVertex(&vertices[1], xmax + L, ymin - L, 1, &half_edges[1]);
    initVertex(&vertices[2], xmax + L, ymax + L, 2, &half_edges[2]);
    initVertex(&vertices[3], xmin - L, ymax + L, 3, &half_edges[3]);
    nf=2;
    initFace(&faces[0], 0, &half_edges[0]);
    initFace(&faces[1], 1, &half_edges[1]);
    nh=6;
    initHalfEdge(&half_edges[0], 0, NULL, &half_edges[4], &half_edges[3], &faces[0], &vertices[0]);
    initHalfEdge(&half_edges[1], 1, NULL, &half_edges[2], &half_edges[5], &faces[1], &vertices[1]);
    initHalfEdge(&half_edges[2], 2, NULL, &half_edges[5], &half_edges[1], &faces[1], &vertices[2]);
    initHalfEdge(&half_edges[3], 3, NULL, &half_edges[0], &half_edges[4], &faces[0], &vertices[3]);
    initHalfEdge(&half_edges[4], 1, &half_edges[5], &half_edges[3], &half_edges[0], &faces[0], &vertices[1]);
    initHalfEdge(&half_edges[5], 2, &half_edges[4], &half_edges[1], &half_edges[2], &faces[1], &vertices[2]);

    // keep track of faces and half-edges to remove
    int *faces_to_remove = (int *)malloc(20 * sizeof(int));
    if (!faces_to_remove) {
        fprintf(stderr, "Memory allocation failed\n");
        free(mesh);
        return NULL;
    }
    int faces_to_remove_count = 0;
    int *half_edges_to_remove = (int *)malloc(100 * sizeof(int));
    if (!half_edges_to_remove) {
        fprintf(stderr, "Memory allocation failed\n");
        free(faces_to_remove);
        free(mesh);
        return NULL;
    }
    int half_edges_to_remove_count = 0;

    int *edges_searching_their_twins = (int *)malloc(100 * sizeof(int));
    if (!edges_searching_their_twins) {
        fprintf(stderr, "Memory allocation failed\n");
        free(half_edges_to_remove);
        free(faces_to_remove);
        free(mesh);
        return NULL;
    }

    for (size_t i = 0; i < n; ++i) {
        initVertex(&vertices[nv], x[i], y[i], nv, &half_edges[0]);
        nv++;
        // parcours des triangles pour savoir pour lesquels le cercle circonscrit contient le point rajouté
        double *d = vertices[nv].coord;
        int cavite[20]; // indices des triangles dans la cavité
        int cavite_count = 0;
        for (size_t j = 0; j < nf; ++j) {
            Face *cur_face = &faces[j];
            double *a = cur_face->half_edge->vertex->coord;
            double *b = cur_face->half_edge->next->vertex->coord;
            double *c = cur_face->half_edge->next->next->vertex->coord;
            if (inCircle(a, b, c, d)) {
                // add to the cavity
                cavite[cavite_count] = j;
                cavite_count++;
            }
        }
        for (size_t j = 0; j < cavite_count; ++j) {
            HalfEdge *v1 = faces[cavite[j]].half_edge;
            HalfEdge *vTwin = v1->twin;
            for (size_t l = 0; l < 3; ++l) {
                int found = 0;
                for (size_t k = 0; k < cavite_count; ++k) {
                    HalfEdge *v2 = faces[cavite[k]].half_edge;
                    if (k != j && (vTwin == v2 || vTwin == v2->next || vTwin == v2->next->next)) {
                        found = 1;
                        break;
                    }
                }
                if (!found) {
                    if (faces_to_remove_count != 0) {
                        int face_idx = faces_to_remove[faces_to_remove_count - 1];
                        initFace(&faces[face_idx], face_idx, v1);
                        faces_to_remove_count--;
                        if (half_edges_to_remove_count != 0) {
                            int he_idx = half_edges_to_remove[half_edges_to_remove_count - 1];
                            initHalfEdge(&half_edges[he_idx], he_idx, NULL, v1, v1->prev, &faces[face_idx], &vertices[nv - 1]);
                            half_edges_to_remove_count--;
                        } else {
                            initHalfEdge(&half_edges[nh], nh, NULL, v1, v1->prev, &faces[face_idx], &vertices[nv - 1]);
                            nh++;
                        }
                    }
                    triangles.append([v1, v2, len(vertices) - 1]);
                }
                v1 = v1->next;
                vTwin = v1->twin;
            }
        for (size_t idx = 0; idx < cavite_count; ++idx) {
            faces_to_remove[faces_to_remove_count] = cavite[idx];
            faces_to_remove_count++;
        }
    }
    return mesh;
}