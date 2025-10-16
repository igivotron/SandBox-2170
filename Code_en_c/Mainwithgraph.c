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

    int *triToCheck = (int *)malloc(100 * sizeof(int)); // contains index of the half-edges for which the twin is in the face we have to check
    if (!triToCheck) {
        fprintf(stderr, "Memory allocation failed\n");
        free(half_edges_to_remove);
        free(faces_to_remove);
        free(mesh);
        return NULL;
    }
    int triToCheck_count = 0;

    for (size_t i = 0; i < n; ++i) {
        initVertex(&vertices[nv], x[i], y[i], nv, &half_edges[0]);
        nv++;
        // parcours des triangles pour savoir pour lesquels le cercle circonscrit contient le point rajouté
        double *d = vertices[nv].coord;
        int cavite = -1; // index of the face that will be the cavity
        for (size_t j = 0; j < nf; ++j) {
            Face *cur_face = &faces[j];
            double *a = cur_face->half_edge->vertex->coord;
            double *b = cur_face->half_edge->next->vertex->coord;
            double *c = cur_face->half_edge->next->next->vertex->coord;
            if (inCircle(a, b, c, d)) {
                cavite = j;
                // add the half-edges to check
                triToCheck[0] = cur_face->half_edge->index;
                triToCheck[1] = cur_face->half_edge->next->index;
                triToCheck[2] = cur_face->half_edge->next->next->index;
                triToCheck_count = 3;
                // mark the face for removal
                faces_to_remove[faces_to_remove_count] = cavite;
                faces_to_remove_count++;
                break;
            }
        }
        while (triToCheck_count > 0) {
            int he_idx = triToCheck[--triToCheck_count];
            HalfEdge *he = &half_edges[he_idx];
            if (he->twin->face == NULL) {
                // this half-edge is on the boundary of the mesh
                continue;
            }
            Face *cur_face = he->twin->face;
            double *a = cur_face->half_edge->vertex->coord;
            double *b = cur_face->half_edge->next->vertex->coord;
            double *c = cur_face->half_edge->next->next->vertex->coord;
            if (inCircle(a, b, c, d)) {
                // mark the face for removal
                faces_to_remove[faces_to_remove_count] = cur_face->index;
                faces_to_remove_count++;

                //update the cavity
                he->prev->next = he->twin->next;
                he->twin->next->prev = he->prev;
                he->next->prev = he->twin->prev;
                he->twin->prev->next = he->next;

                // add the half-edges to check if they are not already in the list
                if (he->twin->next->twin->face != NULL) {
                    int already_in = 0;
                    for (size_t k = 0; k < triToCheck_count; ++k) {
                        if ((&half_edges[triToCheck[k]])->twin->face == he->twin->next->twin->face) {
                            already_in = 1;
                            break;
                        }
                    }
                    if (!already_in) {
                        triToCheck[triToCheck_count] = he->twin->next->index;
                    }
                }
                if (he->twin->next->next->twin->face != NULL) {
                    int already_in = 0;
                    for (size_t k = 0; k < triToCheck_count; ++k) {
                        if ((&half_edges[triToCheck[k]])->twin->face == he->twin->next->next->twin->face) {
                            already_in = 1;
                            break;
                        }
                    }
                    if (!already_in) {
                        triToCheck[triToCheck_count] = he->twin->next->next->index;
                    }
                }
            }
        }

        HalfEdge *he = faces[cavite].half_edge;
        while (1) {
            // for every edge of the cavity, create a new triangle
            int face_idx = nf;
            if (faces_to_remove_count != 0) {
                int face_idx = faces_to_remove[faces_to_remove_count - 1];
                faces_to_remove_count--;
            }
            else {
                nf++;
            }
            initFace(&faces[face_idx], face_idx, he);

            // create the half-edge
            int he_idx1 = nh;
            if (half_edges_to_remove_count != 0) {
                he_idx1 = half_edges_to_remove[half_edges_to_remove_count - 1];
                half_edges_to_remove_count--;
            }
            else {
                nh++;
            }
            initHalfEdge(&half_edges[he_idx1], he_idx1, NULL, he, NULL, &faces[face_idx], &vertices[nv - 1]);

            int he_idx2 = nh;
            if (half_edges_to_remove_count != 0) {
                he_idx2 = half_edges_to_remove[half_edges_to_remove_count - 1];
                half_edges_to_remove_count--;
            }
            else {
                nh++;
            }
            initHalfEdge(&half_edges[he_idx2], he_idx2, NULL, &half_edges[he_idx1], he, &faces[face_idx], he->next->vertex);
            half_edges[he_idx1].prev = &half_edges[he_idx2];
            half_edges[he_idx2].twin = he->next->prev;
            int he_idx = he->index;
            he = he->next;
            half_edges[he_idx].next = &half_edges[he_idx2];
            half_edges[he_idx].prev = &half_edges[he_idx1];
            half_edges[he_idx].face = &faces[face_idx];

            if (he == faces[cavite].half_edge) {
                he->next->twin = half_edges[he_idx].prev;
                half_edges[he_idx].prev->twin = he->next;
                break;
            }
        }
    }
    free(triToCheck);
    free(half_edges_to_remove);
    free(faces_to_remove);
    return mesh;
}