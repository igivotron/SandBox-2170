#include <stdio.h>
#include <stdlib.h>
#include "BowyerWatson.h"
// #include "robust_predicates_c.h"

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

int inCircle(double* a, double* b, double* c, double* d){
    double ax = a[0] - d[0];
    double ay = a[1] - d[1];
    double bx = b[0] - d[0];
    double by = b[1] - d[1];
    double cx = c[0] - d[0];
    double cy = c[1] - d[1];
    double det = (ax * ax + ay * ay) * (bx * cy - by * cx) -
                 (bx * bx + by * by) * (ax * cy - ay * cx) +
                 (cx * cx + cy * cy) * (ax * by - ay * bx);

    return det > 0;
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
    initHalfEdge(&half_edges[4], 4, &half_edges[5], &half_edges[3], &half_edges[0], &faces[0], &vertices[1]);
    initHalfEdge(&half_edges[5], 5, &half_edges[4], &half_edges[1], &half_edges[2], &faces[1], &vertices[3]);

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

    printf("Starting Delaunay triangulation with %zu points\n", n);

    for (size_t i = 0; i < n; ++i) {
        printf("new_point\n");
        initVertex(&vertices[nv], x[i], y[i], nv, &half_edges[0]);
        nv++;
        // parcours des triangles pour savoir pour lesquels le cercle circonscrit contient le point rajouté
        double *d = vertices[nv-1].coord;
        int cavite = -1; // index of the face that will be the cavity
        for (size_t j = 0; j < nf; ++j) {
            Face *cur_face = &faces[j];
            double *a = cur_face->half_edge->vertex->coord;
            double *b = cur_face->half_edge->next->vertex->coord;
            double *c = cur_face->half_edge->next->next->vertex->coord;
            if (inCircle(a, b, c, d)) {
                cavite = j;
                printf("found start of the cavity (%f, %f), (%f, %f), (%f, %f) with point (%f, %f)\n",
                   a[0], a[1], b[0], b[1], c[0], c[1], d[0], d[1]);
                // add the half-edges to check
                if (cur_face->half_edge->twin != NULL) {
                    triToCheck[triToCheck_count++] = cur_face->half_edge->index;
                }
                if (cur_face->half_edge->next->twin != NULL) {
                    triToCheck[triToCheck_count++] = cur_face->half_edge->next->index;
                }
                if (cur_face->half_edge->next->next->twin != NULL) {
                    triToCheck[triToCheck_count++] = cur_face->half_edge->next->next->index;
                }
                break;
            }
        }

        printf("start of the cavity found\n");

        
        while (triToCheck_count > 0) {
            int he_idx = triToCheck[--triToCheck_count];
            HalfEdge *he = &half_edges[he_idx];
            int twin_index = he->twin->index;
            Face *cur_face = he->twin->face;
            double *a = cur_face->half_edge->vertex->coord;
            double *b = cur_face->half_edge->next->vertex->coord;
            double *c = cur_face->half_edge->next->next->vertex->coord;
            // print a b c d
            printf("Checking triangle (%f, %f), (%f, %f), (%f, %f) with point (%f, %f)\n",
                   a[0], a[1], b[0], b[1], c[0], c[1], d[0], d[1]);
            if (inCircle(a, b, c, d)) {
                printf("edge to add to the cavity found\n");
                // mark the face for removal
                faces_to_remove[faces_to_remove_count++] = cur_face->index;
                
                //update the cavity
                he->prev->next = he->twin->next;
                he->twin->next->prev = he->prev;
                he->next->prev = he->twin->prev;
                he->twin->prev->next = he->next;
                // add he to the list of half-edges to remove
                half_edges_to_remove[half_edges_to_remove_count++] = he->index;
                half_edges_to_remove[half_edges_to_remove_count++] = twin_index;
                // add the half-edges to check if they are not already in the list
                if (he->twin->next->twin != NULL) {
                    int already_in = 0;
                    for (size_t k = 0; k < triToCheck_count; ++k) {
                        if ((&half_edges[triToCheck[k]])->twin->face == he->twin->next->twin->face) {
                            already_in = 1;
                            break;
                        }
                    }
                    if (!already_in) {
                        triToCheck[triToCheck_count++] = he->twin->next->index;
                    }
                }
                if (he->twin->next->next->twin != NULL) {
                    int already_in = 0;
                    for (size_t k = 0; k < triToCheck_count; ++k) {
                        if ((&half_edges[triToCheck[k]])->twin->face == he->twin->next->next->twin->face) {
                            already_in = 1;
                            break;
                        }
                    }
                    if (!already_in) {
                        triToCheck[triToCheck_count++] = he->twin->next->next->index;
                    }
                }
                if (he == faces[cavite].half_edge) faces[cavite].half_edge = he->next;
            }
        }
        // printf("cavity created\n");
        // print the vertices of the cavity
        // HalfEdge *he42 = faces[cavite].half_edge;
        // while(1){
        //     printf("(%f, %f)\n", he42->vertex->coord[0], he42->vertex->coord[1]);
        //     he42 = he42->next;
        //     if (he42 == faces[cavite].half_edge) break;
        // }

        // mark the face for removal
        faces_to_remove[faces_to_remove_count++] = cavite;
        HalfEdge *he = faces[cavite].half_edge;
        while (1) {
            printf("coucou\n");
            // for every edge of the cavity, create a new triangle
            int face_idx = nf;
            if (faces_to_remove_count != 0) face_idx = faces_to_remove[--faces_to_remove_count];
            else nf++;
            initFace(&faces[face_idx], face_idx, he);

            // create the half-edge
            int he_idx1 = nh;
            if (half_edges_to_remove_count != 0) he_idx1 = half_edges_to_remove[--half_edges_to_remove_count];
            else nh++;
            initHalfEdge(&half_edges[he_idx1], he_idx1, NULL, he, NULL, &faces[face_idx], &vertices[nv - 1]);

            int he_idx2 = nh;
            if (half_edges_to_remove_count != 0) he_idx2 = half_edges_to_remove[--half_edges_to_remove_count];
            else nh++;
            initHalfEdge(&half_edges[he_idx2], he_idx2, NULL, &half_edges[he_idx1], he, &faces[face_idx], he->next->vertex);
            
            half_edges[he_idx1].prev = &half_edges[he_idx2];
            if (he != faces[cavite].half_edge) {
                half_edges[he_idx1].twin = he->prev->next;
                he->prev->next->twin = &half_edges[he_idx1];
            }
            int he_idx = he->index;
            he = he->next;
            half_edges[he_idx].next = &half_edges[he_idx2];
            half_edges[he_idx].prev = &half_edges[he_idx1];
            half_edges[he_idx].face = &faces[face_idx];

            if (he == faces[cavite].half_edge) {
                he->prev->twin = half_edges[he_idx].next;
                half_edges[he_idx].next->twin = he->prev;
                break;
            }
        }
        printf("nf: %d\n", nf);
        // print the vertices of each face
        for (int f_idx = 0; f_idx < nf; ++f_idx) {
            Face *f = &faces[f_idx];
            HalfEdge *he = f->half_edge;
            printf("Face %d: ", f_idx);
            for (size_t k = 0; k < 3; ++k) {
                printf("(%f, %f) ", he->vertex->coord[0], he->vertex->coord[1]);
                he = he->next;
            }
            printf("\n");
        }
    }

    free(triToCheck);
    free(half_edges_to_remove);
    free(faces_to_remove);
    mesh->vertex_count = nv;
    mesh->face_count   = nf;
    mesh->half_edge_count = nh;
    return mesh;
}

int del2d_py(double *x, double *y, size_t n) {
    TriangularMesh *mesh = del2d(x, y, n);
    if (mesh) {
        printf("Mesh created with %zu vertices, %zu faces, and %zu half-edges.\n",
               mesh->vertex_count, mesh->face_count, mesh->half_edge_count);
        writeMeshToFile(mesh, "output_mesh.txt");
        freeMesh(mesh);
    } else {
        printf("Failed to create mesh.\n");
        return -1;
    }

    return 0;
}
