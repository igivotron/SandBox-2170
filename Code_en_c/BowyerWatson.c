#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "headers/BowyerWatson.h"

double** readPointsFromFile(const char* filename, int* numPoints){
    FILE* file = fopen(filename, "r");
    if (!file) {perror("Failed to open file"); return NULL;}

     if (fscanf(file, "%d", numPoints) != 1 ){
        perror("Failed to read number of points");
        fclose(file);
        return NULL;
     }
     
    double ** vertices = (double**)malloc((*numPoints) * sizeof(double*));
    for (int i = 0; i < (*numPoints); i++) vertices[i] = (double*)malloc(3 * sizeof(double));
    for (int i = 0; i < *numPoints; i++) {
        if (fscanf(file, "%lf %lf", &vertices[i][0], &vertices[i][1]) != 2) {
            perror("Failed to read point coordinates");
            for (int j = 0; j <= i; j++) free(vertices[j]);
            free(vertices);
            fclose(file);
            return NULL;
        }
        vertices[i][2] = 0.0;
    }
    fclose(file);

    return vertices;
}

int freePoints(double** points, int numPoints){
    for (int i = 0; i < numPoints; i++) free(points[i]);
    free(points);
    return 0;
}

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

int isPointInCircumcircle(double* a, double* b, double* c, double* d){
    return incircle(a, b, c, d) > 0;
}

int isTriangleOrientedCCW(double* a, double* b, double* c){
    return orient2d(a, b, c) > 0;
}

void superTriangle(TriangularMesh* mesh , double *x, double *y, size_t n) {
    Vertex *vertices = mesh->vertices;
    HalfEdge *half_edges = mesh->half_edges;
    Face *faces = mesh->faces;

    // vertices=[[x1,y1],[x2,y2],...]
    // triangles=[[4,36,24],[12,97,53],...]
    double xmin = min(x,n);
    double ymin = min(y,n);
    double xmax = max(x,n);
    double ymax = max(y,n);

    // add the four points
    double L = (xmax - xmin + ymax - ymin) / 20;
    mesh->vertex_count=4;
    initVertex(&vertices[0], xmin - L, ymin - L, 0, &half_edges[0]);
    initVertex(&vertices[1], xmax + L, ymin - L, 1, &half_edges[1]);
    initVertex(&vertices[2], xmax + L, ymax + L, 2, &half_edges[2]);
    initVertex(&vertices[3], xmin - L, ymax + L, 3, &half_edges[3]);
    mesh->face_count=2;
    initFace(&faces[0], 0, &half_edges[0]);
    initFace(&faces[1], 1, &half_edges[1]);
    mesh->half_edge_count=6;
    initHalfEdge(&half_edges[0], 0, NULL, &half_edges[4], &half_edges[3], &faces[0], &vertices[0]);
    initHalfEdge(&half_edges[1], 1, NULL, &half_edges[2], &half_edges[5], &faces[1], &vertices[1]);
    initHalfEdge(&half_edges[2], 2, NULL, &half_edges[5], &half_edges[1], &faces[1], &vertices[2]);
    initHalfEdge(&half_edges[3], 3, NULL, &half_edges[0], &half_edges[4], &faces[0], &vertices[3]);
    initHalfEdge(&half_edges[4], 4, &half_edges[5], &half_edges[3], &half_edges[0], &faces[0], &vertices[1]);
    initHalfEdge(&half_edges[5], 5, &half_edges[4], &half_edges[1], &half_edges[2], &faces[1], &vertices[3]);
}

int addPoint(TriangularMesh* mesh, double x, double y) {
    int nv = mesh->vertex_count;
    Vertex *vertices = mesh->vertices;
    int nh = mesh->half_edge_count;
    HalfEdge *half_edges = mesh->half_edges;
    int nf = mesh->face_count;
    Face *faces = mesh->faces;
    
    initVertex(&vertices[nv], x, y, nv, &half_edges[0]);
    nv++;
    // parcours des triangles pour savoir pour lesquels le cercle circonscrit contient le point rajouté
    double *d = vertices[nv-1].coord;
    int cavite = -1; // index of the face that will be the cavity

    for (size_t j = nf-1; j >= 0; --j) {
        Face *cur_face = &faces[j];
        double *a = cur_face->half_edge->vertex->coord;
        double *b = cur_face->half_edge->next->vertex->coord;
        double *c = cur_face->half_edge->next->next->vertex->coord;
    
        if (isPointInCircumcircle(a, b, c, d)) {
            cavite = j;
            // add the half-edges to check
            if (cur_face->half_edge->twin != NULL) {
                mesh->triToCheck[mesh->triToCheck_count++] = cur_face->half_edge->index;
            }
            if (cur_face->half_edge->next->twin != NULL) {
                mesh->triToCheck[mesh->triToCheck_count++] = cur_face->half_edge->next->index;
            }
            if (cur_face->half_edge->next->next->twin != NULL) {
                mesh->triToCheck[mesh->triToCheck_count++] = cur_face->half_edge->next->next->index;
            }
            break;
        }
    }
    
    while (mesh->triToCheck_count > 0) {
        int he_idx = mesh->triToCheck[--mesh->triToCheck_count];
        HalfEdge *he = &half_edges[he_idx];
        int twin_index = he->twin->index;
        Face *cur_face = he->twin->face;
        double *a = cur_face->half_edge->vertex->coord;
        double *b = cur_face->half_edge->next->vertex->coord;
        double *c = cur_face->half_edge->next->next->vertex->coord;

        if (isPointInCircumcircle(a, b, c, d)) {
            // mark the face for removal
            mesh->faces_to_remove[mesh->faces_to_remove_count++] = cur_face->index;
            
            //update the cavity
            he->prev->next = he->twin->next;
            he->twin->next->prev = he->prev;
            he->next->prev = he->twin->prev;
            he->twin->prev->next = he->next;
            // add he to the list of half-edges to remove
            mesh->half_edges_to_remove[mesh->half_edges_to_remove_count++] = he->index;
            mesh->half_edges_to_remove[mesh->half_edges_to_remove_count++] = twin_index;
            // add the half-edges to check if they are not already in the list
            if (he->twin->next->twin != NULL) {
                int already_in = 0;
                for (size_t k = 0; k < mesh->triToCheck_count; ++k) {
                    if ((&half_edges[mesh->triToCheck[k]])->twin->face == he->twin->next->twin->face) {
                        already_in = 1;
                        break;
                    }
                }
                if (!already_in) {
                    mesh->triToCheck[mesh->triToCheck_count++] = he->twin->next->index;
                }
            }
            if (he->twin->next->next->twin != NULL) {
                int already_in = 0;
                for (size_t k = 0; k < mesh->triToCheck_count; ++k) {
                    if ((&half_edges[mesh->triToCheck[k]])->twin->face == he->twin->next->next->twin->face) {
                        already_in = 1;
                        break;
                    }
                }
                if (!already_in) {
                    mesh->triToCheck[mesh->triToCheck_count++] = he->twin->next->next->index;
                }
            }
            if (he == faces[cavite].half_edge) faces[cavite].half_edge = he->next;
        }
    }

    // mark the face for removal
    mesh->faces_to_remove[mesh->faces_to_remove_count++] = cavite;
    HalfEdge *he = faces[cavite].half_edge;
    while (1) {
        // for every edge of the cavity, create a new triangle
        int face_idx = nf;
        if (mesh->faces_to_remove_count != 0) face_idx = mesh->faces_to_remove[--mesh->faces_to_remove_count];
        else nf++;
        initFace(&faces[face_idx], face_idx, he);

        // create the half-edge
        int he_idx1 = nh;
        if (mesh->half_edges_to_remove_count != 0) he_idx1 = mesh->half_edges_to_remove[--mesh->half_edges_to_remove_count];
        else nh++;
        initHalfEdge(&half_edges[he_idx1], he_idx1, NULL, he, NULL, &faces[face_idx], &vertices[nv - 1]);

        int he_idx2 = nh;
        if (mesh->half_edges_to_remove_count != 0) he_idx2 = mesh->half_edges_to_remove[--mesh->half_edges_to_remove_count];
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

    mesh->vertex_count = nv;
    mesh->face_count   = nf;
    mesh->half_edge_count = nh;
    return 0;
}

TriangularMesh *del2d(double *x, double *y, size_t n) {
    TriangularMesh *mesh = allocMesh(n);
    if (!mesh) {
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }

    superTriangle(mesh, x, y, n);

    // keep track of faces and half-edges to remove
    mesh->faces_to_remove = (int *)malloc(2*n*sizeof(int));
    if (!mesh->faces_to_remove) {
        fprintf(stderr, "Memory allocation failed\n");
        free(mesh);
        return NULL;
    }
    mesh->faces_to_remove_count = 0;

    mesh->half_edges_to_remove = (int *)malloc(6*n*sizeof(int));
    if (!mesh->half_edges_to_remove) {
        fprintf(stderr, "Memory allocation failed\n");
        free(mesh->faces_to_remove);
        free(mesh);
        return NULL;
    }
    mesh->half_edges_to_remove_count = 0;

    mesh->triToCheck = (int *)malloc(2*n*sizeof(int)); // contains index of the half-edges for which the twin is in the face we have to check
    if (!mesh->triToCheck) {
        fprintf(stderr, "Memory allocation failed\n");
        free(mesh->half_edges_to_remove);
        free(mesh->faces_to_remove);
        free(mesh);
        return NULL;
    }
    mesh->triToCheck_count = 0;

    printf("Starting Bowyer-Watson with %zu points...\n", n);
    for (size_t i = 0; i < n; ++i) {
        addPoint(mesh, x[i], y[i]);
    }

    free(mesh->triToCheck);
    free(mesh->half_edges_to_remove);
    free(mesh->faces_to_remove);
    return mesh;
}

int del2d_py(double *x, double *y, size_t n) {
    TriangularMesh *mesh = del2d(x, y, n);
    if (mesh) {
        printf("Mesh created with %zu vertices, %zu faces, and %zu half-edges.\n",
               mesh->vertex_count, mesh->face_count, mesh->half_edge_count);
        writeMeshToFile(mesh, "Code_en_c/output_mesh.txt");
        freeMesh(mesh);
    } else {
        printf("Failed to create mesh.\n");
        return -1;
    }

    return 0;
}
