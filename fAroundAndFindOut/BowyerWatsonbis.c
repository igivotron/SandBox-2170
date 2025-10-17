#include "headers/BowyerWatsonbis.h"
#include <stdlib.h>
#include <stdio.h>


double** readPointsFromFile(const char* filename, int* numPoints){
    FILE* file = fopen(filename, "r");
    if (!file) {perror("Failed to open file"); return NULL;}
    fscanf(file, "%d", numPoints);
    double ** vertices = (double**)malloc((*numPoints+4) * sizeof(double*));
    for (int i = 0; i < (*numPoints+4); i++) vertices[i] = (double*)malloc(3 * sizeof(double));
    for (int i = 0; i < *numPoints; i++) {
        fscanf(file, "%lf %lf", &vertices[i][0], &vertices[i][1]);
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

int** findSuperTriangle(double** points, int* numPoints, double L){
    double minX = points[0][0], minY = points[0][1];
    double maxX = points[0][0], maxY = points[0][1];
    for (int i = 1; i < *numPoints; i++){
        if (points[i][0] < minX) minX = points[i][0];
        if (points[i][0] > maxX) maxX = points[i][0];
        if (points[i][1] < minY) minY = points[i][1];
        if (points[i][1] > maxY) maxY = points[i][1];
    }
    for (int i = *numPoints; i < *numPoints + 4; i++) points[i] = (double*)malloc(3 * sizeof(double));
    points[*numPoints][0] = minX - L; points[*numPoints][1] = minY - L; points[*numPoints][2] = 0.0;
    points[*numPoints+1][0] = maxX + L; points[*numPoints+1][1] = minY - L; points[*numPoints+1][2] = 0.0;
    points[*numPoints+2][0] = maxX + L; points[*numPoints+2][1] = maxY + L; points[*numPoints+2][2] = 0.0;
    points[*numPoints+3][0] = minX - L; points[*numPoints+3][1] = maxY + L; points[*numPoints+3][2] = 0.0;

    int** superTriangle = (int**)malloc(2 * sizeof(int*));
    for (int i = 0; i < 2; i++) superTriangle[i] = (int*)malloc(3 * sizeof(int));
    superTriangle[0][0] = *numPoints; superTriangle[0][1] = *numPoints + 1; superTriangle[0][2] = *numPoints + 3;
    superTriangle[1][0] = *numPoints + 1; superTriangle[1][1] = *numPoints + 2; superTriangle[1][2] = *numPoints + 3;

    *numPoints += 4;
    return superTriangle;
}

int freeSuperTriangle(int** superTriangle){
    for (int i = 0; i < 2; i++) free(superTriangle[i]);
    free(superTriangle);
    return 0;
}

// int inCircle(Vertex* v, Face* f){
//     HalfEdge* he1 = f->halfEdge;
//     HalfEdge* he2 = he1->next;
//     HalfEdge* he3 = he2->next;
//     Vertex* v1 = he1->vertex;
//     Vertex* v2 = he2->vertex;
//     Vertex* v3 = he3->vertex;
//     double ax = v1->x - v->x;
//     double ay = v1->y - v->y;
//     double bx = v2->x - v->x;
//     double by = v2->y - v->y;
//     double cx = v3->x - v->x;
//     double cy = v3->y - v->y;
//     double det = (ax * ax + ay * ay) * (bx * cy - by * cx) -
//                  (bx * bx + by * by) * (ax * cy - ay * cx) +
//                  (cx * cx + cy * cy) * (ax * by - ay * bx);

//     return det > 0;
// }

// int triangleOrientation(Face* f){
//     HalfEdge* he1 = f->halfEdge;
//     HalfEdge* he2 = he1->next;
//     HalfEdge* he3 = he2->next;
//     Vertex* v1 = he1->vertex;
//     Vertex* v2 = he2->vertex;
//     Vertex* v3 = he3->vertex;
//     double det = (v2->x - v1->x) * (v3->y - v1->y) - (v3->x - v1->x) * (v2->y - v1->y);
//     return det > 0; // return 1 if counter-clockwise, 0 if clockwise
// }

// int sameEdge(Edge* e1, Edge* e2){
//     return ( (e1->v1 == e2->v1 && e1->v2 == e2->v2) || (e1->v1 == e2->v2 && e1->v2 == e2->v1) );
// }

// void addPointToMesh(TriangularMesh* TD, Vertex* v){
//     // Find cavity
//     Face** cavity = malloc(sizeof(Face*) * TD->numFaces);
//     int cavitySize = 0;
//     for (int i = 0; i < TD->numFaces; i++){
//         if (!TD->faces[i].valid) continue;
//         if (inCircle(v, &TD->faces[i])){
//             cavity[cavitySize++] = &TD->faces[i];
//         }
//     }

//     // Collect boundary edges
//     Edge* boundaryEdges = malloc(sizeof(Edge) * cavitySize * 3);
//     int boundarySize = 0;
//     for (int i = 0; i < cavitySize; i++){
//         HalfEdge* he = cavity[i]->halfEdge;
//         for (int k = 0; k < 3; k++){
//             Edge e = {he->vertex, he->next->vertex};
//             int shared = 0;
//             for (int j = 0; j < cavitySize; j++){
//                 if (i == j) continue;
//                 HalfEdge* he2 = cavity[j]->halfEdge;
//                 for (int m = 0; m < 3; m++){
//                     Edge e2 = {he2->vertex, he2->next->vertex};
//                     if (sameEdge(&e, &e2)){
//                         shared = 1;
//                         break;
//                     }
//                     he2 = he2->next;
//                 }
//                 if (shared) break;
//             }
//             if (!shared){boundaryEdges[boundarySize++] = e;}
//             he = he->next;
//         }
//     }
//     // Remove faces in cavity
//     for (int i = 0; i < cavitySize; i++){
//         cavity[i]->valid = 0;
//         HalfEdge* he = cavity[i]->halfEdge;
//         for (int k = 0; k < 3; k++){
//             he->valid = 0;
//             he = he->next;
//         }
//     }

//     // Create new faces
//     for (int i=0; i<boundarySize; i++){
//         Vertex* v1 = boundaryEdges[i].v1;
//         Vertex* v2 = boundaryEdges[i].v2;

//         // TODO: PUT SIZE OF FACE AND HALFEDGE ARRAYS to allocate more memory if needed
//         // For now, we assume there is enough space mdr

//         HalfEdge* he1 = &TD->halfEdges[TD->numHalfEdges];
//         HalfEdge* he2 = &TD->halfEdges[TD->numHalfEdges + 1];
//         HalfEdge* he3 = &TD->halfEdges[TD->numHalfEdges + 2];
//         he1->index = TD->numHalfEdges;
//         he2->index = TD->numHalfEdges + 1;
//         he3->index = TD->numHalfEdges + 2;
//         he1->valid = 1; he2->valid = 1; he3->valid = 1;
//         he1->vertex = v1; he2->vertex = v2; he3->vertex = v;
//         he1->next = he2; he2->next = he3; he3->next = he1;
//         he1->prev = he3; he2->prev = he1; he3->prev = he2;
//         he1->face = &TD->faces[TD->numFaces];
//         he2->face = &TD->faces[TD->numFaces];
//         he3->face = &TD->faces[TD->numFaces];
//         v1->halfEdge = he1;
//         v2->halfEdge = he2;
//         v->halfEdge = he3;
//         TD->faces[TD->numFaces].index = TD->numFaces;
//         TD->faces[TD->numFaces].halfEdge = he1;
//         TD->faces[TD->numFaces].valid = 1;
//         TD->numHalfEdges += 3;
//         TD->numFaces += 1;

//     }

//     free(cavity);
//     free(boundaryEdges);
// }




int main(){
    const char* filename = "../inputs/100pts";
    int numPoints = 0;
    double** vertices = readPointsFromFile(filename, &numPoints);
    int** superTriangles = findSuperTriangle(vertices, &numPoints, 0.1);

    TriangularMesh* mesh = createTriangularMesh(vertices, numPoints, superTriangles, 2);
    // Cause des problèmes de mémoire environ 10% du temps
    // freePoints(vertices, numPoints);
    // freeSuperTriangle(superTriangles);

    // for (int i = 0; i < numPoints; i++){
    //     addPointToMesh(mesh, &mesh->vertices[i]);
    // }

    saveMeshToOBJ(mesh, "outputbis.obj");
    freeTriangularMesh(mesh);

    return 0;
}