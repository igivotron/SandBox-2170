#include "headers/HalfEdge.h"
#include <stdlib.h>
#include <stdio.h>


TriangularMesh* createTriangularMesh(double** verticies, int numVertices, int** triangles, int numTriangles){
    TriangularMesh* mesh = (TriangularMesh*)malloc(sizeof(TriangularMesh));
    mesh->vertices = (Vertex*) malloc(sizeof(Vertex) * numVertices);
    // mesh->halfEdges = (HalfEdge*) malloc(sizeof(HalfEdge) * numTriangles * 3);
    // mesh->faces = (Face*) malloc(sizeof(Face) * numTriangles);
    mesh->halfEdges = (HalfEdge*) malloc(sizeof(HalfEdge) * numTriangles * 3 * 10000); 
    mesh->faces = (Face*) malloc(sizeof(Face) * numTriangles * 10000);
    mesh->numVertices = numVertices;

    printf("Creating mesh with %d vertices, %d triangles\n", numVertices, numTriangles);
    for (int i =0; i < numVertices; i++) {
        mesh->vertices[i].index = i;
        mesh->vertices[i].x = verticies[i][0];
        mesh->vertices[i].y = verticies[i][1];
        mesh->vertices[i].z = verticies[i][2];
    }

    mesh->numFaces = numTriangles;
    mesh->numHalfEdges = numTriangles * 3;
    int* triangle;
    for (int i=0; i<numTriangles; i++){
        triangle = triangles[i];
        for (int j=0; j<3; j++){
            mesh->halfEdges[i*3 + j].valid = 1;
            mesh->halfEdges[i*3 + j].index = i*3 + j;
            mesh->halfEdges[i*3 + j].vertex = &mesh->vertices[triangle[j]];
            mesh->halfEdges[i*3 + j].face = &mesh->faces[i];
            mesh->halfEdges[i*3 + j].next = &mesh->halfEdges[i*3 + (j+1)%3];
            mesh->halfEdges[i*3 + j].prev = &mesh->halfEdges[i*3 + (j+2)%3];
            mesh->vertices[triangle[j]].halfEdge = &mesh->halfEdges[i*3 + j];
        }
        mesh->faces[i].index = i;
        mesh->faces[i].halfEdge = &mesh->halfEdges[i*3];
        mesh->faces[i].valid = 1;
    }
    // find opposite half-edges
    for (int i=0; i<mesh->numHalfEdges; i++){
        for (int j=0; j<mesh->numHalfEdges; j++){
            mesh->halfEdges[i].opposite = NULL;
            if (i != j 
                && mesh->halfEdges[i].vertex == mesh->halfEdges[j].next->vertex
                && mesh->halfEdges[i].next->vertex == mesh->halfEdges[j].vertex) {
                mesh->halfEdges[i].opposite = &mesh->halfEdges[j];
                mesh->halfEdges[j].opposite = &mesh->halfEdges[i];
                break;
            }
        }
    }
    return mesh;
}

int freeTriangularMesh(TriangularMesh* mesh){
    if (mesh->vertices) free(mesh->vertices);
    if (mesh->halfEdges) free(mesh->halfEdges);
    if (mesh->faces) free(mesh->faces);
    free(mesh);
    return 0;
}

int printMeshInfo(TriangularMesh* mesh){
    if (!mesh) {
        printf("Mesh is NULL\n");
        return -1;
    }
    printf("Mesh Info:\n");
    printf("Number of Vertices: %d\n", mesh->numVertices);
    printf("Number of HalfEdges: %d\n", mesh->numHalfEdges);
    printf("Number of Faces: %d\n", mesh->numFaces);
    return 0;
}

int saveMeshToOBJ(TriangularMesh* mesh, const char* filename) {
    if (!mesh || !filename) {
        return -1;
    }

    FILE* file = fopen(filename, "w");
    if (!file) return -1;

    fprintf(file,"%d\n", mesh->numVertices);
    for (int i = 0; i < mesh->numVertices; i++) {
        fprintf(file, "v %f %f %f\n", mesh->vertices[i].x, mesh->vertices[i].y, mesh->vertices[i].z);
    }

    for (int i = 0; i < mesh->numFaces; i++) {
        if (!mesh->faces[i].valid) continue;
        HalfEdge* he = mesh->faces[i].halfEdge;
        fprintf(file, "f %d %d %d\n", he->vertex->index, he->next->vertex->index, he->next->next->vertex->index);
    }

    for (int i=0; i<mesh->numHalfEdges; i++){
        if (!mesh->halfEdges[i].valid) continue;
        HalfEdge he = mesh->halfEdges[i];
        fprintf(file, "%d: %d to %d, Face %d, Next %d, Prev %d, Opposite %d\n", 
                he.index, 
                he.vertex ? he.vertex->index : -1,
                he.next ? he.next->vertex->index : -1, 
                he.face ? he.face->index : -1, 
                he.next ? he.next->index : -1, 
                he.prev ? he.prev->index : -1, 
                he.opposite ? he.opposite->index : -1);
    }

    fclose(file);
    return 0;
}