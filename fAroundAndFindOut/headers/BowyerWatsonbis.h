#ifndef BOWYERWATSON_H
#define BOWYERWATSON_H

#include "HalfEdge.h"
#include <stdlib.h>
#include <stdio.h>

typedef struct Edge{
    Vertex* v1;
    Vertex* v2;
} Edge;


double ** readPointsFromFile(const char* filename, int* numPoints);
int freePoints(double** points, int numPoints);
int** findSuperTriangle(double** points, int* numPoints, double L);
int freeSuperTriangle(int** superTriangle);

void addPointToMesh(TriangularMesh* TD, Vertex* v);
int inCircle(Vertex* v, Face* f);
int sameEdge(Edge* e1, Edge* e2);
int triangleOrientation(Face* f);

#endif // BOWYERWATSON_H