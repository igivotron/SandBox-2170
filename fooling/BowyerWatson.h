#ifndef BOWYERWATSON_H
#define BOWYERWATSON_H

#include "HalfEdge.h"
#include <stdlib.h>
#include <stdio.h>

double ** readPointsFromFile(const char* filename, int* numPoints);
int freePoints(double** points, int numPoints);
int** findSuperTriangle(double** points, int* numPoints, int L);
int freeSuperTriangle(int** superTriangle);

#endif // BOWYERWATSON_H