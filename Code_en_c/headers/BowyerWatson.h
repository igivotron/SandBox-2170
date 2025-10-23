#ifndef BOWYERWATSON_H
#define BOWYERWATSON_H

#include <stddef.h>
#include "HalfEdge.h"
#include "predicates.h"



/* Utility routines (defined in BowyerWatson.c) */
double min(const double *arr, size_t n);
double max(const double *arr, size_t n);

/* Geometric predicate: returns non-zero if d lies inside the circumcircle of a,b,c */
int isPointInCircumcircle(double* a, double* b, double* c, double* d);
int isTriangleOrientedCCW(double* a, double* b, double* c);

/* Delaunay triangulation entry points */
TriangularMesh *del2d(double *x, double *y, size_t n);
void superTriangle(TriangularMesh* mesh, double *x, double *y, size_t n);
int del2d_py(double *x, double *y, size_t n);
int addPoint(TriangularMesh* mesh, double x, double y);

double **readPointsFromFile(const char* filename, int* numPoints);
int freePoints(double** points, int numPoints);

#endif /* BOWYERWATSON_H */
