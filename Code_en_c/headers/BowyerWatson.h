#ifndef BOWYERWATSON_H
#define BOWYERWATSON_H

#include <stddef.h>
#include "HalfEdge.h"

// #ifdef __cplusplus
// extern "C" {
// #endif

/* Utility routines (defined in BowyerWatson.c) */
double min(const double *arr, size_t n);
double max(const double *arr, size_t n);

/* Geometric predicate: returns non-zero if d lies inside the circumcircle of a,b,c */
int inCircle(double* a, double* b, double* c, double* d);

/* Delaunay triangulation entry points */
TriangularMesh *del2d(double *x, double *y, size_t n);
int del2d_py(double *x, double *y, size_t n);

// #ifdef __cplusplus
// }
// #endif

#endif /* BOWYERWATSON_H */
