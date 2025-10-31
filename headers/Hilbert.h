#ifndef HILBERT_H
#define HILBERT_H
#include "stdio.h"
#include "stdlib.h"

#define DEPTH 12

typedef struct {
    int **hilbert;
    int depth;
} HilbertContext;

void HilbertCoord(double x, double y, double x0, double y0, double xRed, double yRed, double xBlue, double yBlue, int depth, int* bits);
int* sortHilbert(double ** vertices, int N, int depth);
int cmpHilbert(const void* a, const void* b);
#endif
