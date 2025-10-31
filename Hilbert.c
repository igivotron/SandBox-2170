#include "headers/Hilbert.h"
#include "stdlib.h"
#include "stdio.h"


static int **g_hilbertCoords;
static int g_depth;

void HilbertCoord(double x, double y, double x0, double y0, double xRed, double yRed, double xBlue, double yBlue, int depth, int* bits){
    if (depth == 0) return;
    int i;
    double coordRed, coordBlue, temp;

    for (i=0; i < depth; i++){
        coordRed = (x - x0) * xRed + (y - y0) * yRed;
        coordBlue = (x - x0) * xBlue + (y - y0) * yBlue;   
        xRed /= 2;
        yRed /= 2;
        xBlue /= 2;
        yBlue /= 2;
        if (coordRed <= 0 && coordBlue <= 0){ // quadrant 0
            x0 -= (xRed + xBlue);
            y0 -= (yRed + yBlue);

            //SWAP
            temp = xRed;
            xRed = xBlue;
            xBlue = temp;
            temp = yRed;
            yRed = yBlue;
            yBlue = temp;
            bits[i] = 0;

        } else if (coordRed <= 0 && coordBlue >= 0){ // quadrant 1
            x0 += xBlue - xRed;
            y0 += yBlue - yRed;
            bits[i] = 1;
        } else if (coordRed >= 0 && coordBlue >= 0){ // quadrant 2
            x0 += xRed + xBlue;
            y0 += yRed + yBlue;
            bits[i] = 2;
        } else { // quadrant 3
            x0 += xRed - xBlue;
            y0 += yRed - yBlue;

            //SWAP
            temp = xRed;
            xRed = xBlue;
            xBlue = temp;
            temp = yRed;
            yRed = yBlue;
            yBlue = temp;
        
            xBlue = -xBlue;
            yBlue = -yBlue;
            xRed = -xRed;
            yRed = -yRed;
            bits[i] = 3;
        }
    }
    return;
}


int cmpHilbert(const void *a, const void *b) {
    int ia = *(const int *)a;
    int ib = *(const int *)b;

    for (int i = 0; i < g_depth; i++) {
        if (g_hilbertCoords[ia][i] < g_hilbertCoords[ib][i]) return -1;
        if (g_hilbertCoords[ia][i] > g_hilbertCoords[ib][i]) return 1;
    }
    return 0;
}

int* sortHilbert(double ** vertices, int N, int depth){
    // N -= 4; // exclude supertriangle points
    
    int* index = malloc(N * sizeof(int));
    int** hilbertCoords = malloc(N * sizeof(int*));
    double * x = malloc(N * sizeof(double));
    double * y = malloc(N * sizeof(double));

    for (int i = 0; i < N; i++){
        x[i] = vertices[i][0];
        y[i] = vertices[i][1];
    }

    for (int i = 0; i < N; i++){
        hilbertCoords[i] = malloc(depth * sizeof(int));
        HilbertCoord(x[i], y[i], 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, depth, hilbertCoords[i]);
        index[i] = i;
    }

    g_hilbertCoords = hilbertCoords;
    g_depth = depth;

    qsort(index, N, sizeof(int), cmpHilbert);

    for (int i = 0; i < N; i++) free(hilbertCoords[i]);
    free(hilbertCoords);
    free(x);
    free(y);

    return index;
}
