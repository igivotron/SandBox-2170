#include <stdio.h>
#include <stdlib.h>
#include "headers/BowyerWatson.h"
#include "headers/Hilbert.h"
#include "headers/predicates.h"


// TODO: Rajouter le predicat d'orientation des triangles
// Pas suffisemment de mémoire alloué pour les grands nuages de points

int main(void) {
    exactinit();
    const char* input_file = "../inputs/100pts";
    int numPoints = 0;
    double ** points = readPointsFromFile(input_file, &numPoints);
    double *x = (double*)malloc(numPoints * sizeof(double));
    double *y = (double*)malloc(numPoints * sizeof(double));
    int* indexHilbert = sortHilbert(points, numPoints, 9);

    for (int i =0; i < numPoints; i++){
        int idx = indexHilbert[i];
        x[i] = points[idx][0];
        y[i] = points[idx][1];
    }
    TriangularMesh *mesh = del2d(x, y, numPoints);
    writeMeshToFile(mesh, "output_mesh.txt");
    freeMesh(mesh);
    freePoints(points, numPoints);
    free(x);
    free(y);
    free(indexHilbert);
    return 0;

}