#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "headers/BowyerWatson.h"

int main() {
    const size_t n = 10000;
    double *x = malloc(n * sizeof(double));
    double *y = malloc(n * sizeof(double));
    if (!x || !y) {
        fprintf(stderr, "Memory allocation failed.\n");
        return 1;
    }

    // Génération aléatoire de points dans [0,1]x[0,1]
    srand((unsigned)time(NULL));
    for (size_t i = 0; i < n; ++i) {
        x[i] = (double)rand() / RAND_MAX;
        y[i] = (double)rand() / RAND_MAX;
    }

    printf("Démarrage de la triangulation de Delaunay avec %zu points...\n", n);

    clock_t start = clock();
    TriangularMesh *mesh = del2d(x, y, n);
    clock_t end = clock();

    if (mesh) {
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        printf("Triangulation terminée en %.3f secondes.\n", elapsed);
        printf("Maillage : %zu sommets, %zu faces, %zu demi-arêtes.\n",
               mesh->vertex_count, mesh->face_count, mesh->half_edge_count);

        freeMesh(mesh);
    } else {
        printf("Erreur : la triangulation a échoué.\n");
    }

    free(x);
    free(y);
    return 0;
}
