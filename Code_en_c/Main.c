#include <stdio.h>
#include <stdlib.h>
#include "BowyerWatson.h"

int main(void) {
    // Example usage
    double x[] = {0.0, 1.0, 0.5};
    double y[] = {0.0, 0.0, 1.0};


    size_t n = sizeof(x) / sizeof(x[0]);

    TriangularMesh *mesh = del2d(x, y, n);
    if (mesh) {
        printf("Mesh created with %zu vertices, %zu faces, and %zu half-edges.\n",
               mesh->vertex_count, mesh->face_count, mesh->half_edge_count);
        writeMeshToFile(mesh, "output_mesh.txt");
        freeMesh(mesh);
    } else {
        printf("Failed to create mesh.\n");
    }

    return 0;
}