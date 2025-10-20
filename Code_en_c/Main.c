#include <stdio.h>
#include <stdlib.h>
#include "headers/BowyerWatson.h"

int main(void) {
    // Example usage
    double x[] = {0.0, 1.0, 0.25, 0.8, 0.5};
    double y[] = {0.0, 0.1, 1.0, 0.8, 0.6};


    size_t n = sizeof(x) / sizeof(x[0]);

    TriangularMesh *mesh = del2d(x, y, n);
    if (mesh) {
        printf("Mesh created with %zu vertices, %zu faces, and %zu half-edges.\n",
               mesh->vertex_count, mesh->face_count, mesh->half_edge_count);
        // print all half-edges
        HalfEdge *half_edges = mesh->half_edges;
        for (size_t i = 0; i < mesh->half_edge_count; ++i) {
            HalfEdge *he = &half_edges[i];
            printf("Half-edge %zu: (%f, %f) -> (%f, %f)\n",
                   i,
                   he->vertex->coord[0], he->vertex->coord[1],
                   he->next->vertex->coord[0], he->next->vertex->coord[1]);
        }
        writeMeshToFile(mesh, "output_mesh.txt");
        freeMesh(mesh);
    } else {
        printf("Failed to create mesh.\n");
    }

    return 0;
}