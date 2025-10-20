#include "Hilbert.h"
#include "stdlib.h"
#include "stdio.h"

int main(int argc,char *argv[]){ 
    int N = 0;
    char inputFileName[256];
    char outputFileName[256];
    snprintf(inputFileName, sizeof(inputFileName), "../inputs/%s", argv[1]);
    snprintf(outputFileName, sizeof(outputFileName), "../outputs/%s", argv[2]);

    FILE *file = fopen(inputFileName, "r");
    
    if (fscanf(file, "%d", &N) != 1) {
        fclose(file);
        return EXIT_FAILURE;
    }

    double points[N][2];
    for (int i = 0; i < N; i++) {
        if (fscanf(file, "%lf %lf", &points[i][0], &points[i][1]) != 2) {
            fclose(file);
            return EXIT_FAILURE;
        }
    }
    fclose(file);

    int hilbertCoords[N][DEPTH];
    for (int i = 0; i < N; i++) HilbertCoord(points[i][0], points[i][1], 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, DEPTH, hilbertCoords[i]);


    FILE *outputFile = fopen(outputFileName, "w");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < DEPTH; j++) fprintf(outputFile, "%d ", hilbertCoords[i][j]);
        fprintf(outputFile, "\n");
    }
    fclose(outputFile);
    return 0;
}
