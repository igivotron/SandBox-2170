#include "Hilbert.h"
#include "stdlib.h"
#include "stdio.h"

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
// open ../inputs/100pts.txt
int main() {
    int N = 100;
    FILE *file = fopen("../inputs/100pts", "r");
    double points[N][2];
    for (int i = 0; i < N; i++) {
        if (fscanf(file, "%lf %lf", &points[i][0], &points[i][1]) != 2) {
            perror("Error reading point");
            fclose(file);
            return EXIT_FAILURE;
        }
    }
    fclose(file);

    const int depth = 100;
    int hilbertCoords[N][depth];
    for (int i = 0; i < N; i++) HilbertCoord(points[i][0], points[i][1], 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, depth, hilbertCoords[i]);


    FILE *outputFile = fopen("../outputs/h100.txt", "w");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < depth; j++) {
            fprintf(outputFile, "%d ", hilbertCoords[i][j]);
        }
        fprintf(outputFile, "\n");
    }
    fclose(outputFile);

}
