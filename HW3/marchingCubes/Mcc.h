



typedef struct {
    double x, y, z;
} XYZ;

typedef struct {
   XYZ p[3];
} TRIANGLE;

typedef struct {
   XYZ p[8];
   double val[8];
} GRIDCELL;

int marching_cubes_grid(float* field, int nx, int ny, int nz,float iso,float* vertices,int* triangles, int* ntriangles);
int Polygonise(GRIDCELL grid,double isolevel,TRIANGLE *triangles);
XYZ VertexInterp(double isolevel, XYZ p1, XYZ p2, double valp1, double valp2);