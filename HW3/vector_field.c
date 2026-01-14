#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void vector_field(
    int P,
    const double *points,
    const double *normals,
    const double *x, int nx,
    const double *y, int ny,
    const double *z, int nz,
    double sigma,
    float *V
)
{
    double h = x[1] - x[0];

    int r = (int)(3.0 * sigma / h);
    double inv2sigma2 = 1.0 / (2.0 * sigma * sigma);

    double *gx = malloc((2*r+1) * sizeof(double));
    double *gy = malloc((2*r+1) * sizeof(double));
    double *gz = malloc((2*r+1) * sizeof(double));

    for (int p = 0; p < P; ++p)
    {
        double px = points[3*p + 0];
        double py = points[3*p + 1];
        double pz = points[3*p + 2];

        double nxp = normals[3*p + 0];
        double nyp = normals[3*p + 1];
        double nzp = normals[3*p + 2];

        int ix = (int)((px - x[0]) / h);
        int iy = (int)((py - y[0]) / h);
        int iz = (int)((pz - z[0]) / h);

        int i0 = ix - r; if (i0 < 0) i0 = 0;
        int i1 = ix + r; if (i1 >= nx) i1 = nx - 1;

        int j0 = iy - r; if (j0 < 0) j0 = 0;
        int j1 = iy + r; if (j1 >= ny) j1 = ny - 1;

        int k0 = iz - r; if (k0 < 0) k0 = 0;
        int k1 = iz + r; if (k1 >= nz) k1 = nz - 1;

        for (int i = i0; i <= i1; ++i) {
            double dx = x[i] - px;
            gx[i-i0] = exp(-dx*dx * inv2sigma2);
        }
        for (int j = j0; j <= j1; ++j) {
            double dy = y[j] - py;
            gy[j-j0] = exp(-dy*dy * inv2sigma2);
        }
        for (int k = k0; k <= k1; ++k) {
            double dz = z[k] - pz;
            gz[k-k0] = exp(-dz*dz * inv2sigma2);
        }

        for (int i = i0; i <= i1; ++i){
            for (int j = j0; j <= j1; ++j){
                for (int k = k0; k <= k1; ++k) {
                    double w = gx[i-i0] * gy[j-j0] * gz[k-k0];
                    int idx = ((i*ny + j)*nz + k) * 3;
                    V[idx]     += w * nxp;
                    V[idx + 1] += w * nyp;
                    V[idx + 2] += w * nzp;
                }
            }
        }
    }
    free(gx); free(gy); free(gz);
}
