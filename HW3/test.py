import scipy as sp
import numpy as np
from scipy.sparse.linalg import LaplacianNd
import matplotlib.pyplot as plt
import time


def Ig(nx, ny, nz, h):
    return LaplacianNd((nx, ny, nz), boundary_conditions='dirichlet')

def Lio(nx, ny, nz, h):
    nNodes = nx * ny * nz
    A=np.zeros((nNodes, nNodes))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                idx = i*ny*nz + j*nz + k
                A[idx, idx] = -6
                # neighbors in 6 directions
                for offset in [-1, 1]:
                    if 0 <= i + offset < nx:
                        neighbor_idx = (i + offset)*ny*nz + j*nz + k
                        A[idx, neighbor_idx] = 1
                        A[neighbor_idx, idx] = 1
                    if 0 <= j + offset < ny:
                        neighbor_idx = i*ny*nz + (j + offset)*nz + k
                        A[idx, neighbor_idx] = 1
                        A[neighbor_idx, idx] = 1
                    if 0 <= k + offset < nz:
                        neighbor_idx = i*ny*nz + j*nz + (k + offset)
                        A[idx, neighbor_idx] = 1
                        A[neighbor_idx, idx] = 1
    A /= h**2

    return A

nx = 5
ny = 5
nz = 5
h = 1.0

startlio = time.time()
A = Lio(nx, ny, nz, h)
print("Laplacian matrix A shape:", A.shape)
endlio = time.time()

startig = time.time()
B = Ig(nx, ny, nz, h)
print("LaplacianNd matrix B shape:", B.shape)
endig = time.time()

print(f"Lio computation time: {(endlio - startlio)*1e3} ms")
print(f"Ig computation time: {(endig - startig)*1e3} ms")

print("Difference between A and B:", np.linalg.norm(A - B.toarray()))
plt.subplots(1, 2, figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.title("Lio Laplacian")
plt.imshow(A, cmap='viridis')
plt.colorbar()
plt.subplot(1, 2, 2)
plt.title("Ig Laplacian")
plt.imshow(B.toarray(), cmap='viridis')
plt.colorbar()
plt.show()