import numpy as np

#to improve:
# 1. keep the symmetry of the laplacian matrix while imposing boundary conditions
# 2. use conjugate gradient method to solve the sparse linear system

def vector_field(x, y, z, normals, sigma, tree, points):
    """
    x,y,z : 1D arrays defining the grid points
    normals : Nx3 array of normal vectors at data points
    tree : KDTree built from data points
    """
    vec_field=np.zeros((len(x),len(y),len(z), 3))
    for i in range(len(x)):
        for j in range(len(y)):
            for k in range(len(z)):
                p = np.array([x[i], y[j], z[k]])
                indices = tree.query_ball_point(p, r=sigma*3)
                distances = np.linalg.norm(points[indices] - p, axis=1)
                weights = np.exp(- (distances ** 2) / (2 * sigma ** 2))
                vec_field[i, j, k] = np.sum((normals[indices].T * weights), axis=1)
    return vec_field

def solve_poisson(points, normals, sigma, tree, h=None):
    min_point = np.min(points, axis=0)
    max_point = np.max(points, axis=0)
    grid_size = max_point - min_point
    if h is None:
        h=min(grid_size)/20
    if h==0:
        h=0.1
    x = np.arange(min_point[0]-h, max_point[0] + 2*h, h)
    nx=len(x)
    y = np.arange(min_point[1]-h, max_point[1] + 2*h, h)
    ny=len(y)
    z = np.arange(min_point[2]-h, max_point[2] + 2*h, h)
    nz=len(z)
    vec_field=vector_field(x,y,z, normals, sigma, tree, points)
    # create laplacian operator
    nNodes=nx*ny*nz
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
    #create right-hand side
    dVxdx=np.gradient(vec_field[:,:,:,0], h, axis=0)
    dVydy=np.gradient(vec_field[:,:,:,1], h, axis=1)
    dVzdz=np.gradient(vec_field[:,:,:,2], h, axis=2)
    div_field=dVxdx + dVydy + dVzdz
    b = -div_field.flatten()
    #boundary conditions (Dirichlet: phi=0 at boundary)
    for i in range(len(x)):
        for j in range(len(y)):
            for k in range(len(z)):
                if i == 0 or i == len(x)-1 or j == 0 or j == len(y)-1 or k == 0 or k == len(z)-1:
                    idx = i*len(y)*len(z) + j*len(z) + k
                    A[idx, :] = 0
                    A[idx, idx] = 1
                    b[idx] = 0
    #solve Poisson equation
    phi = np.linalg.solve(A, b)
    return x,y,z,phi.reshape((nx, ny, nz))