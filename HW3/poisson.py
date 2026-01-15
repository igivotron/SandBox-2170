import numpy as np
import scipy as sp
from scipy.sparse.linalg import LaplacianNd
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import cg
import time
import os
import ctypes
import pyamg
from scipy.fft import dstn, idstn

lib = ctypes.CDLL(os.path.abspath("./shared_lib/vector_field.dll"))
lib.vector_field.argtypes = [ctypes.c_int, 
                             ctypes.POINTER(ctypes.c_double), 
                             ctypes.POINTER(ctypes.c_double), 
                             ctypes.POINTER(ctypes.c_double), 
                             ctypes.c_int,
                             ctypes.POINTER(ctypes.c_double), 
                             ctypes.c_int,
                             ctypes.POINTER(ctypes.c_double), 
                             ctypes.c_int,
                             ctypes.c_double,
                             ctypes.POINTER(ctypes.c_float)
                            ]
lib.vector_field.restype = None


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

def solve_poisson(points, normals, sigma, tree, N=None):
    min_point = np.min(points, axis=0)
    max_point = np.max(points, axis=0)
    size = np.max(max_point - min_point)*1.15
    middle = (max_point + min_point)/2
    if N is None:
        N=10
    h=size/(N-1)
    x = np.linspace(middle[0]-size/2, middle[0]+size/2, N)
    y = np.linspace(middle[1]-size/2, middle[1]+size/2, N)
    z = np.linspace(middle[2]-size/2, middle[2]+size/2, N)
    x_c=np.ascontiguousarray(x, dtype=np.float64)
    y_c=np.ascontiguousarray(y, dtype=np.float64)
    z_c=np.ascontiguousarray(z, dtype=np.float64)
    
    start_field = time.time()
    vec_field = np.zeros((N, N, N, 3), dtype=np.float32, order='C') 
    points_c = np.ascontiguousarray(points, dtype=np.float64)
    normals_c = np.ascontiguousarray(normals, dtype=np.float64)   
    lib.vector_field(ctypes.c_int(len(points)),
                     points_c.flatten().ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                     normals_c.flatten().ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                     x_c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                     ctypes.c_int(N),
                     y_c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                     ctypes.c_int(N),
                     z_c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                     ctypes.c_int(N),
                     ctypes.c_double(sigma),
                     vec_field.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
    print("\t Vector field computed in", time.time()-start_field, "seconds")
    
    #create right-hand side
    dVxdx=np.gradient(vec_field[:,:,:,0], h, axis=0)
    dVydy=np.gradient(vec_field[:,:,:,1], h, axis=1)
    dVzdz=np.gradient(vec_field[:,:,:,2], h, axis=2)
    div_field=dVxdx + dVydy + dVzdz
    b = (-div_field).ravel()
    
    #solve Poisson equation
    start_solve = time.time()
    rhs = b.reshape((N,N,N))

    rhs_hat = dstn(rhs, type=1)
    kx = np.arange(1, N+1)
    ky = np.arange(1, N+1)
    kz = np.arange(1, N+1)

    denom = (
        (2*np.cos(np.pi*kx/(N+1))-2)[:,None,None] +
        (2*np.cos(np.pi*ky/(N+1))-2)[None,:,None] +
        (2*np.cos(np.pi*kz/(N+1))-2)[None,None,:]
    ) / h**2

    phi = idstn(rhs_hat / denom, type=1)
    print("\t Poisson equation solved in", time.time()-start_solve, "seconds")
    
    return x,y,z,phi.reshape((N, N, N))*h*h