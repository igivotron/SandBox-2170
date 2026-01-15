### IMPORTS ###
from scipy.spatial import KDTree
from scipy.sparse.csgraph import minimum_spanning_tree
import numpy as np
import matplotlib.pyplot as plt
import plyfile as ply
import argparse as ap
from sklearn.neighbors import kneighbors_graph, NearestNeighbors
import time
from scipy.sparse import csr_matrix
import poisson
import ctypes
import skfmm
from plyfile import PlyData, PlyElement
import time

start_all = time.time()

#import C library for marching cubes
lib = ctypes.CDLL("./shared_lib/Mcc.dll")

lib.marching_cubes_grid.argtypes = [
    ctypes.POINTER(ctypes.c_float),
    ctypes.c_int, ctypes.c_int, ctypes.c_int,
    ctypes.c_float,
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]
lib.marching_cubes_grid.restype = ctypes.c_int


### INPUTS ###
# Input Argument 
parser = ap.ArgumentParser()
parser.add_argument('--input', '-i', type=str, default='data/bunny.ply', help='Path to the input PLY file')
parser.add_argument('--output', '-o', type=str, default='output/reconstructed.ply', help='Path to the output PLY file')
parser.add_argument('--k', '-k', type=int, default=10, help='Number of nearest neighbors to find')
parser.add_argument('--N', '-N', type=int, default=400, help='Grid size for Poisson reconstruction')
parser.add_argument('--skfmm','-skfmm', type=int, default=0, help='Use skfmm for distance computation')

args = parser.parse_args()
input_file = args.input
k = args.k
N = args.N
UseSkfmm = args.skfmm

with open(input_file, 'rb') as f: plydata = ply.PlyData.read(f)

elements = plydata['vertex'].data
points = np.array([[elements[i][0], elements[i][1], elements[i][2]] for i in range(len(elements))])
kdtree = KDTree(points)

### FUNCTIONS ###

# Compute normals using PCA
def get_normal(point, points, tree, k):
    _, indices = tree.query(point, k=k+1)
    neighbors = points[indices]
    centroid = np.mean(neighbors, axis=0)
    cov_matrix = np.cov((neighbors - centroid).T)
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
    return eigenvectors[:, np.argmin(eigenvalues)]
    
   
def weights(Graph, normals):
    raws, cols = Graph.nonzero()
    dots = np.sum(normals[raws] * normals[cols], axis=1)
    vals = 2 - np.abs(dots)
    W = csr_matrix((vals, (raws, cols)), shape=Graph.shape)
    return W

# Consistent orientation of normals 
def orientNormals(graph, normals):
    Npoints = normals.shape[0]
    visited = np.zeros(Npoints)
    oriented_normals = normals.copy()
    def dfs(node):
        visited[node] = 1
        neighbors = graph[node].nonzero()[1]
        for neighbor in neighbors:
            if not visited[neighbor]:
                dot = np.dot(oriented_normals[node], oriented_normals[neighbor])
                if dot < 0:
                    oriented_normals[neighbor] = -oriented_normals[neighbor]
                dfs(neighbor)
    
    dfs(0)

    if np.sum(visited) != Npoints: print("Warning: The graph is not fully connected. Some normals may remain unoriented.")
    return oriented_normals

def estimate_sigma(points, k, c=1.5):
    """
    Estime automatiquement σ (à quel point un point influence son voisinage)
    """
    # Crée l'objet kNN
    nbrs = NearestNeighbors(n_neighbors=k+1).fit(points) #k+1 car le point lui-même est inclus
    # Trouve les k+1 plus proches voisins pour chaque point
    dists, _ = nbrs.kneighbors(points)
    # distance au k-ième voisin réel
    dk = dists[:, -1]
    # valeur médiane sur tous les points
    d_median = np.median(dk)
    return c * d_median

start_normal = time.time()
print("Computing normals and orienting them...")
# kNN graph
knnGraph = kneighbors_graph(points, n_neighbors=k, mode='connectivity', include_self=False, n_jobs=-1) # n_jobs=-1 to use all processors
knnGraph = knnGraph.maximum(knnGraph.T)
knnGraph = knnGraph.tocsr()

# Normals
normals = np.array([get_normal(point, points, kdtree, k) for point in points])
normals = normals / np.linalg.norm(normals, axis=1, keepdims=True) 
W = weights(knnGraph, normals)
W = W.maximum(W.T)
W = W.tocsr()

# Minimum Spanning Tree
MST = minimum_spanning_tree(W)
MST = MST.maximum(MST.T)

# Normal orientation
oriented_normals = orientNormals(MST, normals)
oriented_normals /= np.linalg.norm(oriented_normals, axis=1, keepdims=True)
tnorm = time.time() - start_normal
print("Normals computed and oriented in", time.time()-start_normal, "seconds.", "Total time:", time.time()-start_all, "seconds")

# Sigma estimation
sigma = estimate_sigma(points, k)
# print("Sigma:", sigma)

### POISSON SURFACE RECONSTRUCTION ###
print("Starting Poisson surface reconstruction...")
start_poisson = time.time()
x,y,z,chi = poisson.solve_poisson(points, oriented_normals, sigma=sigma, tree=kdtree, N=N)
#refit between 0 and 1
chi-=np.max(chi)
chi/=np.min(chi)

threshold = 0.5 #threshold for the isosurface

if UseSkfmm:
    start_skfmm = time.time()
    chi=skfmm.distance(chi-threshold, dx=(x[1]-x[0], y[1]-y[0], z[1]-z[0]))
    print("\t Distance function computed in", time.time()-start_skfmm, "seconds", "Total time:", time.time()-start_all, "seconds")
    threshold = 0.0
tpoisson = time.time() - start_poisson
print("Poisson surface reconstruction completed in", time.time()-start_poisson, "seconds", "Total time:", time.time()-start_all, "seconds")


##################################################
# save triangles to PLY file
def save_triangles_to_ply(x, y, z, chi, name_file, threshold=0):
    # dimensions
    nx, ny, nz = chi.shape

    # allocation max (5 triangles par cube)
    max_tris = (nx-1)*(ny-1)*(nz-1)*5

    vertices = np.zeros((max_tris*3, 3), dtype=np.float32)
    faces = np.zeros((max_tris, 3), dtype=np.int32)
    ntri = ctypes.c_int()

    # appel C (REMPLACE TOUTE LA BOUCLE PYTHON)
    nverts = lib.marching_cubes_grid(
        chi.astype(np.float32).ravel().ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        nx, ny, nz,
        threshold,
        vertices.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        faces.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.byref(ntri)
    )

    # découpe au bon nombre
    V = vertices[:nverts]
    F = faces[:ntri.value]

    # écriture PLY
    vertices_ply = np.array(
        [(v[0], v[1], v[2]) for v in V],
        dtype=[('x','f4'), ('y','f4'), ('z','f4')]
    )
    faces_ply = np.array(
        [(tuple(face),) for face in F],
        dtype=[('vertex_indices', 'i4', (3,))]
    )

    PlyData([
        PlyElement.describe(vertices_ply, 'vertex'),
        PlyElement.describe(faces_ply, 'face')
    ]).write(name_file)

    print("Saved", len(F), "triangles to", name_file)


start_save = time.time()
print("Saving triangles to ply file...")
save_triangles_to_ply(x, y, z, chi, args.output, threshold)
print("Saved triangles to ply file in", time.time()-start_save, "seconds", "Total time:", time.time()-start_all, "seconds")
tsave = time.time() - start_save

with open("plots/timings.txt", "a") as f:
    f.write(f"{N};{tnorm:.6f};{tpoisson:.6f};{tsave:.6f};{time.time()-start_all:.26f}\n")

# plot_isosurface_marching_cubes(x, y, z, chi, threshold)

