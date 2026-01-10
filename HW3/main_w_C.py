from scipy.spatial import KDTree
from scipy.sparse.csgraph import minimum_spanning_tree
import numpy as np
import matplotlib.pyplot as plt
import plyfile as ply
import argparse as ap
import igl
from sklearn.neighbors import kneighbors_graph, NearestNeighbors
import time
from scipy.sparse import csr_matrix
import poisson
#import marchingCubes.MC as MC
import ctypes
import skfmm
from plyfile import PlyData, PlyElement
import time


# Input Argument 
parser = ap.ArgumentParser()
parser.add_argument('--input', '-i', type=str, default='data/bunny.ply', help='Path to the input PLY file')
parser.add_argument('--k', '-k', type=int, default=10, help='Number of nearest neighbors to find')
args = parser.parse_args()
input_file = args.input
k = args.k

with open(input_file, 'rb') as f:
    plydata = ply.PlyData.read(f)

elements = plydata['vertex'].data
points = np.array([[elements[i][0], elements[i][1], elements[i][2]] for i in range(len(elements))])
kdtree = KDTree(points)

# Compute normals using PCA
def get_normal(point, points, tree, k):
    distances, indices = tree.query(point, k=k+1)
    neighbors = points[indices]
    centroid = np.mean(neighbors, axis=0)
    cov_matrix = np.cov((neighbors - centroid).T)
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
    return eigenvectors[:, np.argmin(eigenvalues)]

normals = np.array([get_normal(point, points, kdtree, k) for point in points])
#Normalize the normals
normals = normals / np.linalg.norm(normals, axis=1, keepdims=True) 
    
   
def weights(Graph, normals):
    raws, cols = Graph.nonzero()
    dots = np.sum(normals[raws] * normals[cols], axis=1)
    vals = 2 - np.abs(dots)
    W = csr_matrix((vals, (raws, cols)), shape=Graph.shape)
    return W

# Consistent orientation of normals 
def preorientNormals(graph, normals):
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

    print("Visited:", np.sum(visited), "out of", Npoints)
    return oriented_normals

# kNN graph
knnGraph = kneighbors_graph(points, n_neighbors=k, mode='connectivity', include_self=False, n_jobs=-1) # n_jobs=-1 to use all processors
# Symmetrize the graph
knnGraph = knnGraph.maximum(knnGraph.T)
knnGraph = knnGraph.tocsr()
# Poids entre les normales
W = weights(knnGraph, normals)
W = W.maximum(W.T)
W = W.tocsr()
# Minimum Spanning Tree
MST = minimum_spanning_tree(W)
MST = MST.maximum(MST.T)
# On redresse les normales
oriented_normals = preorientNormals(MST, normals)
oriented_normals /= np.linalg.norm(oriented_normals, axis=1, keepdims=True)

# from scipy.sparse.csgraph import connected_components
# n_components, labels = connected_components(MST)
# print("Nombre de composantes :", n_components)

def plot_normals(points, normals):
    plt.figure(figsize=(10,10))
    ax = plt.axes(projection='3d')
    ax.quiver(points[:,0], points[:,1], points[:,2],
            oriented_normals[:,0], oriented_normals[:,1], oriented_normals[:,2],
            length=10, normalize=True, color='b', linewidth=0.5)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title('Oriented Normals using MST on kNN Graph')
    plt.show()

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

# how to choose sigma? 
# Mean distance between the neighbors
sigma = estimate_sigma(points, k)
print("Sigma:", sigma)

start = time.time()
print("Solving Poisson equation...")
x,y,z,chi = poisson.solve_poisson(points, oriented_normals, sigma=sigma, tree=kdtree, N=400)
end = time.time()
print("Poisson equation solved in", end-start, "seconds")
#refit between 0 and 1
chi-=np.max(chi)
chi/=np.min(chi)
chi=skfmm.distance(chi-0.5, dx=(x[1]-x[0], y[1]-y[0], z[1]-z[0]))

#plot chi as a grid color points
def plot_chi_grid(x, y, z, chi):
    grid=np.meshgrid(x,y,z)
    plt.figure()
    ax = plt.axes(projection='3d')
    sc = ax.scatter(grid[0], grid[1], grid[2], c=chi.flatten(), cmap='Greys', s=0.01)
    plt.colorbar(sc)
    plt.title('Poisson Solution Values at Points')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.axis('equal')
    plt.show()

threshold = 0 #threshold for the isosurface

#import C library for marching cubes


lib = ctypes.CDLL("./marchingCubes/Mcc.dll")

lib.marching_cubes_grid.argtypes = [
    ctypes.POINTER(ctypes.c_float),
    ctypes.c_int, ctypes.c_int, ctypes.c_int,
    ctypes.c_float,
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]
lib.marching_cubes_grid.restype = ctypes.c_int


##################################################
# save triangles to PLY file
def save_triangles_to_ply(x, y, z, chi, name_file, threshold=0):
    

    start = time.time()

    # dimensions
    nx, ny, nz = chi.shape

    # chargement de la lib C
    lib = ctypes.CDLL("./marchingCubes/Mcc.dll")

    lib.marching_cubes_grid.argtypes = [
        ctypes.POINTER(ctypes.c_float),
        ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.c_float,
        ctypes.POINTER(ctypes.c_float),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.marching_cubes_grid.restype = ctypes.c_int

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

    print("Triangles extracted in", time.time()-start, "seconds")

    # découpe au bon nombre
    V = vertices[:nverts]
    F = faces[:ntri.value]

    # écriture PLY (INCHANGÉE)
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


start = time.time()
save_triangles_to_ply(x, y, z, chi, "test_avec400.ply", threshold)
end = time.time()
print("Saved triangles to ply file in", end-start, "seconds")
##################################################


# plot_isosurface_marching_cubes(x, y, z, chi, threshold)


"""
TODO:
- Homogéniser les normales: KNN graph. Il faut que les normales soient cohérentes entre voisines    DONE
- Utiliser igl pour orienter les normales                                                           IMPOSSIBLE
- Faire Poisson Surface Reconstruction avec les points et les normales
- Utiliser Marching Cubes pour extraire la surface
"""
