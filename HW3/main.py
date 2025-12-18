from scipy.spatial import KDTree
from scipy.sparse.csgraph import minimum_spanning_tree
import numpy as np
import matplotlib.pyplot as plt
import plyfile as ply
import argparse as ap
import igl
from sklearn.neighbors import kneighbors_graph
import time
from scipy.sparse import csr_matrix


parser = ap.ArgumentParser()
parser.add_argument('--input', '-i', type=str, default='data/teapot.ply', help='Path to the input PLY file')
parser.add_argument('--k', '-k', type=int, default=5, help='Number of nearest neighbors to find')
args = parser.parse_args()
input_file = args.input
k = args.k

with open(input_file, 'rb') as f:
    plydata = ply.PlyData.read(f)

elements = plydata['vertex'].data
points = np.array([[elements[i][0], elements[i][1], elements[i][2]] for i in range(len(elements))])
kdtree = KDTree(points)


def get_normal(point, points, tree, k):
    distances, indices = tree.query(point, k=k+1)
    neighbors = points[indices]
    centroid = np.mean(neighbors, axis=0)
    cov_matrix = np.cov((neighbors - centroid).T)
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
    return eigenvectors[:, np.argmin(eigenvalues)]

normals = np.array([get_normal(point, points, kdtree, k) for point in points])
normals = normals / np.linalg.norm(normals, axis=1, keepdims=True)
    
def weights(Graph, normals):
    raws, cols = Graph.nonzero()
    dots = np.sum(normals[raws] * normals[cols], axis=1)
    vals = 1 - np.abs(dots)
    W = csr_matrix((vals, (raws, cols)), shape=Graph.shape)
    return W

def orientNormals(graph, normals):
    Npoints = normals.shape[0]
    visited = np.zeros(Npoints, dtype=bool)
    oriented_normals = normals.copy()
    def dfs(node):
        visited[node] = True
        neighbors = graph[node].nonzero()[1]
        for neighbor in neighbors:
            if not visited[neighbor]:
                dot_product = np.dot(oriented_normals[node], oriented_normals[neighbor])
                if dot_product < 0:
                    oriented_normals[neighbor] = -oriented_normals[neighbor]
                dfs(neighbor)
    for i in range(Npoints):
        if not visited[i]:
            print("Starting DFS at node:", i)
            dfs(i)
    return oriented_normals

# kNN graph
knnGraph = kneighbors_graph(points, n_neighbors=k, mode='connectivity', include_self=False, n_jobs=-1) # n_jobs=-1 to use all processors
# Symmetrize the graph
knnGraph = knnGraph.maximum(knnGraph.T)
# Passe la matrix en CSR pour l'efficacité
knnGraph = knnGraph.tocsr()
# Poids entre les normales
W = weights(knnGraph, normals)
# Minimum Spanning Tree
MST = minimum_spanning_tree(W)
# On redresse les normales
oriented_normals = orientNormals(MST, normals)

plt.figure(figsize=(10,10))
ax = plt.axes(projection='3d')
ax.quiver(points[:,0], points[:,1], points[:,2],
          oriented_normals[:,0], oriented_normals[:,1], oriented_normals[:,2],
          length=5, normalize=True)

ax.quiver(points[:,0], points[:,1], points[:,2],
          normals[:,0], normals[:,1], normals[:,2],
          length=5, color='r', normalize=True)
plt.title('Oriented Normals using MST on kNN Graph')
plt.show()

"""
TODO:
- Homogéniser les normales: KNN graph. Il faut que les normales soient cohérentes entre voisines
- Utiliser igl pour orienter les normales
- Faire Poisson Surface Reconstruction avec les points et les normales
- Utiliser Marching Cubes pour extraire la surface
"""
