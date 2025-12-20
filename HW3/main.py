from scipy.spatial import KDTree
from scipy.sparse.csgraph import minimum_spanning_tree
import numpy as np
import matplotlib.pyplot as plt
import plyfile as ply
import argparse as ap
import igl
from sklearn.neighbors import kneighbors_graph
from scipy.sparse import csr_matrix
from RegularFish import RegularFish


parser = ap.ArgumentParser()
parser.add_argument('--input', '-i', type=str, default='data/teapot.ply', help='Path to the input PLY file')
parser.add_argument('--k', '-k', type=int, default=10, help='Number of nearest neighbors to find')
args = parser.parse_args()
input_file = args.input
k = args.k

def load_point_cloud(input_file):
    with open(input_file, 'rb') as f:
        plydata = ply.PlyData.read(f)

    elements = plydata['vertex'].data
    points = np.array([[elements[i][0], elements[i][1], elements[i][2]] for i in range(len(elements))])
    kdtree = KDTree(points)
    return points, kdtree


def get_normal(point, points, tree, k):
    distances, indices = tree.query(point, k=k+1)
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

    # print("Visited:", np.sum(visited), "out of", Npoints) Normalement il passe partout (Fort Boyard)
    return oriented_normals

points, kdtree = load_point_cloud(input_file)
normals = np.array([get_normal(point, points, kdtree, k) for point in points])
normals = normals / np.linalg.norm(normals, axis=1, keepdims=True)

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

fish = RegularFish(points, kdtree, oriented_normals, h=5)
centers = fish.getCenters()
sol = fish.Solve()
sol /= np.max(np.abs(sol))




# Visualization
plt.figure(figsize=(10,10))
ax = plt.axes(projection='3d')
ax.set_aspect('equal')

# Scatter points
ax.scatter(points[:,0], points[:,1], points[:,2], color='k', s=1)

# Plot oriented normals
ax.quiver(points[:,0], points[:,1], points[:,2],
          oriented_normals[:,0], oriented_normals[:,1], oriented_normals[:,2],
          length=1, normalize=True, color='b', linewidth=0.5)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title('Oriented Normals using MST on kNN Graph')
plt.show()



"""
TODO:
- Homogéniser les normales: KNN graph. Il faut que les normales soient cohérentes entre voisines    DONE
- Utiliser igl pour orienter les normales                                                           IMPOSSIBLE
- Résoudre Poisson
    - Construire la grille régulière
    - Calculer l'approximation finie du Laplacien
    - Gaussian splatting des normales sur la grille
    - Résoudre le système linéaire creux (CG)
- Utiliser Marching Cubes pour extraire la surface
"""
