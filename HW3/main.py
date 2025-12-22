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
import poisson
import marchingCubes.MC as MC
import skfmm
from plyfile import PlyData, PlyElement
import time


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
# remove duplicate points
points = np.unique(points, axis=0)
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

    print("Visited:", np.sum(visited), "out of", Npoints)
    return oriented_normals

    
        

# def orientNormals(normals, points, k):
#     N = len(points)
#     eps = 1e-3 * np.mean(np.linalg.norm(points - np.mean(points, axis=0), axis=1))
#     Pplus = points + eps * normals
#     Pminus = points - eps * normals
#     Wplus = np.zeros(N)
#     Wminus = np.zeros(N)
#     igl.fast_winding_number_for_points(points, points, normals, Pplus, Wplus)
#     igl.fast_winding_number_for_points(points, points, normals, Pminus, Wminus)
#     oriented_normals = normals.copy()
#     for i in range(N):
#         if Wminus[i] > Wplus[i]:
#             oriented_normals[i] = -oriented_normals[i]
#     return oriented_normals

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


# how to choose sigma? density of data points
sigma = (np.max(points, axis=0) - np.min(points, axis=0)) 
sigma = np.cbrt(sigma[0]*sigma[1]*sigma[2]/len(points))
print("Sigma:", sigma)

start = time.time()
print("Solving Poisson equation...")
x,y,z,chi = poisson.solve_poisson(points, oriented_normals, sigma=sigma, tree=kdtree, N=100)
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

##################################################
# save triangles to PLY file
def save_triangles_to_ply(x, y, z, chi, name_file, threshold=0):
    triangles = []
    threshold = 0
    for i in range(len(x)-1):
            for j in range(len(y)-1):
                for k in range(len(z)-1):
                    cube=MC.Cube([i,j,k], 1, [
                        chi[i][j][k],
                        chi[i+1][j][k],
                        chi[i][j+1][k],
                        chi[i+1][j+1][k],
                        chi[i][j][k+1],
                        chi[i+1][j][k+1],
                        chi[i][j+1][k+1],
                        chi[i+1][j+1][k+1],
                    ])
                    tris = cube.getTriangles(threshold=threshold)
                    for tri in tris:
                        triangles.append(tri)
    vertices = []
    vertex_index = {}
    faces= []

    for tri in triangles:
        face = []
        for v in tri:
            v_tuple = (v[0], v[1], v[2])
            if v_tuple not in vertex_index:
                vertex_index[v_tuple] = len(vertices)
                vertices.append(v)
            face.append(vertex_index[v_tuple])
        faces.append(face)
    vertices_ply = np.array(   [(v[0], v[1], v[2]) for v in vertices]  ,  dtype=[('x','f4'), ('y','f4'), ('z','f4')] )
    faces_ply = np.array(   [(face,) for face in faces]  ,  dtype=[('vertex_indices', 'i4', (3,))]     )
    PlyData([ PlyElement.describe(vertices_ply, 'vertex') , PlyElement.describe(faces_ply, 'face')  ]).write(name_file)
start = time.time()
save_triangles_to_ply(x, y, z, chi, "triangles.ply", threshold)
end = time.time()
print("Saved triangles to ply file in", end-start, "seconds")
##################################################

def plot_isosurface_marching_cubes(x, y, z, chi, threshold):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    def draw_cubes(ax, x, y, z, chi, threshold):
        for i in range(len(x)-1):
            for j in range(len(y)-1):
                for k in range(len(z)-1):
                    cube=MC.Cube([i,j,k], 1, [
                        chi[i][j][k],
                        chi[i+1][j][k],
                        chi[i][j+1][k],
                        chi[i+1][j+1][k],
                        chi[i][j][k+1],
                        chi[i+1][j][k+1],
                        chi[i][j+1][k+1],
                        chi[i+1][j+1][k+1],
                    ])
                    cube.draw_without_vertices(ax,alpha=0.3, threshold=threshold)
    draw_cubes(ax, x, y, z, chi, threshold)
    plt.axis('equal')

    # interactive threshold adjustment
    threshold_text = ax.text2D( #afficher le threshold
        0.02, 0.95,
        f"Threshold = {threshold:.2f}",
        transform=ax.transAxes,
        fontsize=12
    )
    scroll_step = (np.max(chi) - np.min(chi)) / 100
    def on_key(event, ax, fig):
        ax.cla()
        global threshold
        threshold += event.step * scroll_step
        threshold = max(-1, min(1, threshold))
        ax.text2D(
            0.02, 0.95,
            f"Threshold = {threshold:.2f}",
            transform=ax.transAxes,
            fontsize=12
        )
        draw_cubes(ax, x, y, z, chi, threshold)
        fig.canvas.draw_idle()
    fig.canvas.mpl_connect(
        'scroll_event',
        lambda event: on_key(event, ax, fig))

    plt.show()

# plot_isosurface_marching_cubes(x, y, z, chi, threshold)


"""
TODO:
- Homogéniser les normales: KNN graph. Il faut que les normales soient cohérentes entre voisines    DONE
- Utiliser igl pour orienter les normales                                                           IMPOSSIBLE
- Faire Poisson Surface Reconstruction avec les points et les normales
- Utiliser Marching Cubes pour extraire la surface
"""
