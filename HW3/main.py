from scipy.spatial import KDTree
import numpy as np
import matplotlib.pyplot as plt
import plyfile as ply
import argparse as ap

parser = ap.ArgumentParser()
parser.add_argument('--input', type=str, default='data/bunny.ply', help='Path to the input PLY file')
parser.add_argument('--k', type=int, default=5, help='Number of nearest neighbors to find')
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
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(points[:, 0], points[:, 1], points[:, 2],
          normals[:, 0], normals[:, 1], normals[:, 2], length=1e-1, normalize=True)
ax.set_title('Point Cloud with Normals')
plt.show()

# plt.figure()
# ax = plt.axes(projection='3d')
# ax.scatter3D(points[:, 0], points[:, 1], points[:, 2], s=1)
# ax.set_title('3D Point Cloud')
# plt.show()