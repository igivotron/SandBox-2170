import scipy as sp
import numpy as np
from scipy.sparse.linalg import LaplacianNd
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
from scipy.sparse.linalg import cg


class RegularFish(object):
    def __init__(self, points, tree, normals, h):
        self.points = points
        self.tree = tree
        self.normals = normals
        self.h = h
        self.dimensions = self.getDimensions()

    def getDimensions(self):
        mult = 1
        x_min = np.min(self.points[:, 0]) - self.h * mult
        x_max = np.max(self.points[:, 0]) + self.h * mult
        y_min = np.min(self.points[:, 1]) - self.h * mult
        y_max = np.max(self.points[:, 1]) + self.h * mult
        z_min = np.min(self.points[:, 2]) - self.h * mult
        z_max = np.max(self.points[:, 2]) + self.h * mult

        self.nx = int( np.ceil((x_max - x_min) / self.h) )
        self.ny = int( np.ceil((y_max - y_min) / self.h) )
        self.nz = int( np.ceil((z_max - z_min) / self.h) )
        return (x_min, x_max, y_min, y_max, z_min, z_max)
    
    def getGrid(self):
        # grid points based on dimensions and h
        x_min, x_max, y_min, y_max, z_min, z_max = self.dimensions
        x = np.linspace(x_min, x_max, self.nx)
        y = np.linspace(y_min, y_max, self.ny)
        z = np.linspace(z_min, z_max, self.nz)
        return np.meshgrid(x, y, z, indexing='ij')

    
    def getCenters(self):
        # Compute cell centers
        x_min, x_max, y_min, y_max, z_min, z_max = self.dimensions
        x_centers = np.linspace(x_min + self.h / 2, x_max - self.h / 2, self.nx - 1)
        y_centers = np.linspace(y_min + self.h / 2, y_max - self.h / 2, self.ny - 1)
        z_centers = np.linspace(z_min + self.h / 2, z_max - self.h / 2, self.nz - 1)
        return np.meshgrid(x_centers, y_centers, z_centers, indexing='ij')
    

    def getA(self):
        # Construct the Laplacian matrix for the grid
        L = LaplacianNd((self.nx-1, self.ny-1, self.nz-1), boundary_conditions='dirichlet')
        L = L.tosparse().tocsr()
        return L / (self.h ** 2)
    

    def getVectorField(self, sigma):
        Xc, Yc, Zc = self.getCenters()
        nx, ny, nz = Xc.shape
        nVertex = nx * ny * nz

        centers = np.vstack((Xc.ravel(), Yc.ravel(), Zc.ravel())).T
        tree = self.tree

        Vx = np.zeros(nVertex)
        Vy = np.zeros(nVertex)
        Vz = np.zeros(nVertex)

        for c_idx, c in enumerate(centers):
            idxs = tree.query_ball_point(c, r=3*sigma)
            if len(idxs) == 0:
                continue

            diff = self.points[idxs] - c
            dist2 = np.sum(diff**2, axis=1)
            weights = np.exp(-dist2 / (2 * sigma**2))

            wsum = np.sum(weights)
            if wsum == 0:
                continue

            Vx[c_idx] = np.sum(weights * self.normals[idxs, 0]) / wsum
            Vy[c_idx] = np.sum(weights * self.normals[idxs, 1]) / wsum
            Vz[c_idx] = np.sum(weights * self.normals[idxs, 2]) / wsum

        Nx = csr_matrix(Vx).T
        Ny = csr_matrix(Vy).T
        Nz = csr_matrix(Vz).T

        return Nx, Ny, Nz
    
    def Solve(self):
        A = self.getA()
        Vx, Vy, Vz = self.getVectorField(sigma=1.0)

        dVVxdx = np.gradient(Vx.toarray().reshape((self.nx-1, self.ny-1, self.nz-1)), self.h, axis=0).ravel()
        dVVydy = np.gradient(Vy.toarray().reshape((self.nx-1, self.ny-1, self.nz-1)), self.h, axis=1).ravel()
        dVVzdz = np.gradient(Vz.toarray().reshape((self.nx-1, self.ny-1, self.nz-1)), self.h, axis=2).ravel()

        b = -(dVVxdx + dVVydy + dVVzdz)

        X, info = cg(A, b)
        print("CG info:", info)
        return X
    

# # sphere example
# num_points = 1000
# phi = np.random.uniform(0, np.pi, num_points)
# theta = np.random.uniform(0, 2 * np.pi, num_points)
# x = np.sin(phi) * np.cos(theta)
# y = np.sin(phi) * np.sin(theta)
# z = np.cos(phi)
# points = np.vstack((x, y, z)).T
# normals = points.copy() 

# tree = sp.spatial.KDTree(points)
# h = 0.1
# fish = RegularFish(points, tree, normals, h)

# centers = fish.getCenters()
# grid = fish.getGrid()

# sol = fish.Solve()

# sol /= np.max(np.abs(sol))

# plt.figure(figsize=(8,8))
# ax = plt.axes(projection='3d')
# ax.set_aspect('equal')
# ax.scatter(points[:,0], points[:,1], points[:,2], color='k', s=1)
# ax.scatter(centers[0].ravel()[sol>0.5], centers[1].ravel()[sol>0.5], centers[2].ravel()[sol>0.5], color='r', s=1)
# plt.show()




