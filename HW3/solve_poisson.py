import numpy as np
import matplotlib.pyplot as plt
import plyfile as ply
import argparse as ap
from scipy.spatial import KDTree
from itertools import combinations, product

class Cell():
    def __init__(self, index, center,h):
        self.index = index
        self.center = center
        self.size=h
        self.field = None

def create_grid(input_points):
    points = input_points.copy()
    min_point = np.min(points, axis=0)
    max_point = np.max(points, axis=0)
    grid_size = max_point - min_point
    h=min(grid_size)/2
    grid=[]
    xs = np.arange(min_point[0], max_point[0] + h, h)
    ys = np.arange(min_point[1], max_point[1] + h, h)
    zs = np.arange(min_point[2], max_point[2] + h, h)
    print(xs,ys,zs)
    for x in xs:
        for y in ys:
            for z in zs:
                center = np.array([x, y, z])
                in_points = find_in_points(center, points, h)
                create_cell(in_points,grid,h,center)
    return grid

def find_in_points(center, points, h):
    mask = np.all((points >= center - h/2) & (points <= center + h/2), axis=1)
    return points[mask]

def create_cell(points,grid,h,center):
    if len(points)<=1:
        grid.append(Cell(len(grid),center,h))
    else:
        offsets = np.array([-h/4, h/4])
        for dx in offsets:
            for dy in offsets:
                for dz in offsets:
                    new_center = center + np.array([dx, dy, dz])
                    in_points = find_in_points(new_center, points, h/2)
                    create_cell(in_points,grid,h/2,new_center)

def vector_field(grid, normals, sigma, tree):
    for p in grid:
        indices = tree.query_ball_point(p.center, r=sigma*3)
        distances = np.linalg.norm(points[indices] - p.center, axis=1)
        weights = np.exp(- (distances ** 2) / (2 * sigma ** 2))
        p.field = np.sum((normals[indices].T * weights), axis=1)
    return

def visualize_vector_field(grid):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for cell in grid:
        if cell.field is not None:
            ax.quiver(cell.center[0], cell.center[1], cell.center[2],
                      cell.field[0], cell.field[1], cell.field[2],
                      length=0.1, color='r')
    plt.show()

def visualise_grid(grid):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(points[:,0], points[:,1], points[:,2], color='b', s=0.3)
    for cell in grid:
        #draw the edges of the cube
        r = [-cell.size/2, cell.size/2]
        for s, e in combinations(np.array(list(product(r, r, r))), 2):
            if np.sum(np.abs(s-e)) == r[1]-r[0]:
                ax.plot3D(*zip(s + cell.center, e + cell.center), color="r")
    # make all axes equal
    ax.set_box_aspect([1,1,1])
    plt.show()

if __name__ == "__main__":
    # For testing, generate random points
    np.random.seed(42)
    n=15
    points = np.random.rand(n, 3)*2-1
    normals = (np.random.rand(n, 3)*2-1)*2
    tree= KDTree(points)

    grid = create_grid(points)
    vector_field(grid, normals, 0.2, tree)

    visualise_grid(grid)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for cell in grid:
        if cell.field is not None:
            ax.quiver(cell.center[0], cell.center[1], cell.center[2],
                      cell.field[0]*2, cell.field[1]*2, cell.field[2]*2,
                      length=0.1, color='r')
    for i in range(n):
        ax.quiver(points[i][0], points[i][1], points[i][2],
                  normals[i][0], normals[i][1], normals[i][2],
                  length=0.1, color='g')
        
    
    plt.show()

    
