import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from .triangleTable import TriangleTable


class Cube(object):

#          6             7
#          +-------------+               +-----6-------+   
#        / |           / |             / |            /|   
#      /   |         /   |          11   7         10   5
#  2 +-----+-------+  3  |         +-----+2------+     |   
#    |   4 +-------+-----+ 5       |     +-----4-+-----+   
#    |   /         |   /           3   8         1   9
#    | /           | /             | /           | /       
#  0 +-------------+ 1             +------0------+   


    def __init__(self, origin, size, verticies):
        self.origin = origin
        self.size = size
        self.vertices = verticies
        self.positions = [
            [origin[0], origin[1], origin[2]],                          # (0 0 0)
            [origin[0] + size, origin[1], origin[2]],                   # (1 0 0)
            [origin[0], origin[1] + size, origin[2]],                   # (0 1 0)
            [origin[0] + size, origin[1] + size, origin[2]],            # (1 1 0)
            [origin[0], origin[1], origin[2] + size],                   # (0 0 1)
            [origin[0] + size, origin[1], origin[2] + size],            # (1 0 1)
            [origin[0], origin[1] + size, origin[2] + size],            # (0 1 1)
            [origin[0] + size, origin[1] + size, origin[2] + size],     # (1 1 1)
        ]
        self.getEdges()
        self.table = TriangleTable.TriangleTable

    def getEdges(self):
        self.edges = [
                        [0,1],
                        [1, 3],
                        [3, 2],
                        [2, 0],
                        [4, 5],
                        [5, 7],
                        [7, 6],
                        [6, 4],
                        [0, 4],
                        [1, 5],
                        [3, 7],
                        [2, 6],
                    ]

    def getBinaryIndex(self,threshold):
        index = 0
        for i in range(8):
            if self.vertices[i] > threshold:
                index |= 1 << i
        return index
    
    def getTriangles(self, threshold=0.5):
        index = self.getBinaryIndex(threshold)
        EdgeIndexes = self.table[index]
        triangles = []
        triangle = []
        i = 0
        for edgeIndex in EdgeIndexes:

            if i%3 == 0 and i != 0:
                triangles.append(triangle)
                triangle = []

            if edgeIndex == -1: return triangles

            # Pas d'interpolation, juste le milieu de l'arête
            v0, v1 = self.edges[edgeIndex]
            mu = (threshold - self.vertices[v0]) / (self.vertices[v1] - self.vertices[v0])
            mid = [
                (self.positions[v0][0] + mu * (self.positions[v1][0] - self.positions[v0][0])),
                (self.positions[v0][1] + mu * (self.positions[v1][1] - self.positions[v0][1])),
                (self.positions[v0][2] + mu * (self.positions[v1][2] - self.positions[v0][2])),
            ]
            triangle.append(mid)
            i += 1
        return triangles
    
    

    def draw(self, ax):
        triangles = self.getTriangles()

        for pos in self.positions:
            color = 'r' if self.vertices[self.positions.index(pos)] > 0.5 else 'b'
            ax.scatter(pos[0], pos[1], pos[2], color=color)

        for triangle in triangles:
            tri = Poly3DCollection([triangle], alpha=0.8, facecolor='r', edgecolor='r')
            ax.add_collection3d(tri)

    def draw_without_vertices(self, ax, alpha=0.5, threshold=0.5):
        triangles = self.getTriangles(threshold)

        for triangle in triangles:
            tri = Poly3DCollection([triangle], alpha=alpha, facecolor='r')
            ax.add_collection3d(tri)

    # Partie Interactive
    def toggleVertex(self, vertexIndex):
        self.vertices[vertexIndex] = 1 - self.vertices[vertexIndex]


    def resetVertices(self):
        self.vertices = [0,0,0,0,0,0,0,0]

    def on_key(self, event, ax, fig):
        if event.key in ['1','2','3','4','5','6', '7', '8']:
            if event.key == 'R': 
                self.resetVertices()
                return
            
            idx = int(event.key) - 1
            self.toggleVertex(idx)
            self.draw(ax)
            fig.canvas.draw_idle()

            
if __name__ == "__main__":
    #exemple with one cube
    cube1 = Cube([0,0,0], 1, [0.7,1,0.5,1,0.2,1,0.3,1])
    triangles = cube1.getTriangles()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    cube1.draw(ax)

    fig.canvas.mpl_connect(
        'key_press_event',
        lambda event: cube1.on_key(event, ax, fig)
    )

    plt.show()

    #exemple with multiple cubes

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')

    # points=np.random.randint(0,2,27).reshape((3,3,3))
    # for i in range(2):
    #     for j in range(2):
    #         for k in range(2):
    #             cube=Cube([i,j,k], 1, [
    #                 points[i][j][k],
    #                 points[i+1][j][k],
    #                 points[i][j+1][k],
    #                 points[i+1][j+1][k],
    #                 points[i][j][k+1],
    #                 points[i+1][j][k+1],
    #                 points[i][j+1][k+1],
    #                 points[i+1][j+1][k+1],
    #             ])
    #             cube.draw(ax)

    # plt.show()
