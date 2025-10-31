import numpy as np



class Point:
    def __init__(self, x, y, j):
        self.x = x
        self.y = y
        self.index = j
        self.cell = None

class Triangle:
    def __init__(self, p1, p2, p3, j):
        self.points = [p1, p2, p3]
        self.center = None
        self.j = j

class Voronoi_cell:
    def __init__(self, center, corners):
        self.center = center
        self.corners = corners
    def add_corner(self, corner):
        self.corners.append(corner)

