from HalfEdge import TriangularMesh, Vertex, Face, Halfedge
import matplotlib.pyplot as plt
import scipy as sp

def import_pts(filename):
    pts = []
    with open(filename, "r") as f:
        lines = f.readlines()
    nr_points = int(lines[0].strip())
    for i in range(1, nr_points + 1):
        x, y = map(float, lines[i].strip().split())
        pts.append((x, y, 0.0))
    return pts, len(pts)

def getBoundingBox(pts, L):
    x_min = min(pts, key=lambda p: p[0])[0]
    x_max = max(pts, key=lambda p: p[0])[0]
    y_min = min(pts, key=lambda p: p[1])[1]
    y_max = max(pts, key=lambda p: p[1])[1]
    return [(x_min-L, y_min-L, 0), (x_max+L, y_min-L, 0), (x_max+L, y_max+L, 0), (x_min-L, y_max+L, 0)]

def initMesh(pts, SuperTriangle):
    mesh = TriangularMesh(pts, SuperTriangle)
    return mesh

def plotMesh(mesh):
    plt.figure()
    xpts = [v.x for v in mesh.vertices]
    ypts = [v.y for v in mesh.vertices]
    for face in mesh.faces:
        he1 = mesh.halfedges[face.halfedge]
        he2 = mesh.halfedges[he1.next]
        he3 = mesh.halfedges[he2.next]
        A = mesh.vertices[he1.vertex]
        B = mesh.vertices[he2.vertex]
        C = mesh.vertices[he3.vertex]
        plt.plot([A.x, B.x], [A.y, B.y], 'b-')
        plt.plot([B.x, C.x], [B.y, C.y], 'b-')
        plt.plot([C.x, A.x], [C.y, A.y], 'b-')
    plt.plot(xpts, ypts, 'ro')
    plt.axis('equal')
    plt.savefig("BowyerWatson.png")


def plot_check(mesh):
    plt.figure()
    xpts = [v.x for v in mesh.vertices]
    ypts = [v.y for v in mesh.vertices]
    plt.plot(xpts, ypts, 'ro')
    tri = sp.spatial.Delaunay([(v.x, v.y) for v in mesh.vertices])
    plt.triplot(xpts, ypts, tri.simplices.copy(), 'g--')
    plt.axis('equal')
    

def inCircle(v, face, mesh):
    he1 = mesh.halfedges[face.halfedge]
    he2 = mesh.halfedges[he1.next]
    he3 = mesh.halfedges[he2.next]
    A = mesh.vertices[he1.vertex]
    B = mesh.vertices[he2.vertex]
    C = mesh.vertices[he3.vertex]
    ax, ay = A.x - v.x, A.y - v.y
    bx, by = B.x - v.x, B.y - v.y
    cx, cy = C.x - v.x, C.y - v.y
    det = (ax * ax + ay * ay) * (bx * cy - by * cx) - (bx * bx + by * by) * (ax * cy - ay * cx) + (cx * cx + cy * cy) * (ax * by - ay * bx)
    return det > 0

def isEdgeShared(he, Cavity, mesh):
    v_start = he.vertex
    v_end = mesh.halfedges[he.next].vertex
    count = 0
    for face in Cavity:
        he_face = mesh.halfedges[face.halfedge]
        for _ in range(3):
            v1 = he_face.vertex
            v2 = mesh.halfedges[he_face.next].vertex
            if (v_start == v2 and v_end == v1 or v_start == v1 and v_end == v2):
                count += 1
                if count > 1:
                    return True
            he_face = mesh.halfedges[he_face.next]
    return False

def BowyerWatson(mesh):
    n_vortices = len(mesh.vertices)-4
    for i in range(n_vortices):
        Cavity = []
        to_remove = []

        for j in range(len(mesh.faces)):
            if inCircle(mesh.vertices[i], mesh.faces[j], mesh):
                Cavity.append(mesh.faces[j])
                to_remove.append(j)
        

        for index in sorted(to_remove, reverse=True):
            del mesh.faces[index]
            
        for j in range(len(Cavity)):
            he = mesh.halfedges[Cavity[j].halfedge]
            for k in range(3):
                if not isEdgeShared(he, Cavity, mesh):
                    new_face = Face(halfedge=len(mesh.halfedges))
                    mesh.faces.append(new_face)
                    he1 = Halfedge(vertex=he.vertex, facet=len(mesh.faces)-1, next=len(mesh.halfedges)+1, opposite=-1)
                    he2 = Halfedge(vertex=mesh.halfedges[he.next].vertex, facet=len(mesh.faces)-1, next=len(mesh.halfedges)+2, opposite=-1)
                    he3 = Halfedge(vertex=i, facet=len(mesh.faces)-1, next=len(mesh.halfedges), opposite=-1)
                    mesh.halfedges.extend([he1, he2, he3])
                    new_face.halfedge = len(mesh.halfedges)-3
                he = mesh.halfedges[he.next]                
    return


if __name__ == "__main__":
    pts, n = import_pts("inputs/10000pts")
    boundary = getBoundingBox(pts, 0.1)
    pts = pts + boundary
    SuperTriangles = [(n, n+1, n+2), (n, n+2, n+3)]
    mesh = initMesh(pts, SuperTriangles)
    BowyerWatson(mesh)
    plotMesh(mesh)
    plot_check(mesh)
    plt.show()
    #TODO: Supprimer les halfedges lorsqu'on surprime une face