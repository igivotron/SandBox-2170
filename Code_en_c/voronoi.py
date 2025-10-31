import matplotlib.pyplot as plt
import numpy as np
from Structure import Triangle, Voronoi_cell, Point
from scipy.spatial import Delaunay
import ctypes
import os

def find_closest_index(vx, vy, xs, ys, tol=1e-6):
    if len(xs) == 0:
        raise ValueError("Liste de points vide")
    arr = np.column_stack((xs, ys))
    dists = np.hypot(arr[:, 0] - vx, arr[:, 1] - vy)
    idx = int(np.argmin(dists))
    # if dists[idx] > tol:
    #     print(f"Warning: closest distance {dists[idx]:.3e} > tol ({tol}) for vertex ({vx:.6f}, {vy:.6f})")
    return idx

def write_points_to_file(x_points, y_points, filename):
    with open(filename, "wb") as f:
        num_points = len(x_points)
        f.write(f"{num_points}\n".encode())
        for x, y in zip(x_points, y_points):
            f.write(f"{x} {y}\n".encode())

def read_triangles_from_file(filename):
    triangles = []
    with open(filename) as f:
        num_triangles = int(f.readline().strip())
        for line in f:
            vals = list(map(float, line.split()))
            if len(vals) == 6:
                triangles.append([(vals[0], vals[1]),
                                  (vals[2], vals[3]),
                                  (vals[4], vals[5])])
    return triangles

def circumcenter(triangle):
    A = triangle.points[0]
    B = triangle.points[1]
    C = triangle.points[2]

    D = 2 * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y))
    if abs(D) < 1e-12:
        return ((A.x + B.x + C.x) / 3.0, (A.y + B.y + C.y) / 3.0)

    Ux = ((A.x**2 + A.y**2) * (B.y - C.y) +
          (B.x**2 + B.y**2) * (C.y - A.y) +
          (C.x**2 + C.y**2) * (A.y - B.y)) / D
    Uy = ((A.x**2 + A.y**2) * (C.x - B.x) +
          (B.x**2 + B.y**2) * (A.x - C.x) +
          (C.x**2 + C.y**2) * (B.x - A.x)) / D

    return (Ux, Uy)

def get_triangle_circumcenters(triangles):
    circumcenters = []
    for triangle in triangles:
        center = circumcenter(triangle)
        circumcenters.append(center)
        triangle.center = center
    return circumcenters

def Delaunay_triangulation(pts):
    write_points_to_file(pts[:, 0], pts[:, 1], input_file)
    res = lib.del2d_py(input_file, output_file)
    if res != 0:
        print(f"Erreur dans l'appel C: code {res}")
    tri = read_triangles_from_file(output_file)
    
    triangles = []
    for j, t in enumerate(tri):
        p1 = Point(t[0][0], t[0][1],find_closest_index(t[0][0], t[0][1], pts[:,0], pts[:,1]))
        p2 = Point(t[1][0], t[1][1],find_closest_index(t[1][0], t[1][1], pts[:,0], pts[:,1]))
        p3 = Point(t[2][0], t[2][1],find_closest_index(t[2][0], t[2][1], pts[:,0], pts[:,1]))
        triangle = Triangle(p1, p2, p3, j)
        center = circumcenter(triangle)
        triangle.center = center
        triangles.append(triangle)
    return triangles

def find_neighbors(triangles):
    """
    Trouve les triangles voisins.
    Retourne un dictionnaire: 
        key: arête 
        val: liste des triangles partageant cette arête. (stocker que l'index des triangles ?)
    """
    edge_dict = {}
    for triangle in triangles:
        points = triangle.points
        edges = [(points[0].index, points[1].index),
                 (points[1].index, points[2].index),
                 (points[2].index, points[0].index)]
        for edge in edges:
            edge_key = tuple(sorted(edge))
            if edge_key in edge_dict:
                edge_dict[edge_key].append(triangle)
            else:
                edge_dict[edge_key] = [triangle]
    for triangle in triangles:
        triangle.neighbors = []
    for neighbors in edge_dict.values():
        if len(neighbors) == 2:
            t1, t2 = neighbors
            t1.neighbors.append(t2)
            t2.neighbors.append(t1)
    return edge_dict


# --- C setup ---
input_file = b"temp_input.txt"
output_file = b"temp_output.txt"
lib = ctypes.CDLL(os.path.abspath("shared_lib/BowyerWatson.dll"))  # Windows
lib.del2d_py.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
lib.del2d_py.restype = ctypes.c_int

# --- Setup figure ---
fig, ax = plt.subplots()
ax.set_aspect("equal")
xmin, xmax = 0, 10
ymin, ymax = 0, 10
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_title("Cliquez pour ajouter des points (Voronoi Cells)")

# structure mutable partagée entre callbacks permet de switcher du mode et ajouter des points
data = {
    'pts': np.array([(5.0, 2.0), (5.0, 8.0), (2.0, 5.0), (8.0, 5.0), (5.0, 4.0), (5.0, 6.0), (9.0, 1.0)]),
    # 'pts': np.array([(2.0, 2.0), (4.0, 2.0), (6.0, 2.0)]),
    'mode': 1  # 0: Delaunay, 1: Voronoi
}

def update(ax, data):
    ax.clear()
    ax.set_aspect("equal")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title("Cliquez pour ajouter des points (Voronoi Cells)")
    triangles = Delaunay_triangulation(data['pts'])
    plot(triangles, data['pts'], ax, data['mode'])
    fig.canvas.draw_idle()


def plot(triangles, pts, ax, mode, extend_length=1000):
    if mode == 0:
        # Dessin des triangles
        for triangle in triangles:
            arr_x = [p.x for p in triangle.points] + [triangle.points[0].x]
            arr_y = [p.y for p in triangle.points] + [triangle.points[0].y]
            ax.plot(arr_x, arr_y, 'b-')
            cx, cy = triangle.center
            ax.plot(cx, cy, 'ro')
    else:
        # Segments de Voronoï (fini + infini)
        edge_dict = find_neighbors(triangles)
        for edge, tris in edge_dict.items():
            if len(tris) == 2:
                t1, t2 = tris
                ax.plot([t1.center[0], t2.center[0]],
                        [t1.center[1], t2.center[1]], 'g-')
            else:
                t = tris[0]
                A = pts[edge[0]]
                B = pts[edge[1]]
                AB = np.array([B[0] - A[0], B[1] - A[1]])
                n = np.array([A[1] - B[1], B[0] - A[0]])
                norm = np.linalg.norm(n)
                if norm < 1e-12:
                    continue
                n = n / norm
                C = np.array(t.center)

                # Sens du vecteur normal
                # On le dirige vers l'extérieur du triangle. Le barycentre est le point moyen du nuage (intérieur au nuage)
                # to_center pointe du centre du triangle vers le barycentre
                # on inverse n si dot (to_center, n) > 0 (cos theta > 0)
                barycenter = np.mean(pts, axis=0)
                to_center = barycenter - C
                if np.dot(n, to_center) > 0:
                    n = -n
                ax.plot([C[0], C[0] + n[0]*extend_length],
                        [C[1], C[1] + n[1]*extend_length],
                        'r--')

    # pts
    ax.plot(pts[:, 0], pts[:, 1], 'ko')


def onclick(event):
    if event.inaxes != ax:
        return
    new_pt = np.array([[event.xdata, event.ydata]])
    data['pts'] = np.vstack([data['pts'], new_pt])
    # data['pts'] = data['pts'][np.lexsort((data['pts'][:,1], data['pts'][:,0]))]
    print(f"Point ajouté : ({event.xdata:.2f}, {event.ydata:.2f})")
    update(ax, data)


def switch_mode(event):
    if event.key == 'r':
        data['mode'] = 1 - data['mode']
        if data['mode'] == 1:
            print("Mode Voronoi Cells")
        else:
            print("Mode Delaunay Triangulation")
        update(ax, data)


cid = fig.canvas.mpl_connect('button_press_event', onclick)
kid = fig.canvas.mpl_connect('key_press_event', switch_mode)

update(ax, data)
plt.show()

# delete temp files
if os.path.exists(input_file.decode()):
    os.remove(input_file.decode())
if os.path.exists(output_file.decode()):
    os.remove(output_file.decode())