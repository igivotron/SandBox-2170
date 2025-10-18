import ctypes
import os
import random

# read points in points_10k.txt
points = []
with open("Code_en_c/points_10k.txt", "r") as f:
    # first line contains the number of points
    first = f.readline().strip()
    n = int(first)
    for i in range(n):
        line = f.readline()
        x, y = map(float, line.split())
        points.append((x, y))

# Load the shared library
lib = ctypes.CDLL(os.path.abspath("Code_en_c/BowyerWatson.dll"))

# Tell ctypes the function signature
lib.del2d_py.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]
lib.del2d_py.restype = ctypes.c_int

# Call the C function
x = (ctypes.c_double * len(points))(*[p[0] for p in points])
y = (ctypes.c_double * len(points))(*[p[1] for p in points])
n = ctypes.c_size_t(len(points))
result = lib.del2d_py(x, y, n)
print(f"C function returned: {result}")


import matplotlib.pyplot as plt

# Lire le fichier de triangles
triangles = []
with open("Code_en_c/output_mesh.txt") as f:
    num_triangles = int(f.readline().strip())
    for line in f:
        vals = list(map(float, line.split()))
        if len(vals) == 6:
            triangles.append([(vals[0], vals[1]),
                              (vals[2], vals[3]),
                              (vals[4], vals[5])])

# Tracé
plt.figure(figsize=(8, 8))
for tri in triangles:
    xs, ys = zip(*tri)
    # fermer le triangle en revenant au premier point
    plt.plot([*xs, xs[0]], [*ys, ys[0]], color='black', linewidth=0.8)
    plt.fill(xs, ys, alpha=0.3)

plt.title("Delaunay Triangulation (from C output)")
plt.xlabel("X")
plt.ylabel("Y")
plt.axis("equal")
plt.show()
