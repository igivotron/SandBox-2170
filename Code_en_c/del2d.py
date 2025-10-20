import numpy as np
import ctypes
import os
import random


# Load the shared library
lib = ctypes.CDLL(os.path.abspath("Code_en_c/BowyerWatson.dll"))
lib2 = ctypes.CDLL(os.path.abspath("Code_en_c/Hilbert.dll"))

# Tell ctypes the function signature
lib.del2d_py.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]
lib.del2d_py.restype = ctypes.c_int

lib2.Hilbert_py.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
lib2.Hilbert_py.restype = ctypes.c_int

# Prepare input for Hilbert function
input_file = b"inputs/10000pts"
hilbert_file = b"outputs/10000pts_hilbert.txt"
Output_file = b"Code_en_c/output_mesh.txt"

# Call the Hilbert C function
print("Hilbert code launch")
result_hilbert = lib2.Hilbert_py(input_file, hilbert_file)
print(f"Hilbert C function returned: {result_hilbert}")

with open(input_file, "r") as f:
    data = f.read().strip().split("\n")

n = int(data[0])
x = np.zeros(n)
y = np.zeros(n)
for i in range(n):
    x[i], y[i] = map(float, data[i + 1].split())

with open(hilbert_file, "r") as f:
    data = f.read().strip().split("\n")

for i in range(n): data[i] = ''.join(data[i].split())

sorted_indices = sorted(range(n), key=lambda i: data[i])
x = x[sorted_indices]
y = y[sorted_indices]

# Call the C function
x=x.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
y=y.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
n=ctypes.c_size_t(n)
print("code launch")
result = lib.del2d_py(x, y, n)
print(f"C function returned: {result}")


import matplotlib.pyplot as plt

# Lire le fichier de triangles
triangles = []
with open(Output_file) as f:
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
