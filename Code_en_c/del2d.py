import numpy as np
import ctypes
import os
import random
import time


with_Hilbert = True

# file names
input_file = b"Code_en_c/input_file/points_100k.txt"
hilbert_file = b"Code_en_c/Hilbert_file/points_100k_hilbert.txt"
Output_file = b"Code_en_c/output_file/output_mesh.txt"

# create the coordinates x, y
if not with_Hilbert:
    # read points in txt
    points = []
    with open(input_file, "r") as f:
        # first line contains the number of points
        first = f.readline().strip()
        n = int(first)
        for i in range(n):
            line = f.readline()
            x, y = map(float, line.split())
            points.append((x, y))
    
    x = (ctypes.c_double * len(points))(*[p[0] for p in points])
    y = (ctypes.c_double * len(points))(*[p[1] for p in points])
    n = ctypes.c_size_t(len(points))

if with_Hilbert:
    lib2 = ctypes.CDLL(os.path.abspath("Code_en_c/shared_lib/Hilbert.dll"))
    lib2.Hilbert_py.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
    lib2.Hilbert_py.restype = ctypes.c_int
    
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


# Load the shared library
lib = ctypes.CDLL(os.path.abspath("Code_en_c/shared_lib/BowyerWatson.dll"))

# Tell ctypes the function signature
lib.del2d_py.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]
lib.del2d_py.restype = ctypes.c_int

# Call the C function
print("code launch")
start_time = time.time()
result = lib.del2d_py(x, y, n)
end_time = time.time()
print(f"C function returned: {result}")
print(f"Execution time of del2d_py: {end_time - start_time:.6f} seconds")


to_plot = False
if to_plot:
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
