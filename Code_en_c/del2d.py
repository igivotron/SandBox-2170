import numpy as np
import ctypes
import os
import random
import time



# file names
input_file = b"input_file/points_1k.txt"
output_file = b"output_file/output_mesh.txt"
save_plot_file = b"output_file/points_1k_plot.pdf"

# Load the shared library
lib = ctypes.CDLL(os.path.abspath("shared_lib/BowyerWatson.so"))   # For Linux/Mac ou wsl
# lib = ctypes.CDLL(os.path.abspath("shared_lib/BowyerWatson.dll"))  # For Windows

# Tell ctypes the function signature
lib.del2d_py.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
lib.del2d_py.restype = ctypes.c_int

# Call the C function
print("code launch")
start_time = time.time()
result = lib.del2d_py(input_file, output_file)
end_time = time.time()
print(f"C function returned: {result}")
print(f"Execution time of del2d_py: {end_time - start_time:.6f} seconds")


to_plot = True
if to_plot:
    import matplotlib.pyplot as plt

    # Lire le fichier de triangles
    triangles = []
    with open("output_file/output_mesh.txt") as f:
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
    plt.savefig(save_plot_file)
    print(f"Graphique enregistré dans {save_plot_file.decode()}")
