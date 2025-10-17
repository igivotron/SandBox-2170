import ctypes
import os

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
lib = ctypes.CDLL(os.path.abspath("Code_en_c/sheee.dll"))

# Tell ctypes the function signature
lib.del2d_py.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]
lib.del2d_py.restype = ctypes.c_int

# Call the C function
x = (ctypes.c_double * len(points))(*[p[0] for p in points])
y = (ctypes.c_double * len(points))(*[p[1] for p in points])
n = ctypes.c_size_t(len(points))
result = lib.del2d_py(x, y, n)
print(f"C function returned: {result}")