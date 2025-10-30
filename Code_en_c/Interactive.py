import matplotlib.pyplot as plt
import numpy as np
import ctypes
import os

input_file = b"temp/temp_input.txt"
output_file = b"temp/temp_output.txt"

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

# --- Initialisation ---
fig, ax = plt.subplots()
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.set_title("Cliquez pour ajouter des points")

# lib = ctypes.CDLL(os.path.abspath("shared_lib/BowyerWatson.so"))   # Linux
lib = ctypes.CDLL(os.path.abspath("shared_lib/BowyerWatsonMACOS.so"))   # Mac
# lib = ctypes.CDLL(os.path.abspath("shared_lib/BowyerWatson.dll"))  # Windows
lib.del2d_py.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
lib.del2d_py.restype = ctypes.c_int

points_x, points_y = [], []
(scatter,) = ax.plot([], [], "ro")
tri_lines = None

def update():
    """Met à jour le graphique après un clic."""
    global tri_lines
    arr_x = np.array(points_x, dtype=np.double)
    arr_y = np.array(points_y, dtype=np.double)

    if tri_lines:
        for t in tri_lines:
            t.remove()

    scatter.set_data(points_x, points_y)

    if len(points_x) >= 3:
        write_points_to_file(points_x, points_y, input_file)
        res = lib.del2d_py(input_file, output_file)
        if res != 0:
            print(f"Erreur dans l'appel C: code {res}")
        tri = read_triangles_from_file(output_file)
        tri_indices = []
        for t in tri:
            i0 = find_closest_index(t[0][0], t[0][1], arr_x, arr_y)
            i1 = find_closest_index(t[1][0], t[1][1], arr_x, arr_y)
            i2 = find_closest_index(t[2][0], t[2][1], arr_x, arr_y)
            tri_indices.append((i0, i1, i2))
        if tri_indices:
            tri_lines = ax.triplot(arr_x, arr_y, np.array(tri_indices), "b-", lw=1)
        else:
            tri_lines = []
    else:
        tri_lines = []

    plt.draw()  # redessine immédiatement

def onclick(event):
    if event.inaxes != ax:
        return
    points_x.append(event.xdata)
    points_y.append(event.ydata)
    print(f"Point ajouté : ({event.xdata:.2f}, {event.ydata:.2f})")
    update()  # mise à jour du plot directement après le clic

fig.canvas.mpl_connect("button_press_event", onclick)

plt.show()
