import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import ctypes
import os

input_file = b"temp_input.txt"
output_file = b"temp_output.txt"

def find_closest_index(vx, vy, xs, ys, tol=1e-6):
    """
    Retourne l'indice i dans xs,ys minimisant la distance euclidienne au point (vx, vy).
    Si la distance minimale dépasse tol, affiche un avertissement (mais renvoie tout de même l'indice).
    """
    if len(xs) == 0:
        raise ValueError("Liste de points vide")
    arr = np.column_stack((xs, ys))
    dists = np.hypot(arr[:, 0] - vx, arr[:, 1] - vy)
    idx = int(np.argmin(dists))
    if dists[idx] > tol:
        # tol peut être ajustée selon la précision souhaitée
        print(f"Warning: closest distance {dists[idx]:.3e} > tol ({tol}) for vertex ({vx:.6f}, {vy:.6f})")
    return idx

def write_points_to_file(x_points, y_points, filename):
    """Écrit les points dans un fichier au format attendu par la fonction C."""
    with open(filename, "wb") as f:
        num_points = len(x_points)
        f.write(f"{num_points}\n".encode())
        for x, y in zip(x_points, y_points):
            f.write(f"{x} {y}\n".encode())

def read_triangles_from_file(filename):
    """Lit les triangles depuis un fichier généré par la fonction C."""
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

# --- Initialisation de la figure ---
fig, ax = plt.subplots()
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.set_title("Cliquez pour ajouter des points (animation + Delaunay)")

# lib = ctypes.CDLL(os.path.abspath("shared_lib/BowyerWatson.so"))   # For Linux/wsl
lib = ctypes.CDLL(os.path.abspath("shared_lib/BowyerWatsonMACOS.so"))   # For Mac
# lib = ctypes.CDLL(os.path.abspath("shared_lib/BowyerWatson.dll"))  # For Windows
lib.del2d_py.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
lib.del2d_py.restype = ctypes.c_int



# --- Données ---
points_x, points_y = [], []

# --- Tracés initiaux ---
(scatter,) = ax.plot([], [], "ro")
tri_lines = None  # objet pour la triangulation

# --- Fonction de mise à jour de l'animation ---
def update(frame):
    global tri_lines
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    arr_x = np.array(points_x, dtype=np.double)
    arr_y = np.array(points_y, dtype=np.double)

    # Nettoyer la triangulation précédente
    if tri_lines:
        for t in tri_lines:
            t.remove()

    # Mettre à jour les points
    scatter.set_data(points_x, points_y)

    # Dessiner la nouvelle triangulation si assez de points
    if len(points_x) >= 3:
        write_points_to_file(points_x, points_y, input_file)
        res = lib.del2d_py(input_file, output_file)
        if res != 0:
            print(f"Erreur dans l'appel C: code {res}")
        tri = read_triangles_from_file(output_file)
        tri_indices = []
        for t in tri:
            i0 = find_closest_index(t[0][0], t[0][1], arr_x, arr_y, tol=1e-6)
            i1 = find_closest_index(t[1][0], t[1][1], arr_x, arr_y, tol=1e-6)
            i2 = find_closest_index(t[2][0], t[2][1], arr_x, arr_y, tol=1e-6)
            tri_indices.append((i0, i1, i2))
        if tri_indices:
            tri_lines = ax.triplot(arr_x, arr_y, np.array(tri_indices), "b-", lw=1)
        else:
            tri_lines = []
    else:
        tri_lines = []

    return scatter,

# --- Fonction de clic ---
def onclick(event):
    if event.inaxes != ax:
        return
    # Ajout d’un point
    points_x.append(event.xdata)
    points_y.append(event.ydata)
    print(f"Point ajouté : ({event.xdata:.2f}, {event.ydata:.2f})")

# --- Connexion du clic ---
fig.canvas.mpl_connect("button_press_event", onclick)

# --- Animation ---
ani = animation.FuncAnimation(fig, update, interval=100, blit=False)

plt.show()
