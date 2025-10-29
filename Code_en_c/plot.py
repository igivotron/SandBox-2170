import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

triangles = []
with open("output_mesh.txt") as f:
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
    plt.plot([*xs, xs[0]], [*ys, ys[0]], color='black', linewidth=0.1)
    plt.fill(xs, ys, alpha=0.3)



test = True

if test:
    with open("../inputs/100pts") as f:
        points = []
        for line in f:
            vals = list(map(float, line.split()))
            if len(vals) == 2:
                points.append((vals[0], vals[1]))

    x= []
    y = []
    for p in points:
        x.append(p[0])
        y.append(p[1])
    
    plt.scatter(x, y, color='red', s=10)
    mtriang = mtri.Triangulation(x, y)
    plt.triplot(mtriang, color='red', linewidth=0.5, alpha=0.5)


plt.title("Delaunay Triangulation (from C output)")
plt.xlabel("X")
plt.ylabel("Y")
plt.axis("equal")
plt.tight_layout()
plt.savefig("output_mesh.png", dpi=600)
plt.show()
