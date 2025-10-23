import matplotlib.pyplot as plt

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
    plt.plot([*xs, xs[0]], [*ys, ys[0]], color='black', linewidth=0.8)
    plt.fill(xs, ys, alpha=0.3)

plt.title("Delaunay Triangulation (from C output)")
plt.xlabel("X")
plt.ylabel("Y")
plt.axis("equal")
plt.show()