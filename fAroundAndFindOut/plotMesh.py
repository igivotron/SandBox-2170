import matplotlib.pyplot as plt
import numpy as np

with open("outputbis.obj", "r") as file:
    lines = file.readlines()

vertices = []
faces = []
for line in lines:
    if line.startswith('v '):
        parts = line.strip().split()
        vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
    elif line.startswith('f '):
        parts = line.strip().split()
        face = [int(part.split('/')[0]) for part in parts[1:]]
        faces.append(face)

vertices = np.array(vertices)
faces = np.array(faces)


fig = plt.figure()
plt.scatter(vertices[:, 0], vertices[:, 1])
for face in faces:
    face_vertices = vertices[face ]
    plt.fill(face_vertices[:, 0], face_vertices[:, 1], edgecolor='k', fill=False)

plt.savefig("mesh_plotbis.png")
plt.show()
