import numpy as np


def write_ply(filename, vertices, faces=None, vertex_colors=None):
    """
    Write a PLY file.

    Parameters:
    - filename: str, path to the output PLY file.
    - vertices: np.ndarray of shape (N, 3), vertex coordinates.
    - faces: np.ndarray of shape (M, 3) or (M, 4), face indices (optional).
    - vertex_colors: np.ndarray of shape (N, 3) or (N, 4), vertex colors (optional).
    """
    num_vertices = vertices.shape[0]
    num_faces = faces.shape[0] if faces is not None else 0

    with open(filename, 'w') as ply_file:
        # Write header
        ply_file.write("ply\n")
        ply_file.write("format ascii 1.0\n")
        ply_file.write(f"element vertex {num_vertices}\n")
        ply_file.write("property float x\n")
        ply_file.write("property float y\n")
        ply_file.write("property float z\n")
        if vertex_colors is not None:
            if vertex_colors.shape[1] == 3:
                ply_file.write("property uchar red\n")
                ply_file.write("property uchar green\n")
                ply_file.write("property uchar blue\n")
            elif vertex_colors.shape[1] == 4:
                ply_file.write("property uchar red\n")
                ply_file.write("property uchar green\n")
                ply_file.write("property uchar blue\n")
                ply_file.write("property uchar alpha\n")
        if num_faces > 0:
            ply_file.write(f"element face {num_faces}\n")
            ply_file.write("property list uchar int vertex_indices\n")
        ply_file.write("end_header\n")

        # Write vertex data
        for i in range(num_vertices):
            vertex_line = f"{vertices[i, 0]} {vertices[i, 1]} {vertices[i, 2]}"
            if vertex_colors is not None:
                color = vertex_colors[i]
                vertex_line += f" {int(color[0])} {int(color[1])} {int(color[2])}"
                if vertex_colors.shape[1] == 4:
                    vertex_line += f" {int(color[3])}"
            ply_file.write(vertex_line + "\n")

        # Write face data
        if num_faces > 0:
            for i in range(num_faces):
                face = faces[i]
                face_line = f"{len(face)} " + " ".join(map(str, face))
                ply_file.write(face_line + "\n")


# # Circle 
# R = 1
# N = 100
# theta = np.linspace(0, 2 * np.pi, N)
# x = R * np.cos(theta)
# y = R * np.sin(theta)
# z = np.zeros(N)
# vertices = np.vstack((x, y, z)).T
# faces = np.array([[i, (i + 1) % N, N] for i in range(N)])
# vertices = np.vstack((vertices, np.array([[0, 0, 0]])))
# write_ply("circle.ply", vertices, faces)


# Sphere
R =1
num_lat = 10
num_lon = 20
vertices = []
faces = []
for i in range(num_lat + 1):
    theta = np.pi * i / num_lat
    for j in range(num_lon):
        phi = 2 * np.pi * j / num_lon
        x = R * np.sin(theta) * np.cos(phi)
        y = R * np.sin(theta) * np.sin(phi)
        z = R * np.cos(theta)
        vertices.append([x, y, z])
vertices = np.array(vertices)
for i in range(num_lat):
    for j in range(num_lon):
        p1 = i * num_lon + j
        p2 = p1 + num_lon
        p3 = p2 + 1 if (j + 1) < num_lon else p2 + 1 - num_lon
        p4 = p1 + 1 if (j + 1) < num_lon else p1 + 1 - num_lon
        if i != 0:
            faces.append([p1, p2, p4])
        if i != (num_lat - 1):
            faces.append([p4, p2, p3])
faces = np.array(faces)
write_ply("sphere.ply", vertices, faces)
