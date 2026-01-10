import open3d as o3d
import polyscope as ps
import numpy as np

ps.init()

mesh = o3d.io.read_triangle_mesh("test_avec400.ply")

vertices = np.asarray(mesh.vertices)
faces = np.asarray(mesh.triangles)

ps.register_surface_mesh("mesh", vertices, faces)

ps.show()