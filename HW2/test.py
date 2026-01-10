import numpy as np
import ctypes
import os

lib = ctypes.CDLL(os.path.abspath("bvh_C.dll"))

class BVHNode(ctypes.Structure):
    _fields_ = [
        ("left", ctypes.c_int),                         # int
        ("right", ctypes.c_int),                        # int
        ("parent", ctypes.c_int),                       # int
        ("bbox", ctypes.POINTER(ctypes.c_double)),      # double*
        ("items", ctypes.POINTER(ctypes.c_int)),        # int*
        ("n_items", ctypes.c_int),                      # int
        ("index", ctypes.c_int),                        # int
        ("state", ctypes.c_int),                        # int
    ]

class bvh_C(ctypes.Structure):
    _fields_ = [
        ("nodes", ctypes.POINTER(BVHNode)),     # BVHNode*
        ("n_nodes", ctypes.c_int),              # int
        ("positions", ctypes.POINTER(ctypes.c_double)),  # double*
        ("radii", ctypes.POINTER(ctypes.c_double)),      # double*
        ("NperLeaf", ctypes.c_int),             # int
        ("root", ctypes.c_int),                 # int
    ]


lib.create_bvh.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int]
lib.create_bvh.restype = ctypes.POINTER(bvh_C)
lib.update.argtypes = [ctypes.POINTER(bvh_C), ctypes.POINTER(BVHNode)]
lib.update.restype = None
lib.free_bvh.argtypes = [ctypes.POINTER(bvh_C)]
lib.free_bvh.restype = None
lib.update_positions.argtypes = [ctypes.POINTER(bvh_C), ctypes.POINTER(ctypes.c_double)]
lib.update_positions.restype = None
lib.build_bvh.argtypes = [ctypes.POINTER(bvh_C)]
lib.build_bvh.restype = None
lib.is_leaf.argtypes = [ctypes.POINTER(BVHNode), ctypes.c_int]
lib.is_leaf.restype = ctypes.c_bool
lib.find_pot_inter.argtypes = [ctypes.POINTER(bvh_C), ctypes.POINTER(ctypes.c_int)]
lib.find_pot_inter.restype = ctypes.c_int


N=6
# pos = np.random.rand(N,3).astype(np.float64).flatten()
pos=np.array([[0, 0, 0],
     [1, 1, 1],
     [2, 2, 2],[2, 2.15, 2],
     [3, 3, 3],
     [4, 4, 4]]).astype(np.float64).flatten()
# radii = np.random.rand(N).astype(np.float64).flatten()
radii=np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1]).astype(np.float64).flatten()
pos_c = pos.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
radii_c = radii.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
bvh_tree = lib.create_bvh(pos_c, radii_c, N, 1)
lib.build_bvh(bvh_tree)

def print_bvh(node_index, depth=0):
    node = bvh_tree.contents.nodes[node_index]
    if lib.is_leaf(ctypes.byref(node), 1):
        bbox = np.ctypeslib.as_array(node.bbox, shape=(2, 3))
        print("  " * depth + f"Node(depth={depth}, bbox=[[{bbox[0][0]:.2f}, {bbox[0][1]:.2f}, {bbox[0][2]:.2f}], [{bbox[1][0]:.2f}, {bbox[1][1]:.2f}, {bbox[1][2]:.2f}]], items={[node.items[i] for i in range(node.n_items)]})")
        return
    bbox = np.ctypeslib.as_array(node.bbox, shape=(2, 3))
    print("  " * depth + f"Node(depth={depth}, bbox=[[{bbox[0][0]:.2f}, {bbox[0][1]:.2f}, {bbox[0][2]:.2f}], [{bbox[1][0]:.2f}, {bbox[1][1]:.2f}, {bbox[1][2]:.2f}]], items={[node.items[i] for i in range(node.n_items)]})")
    print_bvh(node.left, depth + 1)
    print_bvh(node.right, depth + 1)

print_bvh(bvh_tree.contents.root)

#test find_pot_inter
potential_intersections = np.zeros((N*N), dtype=np.int32)
potential_intersections_c = potential_intersections.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
print(lib.find_pot_inter(bvh_tree, potential_intersections_c))
print("Potential intersections (pairs of indices):")
print(potential_intersections)


lib.free_bvh(bvh_tree)

