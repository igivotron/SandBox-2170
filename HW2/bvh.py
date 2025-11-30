import numpy as np
# from simulator import Simulator, Contact, copy_buffer
# from homework import *


class Nodes:
    def __init__(self, item, left=None, right=None, bbox=None, index=None):
        self.left = left
        self.right = right
        self.parent = None
        self.bbox = bbox #[[x_min, y_min, z_min], [x_max, y_max, z_max]]
        self.index = index
        self.item = item #contain the index of the item if leaf
        self.state = False #updated or not
        
    def is_leaf(self):
        return self.left is None and self.right is None

    def update_bbox(self):
        bbox_min = np.min([self.left.bbox[0], self.right.bbox[0]], axis=0)
        bbox_max = np.max([self.left.bbox[1], self.right.bbox[1]], axis=0)
        self.bbox = np.array([bbox_min, bbox_max])
        return

class BVH:
    def __init__(self, positions, radii, NperLeaf=1):
        self.positions = positions
        self.radii = radii
        self.N_items = len(positions)
        self.indexes = np.arange(len(positions))
        self.NperLeaf = NperLeaf
        self.root = Nodes(item=self.indexes)

    def compute_bbox(self, items):
        x_min = np.min(self.positions[items, 0] - self.radii[items])
        x_max = np.max(self.positions[items, 0] + self.radii[items])
        y_min = np.min(self.positions[items, 1] - self.radii[items])
        y_max = np.max(self.positions[items, 1] + self.radii[items])
        z_min = np.min(self.positions[items, 2] - self.radii[items])
        z_max = np.max(self.positions[items, 2] + self.radii[items])
        bbox = np.array([[x_min, y_min, z_min], [x_max, y_max, z_max]])
        return bbox
    def combine_bbox(left_box, right_box):
        x_min = np.min(left_box[0][0], right_box[0][0])
        y_min = np.min(left_box[0][1], right_box[0][1])
        z_min = np.min(left_box[0][2], right_box[0][2])
        x_max = np.min(left_box[1][0], right_box[1][0])
        y_max = np.min(left_box[1][1], right_box[1][1])
        z_max = np.min(left_box[1][2], right_box[1][2])
        return np.array([[x_min, y_min, z_min], [x_max, y_max, z_max]])


    def build_recursive(self, node, k):
        axis = k % 3
        items = node.item
        pts = self.positions[items]

        if len(node.item) <= self.NperLeaf:
            node.bbox = self.compute_bbox(items)
            return
        
        # Find the best split using SAH
        best_cost = float('inf')
        best_split = None
        for i in range(len(pts)):
            left_items = items[pts[:, axis] <= pts[i, axis]]
            right_items = items[pts[:, axis] > pts[i, axis]]

            if len(left_items) == 0 or len(right_items) == 0:
                continue

            left_bbox = self.compute_bbox(left_items)
            right_bbox = self.compute_bbox(right_items)

            left_area = self.surface_area(left_bbox)
            right_area = self.surface_area(right_bbox)

            cost = len(left_items) * left_area + len(right_items) * right_area

            if cost < best_cost:
                best_cost = cost
                best_split = i

        if best_split is None:
            mid = len(items) // 2
            left_items = items[:mid]
            right_items = items[mid:]
        else:
            left_items = items[pts[:, axis] <= pts[best_split, axis]]
            right_items = items[pts[:, axis] > pts[best_split, axis]]
       
        # median = np.median(pts[:, axis])
        # left_indexes = items[pts[:, axis] <= median]
        # right_indexes = items[pts[:, axis] > median]

        # if len(left_indexes) == 0 or len(right_indexes) == 0:
        #     mid = len(items) // 2
        #     left_indexes = items[:mid]
        #     right_indexes = items[mid:]

        node.left = Nodes(item=left_items)
        node.right = Nodes(item=right_items)
        node.left.parent = node
        node.right.parent = node
        node.bbox = self.compute_bbox(items)

        self.build_recursive(node.left, k + 1)
        self.build_recursive(node.right, k + 1)
       
       
    def build(self):
        self.build_recursive(self.root, 0)

    def surface_area(self, bbox):
        d = bbox[1] - bbox[0]
        return 2 * (d[0]*d[1] + d[1]*d[2] + d[2]*d[0])

    def update(self, current=None):

        # On first call, start from the root
        if current is None:
            current = self.root

        # ---------------------------------------------------------
        # 1. If current is a leaf → compute its bbox
        # ---------------------------------------------------------
        if current.is_leaf():
            pos=self.positions[current.item[0]]
            rad=self.radii[current.item[0]]
            current.bbox = np.array([pos - rad, pos + rad])
            return

        # ---------------------------------------------------------
        # 2. Descend until we reach a leaf
        # ---------------------------------------------------------
        self.update(current.left)
        self.update(current.right)

        # ---------------------------------------------------------
        # 3. When both children done → compute internal bbox
        # ---------------------------------------------------------
        current.update_bbox()
        return


        

def print_bvh(node, depth=0):
        if node is None:
            return
        print("  " * depth + f"Node(depth={depth}, bbox=[[{node.bbox[0][0]:.2f}, {node.bbox[0][1]:.2f}, {node.bbox[0][2]:.2f}], [{node.bbox[1][0]:.2f}, {node.bbox[1][1]:.2f}, {node.bbox[1][2]:.2f}]], item={node.item})")
        print_bvh(node.left, depth + 1)
        print_bvh(node.right, depth + 1)



    
if __name__ == "__main__":
    # import matplotlib.pyplot as plt
    ### pts ###
    # pts = np.random.rand(10, 3) * 10
    pts = np.array([[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4], [5, 5, 5],])
    radii = np.array([0.5, 0.5, 0.5, 0.5, 0.5])
    # radii = np.random.rand(100) * 0.5 + 0.1
    bvh = BVH(pts, radii, NperLeaf=1)
    bvh.build()
    
    print_bvh(bvh.root)

    #plot points and bounding boxes
    # def plot_bbox(ax, bbox, color='r'):
    #     x_min, y_min, z_min = bbox[0]
    #     x_max, y_max, z_max = bbox[1]
    #     corners = np.array([[x_min, y_min, z_min],
    #                         [x_max, y_min, z_min],
    #                         [x_max, y_max, z_min],
    #                         [x_min, y_max, z_min],
    #                         [x_min, y_min, z_max],
    #                         [x_max, y_min, z_max],
    #                         [x_max, y_max, z_max],
    #                         [x_min, y_max, z_max]])
    #     edges = [(0, 1), (1, 2), (2, 3), (3, 0),
    #             (4, 5), (5, 6), (6, 7), (7, 4),
    #             (0, 4), (1, 5), (2, 6), (3, 7)]
    #     for edge in edges:
    #         ax.plot3D(*zip(corners[edge[0]], corners[edge[1]]), color=color)



    # def plot_sphere(ax, center, radius, color='g'):
    #     u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    #     x = center[0] + radius * np.cos(u) * np.sin(v)
    #     y = center[1] + radius * np.sin(u) * np.sin(v)
    #     z = center[2] + radius * np.cos(v)
    #     ax.plot_wireframe(x, y, z, color=color, alpha=0.3)


    # def plot_bvh(node):
        # if node is None:
            # return
        # plot_bbox(ax, node.bbox)
        # plot_bvh(node.left)
        # plot_bvh(node.right)
# 

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # # ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], color='b')
    # for i in range(len(pts)):
    #     plot_sphere(ax, pts[i], radii[i], color='b')



    # plot_bvh(bvh.root)
    # plt.savefig('BVH')
    # plt.show()