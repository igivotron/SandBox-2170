import numpy as np
# from simulator import Simulator, Contact, copy_buffer
# from homework import *


#Global constants for cost calculation
C_t = 1.0 # cost of traversal
C_i = 1.0 # cost of intersection
N = 1     # number of primitives per leaf


class Nodes:
    def __init__(self, item, left=None, right=None, bbox=None, index=None):
        self.left = left
        self.right = right
        self.parent = None
        self.bbox = bbox #[[x_min, y_min, z_min], [x_max, y_max, z_max]]
        self.index = index # each node of the bvh has an index (used for flat representation)
        self.item = item #contain the index of the item if leaf
        self.cost = None
        
    def is_leaf(self):
        return self.left is None and self.right is None

    def update_bbox(self):
        bbox_min = np.min([self.left.bbox[0], self.right.bbox[0]], axis=0)
        bbox_max = np.max([self.left.bbox[1], self.right.bbox[1]], axis=0)
        self.bbox = np.array([bbox_min, bbox_max])
        return
    
    def surface_area(self):
        d = self.bbox[1] - self.bbox[0]
        return 2 * (d[0]*d[1] + d[1]*d[2] + d[2]*d[0])
    
    def update_cost(self):
        if self.is_leaf() :
            self.cost = C_t + C_i * N
            return self.cost
        C_L = self.left.update_cost()
        C_R = self.right.update_cost()
        S_L = self.left.surface_area()
        S_R = self.right.surface_area()
        S = self.surface_area()
        self.cost = C_t + (S_L * C_L + S_R * C_R)/S
        return self.cost
        
    def get_cost_quick(self):
        C_L = self.left.cost
        C_R = self.right.cost
        S_L = self.left.surface_area()
        S_R = self.right.surface_area()
        S = self.surface_area()
        return C_t + (S_L * C_L + S_R * C_R)/S
        return self.cost

    def get_cost_recursive(self):
        if self.is_leaf() :
            return C_t + C_i * N
        C_L = self.get_cost_recursive(self.left)
        C_R = self.get_cost_recursive(self.right)
        S_L = self.surface_area(self.left.bbox)
        S_R = self.surface_area(self.right.bbox)
        S = self.surface_area(self.bbox)
        return C_t + (S_L * C_L + S_R * C_R)/S

class BVH:
    def __init__(self, positions, radii, NperLeaf=1):
        self.positions = positions
        self.radii = radii
        self.N_items = len(positions)
        self.indexes = np.arange(len(positions))
        self.NperLeaf = NperLeaf
        self.root = Nodes(item=self.indexes)
        self.nNode = 1

    def compute_bbox(self, items):
        x_min = np.min(self.positions[items, 0] - self.radii[items])
        x_max = np.max(self.positions[items, 0] + self.radii[items])
        y_min = np.min(self.positions[items, 1] - self.radii[items])
        y_max = np.max(self.positions[items, 1] + self.radii[items])
        z_min = np.min(self.positions[items, 2] - self.radii[items])
        z_max = np.max(self.positions[items, 2] + self.radii[items])
        bbox = np.array([[x_min, y_min, z_min], [x_max, y_max, z_max]])
        return bbox
    
    def combine_bbox(self, left_box, right_box):
        bbox_min = np.min([left_box[0], right_box[0]], axis=0)
        bbox_max = np.max([left_box[1], right_box[1]], axis=0)
        return np.array([bbox_min, bbox_max])

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
        best_axis = None
        for axis in range(3):
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
                    best_axis = axis

        if best_split is None:
            mid = len(items) // 2
            left_items = items[:mid]
            right_items = items[mid:]
        else:
            left_items = items[pts[:, best_axis] <= pts[best_split, best_axis]]
            right_items = items[pts[:, best_axis] > pts[best_split, best_axis]]
       
        # median = np.median(pts[:, axis])
        # left_indexes = items[pts[:, axis] <= median]
        # right_indexes = items[pts[:, axis] > median]

        # if len(left_indexes) == 0 or len(right_indexes) == 0:
        #     mid = len(items) // 2
        #     left_indexes = items[:mid]
        #     right_indexes = items[mid:]

        node.left = Nodes(item=left_items,index=self.nNode)
        self.nNode+=1
        node.right = Nodes(item=right_items,index=self.nNode)
        self.nNode+=1
        node.left.parent = node
        node.right.parent = node
        node.bbox = self.compute_bbox(items)

        self.build_recursive(node.left, k + 1)
        self.build_recursive(node.right, k + 1)
           
    def build(self):
        self.nNode = 1
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

    def flat_bvh(self):
        """Return a flat array representation of the BVH for GPU usage."""
        # each node is stored as:
        # [left_index, right_index, min_x, min_y, min_z, max_x, max_y, max_z, item_index]
        flat_bvh = np.zeros((self.nNode, 9), dtype=np.float32) 
        def traverse(node):
            i = node.index
            if node.is_leaf():
                flat_bvh[i] = [-1, -1,
                                 node.bbox[0][0], node.bbox[0][1], node.bbox[0][2],
                                 node.bbox[1][0], node.bbox[1][1], node.bbox[1][2],
                                 node.item[0]]
                return
            flat_bvh[i] = [node.left.index, node.right.index,
                          node.bbox[0][0], node.bbox[0][1], node.bbox[0][2],
                          node.bbox[1][0], node.bbox[1][1], node.bbox[1][2],
                          -1]
            traverse(node.left)
            traverse(node.right)

        traverse(self.root)
        return flat_bvh

    def get_cost(self, node) :
        '''
        C_t = relative cost for a traversal step
        C_i = relative cost for a primitive (sphere) intersection 
        N = primitive per leaf
        '''
        if node.is_leaf() :
            return C_t + C_i * N
        C_L = self.get_cost(node.left)
        C_R = self.get_cost(node.right)
        S_L = self.surface_area(node.left.bbox)
        S_R = self.surface_area(node.right.bbox)
        S = self.surface_area(node.bbox)
        return C_t + (S_L * C_L + S_R * C_R)/S
    
    def rotation(self, node):
        ####################
        #        1
        #       / \
        #      /   \
        #     2     3
        #    / \   / \
        #   4   5 6   7
        
        node1 = node
        node2 = node.left
        node3 = node.right
        node4 = node.left.left
        node5 = node.left.right
        node6 = node.right.left
        node7 = node.right.right

        C0 = node.cost
        C_rot1 = float('inf') ; C_rot2 = float('inf') ; C_rot3= float('inf'); C_rot4= float('inf') 

        #ROTATION 1 (2<->6)
        if node6 is not None :
            S3_= self.surface_area(self.combine_bbox(node2.bbox, node7.bbox))
            C3_ = C_t + (node2.cost*node2.surface_area()+ node7.cost*node7.surface_area())/S3_
            C_rot1 = C_t + (node6.surface_area()*node6.cost+S3_*C3_)/node1.surface_area()

        #ROTATION 2 (2<->7)
        if node7 is not None : 
            S3__= self.surface_area(self.combine_bbox(node2.bbox, node6.bbox))
            C3__ = C_t + (node2.cost*node2.surface_area()+ node6.cost*node6.surface_area())/S3__
            C_rot2 = C_t + (node7.surface_area()*node7.cost+S3__*C3__)/node1.surface_area()

        #ROTATION 3 (3<->4)
        if node4 is not None : 
            S2_ = self.surface_area( self.combine_bbox(node3.bbox, node5.bbox))
            C2_ = C_t + (node3.cost*node3.surface_area()+ node5.cost*node5.surface_area())/S2_
            C_rot3 = C_t + (node4.surface_area()*node4.cost+S2_*C2_)/node1.surface_area()

        #ROTATION 4 (3<->5)
        if node5 is not None : 
            S2__ = self.surface_area(self.combine_bbox(node3.bbox, node4.bbox))
            C2__ = C_t + (node3.cost*node3.surface_area()+ node4.cost*node4.surface_area())/S2__
            C_rot4 = C_t + (node5.surface_area()*node5.cost+S2__*C2__)/node1.surface_area()

        C_min =[C0, C_rot1, C_rot2, C_rot3, C_rot4]
        best = int(np.argmin(C_min))
        # print("rotation costs:",[C0, C_rot1, C_rot2, C_rot3, C_rot4])

        if best == 1 :
            node.left = node6
            node6.parent = node
            node.right.left = node2
            node2.parent = node.right
            node.right.update_bbox()
            # print('rotation 1')

        if best == 2 : 
            node.left = node7
            node7.parent = node
            node.right.right = node2
            node2.parent = node.right
            node.right.update_bbox()
            # print('rotation 2')

        if best == 3 : 
            node.right = node4
            node4.parent = node
            node.left.left = node3 
            node3.parent = node.left
            node.left.update_bbox()
            # print('rotation 3')

        if best == 4 :
            node.right = node5
            node5.parent = node 
            node.left.right = node3 
            node3.parent = node.left
            node.left.update_bbox()
            # print('rotation 4')
        # print("Cost before:", C0)
        # print("Cost after :", C_min[best])
        node.cost=C_min[best]
        return node.cost
        
    
    def optimize_w_rotation(self, node, max_depth, depth=0):
        '''
        For a stable BVH=> height of the tree h = log_2(f)
        with f the number of leafs (N = 2f -1)
        '''
        if depth >= max_depth or node.is_leaf():
            node.update_cost()
            return node.cost

        # go deeper recursively
        C_L=self.optimize_w_rotation(node.left, max_depth, depth + 1)
        C_R=self.optimize_w_rotation(node.right, max_depth, depth + 1)
        S_L = node.left.surface_area()
        S_R = node.right.surface_area()
        S = node.surface_area()
        node.cost = C_t + (S_L * C_L + S_R * C_R)/S
        return self.rotation(node)

def print_bvh(node, depth=0):
        if node is None:
            return
        print("  " * depth + f"Node(index={node.index}, bbox=[[{node.bbox[0][0]:.2f}, {node.bbox[0][1]:.2f}, {node.bbox[0][2]:.2f}], [{node.bbox[1][0]:.2f}, {node.bbox[1][1]:.2f}, {node.bbox[1][2]:.2f}]], item={node.item})")
        print_bvh(node.left, depth + 1)
        print_bvh(node.right, depth + 1)



    
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    ### pts ###
    n=3
    pts = np.random.rand(n, 3) * 10
    radii = np.random.rand(n) * 0.5 + 0.1
    bvh = BVH(pts, radii, NperLeaf=1)
    bvh.build()
    
    print_bvh(bvh.root)
    print(bvh.flat_bvh())

    #plot points and bounding boxes
    def plot_bbox(ax, bbox, color='r'):
        x_min, y_min, z_min = bbox[0]
        x_max, y_max, z_max = bbox[1]
        corners = np.array([[x_min, y_min, z_min],
                            [x_max, y_min, z_min],
                            [x_max, y_max, z_min],
                            [x_min, y_max, z_min],
                            [x_min, y_min, z_max],
                            [x_max, y_min, z_max],
                            [x_max, y_max, z_max],
                            [x_min, y_max, z_max]])
        edges = [(0, 1), (1, 2), (2, 3), (3, 0),
                (4, 5), (5, 6), (6, 7), (7, 4),
                (0, 4), (1, 5), (2, 6), (3, 7)]
        for edge in edges:
            ax.plot3D(*zip(corners[edge[0]], corners[edge[1]]), color=color)



    def plot_sphere(ax, center, radius, color='g'):
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = center[0] + radius * np.cos(u) * np.sin(v)
        y = center[1] + radius * np.sin(u) * np.sin(v)
        z = center[2] + radius * np.cos(v)
        ax.plot_wireframe(x, y, z, color=color, alpha=0.3)


    def plot_bvh(node):
        if node is None:
            return
        plot_bbox(ax, node.bbox)
        plot_bvh(node.left)
        plot_bvh(node.right)


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], color='b')
    for i in range(len(pts)):
        plot_sphere(ax, pts[i], radii[i], color='b')


    print(bvh.nNode)
    plot_bvh(bvh.root)
    plt.savefig('BVH')
    plt.show()