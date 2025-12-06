
#include "bvh.h"
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

BVH* create_bvh(double* positions, double* radii, int n_points, int NperLeaf) {
    BVH* bvh = (BVH*)malloc(sizeof(BVH));
    bvh->positions = malloc(sizeof(double) * n_points * 3);
    memcpy(bvh->positions, positions, sizeof(double) * n_points * 3);
    bvh->radii = malloc(sizeof(double) * n_points);
    memcpy(bvh->radii, radii, sizeof(double) * n_points);
    bvh->NperLeaf = NperLeaf;
    bvh->n_nodes = 0;
    bvh->nodes = malloc(sizeof(BVHNode) * (2 * n_points - 1)); // Max nodes in a binary tree


    // Initialize root node
    int * all_items = malloc(sizeof(int) * n_points);
    for (int i = 0; i < n_points; i++) all_items[i] = i;
    BVHNode root_node = create_node(0, all_items, n_points, -1);
    compute_bbox(&root_node, positions, radii);
    bvh->nodes[bvh->n_nodes++] = root_node;
    bvh->root = root_node.index;
    free(all_items);

    return bvh; 
}

BVHNode create_node(int index, int* items, int n_items, int parent) {
    BVHNode node;
    node.index = index;
    node.bbox = (double*)malloc(sizeof(double) * 6);
    node.n_items = n_items;
    node.parent = parent;
    node.left = -1;
    node.right = -1;

    if (n_items > 0) {
        node.items = (int*)malloc(sizeof(int) * n_items);
        memcpy(node.items, items, sizeof(int) * n_items);
    } else {
        node.items = NULL;
    }

    return node;
}

int is_leaf(BVHNode* node, int NperLeaf) {
    return node->n_items <= NperLeaf;
}

void update_positions(BVH* bvh, double* new_positions, int n_points) {
    memcpy(new_positions, bvh->positions, sizeof(double) * 3 * n_points);
    // bvh->positions = new_positions;
}

double* compute_bbox(BVHNode* node, double* positions, double* radii) {
    // Compute bounding box for the node based on its items positions and radii
    // Take the min and max in each dimension to form the bounding box
    if (node->n_items <= 0) return NULL;
    int* items = node->items;
    int n_items = node->n_items;
    // double* bbox = (double*)malloc(6 * sizeof(double)); // [min_x, min_y, min_z, max_x, max_y, max_z]
    double* bbox = node->bbox;
    double min_x, min_y, min_z;
    double max_x, max_y, max_z;

    for (int i = 0; i < n_items; i++) {
        int idx = items[i];
        double x = positions[3 * idx];
        double y = positions[3 * idx + 1];
        double z = positions[3 * idx + 2];
        double r = radii[idx];

        if (i == 0) {
            min_x = x - r; max_x = x + r;
            min_y = y - r; max_y = y + r;
            min_z = z - r; max_z = z + r;
        } else {
            if (x - r < min_x) min_x = x - r;
            if (x + r > max_x) max_x = x + r;
            if (y - r < min_y) min_y = y - r;
            if (y + r > max_y) max_y = y + r;
            if (z - r < min_z) min_z = z - r;
            if (z + r > max_z) max_z = z + r;
        }
    }
    bbox[0] = min_x; bbox[1] = min_y; bbox[2] = min_z;
    bbox[3] = max_x; bbox[4] = max_y; bbox[5] = max_z;
    return bbox;
}

double surface_area(double* bbox) {
    double dx = bbox[3] - bbox[0];
    double dy = bbox[4] - bbox[1];
    double dz = bbox[5] - bbox[2];
    return 2.0 * (dx * dy + dy * dz + dz * dx);
}

void split_items(BVHNode* node, double* positions, double* radii, int axis, int split_index,
                 int* left_items, int* right_items, int* n_left, int* n_right) {
    *n_left = 0;
    *n_right = 0;
    double threshold = positions[3 * node->items[split_index] + axis];

    for (int i = 0; i < node->n_items; i++) {
        int idx = node->items[i];
        double coord = positions[3*idx + axis];

        if (coord < threshold) left_items[(*n_left)++] = idx;
        else right_items[(*n_right)++] = idx;
    }
}
int best_split_axis(BVHNode* node, double* positions, double* radii, int axis,
                    int* left_items, int* right_items, int* n_left, int* n_right) {
    if (node->n_items == 2){
        left_items[0] = node->items[0];
        right_items[0] = node->items[1];
        *n_left = 1;
        *n_right = 1;
        return 0;
    }
    double best_cost = INFINITY;
    int best_index = -1;
    int n_items = node->n_items;

    int* temp_left = (int*)malloc(sizeof(int) * n_items);
    int* temp_right = (int*)malloc(sizeof(int) * n_items);
    BVHNode temp_node_left = create_node(-1, NULL, 0, -1);
    BVHNode temp_node_right = create_node(-1, NULL, 0, -1);

    if (!temp_left || !temp_right) return -1;

    for (int i = 0; i < n_items; i++) {
        int l = 0, r = 0;
        split_items(node, positions, radii, axis, i, temp_left, temp_right, &l, &r);
        if (l == 0 || r == 0) continue;
        temp_node_left.items = temp_left;
        temp_node_left.n_items = l;
        temp_node_right.items = temp_right;
        temp_node_right.n_items = r;
        compute_bbox(&temp_node_left, positions, radii);
        compute_bbox(&temp_node_right, positions, radii);
        double* left_bbox = temp_node_left.bbox;
        double* right_bbox = temp_node_right.bbox;
        // double* left_bbox = compute_bbox(&(BVHNode){.items = temp_left, .n_items = l}, positions, radii);
        // double* right_bbox = compute_bbox(&(BVHNode){.items = temp_right, .n_items = r}, positions, radii);

        double cost = surface_area(left_bbox) * l + surface_area(right_bbox) * r;
        if (cost < best_cost) {
            best_cost = cost;
            best_index = i;
            *n_left = l;
            *n_right = r;
            memcpy(left_items, temp_left, sizeof(int) * l);
            memcpy(right_items, temp_right, sizeof(int) * r);
        }

        // free(left_bbox);
        // free(right_bbox);
    }

    free_node(&temp_node_left);
    free_node(&temp_node_right);
    // free(temp_left);
    // free(temp_right);
    return best_index;
}

void build_recursion(BVH* bvh, int node_index, int k){
    int axis = k % 3;
    BVHNode* node = &bvh->nodes[node_index];
    if (is_leaf(node, bvh->NperLeaf)){
        // node->bbox = compute_bbox(node, bvh->positions, bvh->radii);
        return;
    }

    // int * items = node->items;
    int n_items = node->n_items;
    double* positions = bvh->positions;
    double* radii = bvh->radii;
    int* left_items = malloc(sizeof(int) * n_items);
    int* right_items = malloc(sizeof(int) * n_items);
    int n_left = 0;
    int n_right = 0;
    int split_index = best_split_axis(node, positions, radii, axis, left_items, right_items, &n_left, &n_right);
    if (n_left == 0 || n_right == 0 || split_index == -1) {
        free(left_items);
        free(right_items);
        return; // Cannot split further
    }

    // double* left_bbox = compute_bbox(&(BVHNode){.items = left_items, .n_items = n_left}, positions, radii);
    // double* right_bbox = compute_bbox(&(BVHNode){.items = right_items, .n_items = n_right}, positions, radii);
    // int left_index = bvh->n_nodes++;
    // int right_index = bvh->n_nodes++;
    // BVHNode left_node = create_node(left_index, left_bbox, left_items, n_left, node_index);
    // BVHNode right_node = create_node(right_index, right_bbox, right_items, n_right, node_index);
    int left_index = bvh->n_nodes++;
    int right_index = bvh->n_nodes++;
    BVHNode left_node = create_node(left_index, left_items, n_left, node_index);
    BVHNode right_node = create_node(right_index, right_items, n_right, node_index);
    compute_bbox(&left_node, positions, radii);
    compute_bbox(&right_node, positions, radii);
    bvh->nodes[left_index] = left_node;
    bvh->nodes[right_index] = right_node;
    node->left = left_index;
    node->right = right_index;
    free(left_items);
    free(right_items);
    build_recursion(bvh, left_index, k + 1);
    build_recursion(bvh, right_index, k + 1);
}

void build_bvh(BVH* bvh) {
    build_recursion(bvh, bvh->root, 0);
}

void update_bbox(BVH* bvh, BVHNode* node){
    double* left_bbox = bvh->nodes[node->left].bbox;
    double* right_bbox = bvh->nodes[node->right].bbox;
    double* bbox = node->bbox;
    bbox[0] = MIN(left_bbox[0], right_bbox[0]);
    bbox[1] = MIN(left_bbox[1], right_bbox[1]);
    bbox[2] = MIN(left_bbox[2], right_bbox[2]);
    bbox[3] = MAX(left_bbox[3], right_bbox[3]);
    bbox[4] = MAX(left_bbox[4], right_bbox[4]);
    bbox[5] = MAX(left_bbox[5], right_bbox[5]);
}

void update(BVH* bvh, BVHNode* current) {
    if (current == NULL) current = &bvh->nodes[bvh->root];
    // 1. If current is a leaf → compute its bbox
    if (is_leaf(current, bvh->NperLeaf)) {
        int idx = current->items[0];
        double radius = bvh->radii[idx];
        current->bbox[0] = bvh->positions[3 * idx] - radius;
        current->bbox[1] = bvh->positions[3 * idx + 1] - radius;
        current->bbox[2] = bvh->positions[3 * idx + 2] - radius;
        current->bbox[3] = bvh->positions[3 * idx] + radius;
        current->bbox[4] = bvh->positions[3 * idx + 1] + radius;
        current->bbox[5] = bvh->positions[3 * idx + 2] + radius;
        return;
    }

    // // 2. Descend until we reach a leaf
    update(bvh, &bvh->nodes[current->left]);
    update(bvh, &bvh->nodes[current->right]);
    update_bbox(bvh, current);
}

int bbox_intersect(double *bb1, double *bb2){
    //box intersection test
    return !(bb1[3] < bb2[0] || bb1[0] > bb2[3] || bb1[4] < bb2[1] || bb1[1] > bb2[4] || bb1[5] < bb2[2] || bb1[2] > bb2[5]);
}

int find_pot_inter(BVH* bvh, int* pot_cont){
    // find potential contacts and store them in pot_cont as pair of indices
    int pcs = 0; // potential contact size
    int* stack; // Ce sont des indices par pair
    int stack_size = 0;
    stack = (int*)malloc(sizeof(int) * 200);
    stack[stack_size++] = bvh->root;
    stack[stack_size++] = bvh->root;
    while (stack_size > 0) {
        int A = stack[--stack_size];
        int B = stack[--stack_size];
        if (A == B) {
            BVHNode* node = &bvh->nodes[A];
            if (!is_leaf(node, bvh->NperLeaf)) {
                stack[stack_size++] = node->left;
                stack[stack_size++] = node->left;
                stack[stack_size++] = node->right;
                stack[stack_size++] = node->right;
                stack[stack_size++] = node->left;
                stack[stack_size++] = node->right;
            }
            continue;
        }
        //then A != B
        BVHNode* nodeA = &bvh->nodes[A];
        BVHNode* nodeB = &bvh->nodes[B];
        if (is_leaf(nodeA, bvh->NperLeaf) && is_leaf(nodeB, bvh->NperLeaf)) {
            //both are leaves, add to potential contacts
            pot_cont[pcs++] = nodeA->items[0];
            pot_cont[pcs++] = nodeB->items[0];
            continue;
        }
        if (bbox_intersect(nodeA->bbox, nodeB->bbox)) {
            //boxes intersect, descend
            if (is_leaf(nodeA, bvh->NperLeaf)) {
                stack[stack_size++] = A;
                stack[stack_size++] = nodeB->left;
                stack[stack_size++] = A;
                stack[stack_size++] = nodeB->right;
            }
            else if (is_leaf(nodeB, bvh->NperLeaf)) {
                stack[stack_size++] = B;
                stack[stack_size++] = nodeA->left;
                stack[stack_size++] = B;
                stack[stack_size++] = nodeA->right;
            }
            else {
                stack[stack_size++] = nodeA->left;
                stack[stack_size++] = nodeB->left;
                stack[stack_size++] = nodeA->left;
                stack[stack_size++] = nodeB->right;
                stack[stack_size++] = nodeA->right;
                stack[stack_size++] = nodeB->left;
                stack[stack_size++] = nodeA->right;
                stack[stack_size++] = nodeB->right;
            }
        }
    }
    free(stack);
    return pcs;
}

void flat_bvh(BVH* bvh, float* flat_bvh){
    // Flatten the BVH into a linear array for GPU usage
    // Each node takes 3 vec4 (12 floats)
    // vec4 0: left index, right index, item + padding
    // vec4 1: bbox min (x,y,z) + padding
    // vec4 2: bbox max (x,y,z) + n_items
    int n_nodes = bvh->n_nodes;
    for (int i = 0; i < n_nodes; i++) {
        BVHNode* node = &bvh->nodes[i];
        int base = i * 12; // 3 vec4 per node

        flat_bvh[base + 0] = (float)node->left;
        flat_bvh[base + 1] = (float)node->right;
        if (is_leaf(node, bvh->NperLeaf)) {
            flat_bvh[base + 2] = (float)node->items[0];
        } else {
            flat_bvh[base + 2] = -1.0f; // Not a leaf
        }
        flat_bvh[base + 3] = 0.0f; // padding

        flat_bvh[base + 4] = node->bbox[0];
        flat_bvh[base + 5] = node->bbox[1];
        flat_bvh[base + 6] = node->bbox[2];
        flat_bvh[base + 7] = 0.0f; // padding

        flat_bvh[base + 8] = node->bbox[3];
        flat_bvh[base + 9] = node->bbox[4];
        flat_bvh[base + 10] = node->bbox[5];
        flat_bvh[base + 11] = 0.0f; // padding
    }
    return;
}

void free_bvh(BVH* bvh) {
    if (bvh) {
        if (bvh->positions) free(bvh->positions);
        if (bvh->radii) free(bvh->radii);
        if (bvh->nodes) {
            for (int i = 0; i < bvh->n_nodes; i++) {
                if (bvh->nodes[i].bbox) free(bvh->nodes[i].bbox);
                // free(bvh->nodes[i].bbox);
                free(bvh->nodes[i].items);
            }
            free(bvh->nodes);
        }
        free(bvh);
    }
}

void free_node(BVHNode* node) {
    if (node) {
        free(node->bbox);
        free(node->items);
        // free(node->index);
    }
}

void printBVH(BVH* bvh, BVHNode* node, int depth) {
    if (node == NULL) return;
    printf("%*sNode(depth=%d, bbox=[[%.2f, %.2f, %.2f], [%.2f, %.2f, %.2f]], n_items=%d)\n", depth * 2, "", depth,
           node->bbox[0], node->bbox[1], node->bbox[2],
           node->bbox[3], node->bbox[4], node->bbox[5],
           node->n_items);
    if (node->left != -1) {printBVH(bvh, &bvh->nodes[node->left], depth + 1);}
    if (node->right != -1) {printBVH(bvh, &bvh->nodes[node->right], depth + 1);}

}

void printBVH2(BVH* bvh) {
    for (int i = 0; i < bvh->n_nodes; i++) {
        BVHNode* node = &bvh->nodes[i];
        printf("Node %d: bbox=[[%.2f, %.2f, %.2f], [%.2f, %.2f, %.2f]], n_items=%d, left=%d, right=%d, parent=%d\n",
               node->index,
               node->bbox[0], node->bbox[1], node->bbox[2],
               node->bbox[3], node->bbox[4], node->bbox[5],
               node->n_items,
               node->left,
               node->right,
               node->parent);
    }
}

int main() {
    // Example usage of BVH
    int n_points = 10;
    double* positions = (double*)malloc(sizeof(double) * n_points * 3);
    double* radii = (double*)malloc(sizeof(double) * n_points);
    for (int i = 0; i < n_points; i++) {
        positions[3*i] = rand() % 100;
        positions[3*i + 1] = rand() % 100;
        positions[3*i + 2] = rand() % 100;
        radii[i] = (rand() % 10) + 1;
    }

    // double positions_array[] = {1.0, 1.0, 1.0,
    //                             2.0, 2.0, 2.0,
    //                             3.0, 3.0, 3.0,
    //                             4.0, 4.0, 4.0,
    //                             5.0, 5.0, 5.0};
    // double radii_array[] = {0.5, 0.5, 0.5, 0.5, 0.5};
    // int n_points = 5;
    // double* positions = (double*)malloc(sizeof(double) * n_points * 3);
    // double* radii = (double*)malloc(sizeof(double) * n_points);
    // memcpy(positions, positions_array, sizeof(double) * n_points * 3);
    // memcpy(radii, radii_array, sizeof(double) * n_points);


    BVH* bvh = create_bvh(positions, radii, n_points, 1);
    build_bvh(bvh);
    printBVH(bvh, &bvh->nodes[bvh->root], 0);
    update(bvh, NULL);
    // Clean up
    free_bvh(bvh);
    free(positions);
    free(radii);
    return 0;
}