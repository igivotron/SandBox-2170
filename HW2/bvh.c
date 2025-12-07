
#include "bvh.h"
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


static inline uint32_t expandBits(uint32_t v) {
    //With v written in bits : b2 b1 b0 
    // expandBits return : (potential zeros here) b2 00 b1 00 b0
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

// Morton code 3D
uint32_t morton3D(double* positions, int i) {
    // returns x0 y0 z0 x1 y1 z1 ... in bits
    uint32_t x =(uint32_t) (positions[3*i] * (1023.0));
    uint32_t y =(uint32_t) (positions[3*i + 1] * (1023.0));
    uint32_t z =(uint32_t) (positions[3*i + 2] * (1023.0));
    return (expandBits(x) << 2) | (expandBits(y) << 1) | expandBits(z);
}

static int Comparaison_Morton(const void *a, const void *b) {
    const Morton_code *Mort_a = (const Morton_code*)a;
    const Morton_code *Mort_b = (const Morton_code*)b;
    if (Mort_a->morton < Mort_b->morton) return -1;
    if (Mort_a->morton > Mort_b->morton) return 1;
    return 0;
}


BVH* create_bvh(double* positions, double* radii, int n_points, int NperLeaf) {
    BVH* bvh = (BVH*)malloc(sizeof(BVH));

    //Morton Code 
    Morton_code* mortonIndices = (Morton_code*)malloc(sizeof(Morton_code) * n_points);
    bvh->morton_codes = mortonIndices;
    for (int i=0; i<n_points; i++){
        mortonIndices[i].index = i;
        mortonIndices[i].morton = morton3D(positions, i); 
    }

    bvh->positions = positions;//malloc(sizeof(double)*n_points*3);
    //memcpy(bvh->positions, sorted_positions, sizeof(double) * n_points * 3);
    bvh->radii = radii;//malloc(sizeof(double) * n_points);
    //memcpy(bvh->radii, sorted_radii, sizeof(double) * n_points);
    bvh->NperLeaf = NperLeaf;
    bvh->n_nodes = 0;
    bvh->nodes = malloc(sizeof(BVHNode) * (2 * n_points - 1)); // Max nodes in a binary tree


    // Initialize root node
    int * all_items = malloc(sizeof(int) * n_points);
    for (int i = 0; i < n_points; i++) all_items[i] = i;
    BVHNode root_node = create_node(0, all_items, n_points, -1, 0, n_points);
    compute_bbox(&root_node, positions, radii);
    bvh->nodes[bvh->n_nodes++] = root_node;
    bvh->root = root_node.index;
    free(all_items);

    return bvh; 
}

BVHNode create_node(int index, int* items, int n_items, int parent, int first, int last) {
    BVHNode node;
    node.index = index;
    node.bbox = (double*)malloc(sizeof(double) * 6);
    node.n_items = n_items;
    node.parent = parent;
    node.left = -1;
    node.right = -1;
    node.first = first;
    node.last = last;

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
    double* positions = bvh->positions;
    for (int i = 0; i < n_points; i++) {
        int orig_index = bvh->morton_codes[i].index;
        positions[3*i] = new_positions[3*orig_index];
        positions[3*i + 1] = new_positions[3*orig_index + 1];
        positions[3*i + 2] = new_positions[3*orig_index + 2];
    }
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

//returns the number of zero before the highest bit : 
// exemple : 0100 and 1000 returns 0 
static inline int count_prefix(uint32_t a, uint32_t b) {
    uint32_t x = a ^ b; // fonction XOR
    if (x == 0) return 32;
    return __builtin_clz(x);  // GCC/Clang builtin
}

int best_split(Morton_code* morton_c, int first, int last) {
    if (last - first == 2) return first + 1;

    uint32_t first_code = morton_c[first].morton;
    uint32_t last_code  = morton_c[last - 1].morton;

    // if identical codes, fallback to middle split
    if (first_code == last_code) {
        return (first + last) >> 1;
    }

    int common = count_prefix(first_code, last_code);

    int lo = first;
    int hi = last - 1;
    while (lo + 1 < hi) {
        int mid = (lo + hi) >> 1;
        int prefix = count_prefix(first_code, morton_c[mid].morton);
        if (prefix > common) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return hi;
}



void build_recursion(BVH* bvh, int node_index, int k){
    BVHNode* node = &bvh->nodes[node_index];
    int first = node->first;
    int last  = node->last;

    int best = best_split(bvh->morton_codes, first, last);
    if (is_leaf(node, bvh->NperLeaf)) {
        return; // Cannot split further
    }
    int left_index  = bvh->n_nodes++;
    int right_index = bvh->n_nodes++;

    bvh->nodes[left_index]  = create_node(left_index, NULL, 0, node_index, first, best);
    bvh->nodes[right_index] = create_node(right_index, NULL, 0, node_index, best, last);

    node->left  = left_index;
    node->right = right_index;
    compute_bbox(&bvh->nodes[left_index], bvh->positions, bvh->radii);
    compute_bbox(&bvh->nodes[right_index], bvh->positions, bvh->radii);
    
    build_recursion(bvh, left_index, k + 1);
    build_recursion(bvh, right_index, k + 1);
}

void build_bvh(BVH* bvh, int N) {
    qsort(bvh->morton_codes, N, sizeof(Morton_code), Comparaison_Morton);
    bvh->nodes[bvh->root].first = 0;
    bvh->nodes[bvh->root].last  = N;
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
    int N = sizeof(positions);
    build_bvh(bvh, N);
    printBVH(bvh, &bvh->nodes[bvh->root], 0);
    update(bvh, NULL);
    // Clean up
    free_bvh(bvh);
    free(positions);
    free(radii);
    return 0;
}