#ifndef BVH_H
#define BVH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>

typedef struct {
    int index;
    uint32_t morton;
} Morton_code;

typedef struct BVHNode {
    int left;   // Index of the left child node
    int right;  // Index of the right child node
    int parent; // Index of the parent node
    double * bbox; // Bounding box: [min_x, min_y, min_z, max_x, max_y, max_z]
    int * items; // Indices of points contained in this node
    int n_items; // Number of items in this node
    int index;   // Index of this node
    int first;
    int last;
} BVHNode;


typedef struct BVH {
    BVHNode* nodes; // Array of BVH nodes
    int n_nodes;     // Number of nodes in the BVH
    double* positions; // Array of point positions (x, y, z)*n_points
    double* radii;    // Array of point radii n_points
    int NperLeaf;    // Maximum number of points per leaf node
    int root;        // Index of the root node
    Morton_code* morton_codes; // Morton codes of the items
} BVH;


static inline uint32_t expandBits(uint32_t v);
uint32_t morton3D(double* positions, int i);
static int Comparaison_Morton(const void *a, const void *b);
BVH* create_bvh(double* positions, double* radii, int n_points, int NperLeaf);
BVHNode create_node(int index, int* items, int n_items, int parent, int first, int last);
int is_leaf(BVHNode* node, int NperLeaf);
double* compute_bbox(BVHNode* node, double* positions, double* radii);
double surface_area(double* bbox);
static inline int count_prefix(uint32_t a, uint32_t b);
int best_split(Morton_code* morton_c, int first, int last);
void build_bvh(BVH* bvh, int N);
void build_recursion(BVH* bvh, int node_index, int k);
void update_bbox(BVH* bvh, BVHNode* node);
void update(BVH* bvh, BVHNode* current);
void update_positions(BVH* bvh, double* new_positions, int n_points);
int bbox_intersect(double *bb1, double *bb2);
int find_pot_inter(BVH* bvh, int* pot_cont);
void flat_bvh(BVH* bvh, float* flat_bvh);

void printBVH(BVH* bvh, BVHNode* node, int depth);
void printBVH2(BVH* bvh);
void free_bvh(BVH* bvh);
void free_node(BVHNode* node);


#endif // BVH_H