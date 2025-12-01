#ifndef BVH_H
#define BVH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


typedef struct BVHNode {
    int left;   // Index of the left child node
    int right;  // Index of the right child node
    int parent; // Index of the parent node
    double * bbox; // Bounding box: [min_x, min_y, min_z, max_x, max_y, max_z]
    int * items; // Indices of points contained in this node
    int n_items; // Number of items in this node
    int index;   // Index of this node
    int state; // State of the node: 0 = not updated, 1 = updated
} BVHNode;


typedef struct BVH {
    BVHNode* nodes; // Array of BVH nodes
    int n_nodes;     // Number of nodes in the BVH
    double* positions; // Array of point positions (x, y, z)*n_points
    double* radii;    // Array of point radii n_points
    int NperLeaf;    // Maximum number of points per leaf node
    int root;        // Index of the root node
} BVH;

BVH* create_bvh(double* positions, double* radii, int n_points, int NperLeaf);
BVHNode create_node(int index, double* bbox, int* items, int n_items, int parent);
int is_leaf(BVHNode* node, int NperLeaf);
double* compute_bbox(BVHNode* node, double* positions, double* radii);
double* combine_bboxes(double* bbox1, double* bbox2);
double surface_area(double* bbox);
void build_bvh(BVH* bvh);
void build_recursion(BVH* bvh, int node_index, int k);
void split_items(BVHNode* node, double* positions, double* radii, int axis, int split_index, int* left_items, int* right_items, int* n_left, int* n_right);
int best_split_axis(BVHNode* node, double* positions, double* radii, int axis, int* left_items, int* right_items, int* n_left, int* n_right);
void update_bbox(BVH* bvh, BVHNode* node);
void update(BVH* bvh, BVHNode* current);
void update_positions(BVH* bvh, double* new_positions);
int bbox_intersect(double *bb1, double *bb2);
int find_pot_inter(BVH* bvh, int* pot_cont);

void printBVH(BVH* bvh, BVHNode* node, int depth);
void free_bvh(BVH* bvh);
void free_node(BVHNode* node);


#endif // BVH_H