struct BVHNode{
  int left;      // index of left child (-1 if leaf)
  int right;     // index of right child (-1 if leaf)
  int parent;    // index of parent node (-1 if root)
  int pad0;      // unused (for memory alignment)
  vec4 bbox_min; // xyz: min coords bbox, w: unused
  vec4 bbox_max; // xyz: max coords bbox, w: unused
  int n_items;   // number of items in the node (1 if leaf, 0 otherwise)
  int item;      // index of the ball if leaf, -1 otherwise
  int index;     // index of this node in the array
  int pad1;      // unused (for memory alignment)

};

layout(set = 2, binding = 0, std430) buffer DataStructure {
  BVHNode bvh_nodes[];
}; 

bool intersect_aabb(Ray ray, vec4 bbox_min, vec4 bbox_max) {
    vec3 inv_dir = 1.0 / ray.direction;
    vec3 t0s = (bbox_min.xyz - ray.origin) * inv_dir;
    vec3 t1s = (bbox_max.xyz - ray.origin) * inv_dir;

    vec3 tmin = min(t0s, t1s);
    vec3 tmax = max(t0s, t1s);

    float t_enter = max(max(tmin.x, tmin.y), tmin.z);
    float t_exit  = min(min(tmax.x, tmax.y), tmax.z);

    return t_exit >= max(t_enter, 0.0);
}

Intersection trace_ray_bvh(Ray ray) {
    Intersection result = EMPTY_INTERSECTION;

    int stack[64]; // assuming a maximum depth of 64
    int sp = 0;

    // On empile le root node
    stack[sp++] = 0; 

    while (sp > 0) {
        int node_idx = stack[--sp];
        BVHNode node = bvh_nodes[node_idx];

        // Test AABB
        if (!intersect_aabb(ray, node.bbox_min, node.bbox_max))
          continue;

        if (node.n_items == 1) {
          // feuille : tester la sphère
          ray_ball_intersection(ray, node.item, result);
        } else {
          BVHNode l = bvh_nodes[node.left];
          BVHNode r = bvh_nodes[node.right];
          r = node.right;
          if (intersect_aabb(ray, l.bbox_min, l.bbox_max))
            stack[sp++] = l.index;
          else if (intersect_aabb(ray, r.bbox_min, r.bbox_max))
            stack[sp++] = r.index;
        }
    }

    return result;
}

Intersection trace_ray(Ray ray) {
  const Plane planes[6] = Plane[](
    Plane(vec3(0.0, +10.0, 0.0), vec3(0.0, -1.0, 0.0),
          Material(vec3(1.0, 1.0, 1.0), ROUGHNESS,
                   vec3(BRIGHTNESS), 0.0)),
    Plane(vec3(0.0, -10.0, 0.0), vec3(0.0, +1.0, 0.0),
          Material(vec3(1.0, 1.0, 1.0), ROUGHNESS,
                   vec3(0.0), 0.0)),
    Plane(vec3(+10.0, 0.0, 0.0), vec3(-1.0, 0.0, 0.0),
          Material(vec3(0.0, 1.0, 0.0), ROUGHNESS,
                   vec3(0.0), 0.0)),
    Plane(vec3(-10.0, 0.0, 0.0), vec3(+1.0, 0.0, 0.0),
          Material(vec3(0.0, 0.0, 1.0), ROUGHNESS,
                   vec3(0.0), 0.0)),
    Plane(vec3(0.0, 0.0, +10.0), vec3(0.0, 0.0, -1.0),
          Material(vec3(1.0, 0.0, 0.0), ROUGHNESS,
                   vec3(0.0), 0.0)), // back
    Plane(vec3(0.0, 0.0, -10.0), vec3(0.0, 0.0, +1.0),
          Material(vec3(1.0, 1.0, 1.0), 0.05,
                   vec3(0.0), 1.0)));



  Intersection result = EMPTY_INTERSECTION;
  Intersection new_result = EMPTY_INTERSECTION;
  for (uint i = 0u; i < 6u; i++) {
    if (ray_plane_intersection(ray, planes[i], new_result))
      merge_intersection(result, new_result);
  }

  /*
   * TODO: For the GPU part of the homework: Replace this code with a function
   * that uses an acceleration data structure.
   */
  merge_intersection(result, trace_ray_bvh(ray));
  // for (uint i = 0u; i < materials.length(); i++) {
  //   ray_ball_intersection(ray, int(i), result);
  // }

  return result;
}

/*
 * The code below is an example of how you may access your data structure (see commented code in
 * homework.py for an example of how to transfer your data from the CPU to the GPU).
 */

/* layout(set = 2, binding = 0, std430) buffer DataStructure { */
/*   vec2 example_data[]; */
/* }; */
