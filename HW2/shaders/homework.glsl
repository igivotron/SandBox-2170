#define STACK_MAX 32u

/* import data */
layout(set = 2, binding = 0, std430) buffer DataStructure {
    vec4 bvh_nodes[];
};

bool ray_bbox_intersection(vec3 inv_dir, vec3 origin, vec3 bbox_min, vec3 bbox_max) {
    vec3 t0s = (bbox_min - origin) * inv_dir;
    vec3 t1s = (bbox_max - origin) * inv_dir;
    vec3 tmin = min(t0s, t1s);
    vec3 tmax = max(t0s, t1s);
    float t_enter = max(max(tmin.x, tmin.y), tmin.z);
    float t_exit  = min(min(tmax.x, tmax.y), tmax.z);
    return t_exit >= max(t_enter, 0.0);
}

Intersection trace_ray_bvh(Ray ray) {
    Intersection result = EMPTY_INTERSECTION;
    vec3 inv_dir = 1.0 / ray.direction;

    uint stack[STACK_MAX];
    uint sp = 0u;
    stack[sp++] = 0u; // root node index 

    while (sp > 0u) {
        uint node_idx = stack[--sp];
        uint base = node_idx*3u;

        vec4 v0 = bvh_nodes[base + 0u];
        vec4 v1 = bvh_nodes[base + 1u];
        vec4 v2 = bvh_nodes[base + 2u];

        int item = int(v0.z);
        if (item != -1) {
            ray_ball_intersection(ray, item, result);
            continue;
        }

        uint left  = uint(v0.x);
        uint right = uint(v0.y);

        vec3 bbox_min = v1.xyz;
        vec3 bbox_max = v2.xyz;

        if (!ray_bbox_intersection(inv_dir, ray.origin, bbox_min, bbox_max)) continue;

        if (sp + 2u >= STACK_MAX) break;
        stack[sp++] = left;
        stack[sp++] = right;
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
