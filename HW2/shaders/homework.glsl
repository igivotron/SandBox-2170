#define get(A,i,j) A[ 1 + (i)*9u + (j) ] 
#define STACK_MAX 32u

/* import data */
layout(set = 2, binding = 0, std430) buffer DataStructure {
  float bvh[];
};

bool ray_bbox_intersection(Ray ray, float x_min, float y_min, float z_min,
                           float x_max, float y_max, float z_max) {
  vec3 bbox_min = vec3(x_min, y_min, z_min);
  vec3 bbox_max = vec3(x_max, y_max, z_max);
  
  vec3 inv_dir = 1.0 / ray.direction;
  vec3 t0s = (bbox_min - ray.origin) * inv_dir;
  vec3 t1s = (bbox_max - ray.origin) * inv_dir;
  vec3 tsmaller = min(t0s, t1s);
  vec3 tbigger = max(t0s, t1s);

  float tmin = max(max(tsmaller.x, tsmaller.y), tsmaller.z);
  float tmax = min(min(tbigger.x, tbigger.y), tbigger.z);

  return tmax >= max(tmin, 0.0);
}

// bool ray_bbox_intersection(Ray ray, float x_min, float y_min, float z_min,
//                            float x_max, float y_max, float z_max) {
//   //rayon de la sqphere antourant la bbox
//   float radius = length(vec3(x_max - x_min, y_max - y_min, z_max - z_min)) * 0.5;
//   vec3 center = vec3((x_min + x_max) * 0.5, (y_min + y_max) * 0.5, (z_min + z_max) * 0.5);
//   vec3 oc = ray.origin - center;
//   float b = dot(oc, ray.direction);
//   float c = dot(oc, oc) - radius * radius;
//   float discriminant = b * b - c;
//   return (discriminant > 0.0);
// }

Intersection trace_ray(Ray ray) {
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

  uint stack[STACK_MAX];

  uint sp = 0u;

  stack[sp++] = 0u;     // root = 0

  while (sp > 0u) {
    uint node = stack[--sp];
    int item = int(get(bvh, node, 8u));

    // leaf
    if (item != -1) {
      ray_ball_intersection(ray, item, result);
      continue;
    }

    uint left  = uint(get(bvh, node, 0u));
    uint right = uint(get(bvh, node, 1u));

    float x0 = get(bvh, node, 2u);
    float y0 = get(bvh, node, 3u);
    float z0 = get(bvh, node, 4u);

    float x1 = get(bvh, node, 5u);
    float y1 = get(bvh, node, 6u);
    float z1 = get(bvh, node, 7u);

    // internal → test bbox
    if (!ray_bbox_intersection(ray, x0, y0, z0, x1, y1, z1))continue;

    // push children on the stack
    stack[sp++] = left;
    stack[sp++] = right;

    if (sp >= STACK_MAX)
      break; // avoid overflow
  }

  // for (int i = 0; i < 4; i++) {
  //   ray_ball_intersection(ray, i, result);
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