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
  for (uint i = 0u; i < materials.length(); i++) {
    ray_ball_intersection(ray, int(i), result);
  }

  return result;
}

/*
 * The code below is an example of how you may access your data structure (see commented code in
 * homework.py for an example of how to transfer your data from the CPU to the GPU).
 */

/* layout(set = 2, binding = 0, std430) buffer DataStructure { */
/*   vec2 example_data[]; */
/* }; */
