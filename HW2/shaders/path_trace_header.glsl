#define WORK_GROUP_SIZE 32
#define NUM_BOUNCES 4
#define BRIGHTNESS 6.0
#define ROUGHNESS 0.5
/* #define SQUARE_LIGHT */
#define LIGHT_RADIUS 9.0

#define MAX_DISTANCE 1000.0

#define M_PI 3.14159265358979323846

layout(local_size_x = WORK_GROUP_SIZE, local_size_y = 1, local_size_z = 1) in;

struct CameraData {
  mat4 inv_view;
  mat4 inv_projection;
};

layout(set=0, binding=0, std140) uniform Camera {
  CameraData camera;
};

layout(set=1, binding=0, std430) buffer Sobol {
  uint SOBOL_DIRECTIONS[];
};

layout(set=1, binding=1, std430) buffer Positions {
  /*
   * WGPU doesn't support the scalar block layout, so use an array of float to
   * reconstruct objects dynamically.
   */
  float positions[];
};

layout(set=1, binding=2, std430) buffer Orientations {
  vec4 orientations[];
};

layout(set=1, binding=3, std430) buffer Radii {
  float radii[];
};

struct Material {
  vec3 albedo;
  float roughness;
  vec3 emittance;
  float metalness;
};

layout(set=1, binding=4, std430) buffer Materials {
  Material materials[];
};

layout(set=1, binding=5, r32f) uniform image2D result;

layout(set=1, binding=6) uniform ResetBlock {
  int reset;
};

struct Ball {
  vec3 pos;
  float radius;
  vec4 rotation;
  Material material;
};

struct Plane {
  vec3 point;
  vec3 normal;
  Material material;
};

struct Ray {
  vec3 origin;
  vec3 direction;
};

struct Intersection {
  float distance;
  vec3 pos;
  vec3 normal;

  Material material;
};

Ball get_ball(int i) {
  return Ball(vec3(positions[3*i], positions[3*i + 1], positions[3*i + 2]), radii[i],
              orientations[i], materials[i]);
}

vec4 quat_multiply(vec4 a, vec4 b) {
  vec3 v = b.w*a.xyz + a.w*b.xyz + cross(a.xyz, b.xyz);
  float w = a.w*b.w - dot(a.xyz, b.xyz);
  return vec4(v, w);
}

vec4 quat_conjugate(vec4 q) {
  return vec4(-q.xyz, q.w);
}

vec3 quat_rotate(vec4 q, vec3 p) {
  vec4 q_conj = quat_conjugate(q);
  vec4 q_p = vec4(p, 0.0);
  return quat_multiply(q, quat_multiply(q_p, q_conj)).xyz;
}

float sobol(uint dim, uint i) {
  uint result = 0u;

  uint j = 0;
  for (; i != 0u; i /= 2) {
    if ((i & 1u) != 0)
      result ^= SOBOL_DIRECTIONS[32*dim + j];
    j++;
  }

  return float(result) * (1.0/float(0xFFFFFFFFu));
}

float rand_value(vec2 p){
  return fract(sin(dot(p, vec2(12.9898,78.233))) * 43758.5453);
}

float sobol_rotation;

void setup_rotation(uint sample_id) {
  vec2 p = vec2(32*float(gl_WorkGroupID.x) + float(gl_LocalInvocationIndex),
                float(gl_WorkGroupID.y));
  sobol_rotation = rand_value(p)*2*M_PI;
  uint i = sample_id * (8*imageSize(result).x*gl_WorkGroupID.y +
                        32*gl_WorkGroupID.x + gl_LocalInvocationIndex);
  sobol_rotation = 2*M_PI*sobol(512, i);
}

#ifdef HAS_BARRIER
shared uint sobol_current_dir[32*3];

void sobol_setup(uint bounce_id) {
  for (uint i = 0; i < 3; i++) {
    for (uint j = 0; j < 32; j += WORK_GROUP_SIZE) {
      uint offset = j + gl_LocalInvocationIndex;
      if (offset < 32)
        sobol_current_dir[32*i + offset] = SOBOL_DIRECTIONS[32*(5*bounce_id + i) + offset];
    }
  }

  /* No memoryBarrierShared instruction in WGPU :( */
  barrier();
}
#else
uint sobol_current_dir[32*3];

void sobol_setup(uint bounce_id) {
  for (uint i = 0; i < 3; i++) {
    for (uint j = 0; j < 32; j++) {
      sobol_current_dir[32 * i + j] = SOBOL_DIRECTIONS[32*(5*bounce_id + i) + j];
    }
  }
}
#endif

float sobol_next(uint dim, uint i) {
  uint result = 0u;
  uint j = 0u;
  for (; i != 0u; i /= 2u) {
    if ((i & 1u) != 0)
      result ^= sobol_current_dir[32*dim + j];
    j++;
  }

  return float(result) * (1.0/float(0xFFFFFFFFu));
}

#define EMPTY_INTERSECTION Intersection(MAX_DISTANCE, vec3(0.0), vec3(0.0), Material(vec3(0.0), 0.0, vec3(0.0), 0.0))

void merge_intersection(inout Intersection a, Intersection b) {
  if (a.distance > b.distance) a = b;
}

bool ray_ball_intersection(Ray ray, int i,
                           inout Intersection intersection) {
  Ball ball = get_ball(i);

  vec3 to_travel = ball.pos - ray.origin;
  float to_pos = dot(ray.direction, to_travel);
  float discr = to_pos*to_pos - dot(to_travel, to_travel) + ball.radius*ball.radius;

  if (discr > 0.0) {
    float t = to_pos - sqrt(discr);
    if (t > 0.0 && t < intersection.distance) {
      intersection.distance = t;
      intersection.pos = ray.origin + t*ray.direction;
      intersection.normal = normalize(intersection.pos - ball.pos);

      vec3 model_pos = quat_rotate(quat_conjugate(ball.rotation),
                                   intersection.pos - ball.pos);
      bool is_rim = abs(model_pos.y) < 0.1 * ball.radius;

      intersection.material = ball.material;
      if (is_rim) intersection.material.albedo = vec3(0.2);

      return true;
    }
  }

  return false;
}

bool ray_plane_intersection(Ray ray, Plane plane,
                            out Intersection intersection) {
  vec3 offset = plane.point - ray.origin;

  float t = dot(offset, plane.normal) / dot(ray.direction, plane.normal);
  if (dot(ray.direction, plane.normal) > 0.0) return false;

  vec3 p = ray.origin + t*ray.direction;

  if (t > 0.0 && max(abs(p.x), max(abs(p.y), abs(p.z))) <= 10.0 + 1e-3) {
    intersection.distance = t;
    intersection.pos = ray.origin + t*ray.direction;
    intersection.normal = plane.normal;
    intersection.material = plane.material;

#ifdef SQUARE_LIGHT
    float r = max(abs(p.x), abs(p.z));
    float r2 = r*r;
#else
    float r2 = dot(p.xz, p.xz);
#endif

    if (r2 > LIGHT_RADIUS*LIGHT_RADIUS)
      intersection.material.emittance = vec3(0.0);
    return true;
  }

  return false;
}
