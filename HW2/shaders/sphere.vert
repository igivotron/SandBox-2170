#version 450 core

struct CameraData {
  mat4 view;
  mat4 projection;
};

layout(set=0, binding=0, std140) uniform Camera {
  CameraData camera;
};

layout(location = 0) in vec3 world_position;
layout(location = 1) in vec3 vertex_normal;

/* Per instance attributes provided by the physics simulation */
layout(location = 2) in vec3 center;
layout(location = 3) in vec4 orientation;
layout(location = 4) in float radius;

/* Material */
layout(location = 5) in vec3 albedo;
layout(location = 6) in float roughness;
layout(location = 7) in vec3 emittance;
layout(location = 8) in float metalness;

layout(location = 0) out vec3 world_normal;
layout(location = 1) out vec3 view_normal;
layout(location = 2) out vec3 world_eye;
layout(location = 3) out vec3 frag_pos;
layout(location = 4) out vec3 color;
layout(location = 5) out float frag_roughness;
layout(location = 6) out float frag_metalness;

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

void main() {
  color = abs(world_position.y) < 0.1 ? vec3(0.0) : albedo;
  frag_roughness = roughness;
  frag_metalness = metalness;

  vec3 pos = center + radius * quat_rotate(orientation, world_position);

  world_normal = quat_rotate(orientation, vertex_normal);
  view_normal = mat3(camera.view) * world_normal;
  world_eye = -camera.view[3].xyz - pos;
  frag_pos = pos;

  gl_Position = camera.projection * (camera.view * vec4(pos, 1.0));
}
