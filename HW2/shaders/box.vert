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

layout(location = 0) out vec3 world_normal;
layout(location = 1) out vec3 view_normal;
layout(location = 2) out vec3 world_eye;
layout(location = 3) out vec3 frag_pos;
layout(location = 4) out vec3 color;
layout(location = 5) out float roughness;
layout(location = 6) out float metalness;

void main() {
  if (vertex_normal.x == -1.0) color = vec3(0.0, 1.0, 0.0);
  if (vertex_normal.x == +1.0) color = vec3(0.0, 0.0, 1.0);
  if (vertex_normal.y == -1.0) color = vec3(1.0, 1.0, 1.0);
  if (vertex_normal.y == +1.0) color = vec3(1.0, 1.0, 1.0);
  if (vertex_normal.z == -1.0) color = vec3(1.0, 0.0, 0.0);
  if (vertex_normal.z == +1.0) color = vec3(1.0, 1.0, 1.0);

  if (vertex_normal.z == +1.0) {
    roughness = 0.05;
    metalness = 1.0;
  }
  else {
    roughness = 0.3;
    metalness = 0.0;
  }

  world_normal = vertex_normal;
  view_normal = mat3(camera.view) * view_normal;
  world_eye = -camera.view[3].xyz - world_position;
  frag_pos = world_position;

  gl_Position = camera.projection * (camera.view * vec4(world_position, 1.0));
}
