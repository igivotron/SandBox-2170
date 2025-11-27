#version 450 core
layout(location=0) in vec2 uv;
layout(set=0, binding=0, r32f) uniform readonly image2D rt_output;

layout(location=0) out vec4 color;

void main() {
  uvec2 size = imageSize(rt_output) / uvec2(4, 1);
  ivec2 coords = ivec2(vec2(size) * vec2(uv.x, 1.0 - uv.y));

  float r = imageLoad(rt_output, ivec2(4, 1) * coords + ivec2(0, 0)).r;
  float g = imageLoad(rt_output, ivec2(4, 1) * coords + ivec2(1, 0)).r;
  float b = imageLoad(rt_output, ivec2(4, 1) * coords + ivec2(2, 0)).r;
  float n = imageLoad(rt_output, ivec2(4, 1) * coords + ivec2(3, 0)).r;

  color = vec4(vec3(r, g, b) / n, 1.0);
}
