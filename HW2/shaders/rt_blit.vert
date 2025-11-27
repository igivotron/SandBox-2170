#version 450 core

layout(location=0) out vec2 uv;

void main() {
  if (gl_VertexIndex == 0) {
    gl_Position = vec4(-1, -1, 0, 1);
    uv = vec2(0.0);
  }
  else if (gl_VertexIndex == 1) {
    gl_Position = vec4(3, -1, 0, 1);
    uv = vec2(2.0, 0.0);
  }
  else if (gl_VertexIndex == 2) {
    gl_Position = vec4(-1, 3, 0, 1);
    uv = vec2(0.0, 2.0);
  }
}
