#version 450 core

layout(location = 0) in vec3 world_normal;
layout(location = 1) in vec3 view_normal;
layout(location = 2) in vec3 world_eye;
layout(location = 3) in vec3 world_pos;
layout(location = 4) in vec3 base_color;
layout(location = 5) in float roughness;
layout(location = 6) in float metalness;

layout(location = 0) out vec4 frag_color;

#define M_PI 3.14159265358979323846

float ggx(vec3 normal, vec3 x, float roughness) {
  float r2 = roughness * roughness;
  float r4 = r2*r2;

  float k = dot(normal, x);
  float k2 = k*k;

  float d = k2*(r4 - 1.0) + 1.0;
  return r4/(M_PI * d*d);
}

vec3 fresnel_schlick(vec3 normal, vec3 view, vec3 r0) {
  return r0 + (vec3(1.0) - r0)*pow(1.0 - clamp(dot(normal, view), 0.0, 1.0), 5.0);
}

float g_cook_torrance(vec3 normal, vec3 half_vector,
                      vec3 incident_vector, vec3 eye_vector) {
  float dot_vh = clamp(dot(eye_vector, half_vector), 1e-6, 1.0);
  float dot_nh = clamp(dot(normal, half_vector), 0.0, 1.0);
  float dot_nv = clamp(dot(normal, eye_vector), 0.0, 1.0);
  float dot_nl = clamp(dot(normal, incident_vector), 0.0, 1.0);

  float x = 2.0*dot_nh/dot_vh;
  return min(1.0, min(x*dot_nl, x*dot_nv));
}

vec3 cook_torrance(vec3 normal,
                    vec3 incident_vector, vec3 eye_vector,
                    float roughness, vec3 f0) {
  eye_vector = normalize(eye_vector);

  float dot_nl = clamp(dot(normal, incident_vector), 0.0, 1.0);
  float dot_nv = clamp(dot(normal, eye_vector), 0.0, 1.0);

  vec3 half_vector = normalize(eye_vector + incident_vector);

  if (dot_nl > 0.0 && dot_nv > 0.0) {
    float D = ggx(normal, half_vector, roughness);
    vec3 F = fresnel_schlick(normal, eye_vector, f0);
    float G = g_cook_torrance(normal, half_vector, incident_vector, eye_vector);

    return (D*F*G)/(4.0*dot_nv*dot_nl);
  }
  else
    return vec3(0.0);
}

void main() {
  const vec3 light_pos = vec3(0.0, 9.5, 0.0);

  vec3 to_light = light_pos - world_pos;
  float light_distance2 = dot(to_light, to_light);
  vec3 light_dir = normalize(to_light);
  float light_intensity = 3200.0 / (4 * M_PI * light_distance2);

  vec3 diffuse = base_color * light_intensity *
    clamp((1 - metalness) * dot(light_dir, normalize(world_normal)), 0.2, 1.0);
  vec3 specular = vec3(light_intensity) *
    cook_torrance(normalize(world_normal), light_dir, world_eye, roughness,
                  mix(vec3(0.04), base_color, metalness));

  frag_color = vec4(diffuse + specular, 1.0);
}
