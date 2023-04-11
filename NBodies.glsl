#version 450 core

layout(std430, binding = 0) buffer pos {
    vec4 positions[]; // x, y, z, t
};

layout(std430, binding = 1) buffer vel {
    vec4 velocities[]; // vx, vy, vz, vt
};

layout(std430, binding = 2) buffer time {
    float dt;
};
layout(local_size_x = 10, local_size_y = 1, local_size_z = 1) in;

const float G = 1.0;
const float rs = 2.0;

void main() {
    // Get easier references to position and velocity
    vec4 position = positions[gl_GlobalInvocationID.x];
    vec4 velocity = velocities[gl_GlobalInvocationID.x];

    // Convert position and velocity to Schwarzschild coordinates
    float r = length(position.xyz);
    float theta = acos(position.z / r);
    float phi = atan(position.y, position.x);
    float vr = dot(position.xyz, velocity.xyz) / r;
    float vt = velocity.w * sin(theta) / r;

    // Calculate the Christoffel symbols
    float Gamma_rrr = -rs / (2.0 * r) * (r - rs);
    float Gamma_rtt = rs * (r - rs) / (2.0 * r * r * r);
    float Gamma_thetathetaphi = (rs - r) * sin(theta) * sin(theta);
    float Gamma_rthetatheta = rs - r;
    float Gamma_phirphi = 1 / r;
    float Gamma_tphi_phi = -sin(theta) * cos(theta);
    float Gamma_thetaphi = 1.0 / tan(theta);
    
    // Calculate the acceleration due to gravity
    float accel_r = -Gamma_rrr * vr * vr - 2.0 * Gamma_rtt * vt * vr;
    float accel_theta = -r * (Gamma_thetathetaphi * vt * vt + Gamma_thetaphi * vr * vt);
    float accel_phi = -r * sin(theta) * sin(theta) * Gamma_tphi_phi * vt * vt;

    // Update velocity and position
    vr += accel_r * position.w;
    vt += accel_theta * position.w;
    position.w += velocity.w * position.w / r * dt; 
    r += vr * position.w;
    theta += vt * position.w / r;
    phi += velocity.w * position.w / (r * r * sin(theta));

    // Convert updated position and velocity to Cartesian coordinates
    float x = r * sin(theta) * cos(phi);
    float y = r * sin(theta) * sin(phi);
    float z = r * cos(theta);
    float vx = vr * sin(theta) * cos(phi) + r * vt * cos(theta) * cos(phi) - velocity.w * sin(theta) * sin(phi) / r;
    float vy = vr * sin(theta) * sin(phi) + r * vt * cos(theta) * sin(phi) + velocity.w * sin(theta) * cos(phi) / r;
    float vz = vr * cos(theta) - r * vt * sin(theta);

    // Save updated position and velocity
    positions[gl_GlobalInvocationID.x] = vec4(x, y, z, position.w);
    velocities[gl_GlobalInvocationID.x] = vec4(vx, vy, vz, velocity.w);
}