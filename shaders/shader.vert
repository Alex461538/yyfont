#version 450

layout (location = 0) in vec2 position;
layout (location = 1) in vec2 UV;

out vec2 iUV;

uniform mat4 MVP;

void main()
{
    gl_Position = MVP * vec4( position, 1.0f, 1.0f);
    iUV = UV;
}
