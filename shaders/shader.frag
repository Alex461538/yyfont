#version 450

in vec2 iUV;

layout (binding = 0) uniform sampler2D glyphs;

out vec4 FragColor;

void main()
{
    vec4 color = texture(glyphs, iUV);
    FragColor = vec4(1.0, 1.0, 0.0, color.r);
}
