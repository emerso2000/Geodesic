#version 450 core
layout (location = 0) in vec2 fragmentTexCoords;
uniform sampler2D basicTexture;
out vec4 finalColor;
void main()
{
    finalColor = texture(basicTexture, fragmentTexCoords);
}