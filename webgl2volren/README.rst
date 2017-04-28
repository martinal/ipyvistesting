
TODO:

- Organize as a package: npm package? widget?

- Test drawElementsInstanced with shaders:

out vec4 color;

void main()
{
    vec3 mesh[8] = {
         1.0,  1.0,  1.0,
         1.0,  0.0,  0.0,
         0.0,  1.0,  0.0,
         0.0,  0.0,  1.0,
        -1.0, -1.0, -1.0,
        -1.0,  0.0,  0.0,
         0.0, -1.0,  0.0,
         0.0,  0.0, -1.0
    };
    vec3 colors[2] = {
        1.0, 0.0. 0.0,
        0.0, 1.0. 0.0,
    }
    color = vec4(colors[gl_instanceId], 1.0);
    gl_Position = vec4(mesh[gl_VertexId], 1.0);
}


in vec4 color;
out vec4 fragColor;

void main()
{
    fragColor = color;
}

