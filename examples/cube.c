#include <pok.h>
#include <stdlib.h>
#include <stdio.h>

#define WIDTH 640
#define HEIGHT 640

#define ATLAS_W 256.f
#define ATLAS_H 256.f

const Pok_Float3 cubePositions[36] = {
    {-0.5f,-0.5f,-0.5f}, {-0.5f,0.5f,-0.5f}, {0.5f,0.5f,-0.5f},   {-0.5f,-0.5f,-0.5f}, {0.5f,0.5f,-0.5f},   {0.5f,-0.5f,-0.5f},
    {0.5f,-0.5f,-0.5f},  {0.5f,0.5f,-0.5f},  {0.5f,0.5f,0.5f},    {0.5f,-0.5f,-0.5f},  {0.5f,0.5f,0.5f},    {0.5f,-0.5f,0.5f},
    {0.5f,-0.5f,0.5f},   {0.5f,0.5f,0.5f},   {-0.5f,0.5f,0.5f},   {0.5f,-0.5f,0.5f},   {-0.5f,0.5f,0.5f},   {-0.5f,-0.5f,0.5f},
    {-0.5f,-0.5f,0.5f},  {-0.5f,0.5f,0.5f},  {-0.5f,0.5f,-0.5f},  {-0.5f,-0.5f,0.5f},  {-0.5f,0.5f,-0.5f},  {-0.5f,-0.5f,-0.5f},
    {-0.5f,0.5f,-0.5f},  {-0.5f,0.5f,0.5f},  {0.5f,0.5f,0.5f},    {-0.5f,0.5f,-0.5f},  {0.5f,0.5f,0.5f},    {0.5f,0.5f,-0.5f},
    {0.5f,-0.5f,0.5f},   {-0.5f,-0.5f,0.5f}, {-0.5f,-0.5f,-0.5f}, {0.5f,-0.5f,0.5f},   {-0.5f,-0.5f,-0.5f}, {0.5f,-0.5f,-0.5f}
};

const Pok_Float2 cubeTexCoords[36] = {
    // Sides.
    {48/ATLAS_W,15/ATLAS_H}, {48/ATLAS_W,0}, {63/ATLAS_W,0}, {48/ATLAS_W,15/ATLAS_H}, {63/ATLAS_W,0}, {63/ATLAS_W,15/ATLAS_H},
    {48/ATLAS_W,15/ATLAS_H}, {48/ATLAS_W,0}, {63/ATLAS_W,0}, {48/ATLAS_W,15/ATLAS_H}, {63/ATLAS_W,0}, {63/ATLAS_W,15/ATLAS_H},
    {48/ATLAS_W,15/ATLAS_H}, {48/ATLAS_W,0}, {63/ATLAS_W,0}, {48/ATLAS_W,15/ATLAS_H}, {63/ATLAS_W,0}, {63/ATLAS_W,15/ATLAS_H},
    {48/ATLAS_W,15/ATLAS_H}, {48/ATLAS_W,0}, {63/ATLAS_W,0}, {48/ATLAS_W,15/ATLAS_H}, {63/ATLAS_W,0}, {63/ATLAS_W,15/ATLAS_H},
    // Grass top.
    {0,15/ATLAS_H}, {0,0}, {15/ATLAS_W,0}, {0,15/ATLAS_H}, {15/ATLAS_W,0}, {15/ATLAS_W,15/ATLAS_H},
    // Dirt bottom.
    {32/ATLAS_W,15/ATLAS_H}, {32/ATLAS_W,0}, {47/ATLAS_W,0}, {32/ATLAS_W,15/ATLAS_H}, {47/ATLAS_W,0}, {47/ATLAS_W,15/ATLAS_H}
};

const Pok_Float3 cubeNormals[6] = {
    { 0,  0,  1},
    {-1,  0,  0},
    { 0,  0, -1},
    { 1,  0,  0},
    { 0,  1,  0},
    { 0, -1,  0}
};

typedef struct {
    Pok_Image* textureAtlas;
    Pok_Float3 normal;
} ShaderData;

Pok_Byte3 PixelShader(Pok_Int2 pixelPos, const Pok_AttributeArray* interpolated, void* userData) {
    ShaderData data = *(ShaderData*)userData;
    Pok_Float2 texCoord = POK_ATTRIB_TO_FLOAT2(interpolated->attributes[0]);

    Pok_Byte* samplAaddress = Pok_ImageSample(data.textureAtlas, texCoord.u, texCoord.v);
    Pok_Byte3 color = {*samplAaddress, *(samplAaddress + 1), *(samplAaddress + 2)};

    return color;
}

int main() {
    Pok_Init();

    Pok_Image* textureAtlas = Pok_ImageLoadBMP("assets/terrain.bmp", POK_R8G8B8A8);
    assert(textureAtlas);

    Pok_ColorBuffer* colorBuffer = Pok_ColorBufferCreate(WIDTH, HEIGHT);
    Pok_DepthBuffer* depthBuffer = Pok_DepthBufferCreate(WIDTH, HEIGHT);

    Pok_Float3 cameraPos = {2, 2, 2};
    Pok_Mat4 projection = Pok_Perspective(HEIGHT / (float)WIDTH, 70 * (M_PI / 180), .1f, 1000.f);

    Pok_Window* window = Pok_WindowCreate("Cube", WIDTH, HEIGHT, NULL);

    while (!Pok_WindowShouldClose(window)) {
        Pok_WindowPollEvents(window);

        Pok_Mat4 view = Pok_LookAt(cameraPos, (Pok_Float3){0, 0, 0}, (Pok_Float3){0, -1, 0});

        Pok_Mat4 model;
        Pok_Mat4InitIdentity(&model);
        model = Pok_Mat4RotateY(model, Pok_GetElapsedTimeMS() * M_PI / 2000);
        model = Pok_Mat4RotateZ(model, Pok_GetElapsedTimeMS() * M_PI / 4000);

        Pok_Mat4 mvp = Pok_Mat4Mul(Pok_Mat4Mul(model, view), projection);

        // Clear buffers.

        Pok_ColorBufferClear(colorBuffer, (Pok_Byte3){0, 0, 0});
        Pok_DepthBufferClear(depthBuffer);

        // Raster cube (12 triangles).
        for (int i = 0; i < 12; i++) {
            Pok_Float4 pos[3];
            Pok_Float2 texCoords[3];
            
            for (int j = 0; j < 3; j++) {
                Pok_Float3 p3 = cubePositions[i * 3 + j];
                pos[j] = (Pok_Float4){p3.x, p3.y, p3.z, 1};
                texCoords[j] = cubeTexCoords[i * 3 + j];

                // Multiply position with MVP matrix, world space -> clip space.
                pos[j] = Pok_Float4MulMat4(pos[j], mvp);
            }

            Pok_Float3 normal = cubeNormals[i / 2];

            ShaderData data;
            data.textureAtlas = textureAtlas;
            data.normal = normal;

            Pok_RenderTriangle3D(colorBuffer, 
                                 depthBuffer, 
                                 pos[0], 
                                 pos[1], 
                                 pos[2], 
                                 POK_ATTRIB_ARRAY(POK_ATTRIB_2(texCoords[0])),
                                 POK_ATTRIB_ARRAY(POK_ATTRIB_2(texCoords[1])),
                                 POK_ATTRIB_ARRAY(POK_ATTRIB_2(texCoords[2])),
                                 1,
                                 PixelShader,
                                 &data);

            Pok_WindowSwapBuffers(window, colorBuffer);
        }
    }
}
