#include <pok.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define WIDTH 640
#define HEIGHT 640

int pauseFlag = 0;
float spinAnimation = 0;

void OnKeyDown(int keyCode, void* userData) {
    // Hit space to pause/play animation.
    if (keyCode == ' ') {
        pauseFlag = !pauseFlag;
    }
}

Pok_Byte3 PixelShader(Pok_Int2 pixel_pos, const Pok_AttributeArray* interpolated, void* userData) {
    Pok_Float3 normal = POK_ATTRIB_TO_FLOAT3(interpolated->attributes[0]);
    
    /*Pok_Byte3 obj_color = {50, 255, 150};

    Pok_Float3 lightDir = *(Pok_Float3*)userData;
    
    float light = Pok_dot(normal, lightDir);
    if (light < 0) {
        light = 0;
    }
    
    Pok_Byte3 color = {obj_color.r * light, obj_color.g * light, obj_color.b * light};*/
    
    Pok_Byte3 color = {
        Pok_Clamp((normal.x + 1) / 2 * 255, 0, 255), 
        Pok_Clamp((normal.y + 1) / 2 * 255, 0, 255), 
        Pok_Clamp((normal.z + 1) / 2 * 255, 0, 255)
    };

    return color;
}

int main() {
    Pok_Mesh* mesh = Pok_MeshLoadOBJ("assets/bullfrog.obj");
    assert(mesh);
    Pok_Image* fontImage = Pok_ImageLoadBMP("assets/font_dos_vga.bmp", POK_R8G8B8A8);
    assert(fontImage);
    Pok_Font* font = Pok_FontLoad(fontImage, "assets/font_dos_vga.txt");
    assert(font);

    Pok_ColorBuffer* colorBuffer = Pok_ColorBufferCreate(WIDTH, HEIGHT);
    Pok_DepthBuffer* depthBuffer = Pok_DepthBufferCreate(WIDTH, HEIGHT);

    Pok_Init();
    Pok_Window* window = Pok_WindowCreate("Spinning Model", WIDTH, HEIGHT, NULL);

    if (!window) {
        printf("Window error.\n");
        return -1;
    }
    
    Pok_WindowSetKeyDownCallback(window, OnKeyDown);

    Pok_Float3 cameraPos = {15, 15, -25};
    Pok_Float3 cameraTarget = {0, 0, 0};
    Pok_Mat4 view = Pok_LookAt(cameraPos, cameraTarget, (Pok_Float3){0, -1, 0});
    Pok_Mat4 projection = Pok_Perspective(HEIGHT / (float)WIDTH, 75.f * (M_PI / 180), 0.1f, 1000.f);

    double lastFrame = Pok_GetElapsedTimeMS();

    const int guiDeltaInterval = 500;
    double guiDeltaLast = lastFrame;
    double guiDeltaValue = 0;
    float guiRenderTime = 0;

    while (!Pok_WindowShouldClose(window)) {
        double currentFrame = Pok_GetElapsedTimeMS();

        Pok_WindowPollEvents(window);

        float delta = currentFrame - lastFrame;
        lastFrame = currentFrame;

        if (!pauseFlag) {
            spinAnimation += delta * M_PI / 1000;
        }

        Pok_DepthBufferClear(depthBuffer);

        // Gradient background.
        /*for (int y = 0; y < HEIGHT; y++) {
            static Pok_Byte3 top = {141, 160, 184};
            static Pok_Byte3 bottom = {20, 20, 20};

            for (int x = 0; x < WIDTH; x++) {
                *Pok_colorBuffer_at(colorBuffer, x, y) = Pok_byte3_lerp(top, bottom, (float)y / HEIGHT);
            }
        }*/
        Pok_ColorBufferClear(colorBuffer, (Pok_Byte3){255, 255, 255});

        Pok_Mat4 model;
        Pok_Mat4InitIdentity(&model);
        
        model = Pok_Mat4Translate(model, (Pok_Float3){-.5, -.5, -.5});
        model = Pok_Mat4RotateX(model, 90.f * (M_PI / 180));
        model = Pok_Mat4RotateY(model, spinAnimation);
        
        Pok_Mat4 mv = Pok_Mat4Mul(model, view);
        Pok_Mat4 mvp = Pok_Mat4Mul(Pok_Mat4Mul(model, view), projection);
        Pok_Mat3 normalMatrix = Pok_Mat4ToMat3(Pok_Mat4Transpose(Pok_Mat4Inverse(model)));

        // Make a copy of the mesh's positions to preserve the original mesh.
        Pok_Float4* positions = Pok_AllocHeap(mesh->positionCount * sizeof(Pok_Float4));
        for (int i = 0; i < mesh->positionCount; i++) {
            positions[i].x = mesh->positions[i].x;
            positions[i].y = mesh->positions[i].y;
            positions[i].z = mesh->positions[i].z;
            positions[i].w = 1;
        }

        int* faceCullFlags = Pok_AllocHeap(mesh->faceCount * sizeof(int));

        // Backface culling.
        for (int i = 0; i < mesh->faceCount; i++) {
            // Pick any one of the triangle's points, put it in model space and check if the triangle should be culled.
            
            Pok_Float4 anyPointOnTriangle = Pok_Float4MulMat4(positions[mesh->faces[i].positionIndices[0]], model);
            
            Pok_Float3 dirTowardsPoint;
            for (int j = 0; j < 3; j++) {
                dirTowardsPoint.at[j] = anyPointOnTriangle.at[j] - cameraPos.at[j];
            }
            dirTowardsPoint = Pok_Normalize(dirTowardsPoint);
            
            Pok_Float3 normal = Pok_Float3MulMat3(mesh->normals[mesh->faces[i].normalIndices[0]], normalMatrix);

            if (Pok_Dot(dirTowardsPoint, normal) >= 0) {
                faceCullFlags[i] = 1;
            }
            else {
                faceCullFlags[i] = 0;
            }
        }

        // Local space to homogeneous clip space.
        for (int i = 0; i < mesh->positionCount; i++) {
            positions[i] = Pok_Float4MulMat4(positions[i], mvp);
        }

        double renderTimeStart = Pok_GetElapsedTimeMS();

        // Render triangles.
        for (int i = 0; i < mesh->faceCount; i++) {
            if (faceCullFlags[i]) {
                continue;
            }

            Pok_Float4 tri[3];
            for (int j = 0; j < 3; j++) {
                tri[j] = positions[mesh->faces[i].positionIndices[j]];
            }

            Pok_Float3 lightDir = {1, 0, 0};
            
            Pok_Float3 n0 = Pok_Float3MulMat3(mesh->normals[mesh->faces[i].normalIndices[0]], normalMatrix);
            Pok_Float3 n1 = Pok_Float3MulMat3(mesh->normals[mesh->faces[i].normalIndices[1]], normalMatrix);
            Pok_Float3 n2 = Pok_Float3MulMat3(mesh->normals[mesh->faces[i].normalIndices[2]], normalMatrix);

            Pok_RenderTriangle3D(colorBuffer, 
                                 depthBuffer, 
                                 tri[0], tri[1], tri[2],
                                 POK_ATTRIB_ARRAY(POK_ATTRIB_3(n0)),
                                 POK_ATTRIB_ARRAY(POK_ATTRIB_3(n1)),
                                 POK_ATTRIB_ARRAY(POK_ATTRIB_3(n2)),
                                 1,
                                 PixelShader,
                                 &lightDir);
        }

        double renderTimeEnd = Pok_GetElapsedTimeMS();

        Pok_Free(faceCullFlags);
        Pok_Free(positions);

        if (currentFrame > (guiDeltaLast + guiDeltaInterval)) {
            guiDeltaLast = currentFrame;
            guiDeltaValue = delta;
            guiRenderTime = renderTimeEnd - renderTimeStart;
        }

        char buf[24];
        snprintf(buf, 24, "Frame:  %.2fms", guiDeltaValue);
        Pok_RenderText(colorBuffer, buf, (Pok_Int2){10, 10}, font, 1);

        snprintf(buf, 24, "Raster: %.2fms", guiRenderTime);
        Pok_RenderText(colorBuffer, buf, (Pok_Int2){10, 24}, font, 1);

        Pok_WindowSwapBuffers(window, colorBuffer);
    }

    Pok_WindowDestroy(window);
    printf("Done.");
}
