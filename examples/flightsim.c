#include <pok.h>
#include <stdio.h>

#define WIDTH 600
#define HEIGHT 600

#define ARROW_DOWN 40
#define ARROW_RIGHT 39
#define ARROW_UP 38
#define ARROW_LEFT 37

typedef struct {
    Pok_Float3 pos;
    Pok_Float3 dir;
    Pok_Float3 up;

    float pitch;
    float yaw;
} Camera;

typedef struct {
    Camera camera;
} StaticShaderData;

void OnKeyDown(int keyCode, void* userData) {
    ((bool*)userData)[keyCode] = true;
}

void OnKeyUp(int keyCode, void* userData) {
    ((bool*)userData)[keyCode] = false;
}

void CameraUpdate(Camera* camera) {
    camera->dir.x = cosf(camera->yaw) * cosf(camera->pitch);
    camera->dir.y = sinf(camera->pitch);
    camera->dir.z = sinf(camera->yaw) * cosf(camera->pitch);
    camera->dir = Pok_Normalize(camera->dir);
}

void CameraInit(Camera* camera, Pok_Float3 pos, float pitch, float yaw) {
    camera->pos = pos;
    camera->pitch = pitch;
    camera->yaw = yaw;
    camera->up = (Pok_Float3){0, 1, 0};
    CameraUpdate(camera);
}

Pok_Byte3 SkyShader(Pok_Int2 pixelPos, const Pok_AttributeArray* attributes, void* userData) {
    StaticShaderData data = *(StaticShaderData*)userData;
    
    Pok_Float3 viewPointLocal = {pixelPos.x / (float)WIDTH * ((data.camera.dir.x + 1) / 2), 
                                 pixelPos.y / (float)HEIGHT * ((data.camera.dir.y + 1) / 2),
                                 ((data.camera.dir.z + 1) / 2)};

    return (Pok_Byte3){255 * viewPointLocal.x, 255 * viewPointLocal.y, 255 * viewPointLocal.z};
} 

Pok_Byte3 TriangleShader(Pok_Int2 pixelPos, const Pok_AttributeArray* attributes, void* userData) {
    return (Pok_Byte3){255, 0, 0};
}

int main() {
    bool keys[255] = {false};

    Pok_Init();
    Pok_Window* window = Pok_WindowCreate("Flight Sim", WIDTH, HEIGHT, keys);
    Pok_WindowSetKeyDownCallback(window, OnKeyDown);
    Pok_WindowSetKeyUpCallback(window, OnKeyUp);

    Pok_ColorBuffer* colorBuffer = Pok_ColorBufferCreate(WIDTH, HEIGHT);
    Pok_DepthBuffer* depthBuffer = Pok_DepthBufferCreate(WIDTH, HEIGHT);

    Camera camera;
    CameraInit(&camera, (Pok_Float3){0, 0, 10}, 0, Pok_Radians(-90));

    Pok_Mat4 projection = Pok_Perspective(HEIGHT / (float)WIDTH, Pok_Radians(70), 0.1f, 1000);

    double last = Pok_GetElapsedTimeMS();

    while (!Pok_WindowShouldClose(window)) {
        double current = Pok_GetElapsedTimeMS();
        float delta = current - last;
        last = current;

        Pok_WindowPollEvents(window);

        double elapsed = Pok_GetElapsedTimeMS();
        
        if (keys[ARROW_UP]) {
            camera.pitch -= delta / 500;
        }
        if (keys[ARROW_DOWN]) {
            camera.pitch += delta / 500;
        }
        if (keys[ARROW_LEFT]) {
            camera.yaw -= delta / 500;
        }
        if (keys[ARROW_RIGHT]) {
            camera.yaw += delta / 500;
        }

        CameraUpdate(&camera);

        Pok_Mat4 model;
        Pok_Mat4InitIdentity(&model);

        Pok_Mat4 view = Pok_LookAt(
            camera.pos, 
            (Pok_Float3){camera.pos.x + camera.dir.x, camera.pos.y + camera.dir.y, camera.pos.z + camera.dir.z}, 
            camera.up);

        Pok_Mat4 mvp = Pok_Mat4Mul(Pok_Mat4Mul(model, view), projection);

        Pok_ColorBufferClear(colorBuffer, (Pok_Byte3){255, 255, 255});
        Pok_DepthBufferClear(depthBuffer);

        // Render sky.

        Pok_Float3 skyVertices[4] = {
            (Pok_Float3){0, 0, 0}, 
            (Pok_Float3){WIDTH - 1, 0, 0}, 
            (Pok_Float3){WIDTH - 1, HEIGHT - 1, 0}, 
            (Pok_Float3){0, HEIGHT - 1, 0}
        };

        StaticShaderData data = {.camera = camera};

        Pok_RenderTriangle3D(colorBuffer, NULL, skyVertices[0], skyVertices[1], skyVertices[2], NULL, NULL, NULL, 0, SkyShader, &data);
        Pok_RenderTriangle3D(colorBuffer, NULL, skyVertices[2], skyVertices[3], skyVertices[0], NULL, NULL, NULL, 0, SkyShader, &data);

        // Render test triangle.

        Pok_Float3 tri[3] = {
            {0, 0, 0},
            {0, 1, 0},
            {1, 0, 0}
        };

        for (int i = 0; i < 3; i++) {
            tri[i] = Pok_Float3MulMat4(tri[i], mvp);

            // Scale from [-1, 1] to viewport size.
            tri[i].x = (tri[i].x + 1) / 2.f * colorBuffer->w;
            tri[i].y = (tri[i].y + 1) / 2.f * colorBuffer->h;
        }

        Pok_RenderTriangle3D(colorBuffer, depthBuffer, tri[0], tri[1], tri[2], NULL, NULL, NULL, 0, TriangleShader, NULL);

        Pok_WindowSwapBuffers(window, colorBuffer);
    }
}