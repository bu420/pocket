#include <pok.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void OnResize(int w, int h, void* userData) {
    Pok_ColorBuffer* buf = (Pok_ColorBuffer*)userData;
    Pok_ColorBufferResize(buf, w, h);
}

void OnKeyDown(int keyCode, void* userData) {
    printf("Key down: %d\n", keyCode);
}

void OnKeyUp(int keyCode, void* userdData) {
    printf("Key up: %d\n", keyCode);
}

int main() {
    Pok_Init();

    Pok_ColorBuffer* colorBuffer = Pok_ColorBufferCreate(400, 400);

    Pok_Window* window = Pok_WindowCreate("Example Window", 400, 400, colorBuffer);

    if (!window) {
        printf("Window error.\n");
        return -1;
    }

    Pok_WindowSetResizeCallback(window, OnResize);
    Pok_WindowSetKeyDownCallback(window, OnKeyDown);
    Pok_WindowSetKeyUpCallback(window, OnKeyUp);

    while (!Pok_WindowShouldClose(window)) {
        Pok_WindowPollEvents(window);

        double time = Pok_GetElapsedTimeMS() / 500;

        for (int x = 0; x < colorBuffer->w; x++) {
            for (int y = 0; y < colorBuffer->h; y++) {
                Pok_Byte3* pixel = Pok_ColorBufferAt(colorBuffer, x, y);

                pixel->r = (sin(time) + 1) / 2 * 255;
                pixel->g = 255 - pixel->r;
                pixel->b = (cos(time) + 1) / 2 * 255;
            }
        }

        Pok_WindowSwapBuffers(window, colorBuffer);
    }

    Pok_WindowDestroy(window);
    printf("Done.\n");
}
