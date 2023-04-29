#ifndef POK_PLATFORM_H
#define POK_PLATFORM_H

#if !(defined(_WIN32) || defined(_WIN64))
#error "Unsupported system. Only Windows is supported."
#endif

#include <pok_core.h>
#include <pok_gfx.h>

typedef struct Pok_Window Pok_Window;

typedef void (*Pok_ResizeCallback)(int w, int h, void* userData);
typedef void (*Pok_KeyDownCallback)(int keyCode, void* userData);
typedef void (*Pok_KeyUpCallback)(int keyCode, void* userData);

void Pok_Init();

Pok_Window* Pok_WindowCreate(char* title, int w, int h, void* userData);
void Pok_WindowDestroy(Pok_Window* window);
void Pok_WindowSetResizeCallback(Pok_Window* window, Pok_ResizeCallback onResize);
void Pok_WindowSetKeyDownCallback(Pok_Window* window, Pok_KeyDownCallback onKeyDown);
void Pok_WindowSetKeyUpCallback(Pok_Window* window, Pok_KeyUpCallback onKeyUp);
bool Pok_WindowShouldClose(Pok_Window* window);
void Pok_WindowPollEvents(Pok_Window* window);
void Pok_WindowSwapBuffers(Pok_Window* window, Pok_ColorBuffer* colorBuffer);

double Pok_GetElapsedTimeMS();
int64_t Pok_GetTicksPerSecond();

void Pok_DebugPrintLastError();

void Pok_Terminate();

#endif
