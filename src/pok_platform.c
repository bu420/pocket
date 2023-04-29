#include "pok_platform.h"

#include <windows.h>
#include <stdlib.h>
#include <stdio.h>

#define PWA_WINDOW_CLASS_NAME "PWA Window Class"

struct Pok_Window {
    HWND hwnd;
    bool should_close;

    Pok_ResizeCallback onResize;
    Pok_KeyDownCallback onKeyDown;
    Pok_KeyUpCallback onKeyUp;

    int bufferWidth;
    int bufferHeight;
    uint32_t* buffer;

    void* userData;
};

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam);

void Pok_Init() {
    WNDCLASS wc = {0};
    wc.lpfnWndProc = WindowProc;
    wc.hInstance = GetModuleHandle(0);
    wc.lpszClassName = PWA_WINDOW_CLASS_NAME;
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.cbWndExtra = sizeof(Pok_Window*);
    RegisterClass(&wc);
}

Pok_Window* Pok_WindowCreate(char* title, int w, int h, void* userData) {
    HWND hwnd = CreateWindowEx(0, PWA_WINDOW_CLASS_NAME, title, WS_OVERLAPPEDWINDOW, 
        CW_USEDEFAULT, CW_USEDEFAULT, w, h, NULL, NULL, GetModuleHandle(0), NULL);

    if (!hwnd) {
        return NULL;
    }

    Pok_Window* window = Pok_AllocHeap(sizeof(Pok_Window));
    window->hwnd = hwnd;
    window->should_close = false;
    window->onResize = NULL;
    window->onKeyDown = NULL;
    window->onKeyUp = NULL;
    window->bufferWidth = 0;
    window->bufferHeight = 0;
    window->buffer = NULL;
    window->userData = userData;
    SetWindowLongPtr(hwnd, 0, (LONG_PTR)window);

    ShowWindow(hwnd, SW_SHOW);

    return window;
}

void Pok_WindowDestroy(Pok_Window* window) {
    if (window->buffer) {
        Pok_Free(window->buffer);
    }

    Pok_Free(window);
    window = NULL;
}

void Pok_WindowSetResizeCallback(Pok_Window* window, Pok_ResizeCallback onResize) {
    window->onResize = onResize;
}

void Pok_WindowSetKeyDownCallback(Pok_Window* window, Pok_KeyDownCallback onKeyDown) {
    window->onKeyDown = onKeyDown;
}

void Pok_WindowSetKeyUpCallback(Pok_Window* window, Pok_KeyUpCallback onKeyUp) {
    window->onKeyUp = onKeyUp;
}

bool Pok_WindowShouldClose(Pok_Window* window) {
    return window->should_close;
}

void Pok_WindowPollEvents(Pok_Window* window) {
    MSG msg;
    while (PeekMessage(&msg, window->hwnd, 0, 0, PM_REMOVE)) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
}

void Pok_WindowSwapBuffers(Pok_Window* window, Pok_ColorBuffer* colorBuffer) {
    bool resize = (colorBuffer->w != window->bufferWidth || colorBuffer->h != window->bufferHeight);
    
    if (resize) {
        Pok_Free(window->buffer);
    }

    if (!window->buffer || resize) {
        window->bufferWidth = colorBuffer->w;
        window->bufferHeight = colorBuffer->h;
        window->buffer = Pok_AllocHeap(colorBuffer->w * colorBuffer->h * sizeof(uint32_t));
    }

    for (int x = 0; x < colorBuffer->w; x++) {
        for (int y = 0; y < colorBuffer->h; y++) {
            Pok_Byte3 color = *Pok_ColorBufferAt(colorBuffer, x, y);
            window->buffer[y * colorBuffer->w + x] = (color.r << 16) | (color.g << 8) | (color.b);
        }
    }

    // Trigger redraw.
    InvalidateRect(window->hwnd, NULL, FALSE);
}

double Pok_GetElapsedTimeMS() {
    LARGE_INTEGER elapsedTime;
    QueryPerformanceCounter(&elapsedTime);
    return (double)(elapsedTime.QuadPart * 1000000 / Pok_GetTicksPerSecond()) / 1000;
}

int64_t Pok_GetTicksPerSecond() {
    LARGE_INTEGER ticksPerSecond;
    QueryPerformanceFrequency(&ticksPerSecond);
    return ticksPerSecond.QuadPart;
}

void Pok_DebugPrintLastError() {
    static int counter = 0;
    
    DWORD id = GetLastError();
    if (id == 0) {
        return;
    }

    LPSTR message = NULL;
    FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, NULL, id, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPSTR)&message, 0, NULL);
    printf("(%d) Win32 Error: %s\n", counter, message);
    LocalFree(message);

    counter++;
}

void Pok_Terminate() {
}

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
    Pok_Window* window = (Pok_Window*)GetWindowLongPtr(hwnd, 0);
    
    switch (uMsg) {
    case WM_CLOSE:
        if (window) {
            window->should_close = true;
        }
        return 0;

    case WM_PAINT:
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hwnd, &ps);
        
        if (window) {
            BITMAPINFO bmi = {0};
            bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
            bmi.bmiHeader.biWidth = window->bufferWidth;
            bmi.bmiHeader.biHeight = -window->bufferHeight;
            bmi.bmiHeader.biPlanes = 1;
            bmi.bmiHeader.biBitCount = 32;
            bmi.bmiHeader.biCompression = BI_RGB;

            StretchDIBits(hdc, 0, 0, window->bufferWidth, window->bufferHeight, 0, 0, window->bufferWidth, window->bufferHeight, window->buffer, &bmi, DIB_RGB_COLORS, SRCCOPY);
        }
        else {
            FillRect(hdc, &ps.rcPaint, (HBRUSH)(COLOR_WINDOW + 1));
        }

        EndPaint(hwnd, &ps);
        return 0;

    case WM_SIZE:
        if (window && window->onResize) {
            window->onResize(LOWORD(lParam), HIWORD(lParam), window->userData);
        }
        return 0;

    case WM_KEYDOWN:
        if (window && window->onKeyDown) {
            window->onKeyDown(wParam, window->userData);
        }
        return 0;

    case WM_KEYUP:
        if (window && window->onKeyUp) {
            window->onKeyUp(wParam, window->userData);
        }
        return 0;
    }

    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}
