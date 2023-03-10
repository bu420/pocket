#include "pwa.h"

#include <windows.h>
#include <stdlib.h>

#define PWA_WINDOW_CLASS_NAME "PWA Window Class"
#define PWA_WINDOW_PROP_NAME "PWA Prop"

struct pwa_window {
    HWND hwnd;
    int should_close;
    pwa_resize_callback on_resize;
    pwa_draw_callback on_draw;
    void* user_data;
};

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam);

void pwa_init() {
    WNDCLASS wc = {0};
    wc.lpfnWndProc = WindowProc;
    wc.hInstance = GetModuleHandle(0);
    wc.lpszClassName = PWA_WINDOW_CLASS_NAME;
    RegisterClass(&wc);
}

pwa_window_t* pwa_window_create(char* title, int w, int h, void* user_data) {
    HWND hwnd = CreateWindowEx(0, PWA_WINDOW_CLASS_NAME, title, WS_OVERLAPPEDWINDOW, 
        CW_USEDEFAULT, CW_USEDEFAULT, w, h, NULL, NULL, GetModuleHandle(0), NULL);

    if (!hwnd) {
        return NULL;
    }

    pwa_window_t* window = malloc(sizeof(pwa_window_t));
    window->hwnd = hwnd;
    window->should_close = 0;
    window->user_data = user_data;
    SetProp(hwnd, PWA_WINDOW_PROP_NAME, window);

    ShowWindow(hwnd, SW_SHOW);

    return window;
}

void pwa_window_destroy(pwa_window_t* window) {
    free(window);
    window = NULL;
}

void pwa_window_set_resize_callback(pwa_window_t* window, pwa_resize_callback on_resize) {
    window->on_resize = on_resize;
}

void pwa_window_set_draw_callback(pwa_window_t* window, pwa_draw_callback on_draw) {
    window->on_draw = on_draw;
}

int pwa_window_should_close(pwa_window_t* window) {
    return window->should_close;
}

void pwa_window_poll_events(pwa_window_t* window) {
    MSG msg;
    while (PeekMessage(&msg, window->hwnd, 0, 0, PM_REMOVE)) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
}

void _pwa_window_close_request(pwa_window_t* window) {
    if (window) {
        window->should_close = 1;
    }
}

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
    pwa_window_t* window = GetProp(hwnd, PWA_WINDOW_PROP_NAME);
    
    switch (uMsg) {
    case WM_CREATE:
        SetTimer(hwnd, 1, 16, NULL);
        return 0;

    case WM_TIMER:
        InvalidateRect(hwnd, NULL, FALSE);
        return 0;

    case WM_CLOSE:
        _pwa_window_close_request(window);
        return 0;

    case WM_PAINT:
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hwnd, &ps);
        
        if (window->on_draw) {
            pwa_pixel_buffer_t buffer = window->on_draw(window->user_data);

            BITMAPINFO bmi = {0};
            bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
            bmi.bmiHeader.biWidth = buffer.w;
            bmi.bmiHeader.biHeight = -buffer.h;
            bmi.bmiHeader.biPlanes = 1;
            bmi.bmiHeader.biBitCount = 32;
            bmi.bmiHeader.biCompression = BI_RGB;

            StretchDIBits(hdc, 0, 0, buffer.w, buffer.h, 0, 0, buffer.w, buffer.h, buffer.pixels, &bmi, DIB_RGB_COLORS, SRCCOPY);
        }
        else {
            FillRect(hdc, &ps.rcPaint, (HBRUSH)(COLOR_WINDOW + 1));
        }

        EndPaint(hwnd, &ps);
        return 0;

    case WM_SIZE:
        if (window->on_resize) {
            window->on_resize(LOWORD(lParam), HIWORD(lParam), window->user_data);
        }
        return 0;
    }

    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}
