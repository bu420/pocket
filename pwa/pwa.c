#include "pwa/pwa.h"

#include <windows.h>
#include <stdlib.h>
#include <stdio.h>

#define PWA_WINDOW_CLASS_NAME "PWA Window Class"

struct pwa_window_t {
    HWND hwnd;
    bool should_close;

    pwa_resize_callback on_resize;
    pwa_key_down_callback on_key_down;
    pwa_key_up_callback on_key_up;

    int buffer_w;
    int buffer_h;
    uint32_t* buffer;

    void* user_data;
};

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam);

void pwa_init() {
    WNDCLASS wc = {0};
    wc.lpfnWndProc = WindowProc;
    wc.hInstance = GetModuleHandle(0);
    wc.lpszClassName = PWA_WINDOW_CLASS_NAME;
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.cbWndExtra = sizeof(pwa_window_t*);
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
    window->should_close = false;
    window->on_resize = NULL;
    window->on_key_down = NULL;
    window->on_key_up = NULL;
    window->buffer_w = 0;
    window->buffer_h = 0;
    window->buffer = NULL;
    window->user_data = user_data;
    SetWindowLongPtr(hwnd, 0, (LONG_PTR)window);

    pwa_print_last_error();

    ShowWindow(hwnd, SW_SHOW);

    return window;
}

void pwa_window_destroy(pwa_window_t* window) {
    if (window->buffer) {
        free(window->buffer);
    }

    free(window);
    window = NULL;
}

void pwa_window_set_resize_callback(pwa_window_t* window, pwa_resize_callback on_resize) {
    window->on_resize = on_resize;
}

void pwa_window_set_key_down_callback(pwa_window_t* window, pwa_key_down_callback on_key_down) {
    window->on_key_down = on_key_down;
}

void pwa_window_set_key_up_callback(pwa_window_t* window, pwa_key_up_callback on_key_up) {
    window->on_key_up = on_key_up;
}

bool pwa_window_should_close(pwa_window_t* window) {
    return window->should_close;
}

void pwa_window_poll_events(pwa_window_t* window) {
    MSG msg;
    while (PeekMessage(&msg, window->hwnd, 0, 0, PM_REMOVE)) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
}

void pwa_window_swap_buffers(pwa_window_t* window, psr_color_buffer_t* color_buffer) {
    bool resize = (color_buffer->w != window->buffer_w || color_buffer->h != window->buffer_h);
    
    if (resize) {
        free(window->buffer);
    }

    if (!window->buffer || resize) {
        window->buffer_w = color_buffer->w;
        window->buffer_h = color_buffer->h;
        window->buffer = malloc(color_buffer->w * color_buffer->h * sizeof(uint32_t));
    }

    for (int x = 0; x < color_buffer->w; x++) {
        for (int y = 0; y < color_buffer->h; y++) {
            psr_byte3_t color = *psr_color_buffer_at(color_buffer, x, y);
            window->buffer[y * color_buffer->w + x] = (color.r << 16) | (color.g << 8) | (color.b);
        }
    }
}

void pwa_window_schedule_redraw(pwa_window_t* window) {
    InvalidateRect(window->hwnd, NULL, FALSE);
}

double pwa_get_elapsed_time_ms() {
    LARGE_INTEGER elapsed_time;
    QueryPerformanceCounter(&elapsed_time);
    return (double)(elapsed_time.QuadPart * 1000000 / pwa_get_ticks_per_second()) / 1000;
}

int64_t pwa_get_ticks_per_second() {
    LARGE_INTEGER ticks_per_second;
    QueryPerformanceFrequency(&ticks_per_second);
    return ticks_per_second.QuadPart;
}

void pwa_print_last_error() {
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

void pwa_terminate() {
}

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
    pwa_window_t* window = (pwa_window_t*)GetWindowLongPtr(hwnd, 0);
    
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
            bmi.bmiHeader.biWidth = window->buffer_w;
            bmi.bmiHeader.biHeight = -window->buffer_h;
            bmi.bmiHeader.biPlanes = 1;
            bmi.bmiHeader.biBitCount = 32;
            bmi.bmiHeader.biCompression = BI_RGB;

            StretchDIBits(hdc, 0, 0, window->buffer_w, window->buffer_h, 0, 0, window->buffer_w, window->buffer_h, window->buffer, &bmi, DIB_RGB_COLORS, SRCCOPY);
        }
        else {
            FillRect(hdc, &ps.rcPaint, (HBRUSH)(COLOR_WINDOW + 1));
        }

        EndPaint(hwnd, &ps);
        return 0;

    case WM_SIZE:
        if (window && window->on_resize) {
            window->on_resize(LOWORD(lParam), HIWORD(lParam), window->user_data);
        }
        return 0;

    case WM_KEYDOWN:
        if (window && window->on_key_down) {
            window->on_key_down(wParam, window->user_data);
        }
        return 0;

    case WM_KEYUP:
        if (window && window->on_key_up) {
            window->on_key_up(wParam, window->user_data);
        }
        return 0;
    }

    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}
