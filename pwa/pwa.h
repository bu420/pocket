#ifndef PWA_H
#define PWA_H

#include <stdint.h>

typedef struct pwa_window pwa_window_t;

typedef struct {
    uint32_t* pixels;
    int w;
    int h;
} pwa_pixel_buffer_t;

typedef void (*pwa_resize_callback)(int w, int h, void* user_data);
typedef pwa_pixel_buffer_t (*pwa_draw_callback)(void* user_data);

void pwa_init();
pwa_window_t* pwa_window_create(char* title, int w, int h, void* user_data);
void pwa_window_destroy(pwa_window_t* window);
void pwa_window_set_resize_callback(pwa_window_t* window, pwa_resize_callback on_resize);
void pwa_window_set_draw_callback(pwa_window_t* window, pwa_draw_callback on_draw);
int pwa_window_should_close(pwa_window_t* window);
void pwa_window_poll_events(pwa_window_t* window);

#endif
