#ifndef PWA_H
#define PWA_H

#include <psr/psr.h>
#include <stdint.h>
#include <stdbool.h>

typedef struct pwa_window_t pwa_window_t;

typedef void (*pwa_resize_callback)(int w, int h, void* user_data);
typedef void (*pwa_key_down_callback)(int key_code, void* user_data);
typedef void (*pwa_key_up_callback)(int key_code, void* user_data);

void pwa_init();

pwa_window_t* pwa_window_create(char* title, int w, int h, void* user_data);
void pwa_window_destroy(pwa_window_t* window);
void pwa_window_set_resize_callback(pwa_window_t* window, pwa_resize_callback on_resize);
void pwa_window_set_key_down_callback(pwa_window_t* window, pwa_key_down_callback on_key_down);
void pwa_window_set_key_up_callback(pwa_window_t* window, pwa_key_up_callback on_key_up);
bool pwa_window_should_close(pwa_window_t* window);
void pwa_window_poll_events(pwa_window_t* window);
void pwa_window_swap_buffers(pwa_window_t* window, psr_color_buffer_t* color_buffer);
void pwa_window_schedule_redraw(pwa_window_t* window);

double pwa_get_elapsed_time_ms();
int64_t pwa_get_ticks_per_second();

void pwa_print_last_error();

void pwa_terminate();

#endif
