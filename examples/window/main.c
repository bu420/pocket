#include <psr/psr.h>
#include <pwa/pwa.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void on_resize(int w, int h, void* user_data) {
    psr_color_buffer_t* buf = (psr_color_buffer_t*)user_data;
    psr_color_buffer_resize(buf, w, h);
}

void on_key_down(int key_code, void* user_data) {
    printf("Key down: %d\n", key_code);
}

void on_key_up(int key_code, void* user_data) {
    printf("Key up: %d\n", key_code);
}

int main() {
    pwa_init();

    psr_color_buffer_t* color_buffer = psr_color_buffer_create(400, 400);

    pwa_window_t* window = pwa_window_create("Example Window", 400, 400, color_buffer);

    if (!window) {
        printf("Window error.\n");
        return -1;
    }

    pwa_window_set_resize_callback(window, on_resize);
    pwa_window_set_key_down_callback(window, on_key_down);
    pwa_window_set_key_up_callback(window, on_key_up);

    while (!pwa_window_should_close(window)) {
        pwa_window_poll_events(window);

        double time = pwa_get_elapsed_time_ms() / 500;

        for (int x = 0; x < color_buffer->w; x++) {
            for (int y = 0; y < color_buffer->h; y++) {
                psr_byte3_t* pixel = psr_color_buffer_at(color_buffer, x, y);

                pixel->r = (sin(time) + 1) / 2 * 255;
                pixel->g = 255 - pixel->r;
                pixel->b = (cos(time) + 1) / 2 * 255;
            }
        }

        pwa_window_swap_buffers(window, color_buffer);
        pwa_window_schedule_redraw(window);
    }

    pwa_window_destroy(window);
    printf("Done.\n");
}
