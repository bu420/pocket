#include <pwa.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

pwa_pixel_buffer_t create_pixel_buffer(int w, int h) {
    pwa_pixel_buffer_t buffer;
    buffer.pixels = malloc(w * h * sizeof(uint32_t));
    buffer.w = w;
    buffer.h = h;

    // Gradient.
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            uint8_t r = (uint8_t)((float)x / w * 255);
            uint8_t g = (uint8_t)((float)y / h * 255);
            uint8_t b = 255 - (uint8_t)((float)y / h * 255);

            buffer.pixels[y * w + x] = (r << 16) | (g << 8) | (b);
        }
    }
    
    return buffer;
}

void on_resize(int w, int h, void* user_data) {
    pwa_pixel_buffer_t* buffer = (pwa_pixel_buffer_t*)user_data;

    free(buffer->pixels);
    *buffer = create_pixel_buffer(w, h);
}

pwa_pixel_buffer_t on_draw(void* user_data) {
    return *(pwa_pixel_buffer_t*)user_data;
}

void on_key_down(int key_code, void* user_data) {
    printf("Key down: %d\n", key_code);
}

void on_key_up(int key_code, void* user_data) {
    printf("Key up: %d\n", key_code);
}

int main() {
    pwa_init();

    pwa_pixel_buffer_t buffer = create_pixel_buffer(400, 400);

    pwa_window_t* window = pwa_window_create("Example Window", 400, 400, &buffer);

    if (!window) {
        printf("Window error.\n");
        return -1;
    }

    pwa_set_resize_callback(window, on_resize);
    pwa_set_draw_callback(window, on_draw);
    pwa_set_key_down_callback(window, on_key_down);
    pwa_set_key_up_callback(window, on_key_up);

    while (!pwa_window_should_close(window)) {
        pwa_poll_events(window);
        pwa_schedule_redraw(window);
    }

    pwa_window_destroy(window);
    printf("Done.\n");
}
