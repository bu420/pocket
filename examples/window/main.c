#include <pwa.h>
#include <stdlib.h>
#include <stdio.h>

pwa_pixel_buffer_t create_pixel_buffer(int w, int h) {
    pwa_pixel_buffer_t buffer;
    buffer.pixels = malloc(w * h * sizeof(uint32_t));
    buffer.w = w;
    buffer.h = h;
    
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h / 2; y++) {
            buffer.pixels[y * w + x] = 0xff0000;
        }
        for (int y = h / 2; y < h; y++) {
            buffer.pixels[y * w + x] = 0x0000ff;
        }
    }
    for (int x = w / 2 - w / 8; x < w / 2 + w / 8; x++) {
        for (int y = 0; y < h; y++) {
            buffer.pixels[y * w + x] = 0x00ff00;
        }
    }
    return buffer;
}

void free_pixel_buffer(pwa_pixel_buffer_t* pixel_buffer) {
    free(pixel_buffer->pixels);
    pixel_buffer->pixels = NULL;
}

void on_resize(int w, int h, void* user_data) {
    pwa_pixel_buffer_t* buffer = (pwa_pixel_buffer_t*)user_data;

    free_pixel_buffer(buffer);
    *buffer = create_pixel_buffer(w, h);
}

pwa_pixel_buffer_t on_draw(void* user_data) {
    return *(pwa_pixel_buffer_t*)user_data;
}

int main() {
    pwa_init();

    pwa_pixel_buffer_t buffer = create_pixel_buffer(400, 400);

    pwa_window_t* window = pwa_window_create("Example Window", 400, 400, &buffer);

    if (!window) {
        printf("Window error.\n");
        return -1;
    }

    pwa_window_set_resize_callback(window, on_resize);
    pwa_window_set_draw_callback(window, on_draw);

    while (!pwa_window_should_close(window)) {
        pwa_window_poll_events(window);
    }

    pwa_window_destroy(window);
    printf("Done.\n");
}
