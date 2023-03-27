#include <psr.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define WIDTH 640
#define HEIGHT 640

int main() {
    psr_image_t image = psr_load_bmp("input.bmp", PSR_R8G8B8);

    psr_color_buffer_t color_buffer;
    psr_color_buffer_init(&color_buffer, image.w, image.h);
    psr_color_buffer_clear(&color_buffer, (psr_byte3_t){0, 0, 0});

    psr_raster_image(&color_buffer, image, 0, 0, 0, 0, image.w, image.h);

    psr_save_bmp("output.bmp", color_buffer);
    printf("Done.");
}
