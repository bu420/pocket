#include <psr.h>
#include <stdlib.h>
#include <stdio.h>

int main() {
    psr_color_buffer_t* canvas = psr_color_buffer_create(100, 100);
    psr_color_buffer_clear(canvas, (psr_byte3_t){0, 0, 0});

    psr_raster_line(canvas, (psr_int2_t){10, 10}, (psr_int2_t){90, 90}, (psr_byte3_t){255, 0, 0}, (psr_byte3_t){0, 255, 0});
    psr_raster_line(canvas, (psr_int2_t){10, 90}, (psr_int2_t){90, 10}, (psr_byte3_t){0, 0, 255}, (psr_byte3_t){255, 255, 0});

    psr_raster_triangle_2d_color(canvas, (psr_int2_t){10, 20}, (psr_int2_t){40, 50}, (psr_int2_t){10, 80}, (psr_byte3_t){255, 0, 0});
    psr_raster_triangle_2d_color(canvas, (psr_int2_t){20, 10}, (psr_int2_t){50, 40}, (psr_int2_t){80, 10}, (psr_byte3_t){0, 255, 0});
    psr_raster_triangle_2d_color(canvas, (psr_int2_t){90, 20}, (psr_int2_t){60, 50}, (psr_int2_t){90, 80}, (psr_byte3_t){0, 0, 255});
    psr_raster_triangle_2d_color(canvas, (psr_int2_t){20, 90}, (psr_int2_t){50, 60}, (psr_int2_t){80, 90}, (psr_byte3_t){255, 255, 0});

    psr_save_bmp("output.bmp", canvas);
    printf("Done.\n");
}
