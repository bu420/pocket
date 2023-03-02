#include <llsr.h>
#include <llsr_io.h>
#include <stdlib.h>
#include <stdio.h>

int main() {
    llsr_color_buffer_t canvas;
    canvas.w = 100;
    canvas.h = 100;
    canvas.data = malloc(100 * 100 * sizeof(llsr_byte3_t));
    llsr_color_buffer_clear(&canvas, (llsr_byte3_t){0, 0, 0});

    llsr_raster_line(&canvas, (llsr_int2_t){10, 10}, (llsr_int2_t){90, 90}, (llsr_byte3_t){255, 0, 0}, (llsr_byte3_t){0, 255, 0});
    llsr_raster_line(&canvas, (llsr_int2_t){10, 90}, (llsr_int2_t){90, 10}, (llsr_byte3_t){0, 0, 255}, (llsr_byte3_t){255, 255, 0});

    llsr_raster_triangle_2d(&canvas, (llsr_int2_t){10, 20}, (llsr_int2_t){40, 50}, (llsr_int2_t){10, 80}, (llsr_byte3_t){255, 0, 0});
    llsr_raster_triangle_2d(&canvas, (llsr_int2_t){20, 10}, (llsr_int2_t){50, 40}, (llsr_int2_t){80, 10}, (llsr_byte3_t){0, 255, 0});
    llsr_raster_triangle_2d(&canvas, (llsr_int2_t){90, 20}, (llsr_int2_t){60, 50}, (llsr_int2_t){90, 80}, (llsr_byte3_t){0, 0, 255});
    llsr_raster_triangle_2d(&canvas, (llsr_int2_t){20, 90}, (llsr_int2_t){50, 60}, (llsr_int2_t){80, 90}, (llsr_byte3_t){255, 255, 0});

    llsr_save_bmp("shapes.bmp", canvas);
    printf("Done.\n");
}
