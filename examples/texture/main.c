#define SRZ_IMPLEMENTATION
#define SRZ_SAVE_AND_LOAD
#include <srz.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define WIDTH 640
#define HEIGHT 640

int main() {
    srz_color_buffer_t color_buffer;
    color_buffer.w = WIDTH;
    color_buffer.h = HEIGHT;
    color_buffer.data = malloc(WIDTH * HEIGHT * sizeof(srz_byte3_t));

    srz_color_buffer_clear(&color_buffer, (srz_byte3_t){200, 200, 200});

    srz_texture_t texture;
    texture.w = 8;
    texture.h = 8;
    texture.data = malloc(8 * 8 * sizeof(srz_byte4_t));

    // Generate checkerboard pattern. 
    for (int x = 0; x < 8; x++) {
        for (int y = 0; y < 8; y++) {
            if ((x + y) % 2 == 0) {
                *srz_texture_at(&texture, x, y) = (srz_byte4_t){0, 0, 0, 255};
            }
            else {
                *srz_texture_at(&texture, x, y) = (srz_byte4_t){255, 255, 255, 255};
            }
        }
    }
    
    srz_int2_t positions[4] = {
        {220, 220},
        {220, 420},
        {420, 420},
        {420, 220}
    };

    srz_float3_t tex_coords[4] = {
        {0, 0, 1},
        {0, 1, 1},
        {1, 1, 1},
        {1, 0, 1}
    };

    srz_raster_triangle_texture(&color_buffer, NULL, positions[2], positions[3], positions[0], tex_coords[0], tex_coords[1], tex_coords[2], texture);
    srz_raster_triangle_texture(&color_buffer, NULL, positions[2], positions[3], positions[0], tex_coords[2], tex_coords[3], tex_coords[0], texture);
    srz_raster_triangle(&color_buffer, NULL, (srz_float3_t){positions[0].x, positions[0].y, 0}, (srz_float3_t){positions[1].x, positions[1].y, 0}, (srz_float3_t){positions[2].x, positions[2].y, 0}, (srz_byte3_t){0, 255, 120});
    srz_raster_triangle(&color_buffer, NULL, (srz_float3_t){positions[2].x, positions[2].y, 0}, (srz_float3_t){positions[3].x, positions[3].y, 0}, (srz_float3_t){positions[0].x, positions[0].y, 0}, (srz_byte3_t){0, 255, 120});

    // Line plotting test.
    srz_int2_t steps[200];
    int count = srz_bresenham_line((srz_int2_t){0, 30}, (srz_int2_t){10, 0}, steps);
    for (int i = 0; i < count; i++) {
        printf("x = %d, y = %d\n", steps[i].x, steps[i].y);
        *srz_color_buffer_at(&color_buffer, steps[i].x, steps[i].y) = (srz_byte3_t){255, 0, 0};
    }

    srz_save_bmp("texture.bmp", color_buffer);
    printf("Done.");
}
