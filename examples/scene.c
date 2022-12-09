#include <stdlib.h>
#include <stdio.h>

#define OSR_IMPLEMENTATION
#define OSR_SAVE_AND_LOAD
#include "osr.h"

#define WIDTH 640
#define HEIGHT 640

void print_matrix(osr_matrix_t m) {
    printf("----------------------------\n");
    for (int y = 0; y < 4; ++y) {
        for (int x = 0; x < 4; ++x) {
            printf("%.2f\t", *osr_matrix_at(&m, x, y));
        }
        printf("\n");
    }
    printf("----------------------------\n\n");
}

void print_float4(osr_float4_t f) {
    printf("----------------------------\n");
    printf("x=%.2f y=%.2f z=%.2f w=%.2f\n", f.x, f.y, f.z, f.w);
    printf("----------------------------\n\n");
}

int main() {
    osr_image_t image;
    void* image_mem = malloc(WIDTH * HEIGHT * sizeof(osr_byte3_t));
    osr_image_init(&image, WIDTH, HEIGHT, image_mem);

    // Gradient background.
    for (int x = 0; x < WIDTH; ++x) {
        for (int y = 0; y < HEIGHT; ++y) {
            osr_byte3_t* pixel = osr_image_at(&image, x, y);

            pixel->r = x * 128 / WIDTH;
            pixel->g = 30;
            pixel->b = (x + y) * 128 / (WIDTH + HEIGHT);
        }
    }

    // Position and color.
    float cube_vertices[36 * 6] = {
        -0.5, -0.5, -0.5, 255, 0, 0,
        0.5, -0.5, -0.5, 255, 0, 0,
        0.5, 0.5, -0.5, 255, 0, 0,
        0.5, 0.5, -0.5, 255, 0, 0,
        -0.5, 0.5, -0.5, 255, 0, 0,
        -0.5, -0.5, -0.5, 255, 0, 0,

        -0.5, -0.5, 0.5, 0, 255, 0,
        0.5, -0.5, 0.5, 0, 255, 0,
        0.5, 0.5, 0.5, 0, 255, 0,
        0.5, 0.5, 0.5, 0, 255, 0,
        -0.5, 0.5, 0.5, 0, 255, 0,
        -0.5, -0.5, 0.5, 0, 255, 0,

        -0.5, 0.5, 0.5, 0, 0, 255,
        -0.5, 0.5, -0.5, 0, 0, 255,
        -0.5, -0.5, -0.5, 0, 0, 255,
        -0.5, -0.5, -0.5, 0, 0, 255,
        -0.5, -0.5, 0.5, 0, 0, 255,
        -0.5, 0.5, 0.5, 0, 0, 255,

        0.5, 0.5, 0.5, 255, 0, 255,
        0.5, 0.5, -0.5, 255, 0, 255,
        0.5, -0.5, -0.5, 255, 0, 255,
        0.5, -0.5, -0.5, 255, 0, 255,
        0.5, -0.5, 0.5, 255, 0, 255,
        0.5, 0.5, 0.5, 255, 0, 255,

        -0.5, -0.5, -0.5, 255, 255, 0,
        0.5, -0.5, -0.5, 255, 255, 0,
        0.5, -0.5, 0.5, 255, 255, 0,
        0.5, -0.5, 0.5, 255, 255, 0,
        -0.5, -0.5, 0.5, 255, 255, 0,
        -0.5, -0.5, -0.5, 255, 255, 0,

        -0.5, 0.5, -0.5, 255, 255, 255,
        0.5, 0.5, -0.5, 255, 255, 255,
        0.5, 0.5, 0.5, 255, 255, 255,
        0.5, 0.5, 0.5, 255, 255, 255,
        -0.5, 0.5, 0.5, 255, 255, 255,
        -0.5, 0.5, -0.5, 255, 255, 255
    };

    osr_matrix_t model;
    osr_matrix_init_identity(&model);
    //model = osr_osr_matrix_translate(model, (osr_float3_t){3, 3, 3});
    //model = osr_osr_matrix_translate(model, (osr_float3_t){WIDTH / 2, HEIGHT / 2, 0});
    //model = osr_matrix_rotate(model, (osr_float3_t){0, 0, OSR_PI / 4});
    //model = osr_matrix_scale(model, (osr_float3_t){2, 2, 2});

    osr_matrix_t view = osr_create_view_matrix((osr_float3_t){0, 0, 8}, (osr_float3_t){0, 0, -1}, (osr_float3_t){0, -1, 0});
    osr_matrix_t projection = osr_create_projection_matrix(HEIGHT / (float)WIDTH, 50.f, 0.1f, 1000.f);

    osr_matrix_t mvp = osr_matrix_mul(osr_matrix_mul(projection, view), model);

    // Multiply cube with matrix.
    for (int i = 0; i < 36; ++i) {
        int j = i * 6;
        osr_float3_t pos = {cube_vertices[j + 0], cube_vertices[j + 1], cube_vertices[j + 2]};
        osr_float3_t result = osr_float3_mul_matrix(pos, mvp);
        //osr_float3_t result = osr_matrix_mul_float3(mvp, pos);

        for (int k = 0; k < 3; ++k) {
            cube_vertices[j + k] = result.values[k];
        }

        // Scale from [-1, 1] to image size.
        cube_vertices[j + 0] = (cube_vertices[j + 0] + 1) / 2.f * WIDTH;
        cube_vertices[j + 1] = (cube_vertices[j + 1] + 1) / 2.f * HEIGHT;
    }

    // Render cube.
    for (int i = 0; i < 12; ++i) {
        int j = i * 18;
        
        osr_int2_t p0 = {osr_roundf(cube_vertices[j]), osr_roundf(cube_vertices[j + 1])};
        osr_int2_t p1 = {osr_roundf(cube_vertices[j + 6]), osr_roundf(cube_vertices[j + 6 + 1])};
        osr_int2_t p2 = {osr_roundf(cube_vertices[j + 12]), osr_roundf(cube_vertices[j + 12 + 1])};

        osr_byte3_t color = {cube_vertices[j + 3], cube_vertices[j + 4], cube_vertices[j + 5]};

        //osr_image_raster_tri(&image, p0, p1, p2, color);
        osr_image_raster_line(&image, p0, p1, color);
        osr_image_raster_line(&image, p1, p2, color);
        osr_image_raster_line(&image, p2, p0, color);
    }

    osr_image_save_bmp("scene.bmp", image);
    free(image_mem);

    printf("Done.");
}
