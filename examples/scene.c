#define SRZ_IMPLEMENTATION
#define SRZ_SAVE_AND_LOAD
#include <srz.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define WIDTH 640
#define HEIGHT 640

int main() {
    srz_image_t image;
    void* image_mem = malloc(WIDTH * HEIGHT * sizeof(srz_byte3_t));
    srz_image_init(&image, WIDTH, HEIGHT, image_mem);

    // Gradient background.
    for (int x = 0; x < WIDTH; ++x) {
        for (int y = 0; y < HEIGHT; ++y) {
            srz_byte3_t* pixel = srz_image_at(&image, x, y);

            pixel->r = x * 128 / WIDTH;
            pixel->g = 0;
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

    srz_float3_t cube_positions[4] = {
        {0, 0, 0},
        {1, 0.25, -0.125},
        {2, 0.5, -0.25},
        {3, 0.75, -0.375}
    };

    srz_matrix_t view = srz_look_at((srz_float3_t){4, 0, 8}, (srz_float3_t){2, 1, 0}, (srz_float3_t){0, 1, 0});
    srz_matrix_t projection = srz_perspective(HEIGHT / (float)WIDTH, 50.f * (SRZ_PI / 180), 0.1f, 1000.f);

    for (int i = 0; i < 4; ++i) {
        srz_matrix_t model;
        srz_matrix_init_identity(&model);
        model = srz_matrix_translate(model, cube_positions[i]);

        srz_matrix_t mvp = srz_matrix_mul(srz_matrix_mul(model, view), projection);

        float vertices_copy[36 * 6];
        memcpy(vertices_copy, cube_vertices, sizeof(float) * 36 * 6);

        // Multiply cube with matrix.
        for (int i = 0; i < 36; ++i) {
            int j = i * 6;
            srz_float3_t pos = {vertices_copy[j + 0], vertices_copy[j + 1], vertices_copy[j + 2]};
            srz_float3_t result = srz_float3_mul_matrix(pos, mvp);

            for (int k = 0; k < 3; ++k) {
                vertices_copy[j + k] = result.values[k];
            }

            // Scale from [-1, 1] to image size.
            vertices_copy[j + 0] = (vertices_copy[j + 0] + 1) / 2.f * WIDTH;
            vertices_copy[j + 1] = (vertices_copy[j + 1] + 1) / 2.f * HEIGHT;
        }

        // Render cube.
        for (int i = 0; i < 12; ++i) {
            int j = i * 18;

            //if (vertices_copy[j + 2] )
            
            srz_int2_t p0 = {srz_roundf(vertices_copy[j]), srz_roundf(vertices_copy[j + 1])};
            srz_int2_t p1 = {srz_roundf(vertices_copy[j + 6]), srz_roundf(vertices_copy[j + 6 + 1])};
            srz_int2_t p2 = {srz_roundf(vertices_copy[j + 12]), srz_roundf(vertices_copy[j + 12 + 1])};

            // Discard triangles behind view.
            //if ()

            srz_byte3_t color = {vertices_copy[j + 3], vertices_copy[j + 4], vertices_copy[j + 5]};

            //srz_image_raster_tri(&image, p0, p1, p2, color);
            srz_image_raster_line(&image, p0, p1, color);
            srz_image_raster_line(&image, p1, p2, color);
            srz_image_raster_line(&image, p2, p0, color);
        }
    }

    srz_image_save_bmp("scene.bmp", image);
    free(image_mem);

    printf("Done.");
}
