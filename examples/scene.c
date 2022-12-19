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

    srz_byte_t cube[36 * 3] = {
        0,0,0,0,1,0,1,1,0,0,0,0,1,1,0,1,0,0,
        1,0,0,1,1,0,1,1,1,1,0,0,1,1,1,1,0,1,
        1,0,1,1,1,1,0,1,1,1,0,1,0,1,1,0,0,1,
        0,0,1,0,1,1,0,1,0,0,0,1,0,1,0,0,0,0,
        0,1,0,0,1,1,1,1,1,0,1,0,1,1,1,1,1,0,
        1,0,1,0,0,1,0,0,0,1,0,1,0,0,0,1,0,0
    };

    char normals[6 * 3] = {
        0, 0, 1,
        -1, 0, 0,
        0, 0, -1,
        1, 0, 0,
        0, -1, 0,
        0, 1, 0
    };

    srz_float3_t cube_positions[4] = {
        {0, 0, 0},
        {1.5, 0.25, 0},
        {3, 0.5, 0},
        {4.5, 0.75, 0}
    };

    srz_float3_t camera_pos = {0, -2, -8};
    srz_matrix_t view = srz_look_at(camera_pos, (srz_float3_t){3, 1, 0}, (srz_float3_t){0, -1, 0});
    srz_matrix_t projection = srz_perspective(HEIGHT / (float)WIDTH, 50.f * (SRZ_PI / 180), 0.1f, 1000.f);

    for (int i = 0; i < 4; ++i) {
        srz_matrix_t model;
        srz_matrix_init_identity(&model);
        model = srz_matrix_translate(model, cube_positions[i]);

        srz_matrix_t mvp = srz_matrix_mul(srz_matrix_mul(model, view), projection);

        float vertices[36 * 3];

        // Vertex multiplication and scale to viewport.
        for (int i = 0; i < 36; ++i) {
            int j = i * 3;

            float x = cube[j + 0];
            float y = cube[j + 1];
            float z = cube[j + 2];

            srz_float3_t vertex = {x, y, z};
            vertex = srz_float3_mul_matrix(vertex, mvp);

            for (int k = 0; k < 3; ++k) {
                vertices[j + k] = vertex.values[k];
            }

            // Scale from [-1, 1] to viewport size.
            vertices[j + 0] = (vertices[j + 0] + 1) / 2.f * WIDTH;
            vertices[j + 1] = (vertices[j + 1] + 1) / 2.f * HEIGHT;
        }

        // Render triangles.
        for (int i = 0; i < 12; ++i) {
            int j = i * 9;

            srz_float3_t v0 = {vertices[j], vertices[j + 1], vertices[j + 2]};
            srz_float3_t v1 = {vertices[j + 3], vertices[j + 3 + 1], vertices[j + 3 + 2]};
            srz_float3_t v2 = {vertices[j + 6], vertices[j + 6 + 1], vertices[j + 6 + 2]};

            int normal_i = srz_floorf(i / 2.f) * 3;
            srz_float3_t normal = {normals[normal_i], normals[normal_i + 1], normals[normal_i + 2]};

            // Discard invisible triangles.
            if (normal.x * (v0.x - camera_pos.x) + normal.y * (v0.y - camera_pos.y) + normal.z * (v0.z - camera_pos.z) < 0) {
                continue;
            }

            srz_float3_t light_dir = srz_normalize((srz_float3_t){0, 1, .5});
            float light = srz_dot(normal, light_dir);
            srz_byte3_t color = {255 * light, 255 * light, 255 * light};
            
            srz_int2_t p0 = {srz_roundf(v0.x), srz_roundf(v0.y)};
            srz_int2_t p1 = {srz_roundf(v1.x), srz_roundf(v1.y)};
            srz_int2_t p2 = {srz_roundf(v2.x), srz_roundf(v2.y)};

            srz_image_raster_tri(&image, p0, p1, p2, color);
            /*srz_image_raster_line(&image, p0, p1, color);
            srz_image_raster_line(&image, p1, p2, color);
            srz_image_raster_line(&image, p2, p0, color);*/
        }
    }

    srz_image_save_bmp("scene.bmp", image);
    free(image_mem);

    printf("Done.");
}
