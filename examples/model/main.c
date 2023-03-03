#include <psr.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define WIDTH 640
#define HEIGHT 640

int main() {
    psr_mesh_t mesh;
    int result = psr_load_obj("bullfrog.obj", &mesh);

    if (result == 0) {
        printf("Failed to load model.\n");
        return -1;
    }
    else if (result == -1) {
        printf("Bad model.");
        return -1;
    }

    psr_color_buffer_t color_buffer;
    color_buffer.w = WIDTH;
    color_buffer.h = HEIGHT;
    color_buffer.data = malloc(WIDTH * HEIGHT * sizeof(psr_byte3_t));

    psr_depth_buffer_t depth_buffer;
    depth_buffer.w = WIDTH;
    depth_buffer.h = HEIGHT;
    depth_buffer.data = malloc(WIDTH * HEIGHT * sizeof(float));
    psr_depth_buffer_clear(&depth_buffer);

    // Gradient background.
    for (int y = 0; y < HEIGHT; y++) {
        psr_byte3_t top = {141, 160, 184};
        psr_byte3_t bottom = {20, 20, 20};

        float amount = (float)y / HEIGHT;
        psr_byte3_t interp = {(bottom.r - top.r) * amount + top.r, (bottom.g - top.g) * amount + top.g, (bottom.b - top.b) * amount + top.b};

        for (int x = 0; x < WIDTH; x++) {
            *psr_color_buffer_at(&color_buffer, x, y) = interp;
        }
    }

    psr_float3_t camera_pos = {15, 15, -25};
    psr_matrix_t view = psr_look_at(camera_pos, (psr_float3_t){0, 0, 0}, (psr_float3_t){0, -1, 0});
    psr_matrix_t projection = psr_perspective(HEIGHT / (float)WIDTH, 75.f * (M_PI / 180), 0.1f, 1000.f);

    psr_matrix_t model;
    psr_matrix_init_identity(&model);
    model = psr_matrix_translate(model, (psr_float3_t){-.5, -.5, -.5});
    model = psr_matrix_rotate_x(model, 90.f * (M_PI / 180));
    model = psr_matrix_rotate_y(model, 180.f * (M_PI / 180));

    psr_matrix_t mvp = psr_matrix_mul(psr_matrix_mul(model, view), projection);

    // Vertex multiplication and scale to viewport.
    for (int i = 0; i < mesh.position_count; i++) {
        mesh.positions[i] = psr_float3_mul_matrix(mesh.positions[i], mvp);

        // Scale from [-1, 1] to viewport size.
        mesh.positions[i].x = (mesh.positions[i].x + 1) / 2.f * WIDTH;
        mesh.positions[i].y = (mesh.positions[i].y + 1) / 2.f * HEIGHT;
    }

    // Render triangles.
    for (int i = 0; i < mesh.face_count; i++) {
        psr_float3_t normal = mesh.normals[mesh.faces[i].normal_indices[0]];
        
        psr_float3_t light_dir = psr_normalize((psr_float3_t){0, 1, .5});
        float light = psr_dot(normal, light_dir);
        psr_byte3_t color = {255 * light, 255 * light, 255 * light};

        psr_int2_t tri[3];
        for (int j = 0; j < 3; j++) {
            tri[j].x = (int)roundf(mesh.positions[mesh.faces[i].position_indices[j]].x);
            tri[j].y = (int)roundf(mesh.positions[mesh.faces[i].position_indices[j]].y);
        }
        psr_raster_triangle_2d_color(&color_buffer, tri[0], tri[1], tri[2], color);
    }

    psr_save_bmp("model.bmp", color_buffer);
    printf("Done.");
}
