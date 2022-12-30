#define SRZ_IMPLEMENTATION
#define SRZ_SAVE_AND_LOAD
#include <srz.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define WIDTH 640
#define HEIGHT 640

int main() {
    srz_mesh_t mesh;
    int result = srz_load_obj("bullfrog.obj", &mesh);

    if (result == 0) {
        printf("Failed to load model.\n");
        return -1;
    }
    else if (result == -1) {
        printf("Bad model.");
        return -1;
    }

    srz_color_buffer_t color_buffer;
    color_buffer.w = WIDTH;
    color_buffer.h = HEIGHT;
    color_buffer.data = malloc(WIDTH * HEIGHT * sizeof(srz_byte3_t));

    srz_depth_buffer_t depth_buffer;
    depth_buffer.w = WIDTH;
    depth_buffer.h = HEIGHT;
    depth_buffer.data = malloc(WIDTH * HEIGHT * sizeof(float));
    srz_depth_buffer_clear(&depth_buffer);

    // Gradient background.
    for (int y = 0; y < HEIGHT; ++y) {
        srz_byte3_t top = {141, 160, 184};
        srz_byte3_t bottom = {20, 20, 20};

        float amount = (float)y / HEIGHT;
        srz_byte3_t interp = {(bottom.r - top.r) * amount + top.r, (bottom.g - top.g) * amount + top.g, (bottom.b - top.b) * amount + top.b};

        for (int x = 0; x < WIDTH; ++x) {
            *srz_color_buffer_at(&color_buffer, x, y) = interp;
        }
    }

    srz_float3_t camera_pos = {15, 15, -25};
    srz_matrix_t view = srz_look_at(camera_pos, (srz_float3_t){0, 0, 0}, (srz_float3_t){0, -1, 0});
    srz_matrix_t projection = srz_perspective(HEIGHT / (float)WIDTH, 75.f * (SRZ_PI / 180), 0.1f, 1000.f);

    srz_matrix_t model;
    srz_matrix_init_identity(&model);
    model = srz_matrix_translate(model, (srz_float3_t){-.5, -.5, -.5});
    model = srz_matrix_rotate_x(model, 90.f * (SRZ_PI / 180));
    model = srz_matrix_rotate_y(model, 180.f * (SRZ_PI / 180));

    srz_matrix_t mvp = srz_matrix_mul(srz_matrix_mul(model, view), projection);

    // Vertex multiplication and scale to viewport.
    for (int i = 0; i < mesh.position_count; ++i) {
        mesh.positions[i] = srz_float3_mul_matrix(mesh.positions[i], mvp);

        // Scale from [-1, 1] to viewport size.
        mesh.positions[i].x = (mesh.positions[i].x + 1) / 2.f * WIDTH;
        mesh.positions[i].y = (mesh.positions[i].y + 1) / 2.f * HEIGHT;
    }

    // Render triangles.
    for (int i = 0; i < mesh.face_count; ++i) {
        srz_float3_t normal = mesh.normals[mesh.faces[i].normal_indices[0]];
        
        srz_float3_t light_dir = srz_normalize((srz_float3_t){0, 1, .5});
        float light = srz_dot(normal, light_dir);
        srz_byte3_t color = {255 * light, 255 * light, 255 * light};

        srz_float3_t tri[3];
        for (int j = 0; j < 3; ++j) {
            tri[j] = mesh.positions[mesh.faces[i].position_indices[j]];
        }
        srz_raster_triangle(&color_buffer, &depth_buffer, tri[0], tri[1], tri[2], color);
    }

    srz_save_bmp("scene.bmp", color_buffer);
    free(color_buffer.data);

    printf("Done.");
}
