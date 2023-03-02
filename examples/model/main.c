#include <llsr.h>
#include <llsr_io.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define WIDTH 640
#define HEIGHT 640

int main() {
    llsr_mesh_t mesh;
    int result = llsr_load_obj("bullfrog.obj", &mesh);

    if (result == 0) {
        printf("Failed to load model.\n");
        return -1;
    }
    else if (result == -1) {
        printf("Bad model.");
        return -1;
    }

    llsr_color_buffer_t color_buffer;
    color_buffer.w = WIDTH;
    color_buffer.h = HEIGHT;
    color_buffer.data = malloc(WIDTH * HEIGHT * sizeof(llsr_byte3_t));

    llsr_depth_buffer_t depth_buffer;
    depth_buffer.w = WIDTH;
    depth_buffer.h = HEIGHT;
    depth_buffer.data = malloc(WIDTH * HEIGHT * sizeof(float));
    llsr_depth_buffer_clear(&depth_buffer);

    // Gradient background.
    for (int y = 0; y < HEIGHT; y++) {
        llsr_byte3_t top = {141, 160, 184};
        llsr_byte3_t bottom = {20, 20, 20};

        float amount = (float)y / HEIGHT;
        llsr_byte3_t interp = {(bottom.r - top.r) * amount + top.r, (bottom.g - top.g) * amount + top.g, (bottom.b - top.b) * amount + top.b};

        for (int x = 0; x < WIDTH; x++) {
            *llsr_color_buffer_at(&color_buffer, x, y) = interp;
        }
    }

    llsr_float3_t camera_pos = {15, 15, -25};
    llsr_matrix_t view = llsr_look_at(camera_pos, (llsr_float3_t){0, 0, 0}, (llsr_float3_t){0, -1, 0});
    llsr_matrix_t projection = llsr_perspective(HEIGHT / (float)WIDTH, 75.f * (LLSR_PI / 180), 0.1f, 1000.f);

    llsr_matrix_t model;
    llsr_matrix_init_identity(&model);
    model = llsr_matrix_translate(model, (llsr_float3_t){-.5, -.5, -.5});
    model = llsr_matrix_rotate_x(model, 90.f * (LLSR_PI / 180));
    model = llsr_matrix_rotate_y(model, 180.f * (LLSR_PI / 180));

    llsr_matrix_t mvp = llsr_matrix_mul(llsr_matrix_mul(model, view), projection);

    // Vertex multiplication and scale to viewport.
    for (int i = 0; i < mesh.position_count; i++) {
        mesh.positions[i] = llsr_float3_mul_matrix(mesh.positions[i], mvp);

        // Scale from [-1, 1] to viewport size.
        mesh.positions[i].x = (mesh.positions[i].x + 1) / 2.f * WIDTH;
        mesh.positions[i].y = (mesh.positions[i].y + 1) / 2.f * HEIGHT;
    }

    // Render triangles.
    for (int i = 0; i < mesh.face_count; i++) {
        llsr_float3_t normal = mesh.normals[mesh.faces[i].normal_indices[0]];
        
        llsr_float3_t light_dir = llsr_normalize((llsr_float3_t){0, 1, .5});
        float light = llsr_dot(normal, light_dir);
        llsr_byte3_t color = {255 * light, 255 * light, 255 * light};

        llsr_int2_t tri[3];
        for (int j = 0; j < 3; j++) {
            tri[j].x = llsr_roundf(mesh.positions[mesh.faces[i].position_indices[j]].x);
            tri[j].y = llsr_roundf(mesh.positions[mesh.faces[i].position_indices[j]].y);
        }
        llsr_raster_triangle_2d(&color_buffer, tri[0], tri[1], tri[2], color);
    }

    llsr_save_bmp("model.bmp", color_buffer);
    printf("Done.");
}
