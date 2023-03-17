#include <psr.h>
#include <pwa.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define WIDTH 640
#define HEIGHT 640

pwa_pixel_buffer_t on_draw(void* user_data) {
    return *(pwa_pixel_buffer_t*)user_data;
}

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
    psr_color_buffer_init(&color_buffer, WIDTH, HEIGHT);

    psr_depth_buffer_t depth_buffer;
    psr_depth_buffer_init(&depth_buffer, WIDTH, HEIGHT);

    pwa_pixel_buffer_t pixel_buffer;
    pixel_buffer.pixels = malloc(WIDTH * HEIGHT * sizeof(uint32_t));
    pixel_buffer.w = WIDTH;
    pixel_buffer.h = HEIGHT;

    pwa_init();
    pwa_window_t* window = pwa_window_create("Spinning Model", WIDTH, HEIGHT, &pixel_buffer);

    if (!window) {
        printf("Window error.\n");
        return -1;
    }
    
    pwa_window_set_draw_callback(window, on_draw);

    while (!pwa_window_should_close(window)) {
        pwa_window_poll_events(window);

        psr_depth_buffer_clear(&depth_buffer);

        // Gradient background.
        for (int y = 0; y < HEIGHT; y++) {
            static psr_byte3_t top = {141, 160, 184};
            static psr_byte3_t bottom = {20, 20, 20};

            for (int x = 0; x < WIDTH; x++) {
                *psr_color_buffer_at(&color_buffer, x, y) = psr_byte3_lerp(top, bottom, (float)y / HEIGHT);
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

        // Copy color buffer into pixel buffer.
        for (int x = 0; x < WIDTH; x++) {
            for (int y = 0; y < HEIGHT; y++) {
                psr_byte3_t color = *psr_color_buffer_at(&color_buffer, x, y);
                pixel_buffer.pixels[y * WIDTH + x] = (color.r << 16) | (color.g << 8) | (color.b);
            }
        }

        pwa_window_schedule_redraw(window);
    }

    pwa_window_destroy(window);

    printf("Done.");
}
