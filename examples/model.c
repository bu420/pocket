#include <psr.h>
#include <pwa.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define WIDTH 640
#define HEIGHT 640

int pause_flag = 0;
float spin_animation = 0;

void on_key_down(int key_code, void* user_data) {
    // Hit space to pause/play animation.
    if (key_code == ' ') {
        pause_flag = !pause_flag;
    }
}

psr_byte3_t pixel_shader(psr_int2_t pixel_pos, const psr_attribute_array_t* interpolated, void* user_data) {
    psr_float3_t normal = PSR_ATTRIB_TO_FLOAT3(interpolated->attributes[0]);
    return (psr_byte3_t){
        fmax(fmin(roundf((normal.x + 1) / 2 * 255), 255), 0), 
        fmax(fmin(roundf((normal.y + 1) / 2 * 255), 255), 0), 
        fmax(fmin(roundf((normal.z + 1) / 2 * 255), 255), 0)};

    /*psr_byte3_t obj_color = {50, 255, 150};

    psr_float3_t light_dir = *(psr_float3_t*)user_data;
    
    float light = psr_dot(normal, light_dir);
    if (light < 0) {
        light = 0;
    }
    
    psr_byte3_t color = {obj_color.r * light, obj_color.g * light, obj_color.b * light};
    
    return color;*/
}

int main() {
    psr_mesh_t* mesh = psr_mesh_load_obj("assets/bullfrog.obj");
    assert(mesh);
    psr_image_t* font_img = psr_image_load_bmp("assets/font_dos_vga.bmp", PSR_R8G8B8A8);
    assert(font_img);
    psr_font_t* font = psr_font_load(font_img, "assets/font_dos_vga.txt");
    assert(font);

    psr_color_buffer_t* color_buffer = psr_color_buffer_create(WIDTH, HEIGHT);
    psr_depth_buffer_t* depth_buffer = psr_depth_buffer_create(WIDTH, HEIGHT);

    pwa_init();
    pwa_window_t* window = pwa_window_create("Spinning Model", WIDTH, HEIGHT, NULL);

    if (!window) {
        printf("Window error.\n");
        return -1;
    }
    
    pwa_window_set_key_down_callback(window, on_key_down);

    psr_float3_t camera_pos = {15, 15, -25};
    psr_float3_t camera_target = {0, 0, 0};
    psr_mat4_t view = psr_look_at(camera_pos, camera_target, (psr_float3_t){0, -1, 0});
    psr_mat4_t projection = psr_perspective(HEIGHT / (float)WIDTH, 75.f * (M_PI / 180), 0.1f, 1000.f);

    double last_frame = pwa_get_elapsed_time_ms();

    const int gui_delta_interval = 500;
    double gui_delta_last = pwa_get_elapsed_time_ms();
    double gui_delta_value = 0;
    float gui_render_time = 0;

    while (!pwa_window_should_close(window)) {
        double current_frame = pwa_get_elapsed_time_ms();

        pwa_window_poll_events(window);

        float delta = current_frame - last_frame;
        last_frame = current_frame;

        if (!pause_flag) {
            spin_animation += delta * M_PI / 1000;
        }

        psr_depth_buffer_clear(depth_buffer);

        // Gradient background.
        /*for (int y = 0; y < HEIGHT; y++) {
            static psr_byte3_t top = {141, 160, 184};
            static psr_byte3_t bottom = {20, 20, 20};

            for (int x = 0; x < WIDTH; x++) {
                *psr_color_buffer_at(color_buffer, x, y) = psr_byte3_lerp(top, bottom, (float)y / HEIGHT);
            }
        }*/
        psr_color_buffer_clear(color_buffer, (psr_byte3_t){255, 255, 255});

        psr_mat4_t model;
        psr_mat4_init_identity(&model);
        model = psr_mat4_translate(model, (psr_float3_t){-.5, -.5, -.5});
        model = psr_mat4_rotate_x(model, 90.f * (M_PI / 180));
        model = psr_mat4_rotate_y(model, spin_animation);
        
        psr_mat4_t mv = psr_mat4_mul(model, view);
        psr_mat4_t mvp = psr_mat4_mul(psr_mat4_mul(model, view), projection);
        psr_mat3_t normal_matrix = psr_mat4_to_mat3(psr_mat4_transpose(psr_mat4_inverse(model)));

        // Make a copy of the mesh's positions to preserve the original mesh.
        psr_float3_t* positions = malloc(mesh->position_count * sizeof(psr_float3_t));
        memcpy(positions, mesh->positions, mesh->position_count * sizeof(psr_float3_t));

        int* face_cull_flags = malloc(mesh->face_count * sizeof(int));

        // Backface culling.
        for (int i = 0; i < mesh->face_count; i++) {
            // Pick any one of the triangle's points, put it in model space and check if the triangle should be culled.
            
            psr_float3_t any_point_on_triangle = psr_float3_mul_mat4(positions[mesh->faces[i].position_indices[0]], model);
            
            psr_float3_t direction_towards_point;
            PSR_SUB(direction_towards_point, any_point_on_triangle, camera_pos, 3);
            direction_towards_point = psr_normalize(direction_towards_point);
            
            psr_float3_t normal = psr_float3_mul_mat3(mesh->normals[mesh->faces[i].normal_indices[0]], normal_matrix);

            if (psr_dot(direction_towards_point, normal) >= 0) {
                face_cull_flags[i] = 1;
            }
            else {
                face_cull_flags[i] = 0;
            }
        }

        // Model/view/projection transform and viewport scaling.
        for (int i = 0; i < mesh->position_count; i++) {
            positions[i] = psr_float3_mul_mat4(positions[i], mvp);

            // Scale from [-1, 1] to viewport size.
            positions[i].x = (positions[i].x + 1) / 2.f * WIDTH;
            positions[i].y = (positions[i].y + 1) / 2.f * HEIGHT;
        }

        double render_time_start = pwa_get_elapsed_time_ms();

        // Render triangles.
        for (int i = 0; i < mesh->face_count; i++) {
            if (face_cull_flags[i]) {
                continue;
            }

            psr_float3_t tri[3];
            for (int j = 0; j < 3; j++) {
                tri[j] = positions[mesh->faces[i].position_indices[j]];
            }

            psr_float3_t light_dir = {1, 0, 0};
            
            psr_float3_t n0 = psr_float3_mul_mat3(mesh->normals[mesh->faces[i].normal_indices[0]], normal_matrix);
            psr_float3_t n1 = psr_float3_mul_mat3(mesh->normals[mesh->faces[i].normal_indices[1]], normal_matrix);
            psr_float3_t n2 = psr_float3_mul_mat3(mesh->normals[mesh->faces[i].normal_indices[2]], normal_matrix);

            psr_raster_triangle_3d(color_buffer, 
                                   depth_buffer, 
                                   tri[0], tri[1], tri[2],
                                   PSR_ATTRIB_ARRAY(PSR_ATTRIB_3(n0)),
                                   PSR_ATTRIB_ARRAY(PSR_ATTRIB_3(n1)),
                                   PSR_ATTRIB_ARRAY(PSR_ATTRIB_3(n2)),
                                   1,
                                   pixel_shader,
                                   &light_dir);
        }

        double render_time_end = pwa_get_elapsed_time_ms();

        free(face_cull_flags);
        free(positions);

        if (current_frame > (gui_delta_last + gui_delta_interval)) {
            gui_delta_last = current_frame;
            gui_delta_value = delta;
            gui_render_time = render_time_end - render_time_start;
        }

        char buf[24];
        snprintf(buf, 24, "Frame:  %.2fms", gui_delta_value);
        psr_raster_text(color_buffer, buf, (psr_int2_t){10, 10}, font, 1);

        snprintf(buf, 24, "Raster: %.2fms", gui_render_time);
        psr_raster_text(color_buffer, buf, (psr_int2_t){10, 24}, font, 1);

        pwa_window_swap_buffers(window, color_buffer);
    }

    pwa_window_destroy(window);
    printf("Done.");
}
