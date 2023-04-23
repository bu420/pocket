#include <psr.h>
#include <pwa.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define WIDTH 640
#define HEIGHT 640

#define ATLAS_W 256.f
#define ATLAS_H 256.f

const psr_float3_t cube_positions[36] = {
    {-0.5f,-0.5f,-0.5f}, {-0.5f,0.5f,-0.5f}, {0.5f,0.5f,-0.5f},   {-0.5f,-0.5f,-0.5f}, {0.5f,0.5f,-0.5f},   {0.5f,-0.5f,-0.5f},
    {0.5f,-0.5f,-0.5f},  {0.5f,0.5f,-0.5f},  {0.5f,0.5f,0.5f},    {0.5f,-0.5f,-0.5f},  {0.5f,0.5f,0.5f},    {0.5f,-0.5f,0.5f},
    {0.5f,-0.5f,0.5f},   {0.5f,0.5f,0.5f},   {-0.5f,0.5f,0.5f},   {0.5f,-0.5f,0.5f},   {-0.5f,0.5f,0.5f},   {-0.5f,-0.5f,0.5f},
    {-0.5f,-0.5f,0.5f},  {-0.5f,0.5f,0.5f},  {-0.5f,0.5f,-0.5f},  {-0.5f,-0.5f,0.5f},  {-0.5f,0.5f,-0.5f},  {-0.5f,-0.5f,-0.5f},
    {-0.5f,0.5f,-0.5f},  {-0.5f,0.5f,0.5f},  {0.5f,0.5f,0.5f},    {-0.5f,0.5f,-0.5f},  {0.5f,0.5f,0.5f},    {0.5f,0.5f,-0.5f},
    {0.5f,-0.5f,0.5f},   {-0.5f,-0.5f,0.5f}, {-0.5f,-0.5f,-0.5f}, {0.5f,-0.5f,0.5f},   {-0.5f,-0.5f,-0.5f}, {0.5f,-0.5f,-0.5f}
};

const psr_float2_t cube_tex_coords[36] = {
    {32/ATLAS_W,15/ATLAS_H}, {32/ATLAS_W,0}, {47/ATLAS_W,0}, {32/ATLAS_W,15/ATLAS_H}, {47/ATLAS_W,0}, {47/ATLAS_W,15/ATLAS_H},
    {32/ATLAS_W,15/ATLAS_H}, {32/ATLAS_W,0}, {47/ATLAS_W,0}, {32/ATLAS_W,15/ATLAS_H}, {47/ATLAS_W,0}, {47/ATLAS_W,15/ATLAS_H},
    {32/ATLAS_W,15/ATLAS_H}, {32/ATLAS_W,0}, {47/ATLAS_W,0}, {32/ATLAS_W,15/ATLAS_H}, {47/ATLAS_W,0}, {47/ATLAS_W,15/ATLAS_H},
    {32/ATLAS_W,15/ATLAS_H}, {32/ATLAS_W,0}, {47/ATLAS_W,0}, {32/ATLAS_W,15/ATLAS_H}, {47/ATLAS_W,0}, {47/ATLAS_W,15/ATLAS_H},
    {32/ATLAS_W,15/ATLAS_H}, {32/ATLAS_W,0}, {47/ATLAS_W,0}, {32/ATLAS_W,15/ATLAS_H}, {47/ATLAS_W,0}, {47/ATLAS_W,15/ATLAS_H},
    {32/ATLAS_W,15/ATLAS_H}, {32/ATLAS_W,0}, {47/ATLAS_W,0}, {32/ATLAS_W,15/ATLAS_H}, {47/ATLAS_W,0}, {47/ATLAS_W,15/ATLAS_H}
};

const psr_float3_t cube_normals[6] = {
    { 0,  0,  1},
    {-1,  0,  0},
    { 0,  0, -1},
    { 1,  0,  0},
    { 0, -1,  0},
    { 0,  1,  0}
};

typedef struct {
    psr_image_t* texture_atlas;
    psr_float3_t normal;
} shader_data_t;

psr_byte3_t pixel_shader(psr_int2_t pixel_pos, const psr_attribute_array_t* interpolated, void* user_data) {
    shader_data_t data = *(shader_data_t*)user_data;
    psr_float2_t tex_coord = PSR_ATTRIB_TO_FLOAT2(interpolated->attributes[0]);

    psr_byte_t* sample_address = psr_image_sample(data.texture_atlas, tex_coord.u, tex_coord.v);
    psr_byte3_t color = {*sample_address, *(sample_address + 1), *(sample_address + 2)};
    return color;
}

int main() {
    pwa_init();

    psr_image_t* texture_atlas = psr_image_load_bmp("assets/minecraft.bmp", PSR_R8G8B8A8);
    assert(texture_atlas);
    printf("%d, %d\n", texture_atlas->w, texture_atlas->h);

    psr_color_buffer_t* color_buffer = psr_color_buffer_create(WIDTH, HEIGHT);
    psr_depth_buffer_t* depth_buffer = psr_depth_buffer_create(WIDTH, HEIGHT);

    psr_float3_t camera_pos = {-2, -2, -2};
    psr_mat4_t projection = psr_perspective(HEIGHT / (float)WIDTH, 70 * (M_PI / 180), .1f, 1000.f);

    pwa_window_t* window = pwa_window_create("Cube", WIDTH, HEIGHT, NULL);

    while (!pwa_window_should_close(window)) {
        pwa_window_poll_events(window);

        // Model, view, projection matrix multiplication.

        psr_mat4_t view = psr_look_at(camera_pos, (psr_float3_t){0, 0, 0}, (psr_float3_t){0, 1, 0});

        psr_mat4_t model;
        psr_mat4_init_identity(&model);
        model = psr_mat4_rotate_y(model, pwa_get_elapsed_time_ms() * M_PI / 2000);
        model = psr_mat4_rotate_z(model, pwa_get_elapsed_time_ms() * M_PI / 4000);

        psr_mat4_t mvp = psr_mat4_mul(psr_mat4_mul(model, view), projection);

        // Clear buffers.

        psr_color_buffer_clear(color_buffer, (psr_byte3_t){0, 0, 0});
        psr_depth_buffer_clear(depth_buffer);

        // Raster cube (12 triangles).
        for (int i = 0; i < 12; i++) {
            psr_float3_t positions[3];
            psr_float2_t tex_coords[3];

            for (int j = 0; j < 3; j++) {
                positions[j] = cube_positions[i * 3 + j];
                tex_coords[j] = cube_tex_coords[i * 3 + j];

                positions[j] = psr_float3_mul_mat4(positions[j], mvp);

                // Scale from [-1, 1] to viewport size.
                positions[j].x = (positions[j].x + 1) / 2.f * color_buffer->w;
                positions[j].y = (positions[j].y + 1) / 2.f * color_buffer->h;
            }

            psr_float3_t normal = cube_normals[i / 2];

            shader_data_t data;
            data.texture_atlas = texture_atlas;
            data.normal = normal;

            psr_raster_triangle_3d(color_buffer, 
                                   depth_buffer, 
                                   positions[0], 
                                   positions[1], 
                                   positions[2], 
                                   PSR_ATTRIB_ARRAY(PSR_ATTRIB_2(tex_coords[0])),
                                   PSR_ATTRIB_ARRAY(PSR_ATTRIB_2(tex_coords[1])),
                                   PSR_ATTRIB_ARRAY(PSR_ATTRIB_2(tex_coords[2])),
                                   1,
                                   pixel_shader,
                                   &data);

            pwa_window_swap_buffers(window, color_buffer);
        }
    }
}
