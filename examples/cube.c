#include <psr.h>
#include <pwa.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define WIDTH 640
#define HEIGHT 640

const char cube_vertex_positions[108] = {
    0,0,0,0,1,0,1,1,0,0,0,0,1,1,0,1,0,0,
    1,0,0,1,1,0,1,1,1,1,0,0,1,1,1,1,0,1,
    1,0,1,1,1,1,0,1,1,1,0,1,0,1,1,0,0,1,
    0,0,1,0,1,1,0,1,0,0,0,1,0,1,0,0,0,0,
    0,1,0,0,1,1,1,1,1,0,1,0,1,1,1,1,1,0,
    1,0,1,0,0,1,0,0,0,1,0,1,0,0,0,1,0,0
};

const char cube_vertex_normals[18] = {
    0, 0, 1,
    -1, 0, 0,
    0, 0, -1,
    1, 0, 0,
    0, -1, 0,
    0, 1, 0
};

typedef struct {
    psr_float3_t normal;
} shader_data_t;

psr_byte3_t pixel_shader(psr_int2_t pixel_pos, const psr_attribute_array_t* interpolated, void* user_data) {
    shader_data_t data = *(shader_data_t*)user_data;
    
    return (psr_byte3_t){
        (data.normal.x + 1) / 2 * 255, 
        (data.normal.y + 1) / 2 * 255, 
        (data.normal.z + 1) / 2 * 255
    };
}

int main(int argc, char** argv) {
    pwa_init();

    psr_color_buffer_t* color_buffer = psr_color_buffer_create(WIDTH, HEIGHT);
    psr_depth_buffer_t* depth_buffer = psr_depth_buffer_create(WIDTH, HEIGHT);

    psr_float3_t camera_pos = {-1, -1, -1};
    psr_mat4_t projection = psr_perspective(HEIGHT / (float)WIDTH, 70 * (M_PI / 180), .1f, 1000.f);

    pwa_window_t* window = pwa_window_create("Cube", WIDTH, HEIGHT, NULL);

    while (!pwa_window_should_close(window)) {
        pwa_window_poll_events(window);

        // Model, view, projection matrix multiplication.

        psr_mat4_t view = psr_look_at(camera_pos, (psr_float3_t){0, 0, 0}, (psr_float3_t){0, 1, 0});

        psr_mat4_t model;
        psr_mat4_init_identity(&model);

        psr_mat4_t mvp = psr_mat4_mul(psr_mat4_mul(model, view), projection);

        // Clear buffers.

        psr_color_buffer_clear(color_buffer, (psr_byte3_t){0, 0, 0});
        psr_depth_buffer_clear(depth_buffer);

        // Raster cube (12 triangles).
        for (int i = 0; i < 12; i++) {
            // Multiply cube vertex position with MVP matrix and scale to viewport.
            psr_float3_t tri[3];
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    tri[j].values[k] = cube_vertex_positions[i * 9 + j * 3 + k];
                }

                tri[j] = psr_float3_mul_mat4(tri[j], mvp);

                // Scale from [-1, 1] to viewport size.
                tri[j].x = (tri[j].x + 1) / 2.f * color_buffer->w;
                tri[j].y = (tri[j].y + 1) / 2.f * color_buffer->h;
            }

            psr_float3_t normal;
            for (int j = 0; j < 3; j++) {
                normal.values[j] = cube_vertex_normals[(i / 2) * 3 + j];
            }

            // Data to pass to pixel shader.
            shader_data_t data;
            data.normal = normal;

            // Raster triangle, no attributes, normal vector is passed directly since
            // it's the same for the whole triangle and doesn't need to be interpolated.
            psr_raster_triangle_3d(color_buffer, 
                                   depth_buffer, 
                                   tri[0], tri[1], tri[2], 
                                   NULL,
                                   NULL,
                                   NULL,
                                   0,
                                   pixel_shader,
                                   &data);

            pwa_window_swap_buffers(window, color_buffer);
        }
    }
}
