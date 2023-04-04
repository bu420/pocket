#include <psr.h>
#include <pwa.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define WIDTH 640
#define HEIGHT 480

static char cube_vertex_positions[108] = {
    0,0,0,0,1,0,1,1,0,0,0,0,1,1,0,1,0,0,
    1,0,0,1,1,0,1,1,1,1,0,0,1,1,1,1,0,1,
    1,0,1,1,1,1,0,1,1,1,0,1,0,1,1,0,0,1,
    0,0,1,0,1,1,0,1,0,0,0,1,0,1,0,0,0,0,
    0,1,0,0,1,1,1,1,1,0,1,0,1,1,1,1,1,0,
    1,0,1,0,0,1,0,0,0,1,0,1,0,0,0,1,0,0
};

static char cube_vertex_normals[18] = {
    0, 0, 1,
    -1, 0, 0,
    0, 0, -1,
    1, 0, 0,
    0, -1, 0,
    0, 1, 0
};

int lvl_idx(int x, int y, int z) {
    return z * 256 + y * 16 + x;
}

pwa_pixel_buffer_t on_draw(void* user_data) {
    return *(pwa_pixel_buffer_t*)user_data;
}

int main(int argc, char** argv) {
    pwa_init();

    psr_color_buffer_t* color_buffer = psr_color_buffer_create(WIDTH, HEIGHT);
    psr_depth_buffer_t* depth_buffer = psr_depth_buffer_create(WIDTH, HEIGHT);

    // List of blocks (16x16x16).
    char* level = malloc(4096 * sizeof(char));
    for (int i = 0; i < 4096; i++) {
        level[i] = 1;
    }

    psr_float3_t pos = {-32, -24, -32};
    psr_mat4_t projection = psr_perspective(HEIGHT / (float)WIDTH, 70 * (M_PI / 180), .1f, 1000.f);

    pwa_pixel_buffer_t buffer;
    buffer.pixels = malloc(WIDTH * HEIGHT * sizeof(uint32_t));
    buffer.w = WIDTH;
    buffer.h = HEIGHT;

    pwa_window_t* window = pwa_window_create("Minecraft", WIDTH, HEIGHT, &buffer);
    pwa_window_set_draw_callback(window, on_draw);

    while (!pwa_window_should_close(window)) {
        pwa_window_poll_events(window);

        psr_mat4_t view = psr_look_at(pos, (psr_float3_t){0, 0, 0}, (psr_float3_t){0, 1, 0});

        psr_byte3_t sky_color = {54, 199, 242};
        psr_color_buffer_clear(color_buffer, sky_color);
        psr_depth_buffer_clear(depth_buffer);

        // Render each block.
        for (int x = 0; x < 16; x++) {
            for (int y = 0; y < 16; y++) {
                for (int z = 0; z < 16; z++) {
                    if (level[lvl_idx(x, y, z)] == 0) {
                        continue;
                    }

                    psr_mat4_t model;
                    psr_mat4_init_identity(&model);
                    model = psr_mat4_translate(model, (psr_float3_t){x, y, z});

                    psr_mat4_t mvp = psr_mat4_mul(psr_mat4_mul(model, view), projection);

                    // 12 triangles in a cube.
                    for (int i = 0; i < 12; i++) {
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

                        psr_float3_t norm;
                        for (int j = 0; j < 3; j++) {
                            norm.values[j] = cube_vertex_normals[(i / 2) * 3 + j];
                        }

                        psr_float3_t light_dir = psr_normalize((psr_float3_t){0, 1, .5});
                        float light = psr_dot(norm, light_dir);
                        psr_byte3_t color = {255 * light, 255 * light, 255 * light};


                        psr_raster_triangle_2d_color(color_buffer, (psr_int2_t){(int)tri[0].x, (int)tri[0].y}, (psr_int2_t){(int)tri[1].x, (int)tri[1].y}, (psr_int2_t){(int)tri[2].x, (int)tri[2].y}, color);
                    }
                }
            }
        }

        // Copy color buffer to window buffer.
        for (int x = 0; x < WIDTH; x++) {
            for (int y = 0; y < HEIGHT; y++) {
                psr_byte3_t color = *psr_color_buffer_at(color_buffer, x, y);
                uint32_t* pixel = &buffer.pixels[y * WIDTH + x];
                
                *pixel = (*pixel & 0xffff00) | (color.b << 0);
                *pixel = (*pixel & 0xff00ff) | (color.g << 8);
                *pixel = (*pixel & 0x00ffff) | (color.r << 16);
            }
        }
    }
}
