#ifndef BOOTLEG_MINECRAFT_H
#define BOOTLEG_MINECRAFT_H

#include <srz.h>

typedef enum {
    BOOTLEG_MINECRAFT_ACTION_DIG            = 1 << 0,
    BOOTLEG_MINECRAFT_ACTION_PLACE          = 1 << 1,
    BOOTLEG_MINECRAFT_ACTION_MOVE_FORWARD   = 1 << 2,
    BOOTLEG_MINECRAFT_ACTION_MOVE_BACKWARD  = 1 << 3,
    BOOTLEG_MINECRAFT_ACTION_MOVE_LEFT      = 1 << 4,
    BOOTLEG_MINECRAFT_ACTION_MOVE_RIGHT     = 1 << 5,
    BOOTLEG_MINECRAFT_ACTION_MOVE_UP        = 1 << 6,
    BOOTLEG_MINECRAFT_ACTION_MOVE_DOWN      = 1 << 7,
    BOOTLEG_MINECRAFT_ACTION_INCREASE_YAW   = 1 << 8,
    BOOTLEG_MINECRAFT_ACTION_DECREASE_YAW   = 1 << 9,
    BOOTLEG_MINECRAFT_ACTION_INCREASE_PITCH = 1 << 10,
    BOOTLEG_MINECRAFT_ACTION_DECREASE_PITCH = 1 << 11,
} bootleg_minecraft_action_t;

typedef enum {
    BOOTLEG_MINECRAFT_BLOCK_NONE,
    BOOTLEG_MINECRAFT_BLOCK_STONE,
    BOOTLEG_MINECRAFT_BLOCK_DIRT,
    BOOTLEG_MINECRAFT_BLOCK_GRASS,
    BOOTLEG_MINECRAFT_BLOCK_LOG,
    BOOTLEG_MINECRAFT_BLOCK_LEAVES,
} bootleg_minecraft_block_t;

typedef struct {
    srz_color_buffer_t color_buffer;
    srz_depth_buffer_t depth_buffer;

    // List of blocks, must be able to contain 4096 (16x16x16) blocks.
    bootleg_minecraft_block_t* level;

    srz_matrix_t projection;

    srz_float3_t pos;
} bootleg_minecraft_t;

void bootleg_minecraft_init(bootleg_minecraft_t* bootleg_minecraft, int w, int h, srz_byte3_t* color_buffer_memory, float* depth_buffer_memory, bootleg_minecraft_block_t* level_4096_blocks_memory);
void bootleg_minecraft_update(bootleg_minecraft_t* bootleg_minecraft, int elapsed_time_ms, int actions);

#ifdef BOOTLEG_MINECRAFT_IMPLEMENTATION
int lvl_idx(int x, int y, int z) {
    return z * 256 + y * 16 + x;
}

void bootleg_minecraft_init(bootleg_minecraft_t* bootleg_minecraft, int w, int h, srz_byte3_t* color_buffer_memory, float* depth_buffer_memory, bootleg_minecraft_block_t* level_4096_blocks_memory) {
    bootleg_minecraft_t* mc = bootleg_minecraft;
    
    mc->color_buffer.w = w;
    mc->color_buffer.h = h;
    mc->color_buffer.data = color_buffer_memory;

    mc->depth_buffer.w = w;
    mc->depth_buffer.h = h;
    mc->depth_buffer.data = depth_buffer_memory;

    mc->level = level_4096_blocks_memory;
    for (int i = 0; i < 4096; i++) {
        mc->level[i] = BOOTLEG_MINECRAFT_BLOCK_STONE;
    }

    mc->pos = (srz_float3_t){-32, -24, -32};

    mc->projection = srz_perspective(h / (float)w, 70 * (SRZ_PI / 180), .1f, 1000.f);
}

void bootleg_minecraft_update(bootleg_minecraft_t* bootleg_minecraft, int elapsed_time_ms, int actions) {
    bootleg_minecraft_t* mc = bootleg_minecraft;

    static srz_byte_t cube_vertex_positions[108] = {
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

    srz_matrix_t view = srz_look_at(mc->pos, (srz_float3_t){0, 0, 0}, (srz_float3_t){0, 1, 0});

    srz_byte3_t sky_color = {54, 199, 242};
    srz_color_buffer_clear(&mc->color_buffer, sky_color);
    srz_depth_buffer_clear(&mc->depth_buffer);

    // Render each block.
    for (int x = 0; x < 16; x++) {
        for (int y = 0; y < 16; y++) {
            for (int z = 0; z < 16; z++) {
                if (mc->level[lvl_idx(x, y, z)] == BOOTLEG_MINECRAFT_BLOCK_NONE) {
                    continue;
                }

                srz_matrix_t model;
                srz_matrix_init_identity(&model);
                model = srz_matrix_translate(model, (srz_float3_t){x, y, z});

                srz_matrix_t mvp = srz_matrix_mul(srz_matrix_mul(model, view), mc->projection);

                // 12 triangles in a cube.
                for (int i = 0; i < 12; i++) {
                    srz_float3_t tri[3];
                    for (int j = 0; j < 3; j++) {
                        for (int k = 0; k < 3; k++) {
                            tri[j].values[k] = cube_vertex_positions[i * 9 + j * 3 + k];
                        }

                        tri[j] = srz_float3_mul_matrix(tri[j], mvp);

                        // Scale from [-1, 1] to viewport size.
                        tri[j].x = (tri[j].x + 1) / 2.f * mc->color_buffer.w;
                        tri[j].y = (tri[j].y + 1) / 2.f * mc->color_buffer.h;
                    }

                    srz_float3_t norm;
                    for (int j = 0; j < 3; j++) {
                        norm.values[j] = cube_vertex_normals[(i / 2) * 3 + j];
                    }

                    srz_float3_t light_dir = srz_normalize((srz_float3_t){0, 1, .5});
                    float light = srz_dot(norm, light_dir);
                    srz_byte3_t color = {255 * light, 255 * light, 255 * light};

                    srz_raster_triangle(&mc->color_buffer, &mc->depth_buffer, tri[0], tri[1], tri[2], color);
                }
            }
        }
    }
}
#endif
#endif
