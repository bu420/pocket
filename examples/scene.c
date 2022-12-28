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

    // Gradient background.
    for (int x = 0; x < WIDTH; ++x) {
        for (int y = 0; y < HEIGHT; ++y) {
            srz_byte3_t* pixel = srz_color_buffer_at(&color_buffer, x, y);

            pixel->r = x * 128 / WIDTH;
            pixel->g = 0;
            pixel->b = (x + y) * 128 / (WIDTH + HEIGHT);
        }
    }

    /*srz_byte_t cube[36 * 3] = {
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

    float dist = 4;
    float ang_x = srz_cosf(60 * (SRZ_PI / 180)) * dist;
    float ang_y = srz_sinf(60 * (SRZ_PI / 180)) * dist;
    srz_float3_t cube_positions[6] = {
        {dist, 0, 0},
        {ang_x, ang_y, 0},
        {-ang_x, ang_y, 0},
        {-dist, 0, 0},
        {-ang_x, -ang_y, 0},
        {ang_x, -ang_y, 0}
    };*/

    srz_float3_t camera_pos = {0, 0, -7};
    srz_matrix_t view = srz_look_at(camera_pos, (srz_float3_t){0, 0, 0}, (srz_float3_t){0, -1, 0});
    srz_matrix_t projection = srz_perspective(HEIGHT / (float)WIDTH, 75.f * (SRZ_PI / 180), 0.1f, 1000.f);

    //for (int cube_i = 0; cube_i < 6; ++cube_i) {
        srz_matrix_t model;
        srz_matrix_init_identity(&model);
        model = srz_matrix_translate(model, (srz_float3_t){-.5, -.5, -.5});
        //model = srz_matrix_translate(model, cube_positions[cube_i]);

        srz_matrix_t mvp = srz_matrix_mul(srz_matrix_mul(model, view), projection);

        //float vertices[36 * 3];

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

            // Discard invisible triangles.
            /*if (normal.x * (v0.x - camera_pos.x) + normal.y * (v0.y - camera_pos.y) + normal.z * (v0.z - camera_pos.z) < 0) {
                continue;
            }*/

            srz_int2_t tri[3];
            for (int j = 0; j < 3; ++j) {
                srz_float3_t pos = mesh.positions[mesh.faces[i].position_indices[j]];
                tri[j].x = srz_roundf(pos.x);
                tri[j].y = srz_roundf(pos.y);
            }

            srz_raster_triangle(&color_buffer, SRZ_NULL, tri[0], tri[1], tri[2], color);
        }
    //}

    srz_save_bmp("scene.bmp", color_buffer);
    free(color_buffer.data);

    printf("Done.");
}
