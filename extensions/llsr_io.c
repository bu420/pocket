#include "llsr_io.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void llsr_save_bmp(char const* filename, llsr_color_buffer_t color_buffer) {
    llsr_byte_t padding[] = {0, 0, 0};
    int padding_size = (4 - color_buffer.w * 3 % 4) % 4;
    int stride = color_buffer.w * 3 + padding_size;
    int file_size = 54 + stride * color_buffer.h;

    llsr_byte_t header[54] = {
        'B', 'M',
        (llsr_byte_t)file_size, (llsr_byte_t)(file_size >> 8), (llsr_byte_t)(file_size >> 16), (llsr_byte_t)(file_size >> 24),
        0, 0, 0, 0,
        54, 0, 0, 0,
        40, 0, 0, 0,
        (llsr_byte_t)color_buffer.w, (llsr_byte_t)(color_buffer.w >> 8), (llsr_byte_t)(color_buffer.w >> 16), (llsr_byte_t)(color_buffer.w >> 24),
        (llsr_byte_t)color_buffer.h, (llsr_byte_t)(color_buffer.h >> 8), (llsr_byte_t)(color_buffer.h >> 16), (llsr_byte_t)(color_buffer.h >> 24),
        1, 0,
        3 * 8, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0
    };

    FILE* file = fopen(filename, "wb");

    fwrite(header, 1, 54, file);

    for (int y = 0; y < color_buffer.h; ++y) {
        for (int x = 0; x < color_buffer.w; ++x) {
            int i = (color_buffer.h - y) * color_buffer.h + x;
            llsr_byte3_t pixel = color_buffer.data[i];
            llsr_byte3_t color = {.r = pixel.b, .g = pixel.g, .b = pixel.r};
            fwrite(&color, 3, 1, file);
        }

        fwrite(padding, 1, padding_size, file);
    }

    fclose(file);
}

int _llsr_split(char* str, char const* delims, char** tokens, int max) {
    int count = 0;
    char* token = strtok(str, delims);
    for (; count < max && token; ++count) {
        tokens[count] = token;
        token = strtok(NULL, delims);
    }
    return count;
}

int llsr_load_obj(char const* filename, llsr_mesh_t* mesh) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        return 0;
    }

    mesh->face_count = 0;
    mesh->position_count = 0;
    mesh->tex_coord_count = 0;
    mesh->normal_count = 0;

    int face_max = 128;
    int position_max = 128;
    int tex_coord_max = 128;
    int normal_max = 128;

    mesh->faces = malloc(face_max * sizeof(llsr_face_t));
    mesh->positions = malloc(position_max * sizeof(llsr_float3_t));
    mesh->tex_coords = malloc(tex_coord_max * sizeof(llsr_float2_t));
    mesh->normals = malloc(normal_max * sizeof(llsr_float3_t));

    char line[256];
    while (fgets(line, 256, file)) {
        char* tokens[5];
        int token_count = _llsr_split(line, " ", tokens, 5);

        if (token_count == 0) {
            continue;
        }

        // Vertex position.
        if (strcmp(tokens[0], "v") == 0) {
            if (mesh->position_count == position_max) {
                mesh->positions = realloc(mesh->positions, (position_max *= 2) * sizeof(llsr_float3_t));
            }
            mesh->positions[mesh->position_count++] = (llsr_float3_t){atof(tokens[1]), atof(tokens[2]), atof(tokens[3])};
        }
        // Vertex texture coordinate.
        else if (strcmp(tokens[0], "vt") == 0) {
            if (mesh->tex_coord_count == tex_coord_max) {
                mesh->tex_coords = realloc(mesh->tex_coords, (tex_coord_max *= 2) * sizeof(llsr_float2_t));
            }
            mesh->tex_coords[mesh->tex_coord_count++] = (llsr_float2_t){atof(tokens[1]), atof(tokens[2])};
        }
        // Vertex normal.
        else if (strcmp(tokens[0], "vn") == 0) {
            if (mesh->normal_count == normal_max) {
                mesh->normals = realloc(mesh->normals, (normal_max *= 2) * sizeof(llsr_float3_t));
            }
            mesh->normals[mesh->normal_count++] = (llsr_float3_t){atof(tokens[1]), atof(tokens[2]), atof(tokens[3])};
        }
        // Face.
        else if (strcmp(tokens[0], "f") == 0) {
            char* face_tokens[4][3];
            int face_token_count = token_count - 1;
            for (int i = 0; i < face_token_count; ++i) {
                // Split into indices.
                if (_llsr_split(tokens[i + 1], "/", face_tokens[i], 3) != 3) {
                    return -1;
                }
            }

            int position_indices[4];
            int tex_coord_indices[4];
            int normal_indices[4];

            for (int i = 0; i < face_token_count; ++i) {
                position_indices[i] = atoi(face_tokens[i][0]) - 1;
                tex_coord_indices[i] = atoi(face_tokens[i][1]) - 1;
                normal_indices[i] = atoi(face_tokens[i][2]) - 1;
            }

            if (face_token_count == 3) {
                llsr_face_t face;
                for (int i = 0; i < 3; ++i) {
                    face.position_indices[i] = position_indices[i];
                    face.tex_coord_indices[i] = tex_coord_indices[i];
                    face.normal_indices[i] = normal_indices[i];
                }

                if (mesh->face_count == face_max) {
                    mesh->faces = realloc(mesh->faces, (face_max *= 2) * sizeof(llsr_face_t));
                }
                mesh->faces[mesh->face_count++] = face;
            }
            else if (face_token_count == 4) {
                llsr_face_t a, b;
                for (int i = 0; i < 3; ++i) {
                    a.position_indices[i] = position_indices[i];
                    a.tex_coord_indices[i] = tex_coord_indices[i];
                    a.normal_indices[i] = normal_indices[i];

                    b.position_indices[i] = position_indices[(i + 2) % 4];
                    b.tex_coord_indices[i] = tex_coord_indices[(i + 2) % 4];
                    b.normal_indices[i] = normal_indices[(i + 2) % 4];
                }

                if (mesh->face_count >= face_max - 1) {
                    mesh->faces = realloc(mesh->faces, (face_max *= 2) * sizeof(llsr_face_t));
                }
                mesh->faces[mesh->face_count++] = a;
                mesh->faces[mesh->face_count++] = b;
            }
            else {
                return -1;
            }
        }
    }
    
    fclose(file);
    return 1;
}

void llsr_mesh_free(llsr_mesh_t* mesh) {
    free(mesh->positions);
    free(mesh->tex_coords);
    free(mesh->normals);
    free(mesh->faces);
}
