#include "psr.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

psr_float3_t psr_normalize(psr_float3_t f) {
    float len = sqrtf(f.x * f.x + f.y * f.y + f.z * f.z);
    psr_float3_t result = {f.x / len, f.y / len, f.z / len};
    return result;
}

psr_float3_t psr_cross(psr_float3_t a, psr_float3_t b) {
    psr_float3_t result = {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
    return result;
}

float psr_dot(psr_float3_t a, psr_float3_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

void psr_swap(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

void psr_swapf(float* a, float* b) {
    float temp = *a;
    *a = *b;
    *b = temp;
}

void psr_int2_swap(psr_int2_t* a, psr_int2_t* b) {
    psr_swap(&a->x, &b->x);
    psr_swap(&a->y, &b->y);
}

void psr_float2_swap(psr_float2_t* a, psr_float2_t* b) {
    psr_swapf(&a->x, &b->x);
    psr_swapf(&a->y, &b->y);
}

void psr_float3_swap(psr_float3_t* a, psr_float3_t* b) {
    psr_swapf(&a->x, &b->x);
    psr_swapf(&a->y, &b->y);
    psr_swapf(&a->z, &b->z);
}

psr_byte3_t psr_byte3_lerp(psr_byte3_t a, psr_byte3_t b, float amount) {
    psr_byte3_t result;
    for (int i = 0; i < 3; i++) {
        result.values[i] = a.values[i] * (1 - amount) + b.values[i] * amount;
    }
    return result;
}

void psr_matrix_init_zero(psr_matrix_t* matrix) {
    for (int x = 0; x < 4; x++) {
        for (int y = 0; y < 4; y++) {
            matrix->m[x][y] = 0;
        }
    }
}

void psr_matrix_init_identity(psr_matrix_t* matrix) {
    psr_matrix_init_zero(matrix);
    for (int i = 0; i < 4; i++) {
        matrix->m[i][i] = 1;
    }
}

psr_matrix_t psr_matrix_mul(psr_matrix_t a, psr_matrix_t b) {
    psr_matrix_t result;
    psr_matrix_init_zero(&result);

    for (int x = 0; x < 4; x++) {
        for (int y = 0; y < 4; y++) {
            for (int i = 0; i < 4; i++) {
                result.m[x][y] += a.m[x][i] * b.m[i][y];
            }
        }
    }
    return result;
}

psr_matrix_t psr_matrix_translate(psr_matrix_t matrix, psr_float3_t f) {
    psr_matrix_t m;
    psr_matrix_init_identity(&m);

    m.m[3][0] = f.x;
    m.m[3][1] = f.y;
    m.m[3][2] = f.z;

    return psr_matrix_mul(matrix, m);
}

psr_matrix_t psr_matrix_rotate_x(psr_matrix_t matrix, float a) {
    float s = sinf(a);
    float c = cosf(a);

    psr_matrix_t x;
    psr_matrix_init_identity(&x);

    x.m[1][1] = c;
    x.m[1][2] = -s;
    x.m[2][1] = s;
    x.m[2][2] = c;

    return psr_matrix_mul(matrix, x);
}

psr_matrix_t psr_matrix_rotate_y(psr_matrix_t matrix, float a) {
    float s = sinf(a);
    float c = cosf(a);

    psr_matrix_t y;
    psr_matrix_init_identity(&y);

    y.m[0][0] = c;
    y.m[0][2] = s;
    y.m[2][0] = -s;
    y.m[2][2] = c;

    return psr_matrix_mul(matrix, y);
}

psr_matrix_t psr_matrix_rotate_z(psr_matrix_t matrix, float a) {
    float s = sinf(a);
    float c = cosf(a);

    psr_matrix_t z;
    psr_matrix_init_identity(&z);

    z.m[0][0] = c;
    z.m[0][1] = -s;
    z.m[1][0] = s;
    z.m[1][1] = c;

    return psr_matrix_mul(matrix, z);
}

psr_matrix_t psr_matrix_scale(psr_matrix_t matrix, psr_float3_t f) {
    psr_matrix_t m;
    psr_matrix_init_identity(&m);

    m.m[0][0] = f.x;
    m.m[1][1] = f.y;
    m.m[2][2] = f.z;

    return psr_matrix_mul(matrix, m);
}

psr_float3_t psr_matrix_mul_float3(psr_matrix_t matrix, psr_float3_t f) {
    psr_float4_t f4;
    for (int x = 0; x < 4; x++) {
        float i = 0;
        for (int y = 0; y < 3; y++) {
            i += matrix.m[x][y] * f.values[y];
        }
        f4.values[x] = i + matrix.m[x][3];
    }
    
    if (f4.w != 0) {
        f4.x /= f4.w;
        f4.y /= f4.w;
        f4.z /= f4.w;
    }

    psr_float3_t result = {f4.x, f4.y, f4.z};
    return result;
}

psr_float4_t psr_matrix_mul_float4(psr_matrix_t matrix, psr_float4_t f) {
    psr_float4_t result;
    for (int x = 0; x < 4; x++) {
        float i = 0;
        for (int y = 0; y < 4; y++) {
            i += matrix.m[x][y] * f.values[y];
        }
        result.values[x] = i;
    }
    return result;
}

psr_float3_t psr_float3_mul_matrix(psr_float3_t f, psr_matrix_t matrix) {
    psr_float4_t f4;
    for (int y = 0; y < 4; y++) {
        float i = 0;
        for (int x = 0; x < 3; x++) {
            i += matrix.m[x][y] * f.values[x];
        }
        f4.values[y] = i + matrix.m[3][y];
    }
    
    if (f4.w != 0) {
        f4.x /= f4.w;
        f4.y /= f4.w;
        f4.z /= f4.w;
    }

    psr_float3_t result = {f4.x, f4.y, f4.z};
    return result;
}

psr_float4_t psr_float4_mul_matrix(psr_float4_t f, psr_matrix_t matrix) {
    psr_float4_t result;
    for (int y = 0; y < 4; y++) {
        float i = 0;
        for (int x = 0; x < 4; x++) {
            i += matrix.m[x][y] * f.values[x];
        }
        result.values[y] = i;
    }
    return result;
}

psr_matrix_t psr_look_at(psr_float3_t pos, psr_float3_t target, psr_float3_t up) {
    psr_float3_t diff = {target.x - pos.x, target.y - pos.y, target.z - pos.z};
    psr_float3_t forward = psr_normalize(diff);
    psr_float3_t right = psr_normalize(psr_cross(forward, up));
    psr_float3_t local_up = psr_normalize(psr_cross(right, forward));

    psr_matrix_t m;
    psr_matrix_init_identity(&m);
    m.m[0][0] = right.x;
    m.m[1][0] = right.y;
    m.m[2][0] = right.z;
    m.m[0][1] = local_up.x;
    m.m[1][1] = local_up.y;
    m.m[2][1] = local_up.z;
    m.m[0][2] = -forward.x;
    m.m[1][2] = -forward.y;
    m.m[2][2] = -forward.z;
    m.m[3][0] = -psr_dot(right, pos);
    m.m[3][1] = -psr_dot(local_up, pos);
    m.m[3][2] = psr_dot(forward, pos);
    return m;
}

psr_matrix_t psr_perspective(float aspect, float fov, float near, float far) {    
    psr_matrix_t m;
    psr_matrix_init_zero(&m);

    float half_tan = tanf(fov / 2);

    m.m[0][0] = 1 / (half_tan * aspect);
    m.m[1][1] = 1 / half_tan;
    m.m[2][2] = -(far + near) / (far - near);
    m.m[2][3] = -1;
    m.m[3][2] = -(2 * far * near) / (far - near);
    return m;
}

psr_byte3_t* psr_color_buffer_at(psr_color_buffer_t* color_buffer, int x, int y) {
    return &color_buffer->data[y * color_buffer->w + x];
}

void psr_color_buffer_clear(psr_color_buffer_t* color_buffer, psr_byte3_t color) {
    for (int i = 0; i < color_buffer->w * color_buffer->h; i++) {
        color_buffer->data[i] = color;
    }
}

float* psr_depth_buffer_at(psr_depth_buffer_t* depth_buffer, int x, int y) {
    return &depth_buffer->data[y * depth_buffer->w + x];
}

void psr_depth_buffer_clear(psr_depth_buffer_t* depth_buffer) {
    for (int i = 0; i < depth_buffer->w * depth_buffer->h; i++) {
        depth_buffer->data[i] = 0;
    }
}

psr_byte4_t* psr_texture_at(psr_texture_t* texture, int x, int y) {
    return &texture->data[y * texture->w + x];
}

psr_byte4_t psr_texture_sample(psr_texture_t texture, float u, float v) {
    return *psr_texture_at(&texture, (int)(u * texture.w), (int)(v * texture.h - 1));
}

psr_bresenham_line_t psr_bresenham_line_create(psr_int2_t start, psr_int2_t end) {
    psr_bresenham_line_t line;
    line.start = start;
    line.end = end;
    line.current = start;
    line.delta = (psr_int2_t){abs(end.x - start.x), -abs(end.y - start.y)};
    line.dir = (psr_int2_t){start.x < end.x ? 1 : -1, start.y < end.y ? 1 : -1};
    line.deviation = line.delta.x + line.delta.y;
    return line;
}

int psr_bresenham_line_step(psr_bresenham_line_t* line) {
    if (line->current.x == line->end.x && line->current.y == line->end.y) {
        return 0;
    }
    if (2 * line->deviation >= line->delta.y) {
        if (line->current.x == line->end.x) {
            return 0;
        }

        line->deviation += line->delta.y;
        line->current.x += line->dir.x;
    }
    if (2 * line->deviation <= line->delta.x) {
        if (line->current.y == line->end.y) {
            return 0;
        }

        line->deviation += line->delta.x;
        line->current.y += line->dir.y;
    }
    return 1;
}

// TODO: find more accurate way.
float _psr_bresenham_line_inverse_lerp(psr_int2_t start, psr_int2_t end, psr_int2_t current) {
    if (end.x - start.x > end.y - start.y) {
        return (current.x - start.x) / (float)(end.x - start.x);
    }
    return (current.y - start.y) / (float)(end.y - start.y);
}

void psr_raster_line(psr_color_buffer_t* color_buffer, psr_int2_t start, psr_int2_t end, psr_byte3_t start_color, psr_byte3_t end_color) {
    psr_bresenham_line_t line = psr_bresenham_line_create(start, end);
    do {
        psr_byte3_t current_color = psr_byte3_lerp(start_color, end_color, _psr_bresenham_line_inverse_lerp(start, end, line.current));
        *psr_color_buffer_at(color_buffer, line.current.x, line.current.y) = current_color;
    } 
    while (psr_bresenham_line_step(&line));
}

void _psr_sort_triangle_vertices_by_height(psr_int2_t* pos0, psr_int2_t* pos1, psr_int2_t* pos2) {
    if (pos0->y > pos1->y) {
        psr_int2_swap(pos0, pos1);
    }
    if (pos0->y > pos2->y) {
        psr_int2_swap(pos0, pos2);
    }
    if (pos1->y > pos2->y) {
        psr_int2_swap(pos1, pos2);
    }
}

int _psr_bresenham_line_step_until_vertical_difference(psr_bresenham_line_t* line) {
    int y = line->current.y;

    while (psr_bresenham_line_step(line)) {
        if (line->current.y != y) {
            return 1;
        }
        y = line->current.y;
    }
    return 0;
}

// Either the top or bottom of the triangle must be flat (the lines must share the same start and end y position).
void _psr_raster_triangle_2d_flat(psr_color_buffer_t* color_buffer, psr_bresenham_line_t line_a, psr_bresenham_line_t line_b, psr_byte3_t color) {
    for (int y = line_a.start.y; y <= line_a.end.y; y++) {
        int start = line_a.current.x;
        int end = line_b.current.x;

        if (start > end) {
            psr_swap(&start, &end);
        }

        for (int x = start; x <= end; x++) {
            // TODO: find way to check bounds.
            if (x >= 0 && x < color_buffer->w && y >= 0 && y < color_buffer->h) {
                *psr_color_buffer_at(color_buffer, x, y) = color;
            }
        }

        _psr_bresenham_line_step_until_vertical_difference(&line_a);
        _psr_bresenham_line_step_until_vertical_difference(&line_b);
    }
}

void psr_raster_triangle_2d_color(psr_color_buffer_t* color_buffer, psr_int2_t pos0, psr_int2_t pos1, psr_int2_t pos2, psr_byte3_t color) {
    _psr_sort_triangle_vertices_by_height(&pos0, &pos1, &pos2);

    // Check if the top of the triangle is flat.
    if (pos0.y == pos1.y) {
        _psr_raster_triangle_2d_flat(color_buffer, psr_bresenham_line_create(pos0, pos2), psr_bresenham_line_create(pos1, pos2), color);
    }
    // Check if the bottom is flat.
    else if (pos1.y == pos2.y) {
        _psr_raster_triangle_2d_flat(color_buffer, psr_bresenham_line_create(pos0, pos1), psr_bresenham_line_create(pos0, pos2), color);
    }
    // Else plit triangle into two smaller triangles.
    else {
        psr_int2_t pos3 = {(int)(pos0.x + ((float)(pos1.y - pos0.y) / (float)(pos2.y - pos0.y)) * (float)(pos2.x - pos0.x)), pos1.y};
        // Top (flat bottom).
        _psr_raster_triangle_2d_flat(color_buffer, psr_bresenham_line_create(pos0, pos1), psr_bresenham_line_create(pos0, pos3), color);
        // Bottom (flat top).
        _psr_raster_triangle_2d_flat(color_buffer, psr_bresenham_line_create(pos3, pos2), psr_bresenham_line_create(pos1, pos2), color);
    }
}

void psr_raster_triangle_2d_callback(psr_color_buffer_t* color_buffer, psr_int2_t pos0, psr_int2_t pos1, psr_int2_t pos2, void (*callback)(psr_int2_t pixel_pos), void* user_data) {
    
}

void psr_raster_triangle_3d(psr_color_buffer_t* color_buffer, psr_depth_buffer_t* depth_buffer, psr_float3_t pos0, psr_float3_t pos1, psr_float3_t pos2, psr_byte3_t color) {
    
}

void psr_save_bmp(char const* filename, psr_color_buffer_t color_buffer) {
    psr_byte_t padding[] = {0, 0, 0};
    int padding_size = (4 - color_buffer.w * 3 % 4) % 4;
    int stride = color_buffer.w * 3 + padding_size;
    int file_size = 54 + stride * color_buffer.h;

    psr_byte_t header[54] = {
        'B', 'M',
        (psr_byte_t)file_size, (psr_byte_t)(file_size >> 8), (psr_byte_t)(file_size >> 16), (psr_byte_t)(file_size >> 24),
        0, 0, 0, 0,
        54, 0, 0, 0,
        40, 0, 0, 0,
        (psr_byte_t)color_buffer.w, (psr_byte_t)(color_buffer.w >> 8), (psr_byte_t)(color_buffer.w >> 16), (psr_byte_t)(color_buffer.w >> 24),
        (psr_byte_t)color_buffer.h, (psr_byte_t)(color_buffer.h >> 8), (psr_byte_t)(color_buffer.h >> 16), (psr_byte_t)(color_buffer.h >> 24),
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
            psr_byte3_t pixel = color_buffer.data[i];
            psr_byte3_t color = {.r = pixel.b, .g = pixel.g, .b = pixel.r};
            fwrite(&color, 3, 1, file);
        }

        fwrite(padding, 1, padding_size, file);
    }

    fclose(file);
}

int _psr_split(char* str, char const* delims, char** tokens, int max) {
    int count = 0;
    char* token = strtok(str, delims);
    for (; count < max && token; ++count) {
        tokens[count] = token;
        token = strtok(NULL, delims);
    }
    return count;
}

int psr_load_obj(char const* filename, psr_mesh_t* mesh) {
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

    mesh->faces = malloc(face_max * sizeof(psr_face_t));
    mesh->positions = malloc(position_max * sizeof(psr_float3_t));
    mesh->tex_coords = malloc(tex_coord_max * sizeof(psr_float2_t));
    mesh->normals = malloc(normal_max * sizeof(psr_float3_t));

    char line[256];
    while (fgets(line, 256, file)) {
        char* tokens[5];
        int token_count = _psr_split(line, " ", tokens, 5);

        if (token_count == 0) {
            continue;
        }

        // Vertex position.
        if (strcmp(tokens[0], "v") == 0) {
            if (mesh->position_count == position_max) {
                mesh->positions = realloc(mesh->positions, (position_max *= 2) * sizeof(psr_float3_t));
            }
            mesh->positions[mesh->position_count++] = (psr_float3_t){atof(tokens[1]), atof(tokens[2]), atof(tokens[3])};
        }
        // Vertex texture coordinate.
        else if (strcmp(tokens[0], "vt") == 0) {
            if (mesh->tex_coord_count == tex_coord_max) {
                mesh->tex_coords = realloc(mesh->tex_coords, (tex_coord_max *= 2) * sizeof(psr_float2_t));
            }
            mesh->tex_coords[mesh->tex_coord_count++] = (psr_float2_t){atof(tokens[1]), atof(tokens[2])};
        }
        // Vertex normal.
        else if (strcmp(tokens[0], "vn") == 0) {
            if (mesh->normal_count == normal_max) {
                mesh->normals = realloc(mesh->normals, (normal_max *= 2) * sizeof(psr_float3_t));
            }
            mesh->normals[mesh->normal_count++] = (psr_float3_t){atof(tokens[1]), atof(tokens[2]), atof(tokens[3])};
        }
        // Face.
        else if (strcmp(tokens[0], "f") == 0) {
            char* face_tokens[4][3];
            int face_token_count = token_count - 1;
            for (int i = 0; i < face_token_count; ++i) {
                // Split into indices.
                if (_psr_split(tokens[i + 1], "/", face_tokens[i], 3) != 3) {
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
                psr_face_t face;
                for (int i = 0; i < 3; ++i) {
                    face.position_indices[i] = position_indices[i];
                    face.tex_coord_indices[i] = tex_coord_indices[i];
                    face.normal_indices[i] = normal_indices[i];
                }

                if (mesh->face_count == face_max) {
                    mesh->faces = realloc(mesh->faces, (face_max *= 2) * sizeof(psr_face_t));
                }
                mesh->faces[mesh->face_count++] = face;
            }
            else if (face_token_count == 4) {
                psr_face_t a, b;
                for (int i = 0; i < 3; ++i) {
                    a.position_indices[i] = position_indices[i];
                    a.tex_coord_indices[i] = tex_coord_indices[i];
                    a.normal_indices[i] = normal_indices[i];

                    b.position_indices[i] = position_indices[(i + 2) % 4];
                    b.tex_coord_indices[i] = tex_coord_indices[(i + 2) % 4];
                    b.normal_indices[i] = normal_indices[(i + 2) % 4];
                }

                if (mesh->face_count >= face_max - 1) {
                    mesh->faces = realloc(mesh->faces, (face_max *= 2) * sizeof(psr_face_t));
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

void psr_mesh_free(psr_mesh_t* mesh) {
    free(mesh->positions);
    free(mesh->tex_coords);
    free(mesh->normals);
    free(mesh->faces);
}
