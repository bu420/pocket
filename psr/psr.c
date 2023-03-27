#include "psr.h"

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define PSR_SWAP(type, a, b) { type _temp = a; a = b; b = _temp; }

#define PSR_SORT_TRIANGLE_VERTICES_BY_HEIGHT(type, a, b, c) \
    if (a.y > b.y)                                        \
        PSR_SWAP(type, a, b);                               \
    if (a.y > c.y)                                        \
        PSR_SWAP(type, a, c);                               \
    if (b.y > c.y)                                        \
        PSR_SWAP(type, b, c);

typedef struct {
    psr_int2_t start;
    psr_int2_t end;
    psr_int2_t current;
    psr_int2_t delta;
    psr_int2_t dir;
    int deviation;
} psr_line_2d_t;

typedef struct {
    struct {
        int x;
        int y;
        float z;
    } start, end, current;
    psr_int2_t delta;
    psr_int2_t dir;
    int deviation;
} psr_line_3d_t;

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

psr_byte3_t psr_byte3_lerp(psr_byte3_t a, psr_byte3_t b, float amount) {
    psr_byte3_t result;
    for (int i = 0; i < 3; i++) {
        result.values[i] = a.values[i] * (1 - amount) + b.values[i] * amount;
    }
    return result;
}

psr_float3_t psr_float3_add(psr_float3_t a, psr_float3_t b) {
    return (psr_float3_t){a.x + b.x, a.y + b.y, a.z + b.z};
}

psr_float3_t psr_float3_sub(psr_float3_t a, psr_float3_t b) {
    return (psr_float3_t){a.x - b.x, a.y - b.y, a.z - b.z};
}

psr_float3_t psr_float3_mul(psr_float3_t a, psr_float3_t b) {
    return (psr_float3_t){a.x * b.x, a.y * b.y, a.z * b.z};
}

psr_float3_t psr_float3_div(psr_float3_t a, psr_float3_t b) {
    return (psr_float3_t){a.x / b.x, a.y / b.y, a.z / b.z};
}

void psr_mat4_init_zero(psr_mat4_t* m) {
    for (int x = 0; x < 4; x++) {
        for (int y = 0; y < 4; y++) {
            m->m[x][y] = 0;
        }
    }
}

void psr_mat4_init_identity(psr_mat4_t* m) {
    psr_mat4_init_zero(m);
    for (int i = 0; i < 4; i++) {
        m->m[i][i] = 1;
    }
}

psr_mat4_t psr_mat4_mul(psr_mat4_t a, psr_mat4_t b) {
    psr_mat4_t result;
    psr_mat4_init_zero(&result);

    for (int x = 0; x < 4; x++) {
        for (int y = 0; y < 4; y++) {
            for (int i = 0; i < 4; i++) {
                result.m[x][y] += a.m[x][i] * b.m[i][y];
            }
        }
    }
    return result;
}

psr_mat4_t psr_mat4_transpose(psr_mat4_t m) {
    psr_mat4_t r;
    for (int x = 0; x < 4; x++) {
        for (int y = 0; y < 4; y++) {
            r.m[y][x] = m.m[x][y];
        }
    }
    return r;
}

psr_mat4_t psr_mat4_inverse(psr_mat4_t m) {
    float A2323 = m.m22 * m.m33 - m.m23 * m.m32;
    float A1323 = m.m21 * m.m33 - m.m23 * m.m31;
    float A1223 = m.m21 * m.m32 - m.m22 * m.m31;
    float A0323 = m.m20 * m.m33 - m.m23 * m.m30;
    float A0223 = m.m20 * m.m32 - m.m22 * m.m30;
    float A0123 = m.m20 * m.m31 - m.m21 * m.m30;
    float A2313 = m.m12 * m.m33 - m.m13 * m.m32;
    float A1313 = m.m11 * m.m33 - m.m13 * m.m31; 
    float A1213 = m.m11 * m.m32 - m.m12 * m.m31;
    float A2312 = m.m12 * m.m23 - m.m13 * m.m22;
    float A1312 = m.m11 * m.m23 - m.m13 * m.m21;
    float A1212 = m.m11 * m.m22 - m.m12 * m.m21;
    float A0313 = m.m10 * m.m33 - m.m13 * m.m30;
    float A0213 = m.m10 * m.m32 - m.m12 * m.m30;
    float A0312 = m.m10 * m.m23 - m.m13 * m.m20;
    float A0212 = m.m10 * m.m22 - m.m12 * m.m20;
    float A0113 = m.m10 * m.m31 - m.m11 * m.m30;
    float A0112 = m.m10 * m.m21 - m.m11 * m.m20;

    float det = 
        m.m00 * (m.m11 * A2323 - m.m12 * A1323 + m.m13 * A1223) 
        - m.m01 * (m.m10 * A2323 - m.m12 * A0323 + m.m13 * A0223) 
        + m.m02 * (m.m10 * A1323 - m.m11 * A0323 + m.m13 * A0123) 
        - m.m03 * (m.m10 * A1223 - m.m11 * A0223 + m.m12 * A0123);
    
    assert(det > 0 && "Matrix is not invertible.");
    det = 1 / det;

    psr_mat4_t result = {
        .m00 = det *   (m.m11 * A2323 - m.m12 * A1323 + m.m13 * A1223),
        .m01 = det * - (m.m01 * A2323 - m.m02 * A1323 + m.m03 * A1223),
        .m02 = det *   (m.m01 * A2313 - m.m02 * A1313 + m.m03 * A1213),
        .m03 = det * - (m.m01 * A2312 - m.m02 * A1312 + m.m03 * A1212),
        .m10 = det * - (m.m10 * A2323 - m.m12 * A0323 + m.m13 * A0223),
        .m11 = det *   (m.m00 * A2323 - m.m02 * A0323 + m.m03 * A0223),
        .m12 = det * - (m.m00 * A2313 - m.m02 * A0313 + m.m03 * A0213),
        .m13 = det *   (m.m00 * A2312 - m.m02 * A0312 + m.m03 * A0212),
        .m20 = det *   (m.m10 * A1323 - m.m11 * A0323 + m.m13 * A0123),
        .m21 = det * - (m.m00 * A1323 - m.m01 * A0323 + m.m03 * A0123),
        .m22 = det *   (m.m00 * A1313 - m.m01 * A0313 + m.m03 * A0113),
        .m23 = det * - (m.m00 * A1312 - m.m01 * A0312 + m.m03 * A0112),
        .m30 = det * - (m.m10 * A1223 - m.m11 * A0223 + m.m12 * A0123),
        .m31 = det *   (m.m00 * A1223 - m.m01 * A0223 + m.m02 * A0123),
        .m32 = det * - (m.m00 * A1213 - m.m01 * A0213 + m.m02 * A0113),
        .m33 = det *   (m.m00 * A1212 - m.m01 * A0212 + m.m02 * A0112)
    };
    return result;
}

psr_mat4_t psr_mat4_translate(psr_mat4_t mat4, psr_float3_t f) {
    psr_mat4_t m;
    psr_mat4_init_identity(&m);

    m.m[3][0] = f.x;
    m.m[3][1] = f.y;
    m.m[3][2] = f.z;

    return psr_mat4_mul(mat4, m);
}

psr_mat4_t psr_mat4_rotate_x(psr_mat4_t m, float a) {
    float s = sinf(a);
    float c = cosf(a);

    psr_mat4_t x;
    psr_mat4_init_identity(&x);

    x.m[1][1] = c;
    x.m[1][2] = -s;
    x.m[2][1] = s;
    x.m[2][2] = c;

    return psr_mat4_mul(m, x);
}

psr_mat4_t psr_mat4_rotate_y(psr_mat4_t m, float a) {
    float s = sinf(a);
    float c = cosf(a);

    psr_mat4_t y;
    psr_mat4_init_identity(&y);

    y.m[0][0] = c;
    y.m[0][2] = s;
    y.m[2][0] = -s;
    y.m[2][2] = c;

    return psr_mat4_mul(m, y);
}

psr_mat4_t psr_mat4_rotate_z(psr_mat4_t m, float a) {
    float s = sinf(a);
    float c = cosf(a);

    psr_mat4_t z;
    psr_mat4_init_identity(&z);

    z.m[0][0] = c;
    z.m[0][1] = -s;
    z.m[1][0] = s;
    z.m[1][1] = c;

    return psr_mat4_mul(m, z);
}

psr_mat4_t psr_mat4_scale(psr_mat4_t mat4, psr_float3_t f) {
    psr_mat4_t m;
    psr_mat4_init_identity(&m);

    m.m[0][0] = f.x;
    m.m[1][1] = f.y;
    m.m[2][2] = f.z;

    return psr_mat4_mul(mat4, m);
}

psr_float3_t psr_mat3_mul_float3(psr_mat3_t m, psr_float3_t f) {
    psr_float3_t result;
    for (int x = 0; x < 3; x++) {
        float i = 0;
        for (int y = 0; y < 3; y++) {
            i += m.m[x][y] * f.values[y];
        }
        result.values[x] = i;
    }
    return result;
}

psr_float3_t psr_float3_mul_mat3(psr_float3_t f, psr_mat3_t m) {
    psr_float3_t result;
    for (int y = 0; y < 3; y++) {
        float i = 0;
        for (int x = 0; x < 3; x++) {
            i += m.m[x][y] * f.values[x];
        }
        result.values[y] = i;
    }
    return result;
}

psr_float3_t psr_mat4_mul_float3(psr_mat4_t m, psr_float3_t f) {
    psr_float4_t f4;
    for (int x = 0; x < 4; x++) {
        float i = 0;
        for (int y = 0; y < 3; y++) {
            i += m.m[x][y] * f.values[y];
        }
        f4.values[x] = i + m.m[x][3];
    }
    
    if (f4.w != 0) {
        f4.x /= f4.w;
        f4.y /= f4.w;
        f4.z /= f4.w;
    }

    psr_float3_t result = {f4.x, f4.y, f4.z};
    return result;
}

psr_float4_t psr_mat4_mul_float4(psr_mat4_t m, psr_float4_t f) {
    psr_float4_t result;
    for (int x = 0; x < 4; x++) {
        float i = 0;
        for (int y = 0; y < 4; y++) {
            i += m.m[x][y] * f.values[y];
        }
        result.values[x] = i;
    }
    return result;
}

psr_float3_t psr_float3_mul_mat4(psr_float3_t f, psr_mat4_t m) {
    psr_float4_t f4;
    for (int y = 0; y < 4; y++) {
        float i = 0;
        for (int x = 0; x < 3; x++) {
            i += m.m[x][y] * f.values[x];
        }
        f4.values[y] = i + m.m[3][y];
    }
    
    if (f4.w != 0) {
        f4.x /= f4.w;
        f4.y /= f4.w;
        f4.z /= f4.w;
    }

    psr_float3_t result = {f4.x, f4.y, f4.z};
    return result;
}

psr_float4_t psr_float4_mul_mat4(psr_float4_t f, psr_mat4_t m) {
    psr_float4_t result;
    for (int y = 0; y < 4; y++) {
        float i = 0;
        for (int x = 0; x < 4; x++) {
            i += m.m[x][y] * f.values[x];
        }
        result.values[y] = i;
    }
    return result;
}

psr_mat3_t psr_mat4_to_mat3(psr_mat4_t m) {
    psr_mat3_t r;
    for (int x = 0; x < 3; x++) {
        for (int y = 0; y < 3; y++) {
            r.m[x][y] = m.m[x][y];
        }
    }
    return r;
}

psr_mat4_t psr_look_at(psr_float3_t pos, psr_float3_t target, psr_float3_t up) {
    psr_float3_t diff = {target.x - pos.x, target.y - pos.y, target.z - pos.z};
    psr_float3_t forward = psr_normalize(diff);
    psr_float3_t right = psr_normalize(psr_cross(forward, up));
    psr_float3_t local_up = psr_normalize(psr_cross(right, forward));

    psr_mat4_t m;
    psr_mat4_init_identity(&m);
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

psr_mat4_t psr_perspective(float aspect, float fov, float near, float far) {    
    psr_mat4_t m;
    psr_mat4_init_zero(&m);

    float half_tan = tanf(fov / 2);

    m.m[0][0] = 1 / (half_tan * aspect);
    m.m[1][1] = 1 / half_tan;
    m.m[2][2] = -(far + near) / (far - near);
    m.m[2][3] = -1;
    m.m[3][2] = -(2 * far * near) / (far - near);
    return m;
}

void psr_color_buffer_init(psr_color_buffer_t* color_buffer, int width, int height) {
    color_buffer->w = width;
    color_buffer->h = height;
    color_buffer->data = malloc(width * height * sizeof(psr_byte3_t));
}

void psr_color_buffer_free(psr_color_buffer_t* color_buffer) {
    free(color_buffer->data);
    color_buffer->data = NULL;
}

psr_byte3_t* psr_color_buffer_at(psr_color_buffer_t* color_buffer, int x, int y) {
    return &color_buffer->data[y * color_buffer->w + x];
}

void psr_color_buffer_clear(psr_color_buffer_t* color_buffer, psr_byte3_t color) {
    for (int i = 0; i < color_buffer->w * color_buffer->h; i++) {
        color_buffer->data[i] = color;
    }
}

void psr_depth_buffer_init(psr_depth_buffer_t* depth_buffer, int width, int height) {
    depth_buffer->w = width;
    depth_buffer->h = height;
    depth_buffer->data = malloc(width * height * sizeof(float));
}

void psr_depth_buffer_free(psr_depth_buffer_t* depth_buffer) {
    free(depth_buffer->data);
    depth_buffer->data = NULL;
}

float* psr_depth_buffer_at(psr_depth_buffer_t* depth_buffer, int x, int y) {
    return &depth_buffer->data[y * depth_buffer->w + x];
}

void psr_depth_buffer_clear(psr_depth_buffer_t* depth_buffer) {
    for (int i = 0; i < depth_buffer->w * depth_buffer->h; i++) {
        depth_buffer->data[i] = 0;
    }
}

// Bytes per pixel.
int _psr_pixel_size(psr_color_depth_t color_depth) {
    switch (color_depth) {
    case PSR_R5G6B5:
    case PSR_A1R5G5B5:
        return 2;
    case PSR_R8G8B8:
        return 3;
    case PSR_A8R8G8B8:
        return 4;
    default:
        assert(!"Not implemented.");
    }
}

int _psr_channels(psr_color_depth_t color_depth) {
    switch (color_depth) {
    case PSR_R5G6B5:
    case PSR_R8G8B8:
        return 3;
    case PSR_A1R5G5B5:
    case PSR_A8R8G8B8:
        return 4;
    default:
        assert(!"Not implemented.");
    }
}

void psr_image_init(psr_image_t* image, psr_color_depth_t color_depth, int w, int h) {
    image->w = w;
    image->h = h;
    image->color_depth = color_depth;
    image->data = malloc(w * h * _psr_pixel_size(color_depth));
}

void psr_image_free(psr_image_t* image) {
    free(image->data);
    image->data = NULL;
}

psr_byte_t* psr_image_at(psr_image_t* image, int x, int y) {
    return &image->data[(y * image->w + x) * _psr_pixel_size(image->color_depth)];
}

psr_line_2d_t _psr_line_2d_create(psr_int2_t start, psr_int2_t end) {
    psr_line_2d_t line;
    line.start = start;
    line.end = end;
    line.current = start;
    line.delta = (psr_int2_t){abs(end.x - start.x), -abs(end.y - start.y)};
    line.dir = (psr_int2_t){start.x < end.x ? 1 : -1, start.y < end.y ? 1 : -1};
    line.deviation = line.delta.x + line.delta.y;
    return line;
}

int _psr_line_2d_step(psr_line_2d_t* line) {
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
float _psr_line_2d_inverse_lerp(psr_int2_t start, psr_int2_t end, psr_int2_t current) {
    if (end.x - start.x > end.y - start.y) {
        return (current.x - start.x) / (float)(end.x - start.x);
    }
    return (current.y - start.y) / (float)(end.y - start.y);
}

void psr_raster_line(psr_color_buffer_t* color_buffer, psr_int2_t start, psr_int2_t end, psr_byte3_t start_color, psr_byte3_t end_color) {
    psr_line_2d_t line = _psr_line_2d_create(start, end);
    do {
        psr_byte3_t current_color = psr_byte3_lerp(start_color, end_color, _psr_line_2d_inverse_lerp(start, end, line.current));
        *psr_color_buffer_at(color_buffer, line.current.x, line.current.y) = current_color;
    } 
    while (_psr_line_2d_step(&line));
}

int _psr_line_2d_step_until_vertical_difference(psr_line_2d_t* line) {
    int y = line->current.y;

    while (_psr_line_2d_step(line)) {
        if (line->current.y != y) {
            return 1;
        }
        y = line->current.y;
    }
    return 0;
}

// Flat top or flat bottom.
void _psr_raster_triangle_2d_flat(psr_color_buffer_t* color_buffer, psr_line_2d_t line_a, psr_line_2d_t line_b, psr_byte3_t color) {
    if (line_a.start.x > line_b.start.x || line_a.end.x > line_b.end.x) {
        PSR_SWAP(psr_line_2d_t, line_a, line_b);
    }
    
    for (int y = line_a.start.y; y <= line_a.end.y; y++) {
        int start = line_a.current.x;
        int end = line_b.current.x;

        for (int x = start; x <= end; x++) {
            *psr_color_buffer_at(color_buffer, x, y) = color;
        }

        _psr_line_2d_step_until_vertical_difference(&line_a);
        _psr_line_2d_step_until_vertical_difference(&line_b);
    }
}

void psr_raster_triangle_2d_color(psr_color_buffer_t* color_buffer, psr_int2_t pos0, psr_int2_t pos1, psr_int2_t pos2, psr_byte3_t color) {
    // HACK: discard triangle if any of it's vertices are outside view.
    if (pos0.x < 0 || pos0.x >= color_buffer->w || pos0.y < 0 || pos0.y >= color_buffer->h ||
        pos1.x < 0 || pos1.x >= color_buffer->w || pos1.y < 0 || pos1.y >= color_buffer->h ||
        pos2.x < 0 || pos2.x >= color_buffer->w || pos2.y < 0 || pos2.y >= color_buffer->h) {
        return;
    }
    
    PSR_SORT_TRIANGLE_VERTICES_BY_HEIGHT(psr_int2_t, pos0, pos1, pos2);

    // Check if the top of the triangle is flat.
    if (pos0.y == pos1.y) {
        _psr_raster_triangle_2d_flat(color_buffer, _psr_line_2d_create(pos0, pos2), _psr_line_2d_create(pos1, pos2), color);
    }
    // Check if the bottom is flat.
    else if (pos1.y == pos2.y) {
        _psr_raster_triangle_2d_flat(color_buffer, _psr_line_2d_create(pos0, pos1), _psr_line_2d_create(pos0, pos2), color);
    }
    // Else plit triangle into two smaller triangles.
    else {
        psr_int2_t pos3 = {(int)(pos0.x + ((float)(pos1.y - pos0.y) / (float)(pos2.y - pos0.y)) * (float)(pos2.x - pos0.x)), pos1.y};
        // Top (flat bottom).
        _psr_raster_triangle_2d_flat(color_buffer, _psr_line_2d_create(pos0, pos1), _psr_line_2d_create(pos0, pos3), color);
        // Bottom (flat top).
        _psr_raster_triangle_2d_flat(color_buffer, _psr_line_2d_create(pos3, pos2), _psr_line_2d_create(pos1, pos2), color);
    }
}

void psr_raster_triangle_2d_callback(psr_color_buffer_t* color_buffer, psr_int2_t pos0, psr_int2_t pos1, psr_int2_t pos2, void (*callback)(psr_int2_t pixel_pos), void* user_data) {
    
}

// Flat top or flat bottom.
void _psr_raster_triangle_3d_flat(psr_color_buffer_t* color_buffer, psr_depth_buffer_t* depth_buffer, psr_line_2d_t line_a, psr_line_2d_t line_b, psr_byte3_t color) {
    if (line_a.start.x > line_b.start.x || line_a.end.x > line_b.end.x) {
        PSR_SWAP(psr_line_2d_t, line_a, line_b);
    }
    
    for (int y = line_a.start.y; y <= line_a.end.y; y++) {
        int start = line_a.current.x;
        int end = line_b.current.x;

        for (int x = start; x <= end; x++) {
            *psr_color_buffer_at(color_buffer, x, y) = color;
        }

        _psr_line_2d_step_until_vertical_difference(&line_a);
        _psr_line_2d_step_until_vertical_difference(&line_b);
    }
}

void psr_raster_triangle_3d(psr_color_buffer_t* color_buffer, psr_depth_buffer_t* depth_buffer, psr_float3_t pos0, psr_float3_t pos1, psr_float3_t pos2, psr_byte3_t color) {
    // HACK: discard triangle if any of it's vertices are outside view.
    if (pos0.x < 0 || pos0.x >= color_buffer->w || pos0.y < 0 || pos0.y >= color_buffer->h ||
        pos1.x < 0 || pos1.x >= color_buffer->w || pos1.y < 0 || pos1.y >= color_buffer->h ||
        pos2.x < 0 || pos2.x >= color_buffer->w || pos2.y < 0 || pos2.y >= color_buffer->h) {
        return;
    }
    
    PSR_SORT_TRIANGLE_VERTICES_BY_HEIGHT(psr_float3_t, pos0, pos1, pos2);

    psr_int2_t xy0 = {(int)roundf(pos0.x), (int)roundf(pos0.y)};
    float z0 = pos0.z;
    psr_int2_t xy1 = {(int)roundf(pos1.x), (int)roundf(pos1.y)};
    float z1 = pos1.z;
    psr_int2_t xy2 = {(int)roundf(pos2.x), (int)roundf(pos2.y)};
    float z2 = pos2.z;

    // Check if the top of the triangle is flat.
    if (xy0.y == xy1.y) {
        _psr_raster_triangle_2d_flat(color_buffer, _psr_line_2d_create(xy0, xy2), _psr_line_2d_create(xy1, xy2), color);
    }
    // Check if the bottom is flat.
    else if (xy1.y == xy2.y) {
        _psr_raster_triangle_2d_flat(color_buffer, _psr_line_2d_create(xy0, xy1), _psr_line_2d_create(xy0, xy2), color);
    }
    // Else plit triangle into two smaller triangles.
    else {
        psr_int2_t xy3 = {(int)(xy0.x + ((float)(xy1.y - xy0.y) / (float)(xy2.y - xy0.y)) * (float)(xy2.x - xy0.x)), xy1.y};
        // Top (flat bottom).
        _psr_raster_triangle_2d_flat(color_buffer, _psr_line_2d_create(xy0, xy1), _psr_line_2d_create(xy0, xy3), color);
        // Bottom (flat top).
        _psr_raster_triangle_2d_flat(color_buffer, _psr_line_2d_create(xy3, xy2), _psr_line_2d_create(xy1, xy2), color);
    }
}

void psr_raster_image(psr_color_buffer_t* color_buffer, psr_image_t image, int x, int y, int sx, int sy, int sw, int sh) {
    for (int _x = 0; _x < sw; _x++) {
        for (int _y = 0; _y < sh; _y++) {
            psr_byte3_t* dst = psr_color_buffer_at(color_buffer, _x + x, _y + y);
            void* src = psr_image_at(&image, _x + sx, _y + sy);

            switch (image.color_depth) {
            case PSR_R8G8B8:
                for (int i = 0; i < 3; i++) {
                    dst->values[i] = *((psr_byte_t*)src + i);
                }
                break;
            case PSR_R5G6B5:
            case PSR_A1R5G5B5:
            case PSR_A8R8G8B8:
            default:
                assert(!"Not implemented.");
            }
        }
    }
}

psr_image_t psr_load_bmp(char* path, psr_color_depth_t color_depth) {
    FILE* file = fopen(path, "rb");

    psr_byte_t header[54];
    fread(header, 1, 54, file);

    int file_size = (header[2]) | (header[3] << 8) | (header[4] << 16) | (header[5] << 24);
    int width = (header[18]) | (header[19] << 8) | (header[20] << 16) | (header[21] << 24);
    int height = (header[22]) | (header[23] << 8) | (header[24] << 16) | (header[25] << 24);
    int pixel_size = header[28] / 8;
    int padding_size = (4 - width * 3 % 4) % 4;

    assert(pixel_size == _psr_pixel_size(color_depth) && "Requested color depth does not match image color depth.");

    psr_image_t image;
    psr_image_init(&image, color_depth, width, height);

    psr_byte_t* pixel_data = malloc(file_size - 54);
    fread(pixel_data, 1, file_size - 54, file);
    fclose(file);

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int i = ((height - y) * (width + padding_size) + x) * pixel_size;

            switch (color_depth) {
            case PSR_R8G8B8:
            case PSR_A8R8G8B8:
                for (int j = 0; j < _psr_channels(color_depth); j++) {
                    *(psr_image_at(&image, x, y) + j) = *(pixel_data + i + (_psr_channels(color_depth) - 1 - j));
                }
                break;
            case PSR_R5G6B5:
            case PSR_A1R5G5B5:
            default:
                assert(!"Not implemented.");
            }
        }
    }
    free(pixel_data);

    return image;
}

void psr_save_bmp(char* path, psr_color_buffer_t color_buffer) {
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

    FILE* file = fopen(path, "wb");

    fwrite(header, 1, 54, file);

    for (int y = 0; y < color_buffer.h; y++) {
        for (int x = 0; x < color_buffer.w; x++) {
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

int psr_load_obj(char* path, psr_mesh_t* mesh) {
    FILE* file = fopen(path, "r");
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
    mesh->positions = NULL;
    mesh->tex_coords = NULL;
    mesh->normals = NULL;
    mesh->faces = NULL;
}
