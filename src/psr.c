#include "psr.h"

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>

typedef struct {
    psr_vertex_t start;
    psr_vertex_t end;
    psr_vertex_t* current;

    // Do not use this to plot pixels, use _psr_line_t::current instead.
    psr_float2_t _current_pos;

    int i;
    int steps;

    psr_float3_t increment;
    psr_float4_t* attribute_increments;
} _psr_line_t;

typedef struct {
    int id;
    psr_rect_t src;
    psr_int2_t offset;
    int x_advance;
} _psr_char_info_t;

struct psr_font_t {
    psr_image_t* image;
    _psr_char_info_t* char_infos;
    int char_info_count;
    int size;
};

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

float psr_lerp(float a, float b, float amount) {
    return a * (1 - amount) + (b * amount);
}

psr_byte3_t psr_byte3_lerp(psr_byte3_t a, psr_byte3_t b, float amount) {
    psr_byte3_t result;
    for (int i = 0; i < 3; i++) {
        result.values[i] = a.values[i] * (1 - amount) + b.values[i] * amount;
    }
    return result;
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

psr_color_buffer_t* psr_color_buffer_create(int w, int h) {
    psr_color_buffer_t* buf = malloc(sizeof(psr_color_buffer_t));
    buf->w = w;
    buf->h = h;
    buf->data = malloc(w * h * sizeof(psr_byte3_t));
    return buf;
}

void psr_color_buffer_destroy(psr_color_buffer_t* color_buffer) {
    free(color_buffer->data);
    free(color_buffer);
    color_buffer = NULL;
}

void psr_color_buffer_resize(psr_color_buffer_t* color_buffer, int w, int h) {
    color_buffer->w = w;
    color_buffer->h = h;
    free(color_buffer->data);
    color_buffer->data = malloc(w * h * sizeof(psr_byte3_t));
}

psr_byte3_t* psr_color_buffer_at(psr_color_buffer_t* color_buffer, int x, int y) {
    return &color_buffer->data[y * color_buffer->w + x];
}

void psr_color_buffer_clear(psr_color_buffer_t* color_buffer, psr_byte3_t color) {
    for (int i = 0; i < color_buffer->w * color_buffer->h; i++) {
        color_buffer->data[i] = color;
    }
}

psr_depth_buffer_t* psr_depth_buffer_create(int w, int h) {
    psr_depth_buffer_t* buf = malloc(sizeof(psr_depth_buffer_t));
    buf->w = w;
    buf->h = h;
    buf->data = malloc(w * h * sizeof(float));
    return buf;
}

void psr_depth_buffer_destroy(psr_depth_buffer_t* depth_buffer) {
    free(depth_buffer->data);
    free(depth_buffer);
    depth_buffer = NULL;
}

void psr_depth_buffer_resize(psr_depth_buffer_t* depth_buffer, int w, int h) {
    depth_buffer->w = w;
    depth_buffer->h = h;
    free(depth_buffer->data);
    depth_buffer->data = malloc(w * h * sizeof(float));
}

float* psr_depth_buffer_at(psr_depth_buffer_t* depth_buffer, int x, int y) {
    return &depth_buffer->data[y * depth_buffer->w + x];
}

void psr_depth_buffer_clear(psr_depth_buffer_t* depth_buffer) {
    for (int i = 0; i < depth_buffer->w * depth_buffer->h; i++) {
        depth_buffer->data[i] = 1;
    }
}

// Bytes per pixel.
int _psr_pixel_size(psr_color_depth_t color_depth) {
    switch (color_depth) {
    case PSR_R5G6B5:
    case PSR_R5G5B5A1:
        return 2;
    case PSR_R8G8B8:
        return 3;
    case PSR_R8G8B8A8:
        return 4;
    default:
        assert(!"Not implemented.");
        return -1;
    }
}

int _psr_channels(psr_color_depth_t color_depth) {
    switch (color_depth) {
    case PSR_R5G6B5:
    case PSR_R8G8B8:
        return 3;
    case PSR_R5G5B5A1:
    case PSR_R8G8B8A8:
        return 4;
    default:
        assert(!"Not implemented.");
        return -1;
    }
}

psr_image_t* psr_image_create(psr_color_depth_t color_depth, int w, int h) {
    psr_image_t* img = malloc(sizeof(psr_image_t));
    img->w = w;
    img->h = h;
    img->color_depth = color_depth;
    img->data = malloc(w * h * _psr_pixel_size(color_depth));
    return img;
}

void psr_image_destroy(psr_image_t* image) {
    free(image->data);
    free(image);
    image = NULL;
}

psr_byte_t* psr_image_at(psr_image_t* image, int x, int y) {
    return &image->data[(y * image->w + x) * _psr_pixel_size(image->color_depth)];
}

// Remember to free return value.
psr_vertex_t* _psr_vertex_copy(psr_vertex_t* src) {
    psr_vertex_t* v = malloc(sizeof(psr_vertex_t));

    // Copy value fields.
    memcpy(v, src, sizeof(psr_vertex_t));

    // Allocate memory for and copy attributes.
    v->attributes = malloc(src->attribute_count * sizeof(psr_vertex_attribute_t));
    memcpy(v->attributes, src->attributes, src->attribute_count * sizeof(psr_vertex_attribute_t));
    
    return v;
}

_psr_line_t* _psr_line_create(psr_vertex_t start, psr_vertex_t end) {    
    _psr_line_t* line = malloc(sizeof(_psr_line_t));
    
    line->start = start;
    line->end = end;

    line->current = _psr_vertex_copy(&start);

    line->_current_pos.x = start._x;
    line->_current_pos.y = start._y;

    psr_int2_t delta = {end.x - start.x, end.y - start.y};

    line->i = 0;
    line->steps = (abs(delta.x) >= abs(delta.y)) ? abs(delta.x) : abs(delta.y);
    line->increment.x = delta.x / (float)line->steps;
    line->increment.y = delta.y / (float)line->steps;

    return line;
}

void _psr_line_destroy(_psr_line_t* line) {
    free(line->current->attributes);
    free(line);
    line = NULL;
}

bool _psr_line_step(_psr_line_t* line) {
    if (line->i++ <= line->steps) {
        // X, Y, Z increment.

        line->_current_pos.x += line->increment.x;
        line->_current_pos.y += line->increment.y;
        line->current->z += line->increment.z;

        line->current->x = (int)roundf(line->_current_pos.x);
        line->current->y = (int)roundf(line->_current_pos.y);

        // Vertex attributes increment.

        for (int i = 0; i < line->current->attribute_count; i++) {
            PSR_ADD(
                line->current->attributes[i].data, 
                line->current->attributes[i].data, 
                line->attribute_increments[i], 
                line->current->attributes[i].float_count);            
        }

        return true;
    }
    return false;
}

bool _psr_line_step_until_vertical_difference(_psr_line_t* line) {
    int y = line->current->y;

    while (_psr_line_step(line)) {
        if (line->current->y != y) {
            return true;
        }
        y = line->current->y;
    }
    return false;
}

psr_vertex_t psr_vertex_create(psr_float3_t pos, psr_vertex_attribute_t* attributes, int attribute_count) {
    psr_vertex_t v;
    v.x = (int)roundf(pos.x);
    v.y = (int)roundf(pos.y);
    v._x = pos.x;
    v._y = pos.y;
    v.z = pos.z;
    v.attributes = attributes;
    v.attribute_count = attribute_count;
    return v;
}

// Remember to free return value.
psr_vertex_t* _psr_vertex_lerp(const psr_vertex_t* a, const psr_vertex_t* b, float amount) {
    assert(a->attribute_count == b->attribute_count);

    psr_vertex_t* v = malloc(sizeof(psr_vertex_t));
    v->attributes = malloc(a->attribute_count * sizeof(psr_vertex_attribute_t));
    v->attribute_count = a->attribute_count;
    v->x = psr_lerp(a->x, b->x, amount);
    v->y = psr_lerp(a->y, b->y, amount);
    v->z = psr_lerp(a->z, b->z, amount);
    v->_x = psr_lerp(a->_x, b->_x, amount);
    v->_y = psr_lerp(a->_y, b->_y, amount);

    for (int i = 0; i < a->attribute_count; i++) {
        for (int j = 0; j < a->attributes[i].float_count; j++) {
            v->attributes[i].data.values[j] = 
                psr_lerp(
                    a->attributes[i].data.values[j], 
                    b->attributes[i].data.values[j], 
                    amount);
        }
    }
    return v;
}

void _psr_vertex_destroy(psr_vertex_t* v) {
    free(v->attributes);
    free(v);
    v = NULL;
}

void _psr_raster_triangle_3d_flat_top_or_flat_bottom(psr_color_buffer_t* color_buffer, psr_depth_buffer_t* depth_buffer, _psr_line_t* line_a, _psr_line_t* line_b, psr_pixel_shader_callback pixel_shader, void* user_data) {
    if (line_a->start.x > line_b->start.x || line_a->end.x > line_b->end.x) {
        PSR_SWAP(_psr_line_t*, line_a, line_b);
    }
    
    for (int y = line_a->start.y; y <= line_a->end.y; y++) {
        psr_vertex_t start = *line_a->current;
        psr_vertex_t end = *line_b->current;

        for (int x = start.x; x <= end.x; x++) {
            float how_far_along_are_we = (x - start._x) / (end._x - start._x + 1);
            psr_vertex_t* interpolated = _psr_vertex_lerp(&start, &end, how_far_along_are_we);
            
            psr_byte3_t color = pixel_shader(interpolated->attributes, user_data);
            
            if (depth_buffer) {
                if (interpolated->z < *psr_depth_buffer_at(depth_buffer, x, y)) {
                    *psr_depth_buffer_at(depth_buffer, x, y) = interpolated->z;
                    *psr_color_buffer_at(color_buffer, x, y) = color;
                }
            }
            else {
                *psr_color_buffer_at(color_buffer, x, y) = color;
            }

            _psr_vertex_destroy(interpolated);
        }

        _psr_line_step_until_vertical_difference(line_a);
        _psr_line_step_until_vertical_difference(line_b);
    }
}

void psr_raster_triangle_3d(psr_color_buffer_t* color_buffer, psr_depth_buffer_t* depth_buffer, psr_vertex_t v0, psr_vertex_t v1, psr_vertex_t v2, psr_pixel_shader_callback pixel_shader, void* user_data) {
    // Sort vertices by height.
    if (v0.y > v1.y) {
        PSR_SWAP(psr_vertex_t, v0, v1);
    }
    if (v0.y > v2.y) {
        PSR_SWAP(psr_vertex_t, v0, v2);
    }
    if (v1.y > v2.y) {
        PSR_SWAP(psr_vertex_t, v1, v2);
    }

    // HACK: discard triangle if outside view.
    if (v0.x < 0 || v0.x >= color_buffer->w || v0.y < 0 || v0.y >= color_buffer->h ||
        v1.x < 0 || v1.x >= color_buffer->w || v1.y < 0 || v1.y >= color_buffer->h ||
        v2.x < 0 || v2.x >= color_buffer->w || v2.y < 0 || v2.y >= color_buffer->h) {
        return;
    }

    // Check if the top of the triangle is flat.
    if (v0.y == v1.y) {
        _psr_line_t* line_a = _psr_line_create(v0, v2);
        _psr_line_t* line_b = _psr_line_create(v1, v2);

        _psr_raster_triangle_3d_flat_top_or_flat_bottom(color_buffer, depth_buffer, line_a, line_b, pixel_shader, user_data);
        
        _psr_line_destroy(line_b);
        _psr_line_destroy(line_a);
    }
    // Check if the bottom is flat.
    else if (v1.y == v2.y) {
        _psr_line_t* line_a = _psr_line_create(v0, v1);
        _psr_line_t* line_b = _psr_line_create(v0, v2);

        _psr_raster_triangle_3d_flat_top_or_flat_bottom(color_buffer, depth_buffer, line_a, line_b, pixel_shader, user_data);

        _psr_line_destroy(line_b);
        _psr_line_destroy(line_a);
    }
    // Else split into two smaller triangles.
    else {
        // 0 = v0._y, 1 = v2._y
        float alpha = (v1._y - v0._y) / (v2._y - v0._y);

        psr_vertex_t* v3 = _psr_vertex_lerp(&v0, &v2, alpha);

        _psr_line_t* line_0_1 = _psr_line_create(v0, v1);
        _psr_line_t* line_0_3 = _psr_line_create(v0, *v3);
        _psr_line_t* line_3_2 = _psr_line_create(*v3, v2);
        _psr_line_t* line_1_2 = _psr_line_create(v1, v2);

        // Top (flat bottom).
        _psr_raster_triangle_3d_flat_top_or_flat_bottom(color_buffer, depth_buffer, line_0_1, line_0_3, pixel_shader, user_data);
        // Bottom (flat top).
        _psr_raster_triangle_3d_flat_top_or_flat_bottom(color_buffer, depth_buffer, line_3_2, line_1_2, pixel_shader, user_data);

        _psr_line_destroy(line_1_2);
        _psr_line_destroy(line_3_2);
        _psr_line_destroy(line_0_3);
        _psr_line_destroy(line_0_1);

        _psr_vertex_destroy(v3);
    }
}

void psr_raster_image(psr_color_buffer_t* color_buffer, psr_image_t* image, psr_rect_t src, psr_rect_t dst) {    
    psr_float2_t scale_factor = {src.w / (float)dst.w, src.h / (float)dst.h};
    
    for (int x = 0; x < dst.w; x++) {
        for (int y = 0; y < dst.h; y++) {
            psr_int2_t sample = {src.x + (int)floorf(x * scale_factor.x), src.y + (int)floorf(y * scale_factor.y)};

            psr_byte_t* src_pixel = psr_image_at(image, sample.x, sample.y);
            psr_byte3_t* dst_pixel = psr_color_buffer_at(color_buffer, dst.x + x, dst.y + y);

            switch (image->color_depth) {
            case PSR_R8G8B8:
                for (int i = 0; i < 3; i++) {
                    dst_pixel->values[i] = *(src_pixel + i);
                }
                break;
            case PSR_R8G8B8A8:
                // HACK: skip pixel if not fully opaque.
                if (*(src_pixel + 3) < 255) {
                    break;
                }
                for (int i = 0; i < 3; i++) {
                    dst_pixel->values[i] = *(src_pixel + i);
                }
                break;
            default:
                assert(!"Not implemented.");
            }
        }
    }
}

void psr_raster_text(psr_color_buffer_t* color_buffer, char* text, psr_int2_t pos, psr_font_t* font, int scale) {       
    for (size_t i = 0; i < strlen(text); i++) {
        _psr_char_info_t info = font->char_infos[(int)text[i]];

        psr_rect_t dst = {pos.x + info.offset.x * scale, pos.y + info.offset.y * scale, info.src.w * scale, info.src.h * scale};
        psr_raster_image(color_buffer, font->image, info.src, dst);

        pos.x += info.x_advance * scale;
    }
}

psr_image_t* psr_image_load_bmp(char* path, psr_color_depth_t color_depth) {
    FILE* file = fopen(path, "rb");

    if (!file) {
        return NULL;
    }

    psr_byte_t header[54];
    fread(header, 1, 54, file);

    int file_size = (header[2]) | (header[3] << 8) | (header[4] << 16) | (header[5] << 24);
    int header_size = (header[10]) | (header[11] << 8) | (header[12] << 16) | (header[13] << 24);
    int width = (header[18]) | (header[19] << 8) | (header[20] << 16) | (header[21] << 24);
    int height = (header[22]) | (header[23] << 8) | (header[24] << 16) | (header[25] << 24);
    int pixel_size = header[28] / 8;
    int padding_size = (4 - width * pixel_size % 4) % 4;

    assert(pixel_size == _psr_pixel_size(color_depth) && "Incompatible color depths.");

    psr_image_t* image = psr_image_create(color_depth, width, height);

    psr_byte_t* pixel_data = malloc(file_size - header_size);
    fseek(file, header_size, SEEK_SET);
    fread(pixel_data, 1, file_size - header_size, file);
    fclose(file);

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int i = ((height - 1 - y) * (width + padding_size) + x) * pixel_size;

            switch (color_depth) {
            case PSR_R8G8B8:
                for (int j = 0; j < 3; j++) {
                    *(psr_image_at(image, x, y) + j) = *(pixel_data + i + (2 - j));
                }
                break;
            case PSR_R8G8B8A8:
                // Alpha.
                *(psr_image_at(image, x, y) + 3) = *(pixel_data + i + 3);
                // RGB.
                for (int j = 0; j < 3; j++) {
                    *(psr_image_at(image, x, y) + j) = *(pixel_data + i + (2 - j));
                }
                break;
            default:
                assert(!"Not implemented.");
            }
        }
    }
    free(pixel_data);

    return image;
}

void psr_save_bmp(char* path, psr_color_buffer_t* color_buffer) {
    psr_color_buffer_t* buf = color_buffer;
    
    psr_byte_t padding[] = {0, 0, 0};
    int padding_size = (4 - buf->w * 3 % 4) % 4;
    int stride = buf->w * 3 + padding_size;
    int file_size = 54 + stride * buf->h;

    psr_byte_t header[54] = {
        'B', 'M',
        (psr_byte_t)file_size, (psr_byte_t)(file_size >> 8), (psr_byte_t)(file_size >> 16), (psr_byte_t)(file_size >> 24),
        0, 0, 0, 0,
        54, 0, 0, 0,
        40, 0, 0, 0,
        (psr_byte_t)buf->w, (psr_byte_t)(buf->w >> 8), (psr_byte_t)(buf->w >> 16), (psr_byte_t)(buf->w >> 24),
        (psr_byte_t)buf->h, (psr_byte_t)(buf->h >> 8), (psr_byte_t)(buf->h >> 16), (psr_byte_t)(buf->h >> 24),
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

    for (int y = 0; y < buf->h; y++) {
        for (int x = 0; x < buf->w; x++) {
            int i = (buf->h - 1 - y) * buf->h + x;
            psr_byte3_t pixel = buf->data[i];
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

psr_mesh_t* psr_mesh_load_obj(char* path) {
    psr_mesh_t* mesh = malloc(sizeof(psr_mesh_t));

    FILE* file = fopen(path, "r");
    if (!file) {
        return NULL;
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
                    return NULL;
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
                return NULL;
            }
        }
    }
    fclose(file);
    return mesh;
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

// NOTE: leaks a lot of memory, but the functions that call this one
// are typically only called a few times in the beginning of the
// apps lifetime so that's why it hasn't been fixed.
char** _psr_split2(char* str, char* delims) {
    int count = 2;
    char** tokens = malloc(count * sizeof(char*));
    char* copy = strdup(str);

    char* token = strtok(copy, delims);
    int i = 0;
    while (token != NULL) {
        tokens[i] = strdup(token);
        token = strtok(NULL, delims);

        if (++i == count) {
            tokens = realloc(tokens, (count *= 2) * sizeof(char*));
        }
    }
    tokens[i] = NULL;
    return tokens;
}

int _psr_parse_int(char* str) {
    char** tokens = _psr_split2(str, "=");
    return atoi(tokens[1]);
}

psr_font_t* psr_font_load(psr_image_t* image, char* info_path) {  
    psr_font_t* font = malloc(sizeof(psr_font_t));
    
    font->image = image;
    font->char_info_count = 256;
    font->char_infos = malloc(256 * sizeof(_psr_char_info_t));
    font->size = -1;

    for (int i = 0; i < 256; i++) {
        font->char_infos[i].id = 0;
    }

    // Read file content into string.
    char* content;
    FILE* file = fopen(info_path, "rb");
    if (!file) {
        return NULL;
    }
    fseek(file, 0, SEEK_END);
    int size = ftell(file);
    fseek(file, 0, SEEK_SET);
    content = malloc(size + 1);
    fread(content, 1, size, file);
    content[size] = '\0';
    fclose(file);

    // Split into lines.
    char** lines = _psr_split2(content, "\n");
    for (char** line = lines; *line; line++) {
        // Split each line into words.
        char** words = _psr_split2(*line, " ");

        if (strcmp(words[0], "info") == 0) {
            font->size = _psr_parse_int(words[2]);
        }
        else if (strcmp(words[0], "char") == 0) {
            int id = _psr_parse_int(words[1]);

            if (id < 0 || id >= 256) {
                continue;
            }

            font->char_infos[id].id = id;
            font->char_infos[id].src.x = _psr_parse_int(words[2]);
            font->char_infos[id].src.y = _psr_parse_int(words[3]);
            font->char_infos[id].src.w = _psr_parse_int(words[4]);
            font->char_infos[id].src.h = _psr_parse_int(words[5]);
            font->char_infos[id].offset.x = _psr_parse_int(words[6]);
            font->char_infos[id].offset.y = _psr_parse_int(words[7]);
            font->char_infos[id].x_advance = _psr_parse_int(words[8]);
        }
    }
    return font;
}

void psr_font_destroy(psr_font_t* font) {
    free(font->char_infos);
    free(font);
    font = NULL;
}
