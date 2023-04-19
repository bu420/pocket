#include "psr.h"

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>

#define _PSR_ABS(x) ((x) < 0 ? -(x) : (x))

typedef enum {
    PSR_LINE_INCREMENT_ANY,
    PSR_LINE_INCREMENT_VERTICAL,
    PSR_LINE_INCREMENT_HORIZONTAL
} _psr_line_increment_type_t;

typedef struct {
    int attribute_count;

    psr_attribute_array_t start;
    psr_attribute_array_t end;
    psr_attribute_array_t current;

    _psr_line_increment_type_t increment_type;
    psr_float4_t increments[PSR_MAX_ATTRIBUTES];

    int i;
    float steps;

    // Floating point positions are stored in the attribute arrays at index 0.
    psr_int2_t start_pixel_pos;
    psr_int2_t end_pixel_pos;
    psr_int2_t current_pixel_pos;
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

void* _psr_malloc(int size) {
    void* mem = malloc(size);
    assert(mem && "Memory allocation failed.");
    return mem;
}

void _psr_free(void* mem) {
    free(mem);
}

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
        result.values[i] = (psr_byte_t)roundf(psr_lerp(a.values[i], b.values[i], amount));
    }
    return result;
}

psr_float3_t psr_float3_lerp(psr_float3_t a, psr_float3_t b, float amount) {
    psr_float3_t result;
    for (int i = 0; i < 3; i++) {
        result.values[i] = psr_lerp(a.values[i], b.values[i], amount);
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
    psr_color_buffer_t* buf = _psr_malloc(sizeof(psr_color_buffer_t));
    buf->w = w;
    buf->h = h;
    buf->data = _psr_malloc(w * h * sizeof(psr_byte3_t));
    return buf;
}

void psr_color_buffer_destroy(psr_color_buffer_t* color_buffer) {
    _psr_free(color_buffer->data);
    _psr_free(color_buffer);
    color_buffer = NULL;
}

void psr_color_buffer_resize(psr_color_buffer_t* color_buffer, int w, int h) {
    color_buffer->w = w;
    color_buffer->h = h;
    _psr_free(color_buffer->data);
    color_buffer->data = _psr_malloc(w * h * sizeof(psr_byte3_t));
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
    psr_depth_buffer_t* buf = _psr_malloc(sizeof(psr_depth_buffer_t));
    buf->w = w;
    buf->h = h;
    buf->data = _psr_malloc(w * h * sizeof(float));
    return buf;
}

void psr_depth_buffer_destroy(psr_depth_buffer_t* depth_buffer) {
    _psr_free(depth_buffer->data);
    _psr_free(depth_buffer);
    depth_buffer = NULL;
}

void psr_depth_buffer_resize(psr_depth_buffer_t* depth_buffer, int w, int h) {
    depth_buffer->w = w;
    depth_buffer->h = h;
    _psr_free(depth_buffer->data);
    depth_buffer->data = _psr_malloc(w * h * sizeof(float));
}

float* psr_depth_buffer_at(psr_depth_buffer_t* depth_buffer, int x, int y) {
    return &depth_buffer->data[y * depth_buffer->w + x];
}

void psr_depth_buffer_clear(psr_depth_buffer_t* depth_buffer) {
    for (int i = 0; i < depth_buffer->w * depth_buffer->h; i++) {
        depth_buffer->data[i] = 1;
    }
}

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
    psr_image_t* img = _psr_malloc(sizeof(psr_image_t));
    img->w = w;
    img->h = h;
    img->color_depth = color_depth;
    img->data = _psr_malloc(w * h * _psr_pixel_size(color_depth));
    return img;
}

void psr_image_destroy(psr_image_t* image) {
    _psr_free(image->data);
    _psr_free(image);
    image = NULL;
}

psr_byte_t* psr_image_at(psr_image_t* image, int x, int y) {
    return &image->data[(y * image->w + x) * _psr_pixel_size(image->color_depth)];
}

// "start" and "end" must each contain a 3d position at index 0.
void _psr_line_create_3d(_psr_line_t* out, 
                         const psr_attribute_array_t* start, 
                         const psr_attribute_array_t* end, 
                         int attribute_count,
                         _psr_line_increment_type_t increment_type) {
    assert(attribute_count >= 1 && attribute_count <= PSR_MAX_ATTRIBUTES);
    out->attribute_count = attribute_count;

    assert(start->attributes[0].nr_of_floats == 3 && end->attributes[0].nr_of_floats == 3);
    psr_float3_t start_pos = PSR_ATTRIB_TO_FLOAT3(start->attributes[0]);
    psr_float3_t end_pos = PSR_ATTRIB_TO_FLOAT3(end->attributes[0]);

    out->i = 0;
    out->increment_type = increment_type;

    psr_float2_t pos_delta = {end_pos.x - start_pos.x, end_pos.y - start_pos.y};

    switch (increment_type) {
    case PSR_LINE_INCREMENT_ANY:        out->steps = (_PSR_ABS(pos_delta.x) >= _PSR_ABS(pos_delta.y)) ? _PSR_ABS(pos_delta.x) : _PSR_ABS(pos_delta.y); break;
    case PSR_LINE_INCREMENT_VERTICAL:   out->steps = _PSR_ABS(pos_delta.y); break;
    case PSR_LINE_INCREMENT_HORIZONTAL: out->steps = _PSR_ABS(pos_delta.x); break;
    default: assert(!"Non-reachable entry reached.");
    }

    memcpy(&out->start, start, sizeof(psr_attribute_array_t));
    memcpy(&out->end, end, sizeof(psr_attribute_array_t));
    memcpy(&out->current, &out->start, sizeof(psr_attribute_array_t));

    // Setup attribute increments.
    if (out->steps < 1) {
        for (int i = 0; i < attribute_count; i++) {
            out->increments[i] = out->start.attributes[i].data;
        }
    }
    else {
        for (int i = 0; i < attribute_count; i++) {
            for (int j = 0; j < out->current.attributes[i].nr_of_floats; j++) {
                float delta = end->attributes[i].data.values[j] - start->attributes[i].data.values[j];
                out->increments[i].values[j] = delta / out->steps;
            }
        }
    }

    out->start_pixel_pos.x = (int)roundf(start_pos.x);
    out->start_pixel_pos.y = (int)roundf(start_pos.y);
    out->end_pixel_pos.x = (int)roundf(end_pos.x);
    out->end_pixel_pos.y = (int)roundf(end_pos.y);
    out->current_pixel_pos = out->start_pixel_pos;
}

void _psr_attributes_interpolate(psr_attribute_array_t* out, const psr_attribute_array_t* a, const psr_attribute_array_t* b, int attribute_count, float amount) {
    assert(attribute_count <= PSR_MAX_ATTRIBUTES);
    
    for (int i = 0; i < attribute_count; i++) {
        assert(a->attributes[i].nr_of_floats == b->attributes[i].nr_of_floats);

        for (int j = 0; j < a->attributes[i].nr_of_floats; j++) {
            out->attributes[i].data.values[j] = psr_lerp(a->attributes[i].data.values[j], b->attributes[i].data.values[j], amount);
            out->attributes[i].nr_of_floats = a->attributes[i].nr_of_floats;
        }
    }
}

bool _psr_line_step(_psr_line_t* line) {
    if (++line->i <= line->steps) {
        for (int i = 0; i < line->attribute_count; i++) {
            PSR_ADD(
                line->current.attributes[i].data, 
                line->current.attributes[i].data,
                line->increments[i],
                line->current.attributes[i].nr_of_floats);
        }

        line->current_pixel_pos.x = (int)roundf(line->current.attributes[0].data.x);
        line->current_pixel_pos.y = (int)roundf(line->current.attributes[0].data.y);

        return true;
    }
    return false;
}

void _psr_raster_triangle_3d_flat_top_or_flat_bottom(psr_color_buffer_t* color_buffer, 
                                                     psr_depth_buffer_t* depth_buffer,
                                                     const psr_attribute_array_t* line_a_start_attribs,
                                                     const psr_attribute_array_t* line_a_end_attribs,
                                                     const psr_attribute_array_t* line_b_start_attribs,
                                                     const psr_attribute_array_t* line_b_end_attribs,
                                                     int attribute_count,
                                                     psr_pixel_shader_callback pixel_shader, 
                                                     void* user_data) {
    _psr_line_t line_a;
    _psr_line_create_3d(&line_a, line_a_start_attribs, line_a_end_attribs, attribute_count, PSR_LINE_INCREMENT_VERTICAL);
    _psr_line_t line_b;
    _psr_line_create_3d(&line_b, line_b_start_attribs, line_b_end_attribs, attribute_count, PSR_LINE_INCREMENT_VERTICAL);
    
    if (line_a.start_pixel_pos.x > line_b.start_pixel_pos.x || line_a.end_pixel_pos.x > line_b.end_pixel_pos.x) {
        PSR_SWAP(_psr_line_t, line_a, line_b);
    }
    
    do {
        assert(line_a.current_pixel_pos.y == line_b.current_pixel_pos.y);

        // Create a horizontal line between line a and line b.
        _psr_line_t line_x;
        _psr_line_create_3d(&line_x, &line_a.current, &line_b.current, line_a.attribute_count, PSR_LINE_INCREMENT_HORIZONTAL);

        do {
            int x = line_x.current_pixel_pos.x;
            int y = line_a.current_pixel_pos.y;
            float z = line_x.current.attributes[0].data.z;

            // Remove first attribute (position) from array.
            psr_attribute_array_t interpolated;
            for (int i = 0; i < line_x.attribute_count; i++) {
                interpolated.attributes[i] = line_x.current.attributes[i + 1];
            }
            
            psr_byte3_t color = pixel_shader(line_x.current_pixel_pos, &interpolated, user_data);
            
            if (depth_buffer) {
                if (z < *psr_depth_buffer_at(depth_buffer, x, y)) {
                    *psr_depth_buffer_at(depth_buffer, x, y) = z;
                    *psr_color_buffer_at(color_buffer, x, y) = color;
                }
            }
            else {
                *psr_color_buffer_at(color_buffer, x, y) = color;
            }
        }
        while (_psr_line_step(&line_x));
    }
    while (_psr_line_step(&line_a) && _psr_line_step(&line_b));
}

void psr_raster_triangle_3d(psr_color_buffer_t* color_buffer, 
                            psr_depth_buffer_t* depth_buffer, 
                            psr_float3_t p0, 
                            psr_float3_t p1, 
                            psr_float3_t p2, 
                            const psr_attribute_array_t* attributes0,
                            const psr_attribute_array_t* attributes1,
                            const psr_attribute_array_t* attributes2,
                            int attribute_count,
                            psr_pixel_shader_callback pixel_shader, 
                            void* user_data) {
    p0.x = floorf(p0.x);
    p0.y = floorf(p0.y);
    p1.x = floorf(p1.x);
    p1.y = floorf(p1.y);
    p2.x = floorf(p2.x);
    p2.y = floorf(p2.y);

    // HACK: discard triangle if outside view.
    if (p0.x < 0 || p0.x >= color_buffer->w || p0.y < 0 || p0.y >= color_buffer->h ||
        p1.x < 0 || p1.x >= color_buffer->w || p1.y < 0 || p1.y >= color_buffer->h ||
        p2.x < 0 || p2.x >= color_buffer->w || p2.y < 0 || p2.y >= color_buffer->h) {
        return;
    }
    
    // Sort vertices by height.
    if (p0.y > p1.y) {
        PSR_SWAP(psr_float3_t, p0, p1);
        PSR_SWAP(const psr_attribute_array_t*, attributes0, attributes1);
    }
    if (p0.y > p2.y) {
        PSR_SWAP(psr_float3_t, p0, p2);
        PSR_SWAP(const psr_attribute_array_t*, attributes0, attributes2);
    }
    if (p1.y > p2.y) {
        PSR_SWAP(psr_float3_t, p1, p2);
        PSR_SWAP(const psr_attribute_array_t*, attributes1, attributes2);
    }

    // Create new attribute arrays with position at index 0.

    psr_attribute_array_t a0 = {.attributes[0] = PSR_ATTRIB_3(p0)};
    psr_attribute_array_t a1 = {.attributes[0] = PSR_ATTRIB_3(p1)};
    psr_attribute_array_t a2 = {.attributes[0] = PSR_ATTRIB_3(p2)};

    assert(attribute_count <= PSR_MAX_ATTRIBUTES - 1);
    for (int i = 0; i < attribute_count; i++) {
        a0.attributes[i + 1] = attributes0->attributes[i];
        a1.attributes[i + 1] = attributes1->attributes[i];
        a2.attributes[i + 1] = attributes2->attributes[i];
    }

    attribute_count += 1;

    // Check if the top of the triangle is flat.
    if (p0.y == p1.y) {
        _psr_raster_triangle_3d_flat_top_or_flat_bottom(color_buffer, 
                                                        depth_buffer, 
                                                        &a0, &a2,
                                                        &a1, &a2,
                                                        attribute_count,
                                                        pixel_shader,
                                                        user_data);
    }
    // Check if the bottom is flat.
    else if (p1.y == p2.y) {
        _psr_raster_triangle_3d_flat_top_or_flat_bottom(color_buffer, 
                                                        depth_buffer, 
                                                        &a0, &a1,
                                                        &a0, &a2,
                                                        attribute_count,
                                                        pixel_shader, 
                                                        user_data);
    }
    // Else split into two smaller triangles.
    else {
        psr_attribute_array_t a3;
        _psr_attributes_interpolate(&a3, &a0, &a2, attribute_count, (p1.y - p0.y) / (p2.y - p0.y));

        // Top (flat bottom).
        _psr_raster_triangle_3d_flat_top_or_flat_bottom(color_buffer, 
                                                        depth_buffer, 
                                                        &a0, &a1,
                                                        &a0, &a3,
                                                        attribute_count,
                                                        pixel_shader, 
                                                        user_data);
        
        // Bottom (flat top).
        _psr_raster_triangle_3d_flat_top_or_flat_bottom(color_buffer, 
                                                        depth_buffer,  
                                                        &a3, &a2,
                                                        &a1, &a2,
                                                        attribute_count,
                                                        pixel_shader,
                                                        user_data);
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

psr_image_t* psr_image_load_bmp(const char* path, psr_color_depth_t color_depth) {
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

    psr_byte_t* pixel_data = _psr_malloc(file_size - header_size);
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
    _psr_free(pixel_data);

    return image;
}

void psr_save_bmp(const char* path, const psr_color_buffer_t* color_buffer) {
    const psr_color_buffer_t* buf = color_buffer;
    
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

// Loads an entire file into a string. Remember to free return value.
char* _psr_load_file_binary_blob(const char* path) {
    FILE* file = fopen(path, "rb");
    if (!file) {
        return NULL;
    }
    fseek(file, 0, SEEK_END);
    int size = ftell(file);
    fseek(file, 0, SEEK_SET);
    char* content = _psr_malloc(size + 1);
    fread(content, 1, size, file);
    content[size] = '\0';
    fclose(file);
    return content;
}

psr_mesh_t* psr_mesh_load_obj(const char* path) {
    psr_mesh_t* mesh = _psr_malloc(sizeof(psr_mesh_t));

    mesh->face_count        = 0;
    mesh->position_count    = 0;
    mesh->tex_coord_count   = 0;
    mesh->normal_count      = 0;

    int face_max        = 128;
    int position_max    = 128;
    int tex_coord_max   = 128;
    int normal_max      = 128;

    mesh->faces         = _psr_malloc(face_max      * sizeof(psr_face_t));
    mesh->positions     = _psr_malloc(position_max  * sizeof(psr_float3_t));
    mesh->tex_coords    = _psr_malloc(tex_coord_max * sizeof(psr_float2_t));
    mesh->normals       = _psr_malloc(normal_max    * sizeof(psr_float3_t));

    // Load OBJ file.
    char* file_content = _psr_load_file_binary_blob(path);
    if (!file_content) {
        return NULL;
    }

    // Split into lines.
    int line_count;
    char** lines = psr_str_split(file_content, "\n", false, &line_count);

    for (int i = 0; i < line_count; i++) {
        // Split each line into tokens.
        int token_count;
        char** tokens = psr_str_split(lines[i], " ", false, &token_count);

        if (token_count == 0) {
            continue;
        }

        // Check the first token of each line.
        // The first token of relevant lines contains either 'v' for 
        // vertex position, 'vt' for vertex texture coordinate, 
        // 'vn' for vertex normal or 'f' for face.

#define _PUSH_BACK(type, nr_elements, count, max, array)    \
    assert(token_count >= 1 + nr_elements && "Bad model."); \
    if (count == max)                                       \
        array = realloc(array, (max *= 2) * sizeof(type));  \
    for (int _i = 0; _i < nr_elements; _i++)                \
        array[count].values[_i] = atof(tokens[_i + 1]);     \
    count++;

        if (strcmp(tokens[0], "v") == 0) {
            _PUSH_BACK(psr_float3_t, 3, mesh->position_count, position_max, mesh->positions);
        }
        else if (strcmp(tokens[0], "vt") == 0) {
            _PUSH_BACK(psr_float2_t, 2, mesh->tex_coord_count, tex_coord_max, mesh->tex_coords);
        }
        else if (strcmp(tokens[0], "vn") == 0) {
            _PUSH_BACK(psr_float3_t, 3, mesh->normal_count, normal_max, mesh->normals);
        }
        else if (strcmp(tokens[0], "f") == 0) {
            // Face lines contain 3 OR 4 additional tokens that describe 
            // which positions, texture coordinates and normals to use
            // for each face of the model.

            // Each face token contains 3 indices separated by '/'.
            // These indices point to which vertex position,
            // texture coordinate and normal to use.
            // For example, a face line with 4 additional tokens
            // might look like this:
            // f 2667/4489/9 642/718/10 2670/4490/11 6293/2381/12 

            int face_token_count = token_count - 1 > 4 ? 4 : token_count - 1;
            assert(face_token_count == 3 || face_token_count == 4);

            // Split each token into indices.

            int v_indices[4];
            int vt_indices[4];
            int vn_indices[4];
            
            for (int j = 0; j < face_token_count; j++) {
                int index_count;
                char** indices = psr_str_split(tokens[j + 1], "/", true, &index_count);
                assert(index_count == 3);

                // OBJ indices start at 1 and we want them to start at 0.
                v_indices[j] = atoi(indices[0]) - 1;
                vt_indices[j] = atoi(indices[1]) - 1;
                vn_indices[j] = atoi(indices[2]) - 1;

                psr_str_split_free(indices);
            }

            // Single face.
            if (face_token_count == 3) {
                psr_face_t face;

                for (int j = 0; j < 3; j++) {
                    face.position_indices[j] = v_indices[j];
                    face.tex_coord_indices[j] = vt_indices[j];
                    face.normal_indices[j] = vn_indices[j];
                }

                if (mesh->face_count == face_max) {
                    mesh->faces = realloc(mesh->faces, (face_max *= 2) * sizeof(psr_face_t));
                }
                mesh->faces[mesh->face_count++] = face;
            }
            // 2 faces.
            // A face could theoretically contain 4 vertices but we keep
            // it simple and break it up into 2 faces with 3 vertices each.
            else if (face_token_count == 4) {
                psr_face_t face_0;
                psr_face_t face_1;

                for (int j = 0; j < 3; j++) {
                    face_0.position_indices[j] = v_indices[j];
                    face_0.tex_coord_indices[j] = vt_indices[j];
                    face_0.normal_indices[j] = vn_indices[j];

                    face_1.position_indices[j] = v_indices[(j + 2) % 4];
                    face_1.tex_coord_indices[j] = vt_indices[(j + 2) % 4];
                    face_1.normal_indices[j] = vn_indices[(j + 2) % 4];
                }

                if (mesh->face_count >= face_max - 1) {
                    mesh->faces = realloc(mesh->faces, (face_max *= 2) * sizeof(psr_face_t));
                }
                mesh->faces[mesh->face_count++] = face_0;
                mesh->faces[mesh->face_count++] = face_1;
            }
        }

        psr_str_split_free(tokens);
    }

    psr_str_split_free(lines);
    _psr_free(file_content);

    return mesh;
}

void psr_mesh_free(psr_mesh_t* mesh) {
    _psr_free(mesh->positions);
    _psr_free(mesh->tex_coords);
    _psr_free(mesh->normals);
    _psr_free(mesh->faces);
    mesh->positions = NULL;
    mesh->tex_coords = NULL;
    mesh->normals = NULL;
    mesh->faces = NULL;
    _psr_free(mesh);
}

int _psr_parse_int(const char* str) {
    int count;
    char** tokens = psr_str_split(str, "=", false, &count);
    assert(count == 2);
    int parsed = atoi(tokens[1]);
    psr_str_split_free(tokens);
    return parsed;
}

psr_font_t* psr_font_load(psr_image_t* image, char* info_path) { 
    psr_font_t* font = _psr_malloc(sizeof(psr_font_t));
    
    font->image = image;
    font->char_info_count = 256;
    font->char_infos = _psr_malloc(256 * sizeof(_psr_char_info_t));
    font->size = -1;

    for (int i = 0; i < 256; i++) {
        font->char_infos[i].id = 0;
    }

    char* content = _psr_load_file_binary_blob(info_path);

    // Split into lines.
    int line_count;
    char** lines = psr_str_split(content, "\n", false, &line_count);

    for (int i = 0; i < line_count; i++) {
        // Split each line into words.
        int word_count;
        char** words = psr_str_split(lines[i], " ", false, &word_count);

        if (word_count == 0) {
            continue;
        }

        // Search for lines that start with either 'info' or 'char'.

        if (strcmp(words[0], "info") == 0) {
            font->size = _psr_parse_int(words[2]);
        }
        else if (strcmp(words[0], "char") == 0) {
            assert(word_count >= 9);

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

        psr_str_split_free(words);
    }

    psr_str_split_free(lines);
    _psr_free(content);
 
    return font;
}

void psr_font_destroy(psr_font_t* font) {
    _psr_free(font->char_infos);
    _psr_free(font);
}

char** psr_str_split(const char* str, const char* delim, bool count_consecutive_delimiters, int* out_count) {
    const int str_len = strlen(str);
	const int delim_len = strlen(delim);

	char** tokens = _psr_malloc((str_len / delim_len + 1) * sizeof(char*));
	int count = 0;

	const char* start = str;

	while (true) {
		// Find start of next delimiter (end of current token).
		const char* end = strstr(start, delim);

		// If there are more tokens.
		if (end) {
			int len = end - start;

            // Only copy token into array if delimiter is not consecutive
            // of if it is consecutive but flag is set to OK. 
            if (len > 0 || count_consecutive_delimiters) {
                char* token = _psr_malloc(len + 1);
                memcpy(token, start, len);
                token[len] = '\0';

                tokens[count++] = token;
            }
		}
		// If this is the last token.
		else if (start < str + str_len) {
			tokens[count++] = _strdup(start);
		}

		// Break if this was the last token.
		if (!end) {
			break;
		}

		// Advance starting point.
		start = end + delim_len;
	}
	tokens[count] = NULL;

	if (out_count) {
		*out_count = count;
	}

	return tokens;
}

void psr_str_split_free(char** tokens) {
	for (char** ptr = tokens; *ptr; ptr++) {
		_psr_free(*ptr);
	}
	_psr_free(tokens);
	tokens = NULL;
}
