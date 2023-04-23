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
    psr_int2_t start;
    psr_int2_t end;
    psr_int2_t current;
    
    psr_attribute_array_t attributes;
    psr_float4_t increments[PSR_MAX_ATTRIBUTES];
    int attribute_count;

    int steps;
    _psr_line_increment_type_t increment_type;

    // Data used with the bresenham algorithm.
    psr_int2_t delta;
    psr_int2_t dir;
    int error;
} _psr_line_t;

typedef struct {
    psr_float3_t actual_pos;
    psr_int2_t pixel_pos;
} _psr_3d_pos_t;

typedef struct {
    _psr_line_t line;

    psr_float3_t actual_start;
    psr_float3_t actual_end;
    psr_float3_t actual_current;

    psr_float3_t pos_increment;
} _psr_line_3d_t;

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

float psr_clamp(float f, float min, float max) {
    float t = f < min ? min : f;
    return t > max ? max : t;
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

psr_byte_t* psr_image_sample(psr_image_t* image, float u, float v) {
    return psr_image_at(image, u * image->w, v * image->h);
}

void _psr_line_init(_psr_line_t* line, 
                    psr_int2_t start,
                    psr_int2_t end,
                    const psr_attribute_array_t* start_attributes, 
                    const psr_attribute_array_t* end_attributes, 
                    int attribute_count,
                    _psr_line_increment_type_t increment_type) {
    line->start = start;
    line->end = end;
    line->current = start;
    
    line->attribute_count = attribute_count;
    line->increment_type = increment_type;

    if (start_attributes) {
        memcpy(&line->attributes, start_attributes, sizeof(psr_attribute_array_t));
    }

    // Calculate number of steps.

    psr_int2_t delta_abs = {_PSR_ABS(end.x - start.x), _PSR_ABS(end.y - start.y)};

    switch (increment_type) {
        case PSR_LINE_INCREMENT_ANY:        line->steps = (delta_abs.x >= delta_abs.y) ? delta_abs.x : delta_abs.y; break;
        case PSR_LINE_INCREMENT_VERTICAL:   line->steps = delta_abs.y; break;
        case PSR_LINE_INCREMENT_HORIZONTAL: line->steps = delta_abs.x; break;
        default:                            assert(!"Unhandled enum.");
    }

    // Setup increments.

    if (line->steps > 0) {
        for (int i = 0; i < attribute_count; i++) {
            for (int j = 0; j < line->attributes.attributes[i].nr_of_floats; j++) {
                float delta = end_attributes->attributes[i].data.values[j] - start_attributes->attributes[i].data.values[j];
                line->increments[i].values[j] = delta / line->steps;
            }
        }
    }

    // Setup bresenham members.

    line->delta = (psr_int2_t){delta_abs.x, -delta_abs.y};
    line->dir = (psr_int2_t){start.x < end.x ? 1 : -1, start.y < end.y ? 1 : -1};
    line->error = line->delta.x + line->delta.y;
}

void _psr_line_3d_init(_psr_line_3d_t* line,
                       _psr_3d_pos_t start,
                       _psr_3d_pos_t end,
                       const psr_attribute_array_t* start_attributes,
                       const psr_attribute_array_t* end_attributes,
                       int attribute_count,
                       _psr_line_increment_type_t increment_type) {
    _psr_line_init(&line->line, start.pixel_pos, end.pixel_pos, start_attributes, end_attributes, attribute_count, increment_type);

    line->actual_start = start.actual_pos;
    line->actual_end = end.actual_pos;
    line->actual_current = start.actual_pos;

    // Setup position increment.

    for (int i = 0; i < 3; i++) {
        line->pos_increment.values[i] = (end.actual_pos.values[i] - start.actual_pos.values[i]) / line->line.steps;
    }
}

void _psr_attribute_array_interpolate(psr_attribute_array_t* out, const psr_attribute_array_t* a, const psr_attribute_array_t* b, int attribute_count, float amount) {
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
    psr_int2_t prev = line->current;

    while (true) {
        if (line->current.x == line->end.x && line->current.y == line->end.y) {
            return false;
        }

        // Bresenham algorithm.
        int e2 = line->error << 1;
        if (e2 >= line->delta.y) {
            line->error += line->delta.y;
            line->current.x += line->dir.x;
        }
        if (e2 <= line->delta.x) {
            line->error += line->delta.x;
            line->current.y += line->dir.y;
        }

        // Step until difference in specified coordinate.
        // If for example the type of increment is set to vertical,
        // we must step until there is a difference in Y.

        if (line->increment_type == PSR_LINE_INCREMENT_HORIZONTAL) {
            if (prev.x != line->current.x) {
                break;
            }
        }
        else if (line->increment_type == PSR_LINE_INCREMENT_VERTICAL) {
            if (prev.y != line->current.y) {
                break;
            }
        }
        else {
            break;
        }
    }

    // Increment attributes.
    for (int i = 0; i < line->attribute_count; i++) {
        for (int j = 0; j < line->attributes.attributes[i].nr_of_floats; j++) {
            line->attributes.attributes[i].data.values[j] += line->increments[i].values[j];
        }
    }

    return true;
}

bool _psr_line_3d_step(_psr_line_3d_t* line) {
    if (_psr_line_step(&line->line)) {
        PSR_ADD(line->actual_current, line->actual_current, line->pos_increment, 3);
        return true;
    }
    return false;
}

void _psr_raster_triangle_3d_between_two_vertical_lines(psr_color_buffer_t* color_buffer, 
                                                        psr_depth_buffer_t* depth_buffer,
                                                        _psr_line_3d_t line_a,
                                                        _psr_line_3d_t line_b,
                                                        psr_pixel_shader_callback pixel_shader, 
                                                        void* user_data) {    
    if (line_a.line.start.x > line_b.line.start.x || line_a.line.end.x > line_b.line.end.x) {
        PSR_SWAP(_psr_line_3d_t, line_a, line_b);
    }
    
    do {
        assert(line_a.line.current.y == line_b.line.current.y);

        // Create a horizontal line between line a and line b.
        _psr_line_3d_t line_x;
        _psr_line_3d_init(&line_x, 
                          (_psr_3d_pos_t){.pixel_pos = line_a.line.current, .actual_pos = line_a.actual_current}, 
                          (_psr_3d_pos_t){.pixel_pos = line_b.line.current, .actual_pos = line_b.actual_current}, 
                          &line_a.line.attributes, 
                          &line_b.line.attributes, 
                          line_a.line.attribute_count, 
                          PSR_LINE_INCREMENT_HORIZONTAL);

        do {
            int x = line_x.line.current.x;
            int y = line_x.line.current.y;
            
            // Update depth buffer.
            if (depth_buffer) {
                assert(x >= 0 && x < depth_buffer->w && y >= 0 && y < depth_buffer->h);

                float z = line_x.actual_current.z;

                if (z < *psr_depth_buffer_at(depth_buffer, x, y)) {
                    *psr_depth_buffer_at(depth_buffer, x, y) = z;
                }
                // If pixel is invisible, skip color buffer update.
                else {
                    continue;
                }
            }

            // Update color buffer.
            if (color_buffer) {
                assert(x >= 0 && x < color_buffer->w && y >= 0 && y < color_buffer->h);

                // Invoke pixel shader.
                psr_byte3_t color = 
                    pixel_shader(line_x.line.current, &line_x.line.attributes, user_data);

                *psr_color_buffer_at(color_buffer, x, y) = color;
            }
        }
        while (_psr_line_3d_step(&line_x));
    }
    while (_psr_line_3d_step(&line_a) && _psr_line_3d_step(&line_b));
}

void _psr_raster_triangle_3d_flat_top(psr_color_buffer_t* color_buffer, 
                                      psr_depth_buffer_t* depth_buffer,
                                      _psr_3d_pos_t top_left_pos,
                                      _psr_3d_pos_t top_right_pos,
                                      _psr_3d_pos_t bottom_pos,
                                      const psr_attribute_array_t* top_left_attribs,
                                      const psr_attribute_array_t* top_right_attribs,
                                      const psr_attribute_array_t* bottom_attribs,
                                      int attribute_count,
                                      psr_pixel_shader_callback pixel_shader, 
                                      void* user_data) {
    assert(top_left_pos.pixel_pos.y == top_right_pos.pixel_pos.y);
    
    _psr_line_3d_t line_a;
    _psr_line_3d_init(&line_a, top_left_pos, bottom_pos, top_left_attribs, bottom_attribs, 
        attribute_count, PSR_LINE_INCREMENT_VERTICAL);
    
    _psr_line_3d_t line_b;
    _psr_line_3d_init(&line_b, top_right_pos, bottom_pos, top_right_attribs, bottom_attribs, 
        attribute_count, PSR_LINE_INCREMENT_VERTICAL);

    _psr_raster_triangle_3d_between_two_vertical_lines(color_buffer, depth_buffer, 
        line_a, line_b, pixel_shader, user_data);
}

void _psr_raster_triangle_3d_flat_bottom(psr_color_buffer_t* color_buffer, 
                                         psr_depth_buffer_t* depth_buffer,
                                         _psr_3d_pos_t top_pos,
                                         _psr_3d_pos_t bottom_left_pos,
                                         _psr_3d_pos_t bottom_right_pos,
                                         const psr_attribute_array_t* top_attribs,
                                         const psr_attribute_array_t* bottom_left_attribs,
                                         const psr_attribute_array_t* bottom_right_attribs,
                                         int attribute_count,
                                         psr_pixel_shader_callback pixel_shader, 
                                         void* user_data) {
    assert(bottom_left_pos.pixel_pos.y == bottom_right_pos.pixel_pos.y);
    
    _psr_line_3d_t line_a;
    _psr_line_3d_init(&line_a, top_pos, bottom_left_pos, top_attribs, bottom_left_attribs, 
        attribute_count, PSR_LINE_INCREMENT_VERTICAL);
    
    _psr_line_3d_t line_b;
    _psr_line_3d_init(&line_b, top_pos, bottom_right_pos, top_attribs, bottom_right_attribs, 
        attribute_count, PSR_LINE_INCREMENT_VERTICAL);

    _psr_raster_triangle_3d_between_two_vertical_lines(color_buffer, depth_buffer, 
        line_a, line_b, pixel_shader, user_data);
}

void psr_raster_triangle_3d(psr_color_buffer_t* color_buffer, 
                            psr_depth_buffer_t* depth_buffer, 
                            psr_float3_t pos0, 
                            psr_float3_t pos1, 
                            psr_float3_t pos2, 
                            const psr_attribute_array_t* attributes0,
                            const psr_attribute_array_t* attributes1,
                            const psr_attribute_array_t* attributes2,
                            int attribute_count,
                            psr_pixel_shader_callback pixel_shader, 
                            void* user_data) {
    // Floor X and Y, otherwise there's missing pixel artifacts.
    for (int i = 0; i < 2; i++) {
        pos0.values[i] = floorf(pos0.values[i]);
        pos1.values[i] = floorf(pos1.values[i]);
        pos2.values[i] = floorf(pos2.values[i]);
    }
    
    _psr_3d_pos_t p0 = {.actual_pos = pos0, .pixel_pos = {pos0.x, pos0.y}};
    _psr_3d_pos_t p1 = {.actual_pos = pos1, .pixel_pos = {pos1.x, pos1.y}};
    _psr_3d_pos_t p2 = {.actual_pos = pos2, .pixel_pos = {pos2.x, pos2.y}};

    // HACK: discard triangle if outside view.
    if (p0.pixel_pos.x < 0 || p0.pixel_pos.x >= color_buffer->w || p0.pixel_pos.y < 0 || p0.pixel_pos.y >= color_buffer->h ||
        p1.pixel_pos.x < 0 || p1.pixel_pos.x >= color_buffer->w || p1.pixel_pos.y < 0 || p1.pixel_pos.y >= color_buffer->h ||
        p2.pixel_pos.x < 0 || p2.pixel_pos.x >= color_buffer->w || p2.pixel_pos.y < 0 || p2.pixel_pos.y >= color_buffer->h) {
        return;
    }
    
    // Sort vertices by height.
    if (p0.pixel_pos.y > p1.pixel_pos.y) {
        PSR_SWAP(_psr_3d_pos_t, p0, p1);
        PSR_SWAP(const psr_attribute_array_t*, attributes0, attributes1);
    }
    if (p0.pixel_pos.y > p2.pixel_pos.y) {
        PSR_SWAP(_psr_3d_pos_t, p0, p2);
        PSR_SWAP(const psr_attribute_array_t*, attributes0, attributes2);
    }
    if (p1.pixel_pos.y > p2.pixel_pos.y) {
        PSR_SWAP(_psr_3d_pos_t, p1, p2);
        PSR_SWAP(const psr_attribute_array_t*, attributes1, attributes2);
    }

    // Check if the top of the triangle is flat.
    if (p0.pixel_pos.y == p1.pixel_pos.y) {
        _psr_raster_triangle_3d_flat_top(color_buffer, depth_buffer, 
            p0, p1, p2, attributes0, attributes1, attributes2, attribute_count,
            pixel_shader, user_data);
    }
    // Check if the bottom is flat.
    else if (p1.pixel_pos.y == p2.pixel_pos.y) {
        _psr_raster_triangle_3d_flat_bottom(color_buffer, depth_buffer,
            p0, p1, p2, attributes0, attributes1, attributes2, attribute_count, 
            pixel_shader, user_data);
    }
    // Else split into two smaller triangles.
    else {
        float alpha_split = 
            (p1.actual_pos.y - p0.actual_pos.y) / (p2.actual_pos.y - p0.actual_pos.y);

        _psr_3d_pos_t p3;
        p3.actual_pos.y = p1.actual_pos.y;
        p3.actual_pos.x = psr_lerp(p0.actual_pos.x, p2.actual_pos.x, alpha_split);
        p3.actual_pos.z = psr_lerp(p0.actual_pos.z, p2.actual_pos.z, alpha_split);
        p3.pixel_pos = (psr_int2_t){(int)roundf(p3.actual_pos.x), (int)roundf(p3.actual_pos.y)};

        psr_attribute_array_t attributes3;
        _psr_attribute_array_interpolate(&attributes3, attributes0, attributes2, attribute_count, alpha_split);

        // Bottom (flat top).
        _psr_raster_triangle_3d_flat_top(color_buffer, depth_buffer,
            p1, p3, p2, attributes1, &attributes3, attributes2, attribute_count, 
            pixel_shader, user_data);
        
        // Top (flat bottom).
        _psr_raster_triangle_3d_flat_bottom(color_buffer, depth_buffer,
            p0, p1, p3, attributes0, attributes1, &attributes3, attribute_count, 
            pixel_shader, user_data);
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
