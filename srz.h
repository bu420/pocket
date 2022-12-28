/*
software-rasterizer
https://github.com/bu420/software-rasterizer

Define 'SRZ_IMPLEMENTATION' once.
Define 'SRZ_SAVE_AND_LOAD' to enable saving to and loading from disk.

Example:
#define SRZ_IMPLEMENTATION
#define SRZ_SAVE_AND_LOAD
#include <srz.h>

Credits to 'el nora' for sine and cosine implementations.
*/

#pragma once

#define SRZ_PI 3.1415926f
#ifdef __cplusplus
#define SRZ_NULL nullptr
#else
#define SRZ_NULL ((void*)0)
#endif

typedef unsigned char srz_byte_t;

typedef union {
    struct {
        srz_byte_t x, y, z;
    };
    struct {
        srz_byte_t r, g, b;
    };
    srz_byte_t values[3];
} srz_byte3_t;

typedef union {
    struct {
        srz_byte_t x, y, z, w;
    };
    struct {
        srz_byte_t r, g, b, a;
    };
    srz_byte_t values[4];
} srz_byte4_t;

typedef union {
    struct {
        int x, y;
    };
    struct {
        int r, g;
    };
    int values[2];
} srz_int2_t;

typedef union {
    struct {
        float x, y;
    };
    struct {
        float r, g;
    };
    float values[2];
} srz_float2_t;

typedef union {
    struct {
        float x, y, z;
    };
    struct {
        float r, g, b;
    };
    float values[3];
} srz_float3_t;

typedef union {
    struct {
        float x, y, z, w;
    };
    struct {
        float r, g, b, a;
    };
    float values[4];
} srz_float4_t;

typedef struct {
    float m[4][4];
} srz_matrix_t;

typedef struct {
    int w, h;
    srz_byte3_t* data;
} srz_color_buffer_t;

typedef struct {
    int w, h;
    float* data;
} srz_depth_buffer_t;

typedef struct {
    int w, h;
    srz_byte4_t* data;
} srz_texture_t;

#ifdef SRZ_SAVE_AND_LOAD
// Vertex, texture coordinate and normal indices.
typedef struct {
    int position_indices[3];
    // -1 means there is no texture coordinate associated.
    int tex_coord_indices[3];
    int normal_indices[3];
} srz_face_t;

typedef struct {
    srz_float3_t* positions;
    srz_float2_t* tex_coords;
    srz_float3_t* normals;
    srz_face_t* faces;
    int position_count;
    int tex_coord_count;
    int normal_count;
    int face_count;
} srz_mesh_t;
#endif

int srz_max(int a, int b);
int srz_min(int a, int b);
int srz_max3(int a, int b, int c);
int srz_min3(int a, int b, int c);
int srz_abs(int n);
float srz_maxf(float a, float b);
float srz_minf(float a, float b);
float srz_absf(float n);
float srz_floorf(float n);
float srz_ceilf(float n);
float srz_roundf(float n);
float srz_sqrtf(float n);
float srz_fmodf(float a, float b);
float srz_sinf(float x);
float srz_cosf(float x);
float srz_tanf(float n);
srz_float3_t srz_normalize(srz_float3_t f);
srz_float3_t srz_cross(srz_float3_t a, srz_float3_t b);
float srz_dot(srz_float3_t a, srz_float3_t b);

void srz_matrix_init_zero(srz_matrix_t* matrix);
void srz_matrix_init_identity(srz_matrix_t* matrix);
srz_matrix_t srz_matrix_mul(srz_matrix_t a, srz_matrix_t b);
srz_matrix_t srz_matrix_translate(srz_matrix_t matrix, srz_float3_t f);
srz_matrix_t srz_matrix_rotate_x(srz_matrix_t matrix, float a);
srz_matrix_t srz_matrix_rotate_y(srz_matrix_t matrix, float a);
srz_matrix_t srz_matrix_rotate_z(srz_matrix_t matrix, float a);
srz_matrix_t srz_matrix_scale(srz_matrix_t matrix, srz_float3_t f);
srz_float3_t srz_matrix_mul_float3(srz_matrix_t matrix, srz_float3_t f);
srz_float4_t srz_matrix_mul_float4(srz_matrix_t matrix, srz_float4_t f);
srz_float3_t srz_float3_mul_matrix(srz_float3_t f, srz_matrix_t matrix);
srz_float4_t srz_float4_mul_matrix(srz_float4_t f, srz_matrix_t matrix);
srz_matrix_t srz_look_at(srz_float3_t pos, srz_float3_t target, srz_float3_t up);
srz_matrix_t srz_perspective(float aspect, float fov, float near, float far);

srz_byte3_t* srz_color_buffer_at(srz_color_buffer_t* color_buffer, int x, int y);
void srz_color_buffer_clear(srz_color_buffer_t* color_buffer, srz_byte3_t color);
float* srz_depth_buffer_at(srz_depth_buffer_t* depth_buffer, int x, int y);
void srz_depth_buffer_clear(srz_depth_buffer_t* depth_buffer);
srz_byte4_t* srz_texture_at(srz_texture_t* texture, int x, int y);
// @param depth_buffer pass NULL if depth buffering should not be used.
void srz_raster_line(srz_color_buffer_t* color_buffer, srz_depth_buffer_t* depth_buffer, srz_int2_t a, srz_int2_t b, srz_byte3_t color);
// @param depth_buffer pass NULL if depth buffering should not be used.
void srz_raster_triangle(srz_color_buffer_t* color_buffer, srz_depth_buffer_t* depth_buffer, srz_int2_t a, srz_int2_t b, srz_int2_t c, srz_byte3_t color);
#ifdef SRZ_SAVE_AND_LOAD
void srz_save_bmp(char const* filename, srz_color_buffer_t color_buffer);
// Simple OBJ loader that WILL break on anything complex format-wise.
// @return 0 on failure to open file and -1 on bad model.
int srz_load_obj(char const* filename, srz_mesh_t* mesh);
void srz_mesh_free(srz_mesh_t* mesh);
#endif

#ifdef SRZ_IMPLEMENTATION
#ifdef SRZ_SAVE_AND_LOAD
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#endif

int srz_max(int a, int b) {
    return a > b ? a : b;
}

int srz_min(int a, int b) {
    return a < b ? a : b;
}

int srz_max3(int a, int b, int c) {
    return srz_max(srz_max(a, b), c);
}

int srz_min3(int a, int b, int c) {
    return srz_min(srz_min(a, b), c);
}

int srz_abs(int n) {
    return n >= 0 ? n : -n;
}

float srz_maxf(float a, float b) {
    return a > b ? a : b;
}

float srz_minf(float a, float b) {
    return a < b ? a : b;
}

float srz_absf(float n) {
    return n >= 0 ? n : -n;
}

float srz_floorf(float n) {
    return n >= 0 ? (int)n : (int)n - 1;
}

float srz_ceilf(float n) {
    return -srz_floorf(-n);
}

float srz_roundf(float n) {
    return (n >= 0) ? srz_floorf(n + .5f) : srz_ceilf(n - .5f);
}

float srz_sqrtf(float n) {
    float lo = srz_minf(1, n);
    float hi = srz_maxf(1, n);
    float mid;

    while (100 * lo * lo < n) {
        lo *= 10;
    }
    while (0.01 * hi * hi > n) {
        hi *= 0.1;
    }

    for (int i = 0; i < 100; ++i) {
        mid = (lo + hi) / 2;
        if (mid * mid == n) {
            return mid;
        }

        if (mid * mid > n) {
            hi = mid;
        }
        else {
            lo = mid;
        }
    }
    return mid;
}

float srz_fmodf(float a, float b) {
    return a - (srz_roundf(a / b) * b);
}

float _srz_sf(float x) {
    float a3 = -53272705 / 360869676.0f;
    float a5 = 38518909 / 7217393520.0f;
    float a7 = -269197963 / 3940696861920.0f;
    float a9 = 4585922449 / 15605159573203200.0f;
    float b2 = 2290747 / 120289892.0f;
    float b4 = 1281433 / 7217393520.0f;
    float b6 = 560401 / 562956694560.0f;
    float b8 = 1029037 / 346781323848960.0f;
    float t2 = x * x;
    float t4 = t2 * t2;
    float t6 = t2 * t4;
    float t8 = t4 * t4;
    float n = (1 + a3 * t2 + a5 * t4 + a7 * t6 + a9 * t8) * x;
    float d = 1 + b2 * t2 + b4 * t4 + b6 * t6 + b8 * t8;
    return n / d;
}

float _srz_cf(float x) {
    float a2 = -260735 / 545628.0f;
    float a4 = 4375409 / 141863280.0f;
    float a6 = -7696415 / 13108167072.0f;
    float a8 = 80737373 / 23594700729600.0f;
    float b2 = 12079 / 545628.0f;
    float b4 = 34709 / 141863280.0f;
    float b6 = 109247 / 65540835360.0f;
    float b8 = 11321 / 1814976979200.0f;
    float t2 = x * x;
    float t4 = t2 * t2;
    float t6 = t2 * t4;
    float t8 = t4 * t4;
    float n = 1 + a2 * t2 + a4 * t4 + a6 * t6 + a8 * t8;
    float d = 1 + b2 * t2 + b4 * t4 + b6 * t6 + b8 * t8;
    return n / d;
}

float srz_sinf(float x) {
    int i = x / SRZ_PI / 2;
    float a = x - 2 * SRZ_PI * i;
    if (0 <= a && a < SRZ_PI / 4) {
        return _srz_sf(a);
    }
    if (SRZ_PI / 4 <= a && a < SRZ_PI / 2) {
        return _srz_cf(SRZ_PI / 2 - a);
    }
    if (SRZ_PI / 2 <= a && a < 3 * SRZ_PI / 4) {
        return _srz_cf(a - SRZ_PI / 2);
    }
    if (3 * SRZ_PI / 4 <= a && a < SRZ_PI) {
        return _srz_sf(SRZ_PI - a);
    }
    return -srz_sinf(a - SRZ_PI);
}

float srz_cosf(float x) {
    return srz_sinf(x + SRZ_PI / 2);
}

float srz_tanf(float n) {
    return srz_sinf(n) / srz_cosf(n);
}

srz_float3_t srz_normalize(srz_float3_t f) {
    float len = srz_sqrtf(f.x * f.x + f.y * f.y + f.z * f.z);
    srz_float3_t result = {f.x / len, f.y / len, f.z / len};
    return result;
}

srz_float3_t srz_cross(srz_float3_t a, srz_float3_t b) {
    srz_float3_t result = {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
    return result;
}

float srz_dot(srz_float3_t a, srz_float3_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

void srz_matrix_init_zero(srz_matrix_t* matrix) {
    for (int x = 0; x < 4; ++x) {
        for (int y = 0; y < 4; ++y) {
            matrix->m[x][y] = 0;
        }
    }
}

void srz_matrix_init_identity(srz_matrix_t* matrix) {
    srz_matrix_init_zero(matrix);
    for (int i = 0; i < 4; ++i) {
        matrix->m[i][i] = 1;
    }
}

srz_matrix_t srz_matrix_mul(srz_matrix_t a, srz_matrix_t b) {
    srz_matrix_t result;
    srz_matrix_init_zero(&result);

    for (int x = 0; x < 4; ++x) {
        for (int y = 0; y < 4; ++y) {
            for (int i = 0; i < 4; ++i) {
                result.m[x][y] += a.m[x][i] * b.m[i][y];
            }
        }
    }
    return result;
}

srz_matrix_t srz_matrix_translate(srz_matrix_t matrix, srz_float3_t f) {
    srz_matrix_t m;
    srz_matrix_init_identity(&m);

    m.m[3][0] = f.x;
    m.m[3][1] = f.y;
    m.m[3][2] = f.z;

    return srz_matrix_mul(matrix, m);
}

srz_matrix_t srz_matrix_rotate_x(srz_matrix_t matrix, float a) {
    float s = srz_sinf(a);
    float c = srz_cosf(a);

    srz_matrix_t x;
    srz_matrix_init_identity(&x);

    x.m[1][1] = c;
    x.m[1][2] = -s;
    x.m[2][1] = s;
    x.m[2][2] = c;

    return srz_matrix_mul(matrix, x);
}

srz_matrix_t srz_matrix_rotate_y(srz_matrix_t matrix, float a) {
    float s = srz_sinf(a);
    float c = srz_cosf(a);

    srz_matrix_t y;
    srz_matrix_init_identity(&y);

    y.m[0][0] = c;
    y.m[0][2] = s;
    y.m[2][0] = -s;
    y.m[2][2] = c;

    return srz_matrix_mul(matrix, y);
}

srz_matrix_t srz_matrix_rotate_z(srz_matrix_t matrix, float a) {
    float s = srz_sinf(a);
    float c = srz_cosf(a);

    srz_matrix_t z;
    srz_matrix_init_identity(&z);

    z.m[0][0] = c;
    z.m[0][1] = -s;
    z.m[1][0] = s;
    z.m[1][1] = c;

    return srz_matrix_mul(matrix, z);
}

srz_matrix_t srz_matrix_scale(srz_matrix_t matrix, srz_float3_t f) {
    srz_matrix_t m;
    srz_matrix_init_identity(&m);

    m.m[0][0] = f.x;
    m.m[1][1] = f.y;
    m.m[2][2] = f.z;

    return srz_matrix_mul(matrix, m);
}

srz_float3_t srz_matrix_mul_float3(srz_matrix_t matrix, srz_float3_t f) {
    srz_float4_t f4;
    for (int x = 0; x < 4; ++x) {
        float i = 0;
        for (int y = 0; y < 3; ++y) {
            i += matrix.m[x][y] * f.values[y];
        }
        f4.values[x] = i + matrix.m[x][3];
    }
    
    if (f4.w != 0) {
        f4.x /= f4.w;
        f4.y /= f4.w;
        f4.z /= f4.w;
    }

    srz_float3_t result = {f4.x, f4.y, f4.z};
    return result;
}

srz_float4_t srz_matrix_mul_float4(srz_matrix_t matrix, srz_float4_t f) {
    srz_float4_t result;
    for (int x = 0; x < 4; ++x) {
        float i = 0;
        for (int y = 0; y < 4; ++y) {
            i += matrix.m[x][y] * f.values[y];
        }
        result.values[x] = i;
    }
    return result;
}

srz_float3_t srz_float3_mul_matrix(srz_float3_t f, srz_matrix_t matrix) {
    srz_float4_t f4;
    for (int y = 0; y < 4; ++y) {
        float i = 0;
        for (int x = 0; x < 3; ++x) {
            i += matrix.m[x][y] * f.values[x];
        }
        f4.values[y] = i + matrix.m[3][y];
    }
    
    if (f4.w != 0) {
        f4.x /= f4.w;
        f4.y /= f4.w;
        f4.z /= f4.w;
    }

    srz_float3_t result = {f4.x, f4.y, f4.z};
    return result;
}

srz_float4_t srz_float4_mul_matrix(srz_float4_t f, srz_matrix_t matrix) {
    srz_float4_t result;
    for (int y = 0; y < 4; ++y) {
        float i = 0;
        for (int x = 0; x < 4; ++x) {
            i += matrix.m[x][y] * f.values[x];
        }
        result.values[y] = i;
    }
    return result;
}

srz_matrix_t srz_look_at(srz_float3_t pos, srz_float3_t target, srz_float3_t up) {
    srz_float3_t diff = {target.x - pos.x, target.y - pos.y, target.z - pos.z};
    srz_float3_t forward = srz_normalize(diff);
    srz_float3_t right = srz_normalize(srz_cross(forward, up));
    srz_float3_t local_up = srz_normalize(srz_cross(right, forward));

    srz_matrix_t m;
    srz_matrix_init_identity(&m);
    m.m[0][0] = right.x;
    m.m[1][0] = right.y;
    m.m[2][0] = right.z;
    m.m[0][1] = local_up.x;
    m.m[1][1] = local_up.y;
    m.m[2][1] = local_up.z;
    m.m[0][2] = -forward.x;
    m.m[1][2] = -forward.y;
    m.m[2][2] = -forward.z;
    m.m[3][0] = -srz_dot(right, pos);
    m.m[3][1] = -srz_dot(local_up, pos);
    m.m[3][2] = srz_dot(forward, pos);
    return m;
}

srz_matrix_t srz_perspective(float aspect, float fov, float near, float far) {    
    srz_matrix_t m;
    srz_matrix_init_zero(&m);

    float half_tan = srz_tanf(fov / 2);

    m.m[0][0] = 1 / (half_tan * aspect);
    m.m[1][1] = 1 / half_tan;
    m.m[2][2] = -(far + near) / (far - near);
    m.m[2][3] = -1;
    m.m[3][2] = -(2 * far * near) / (far - near);
    return m;
}

srz_byte3_t* srz_color_buffer_at(srz_color_buffer_t* color_buffer, int x, int y) {
    return &color_buffer->data[y * color_buffer->w + x];
}

void srz_color_buffer_clear(srz_color_buffer_t* color_buffer, srz_byte3_t color) {
    for (int i = 0; i < color_buffer->w * color_buffer->h; ++i) {
        color_buffer->data[i] = color;
    }
}

float* srz_depth_buffer_at(srz_depth_buffer_t* depth_buffer, int x, int y) {
    return &depth_buffer->data[y * depth_buffer->w + x];
}

void srz_depth_buffer_clear(srz_depth_buffer_t* depth_buffer) {
    for (int i = 0; i < depth_buffer->w * depth_buffer->h; ++i) {
        depth_buffer->data[i] = 0;
    }
}

srz_byte4_t* srz_texture_at(srz_texture_t* texture, int x, int y) {
    return &texture->data[y * texture->w + x];
}

void srz_raster_line(srz_color_buffer_t* color_buffer, srz_depth_buffer_t* depth_buffer, srz_int2_t a, srz_int2_t b, srz_byte3_t color) {
    if (!depth_buffer) {
        // No depth buffer used.
    }
    
    srz_int2_t diff = {b.x - a.x, b.y - a.y};
    srz_int2_t diff_abs = {srz_abs(diff.x), srz_abs(diff.y)};
    srz_int2_t p = {2 * diff_abs.y - diff_abs.x, 2 * diff_abs.x - diff_abs.y};
    srz_int2_t pos;
    srz_int2_t e;

    if (diff_abs.y <= diff_abs.x) {
        if (diff.x >= 0) {
            pos.x = a.x;
            pos.y = a.y;
            e.x = b.x;
        }
        else {
            pos.x = b.x;
            pos.y = b.y;
            e.x = a.x;
        }

        int idx = pos.y * color_buffer->w + pos.x;
        if (idx >= 0 && idx < color_buffer->w * color_buffer->h) {
            *srz_color_buffer_at(color_buffer, pos.x, pos.y) = color;
        }

        for (int i = 0; pos.x < e.x; i++) {
            pos.x++;
            if (p.x < 0) {
                p.x += 2 * diff_abs.y;
            }
            else {
                if ((diff.x < 0 && diff.y < 0) || (diff.x > 0 && diff.y > 0)) {
                    pos.y++;
                }
                else {
                    pos.y--;
                }
                p.x += 2 * (diff_abs.y - diff_abs.x);
            }

            int idx = pos.y * color_buffer->w + pos.x;
            if (idx >= 0 && idx < color_buffer->w * color_buffer->h) {
                *srz_color_buffer_at(color_buffer, pos.x, pos.y) = color;
            }
        }
    }
    else {
        if (diff.y >= 0) {
            pos.x = a.x;
            pos.y = a.y;
            e.y = b.y;
        }
        else {
            pos.x = b.x;
            pos.y = b.y;
            e.y = a.y;
        }

        int idx = pos.y * color_buffer->w + pos.x;
        if (idx >= 0 && idx < color_buffer->w * color_buffer->h) {
            *srz_color_buffer_at(color_buffer, pos.x, pos.y) = color;
        }

        for (int i = 0; pos.y < e.y; i++) {
            pos.y++;
            if (p.y <= 0) {
                p.y += 2 * diff_abs.x;
            }
            else {
                if ((diff.x < 0 && diff.y < 0) || (diff.x > 0 && diff.y > 0)) {
                    pos.x++;
                }
                else {
                    pos.x--;
                }
                p.y += 2 * (diff_abs.x - diff_abs.y);
            }

            int idx = pos.y * color_buffer->w + pos.x;
            if (idx >= 0 && idx < color_buffer->w * color_buffer->h) {
                *srz_color_buffer_at(color_buffer, pos.x, pos.y) = color;
            }
        }
    }
}

int _srz_t(srz_int2_t a, srz_int2_t b, srz_int2_t c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

void srz_raster_triangle(srz_color_buffer_t* color_buffer, srz_depth_buffer_t* depth_buffer, srz_int2_t a, srz_int2_t b, srz_int2_t c, srz_byte3_t color) {
    if (!depth_buffer) {
        // No depth buffer used.
    }
    
    // Bounding box.
    int min_x = srz_min3(a.x, b.x, c.x);
    int min_y = srz_min3(a.y, b.y, c.y);
    int max_x = srz_max3(a.x, b.x, c.x);
    int max_y = srz_max3(a.y, b.y, c.y);

    // Clip.
    min_x = srz_max(min_x, 0);
    min_y = srz_max(min_y, 0);
    max_x = srz_min(max_x, color_buffer->w - 1);
    max_y = srz_min(max_y, color_buffer->h - 1);

    // Rasterize.
    srz_int2_t p;
    for (p.x = min_x; p.x <= max_x; ++p.x) {
        for (p.y = min_y; p.y <= max_y; ++p.y) {
            int t0 = _srz_t(b, c, p);
            int t1 = _srz_t(c, a, p);
            int t2 = _srz_t(a, b, p);

            if (t0 >= 0 && t1 >= 0 && t2 >= 0) {
                *srz_color_buffer_at(color_buffer, p.x, p.y) = color;
            }
        }
    }
}

#ifdef SRZ_SAVE_AND_LOAD
void srz_save_bmp(char const* filename, srz_color_buffer_t color_buffer) {
    srz_byte_t padding[] = {0, 0, 0};
    int padding_size = (4 - color_buffer.w * 3 % 4) % 4;
    int stride = color_buffer.w * 3 + padding_size;
    int file_size = 54 + stride * color_buffer.h;

    srz_byte_t header[54] = {
        'B', 'M',
        (srz_byte_t)file_size, (srz_byte_t)(file_size >> 8), (srz_byte_t)(file_size >> 16), (srz_byte_t)(file_size >> 24),
        0, 0, 0, 0,
        54, 0, 0, 0,
        40, 0, 0, 0,
        (srz_byte_t)color_buffer.w, (srz_byte_t)(color_buffer.w >> 8), (srz_byte_t)(color_buffer.w >> 16), (srz_byte_t)(color_buffer.w >> 24),
        (srz_byte_t)color_buffer.h, (srz_byte_t)(color_buffer.h >> 8), (srz_byte_t)(color_buffer.h >> 16), (srz_byte_t)(color_buffer.h >> 24),
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
            srz_byte3_t pixel = color_buffer.data[i];

            // RBG to BGR.
            srz_byte3_t color = {.r = pixel.b, .g = pixel.g, .b = pixel.r};

            fwrite(&color, 3, 1, file);
        }

        fwrite(padding, 1, padding_size, file);
    }

    fclose(file);
}

int _srz_split(char* str, char const* delims, char** tokens, int max) {
    int count = 0;
    char* token = strtok(str, delims);
    for (; count < max && token; ++count) {
        tokens[count] = token;
        token = strtok(SRZ_NULL, delims);
    }
    return count;
}

int srz_load_obj(char const* filename, srz_mesh_t* mesh) {
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

    mesh->faces = malloc(face_max * sizeof(srz_face_t));
    mesh->positions = malloc(position_max * sizeof(srz_float3_t));
    mesh->tex_coords = malloc(tex_coord_max * sizeof(srz_float2_t));
    mesh->normals = malloc(normal_max * sizeof(srz_float3_t));

    char line[256];
    while (fgets(line, 256, file)) {
        char* tokens[5];
        int token_count = _srz_split(line, " ", tokens, 5);

        if (token_count == 0) {
            continue;
        }

        // Vertex position.
        if (strcmp(tokens[0], "v") == 0) {
            if (mesh->position_count == position_max) {
                mesh->positions = realloc(mesh->positions, (position_max *= 2) * sizeof(srz_float3_t));
            }
            mesh->positions[mesh->position_count++] = (srz_float3_t){atof(tokens[1]), atof(tokens[2]), atof(tokens[3])};
        }
        // Vertex texture coordinate.
        else if (strcmp(tokens[0], "vt") == 0) {
            if (mesh->tex_coord_count == tex_coord_max) {
                mesh->tex_coords = realloc(mesh->tex_coords, (tex_coord_max *= 2) * sizeof(srz_float2_t));
            }
            mesh->tex_coords[mesh->tex_coord_count++] = (srz_float2_t){atof(tokens[1]), atof(tokens[2])};
        }
        // Vertex normal.
        else if (strcmp(tokens[0], "vn") == 0) {
            if (mesh->normal_count == normal_max) {
                mesh->normals = realloc(mesh->normals, (normal_max *= 2) * sizeof(srz_float3_t));
            }
            mesh->normals[mesh->normal_count++] = (srz_float3_t){atof(tokens[1]), atof(tokens[2]), atof(tokens[3])};
        }
        // Face.
        else if (strcmp(tokens[0], "f") == 0) {
            char* face_tokens[4][3];
            int face_token_count = token_count - 1;
            for (int i = 0; i < face_token_count; ++i) {
                // Split into indices.
                if (_srz_split(tokens[i + 1], "/", face_tokens[i], 3) != 3) {
                    return -1;
                }
            }

            int position_indices[4];
            int tex_coord_indices[4];
            int normal_indices[4];

            for (int i = 0; i < face_token_count; ++i) {
                position_indices[i] = atoi(face_tokens[i][0]);
                tex_coord_indices[i] = atoi(face_tokens[i][1]);
                normal_indices[i] = atoi(face_tokens[i][2]);
            }

            if (face_token_count == 3) {
                srz_face_t face;
                for (int i = 0; i < 3; ++i) {
                    face.position_indices[i] = position_indices[i];
                    face.tex_coord_indices[i] = tex_coord_indices[i];
                    face.normal_indices[i] = normal_indices[i];
                }

                if (mesh->face_count == face_max) {
                    mesh->faces = realloc(mesh->faces, (face_max *= 2) * sizeof(srz_face_t));
                }
                mesh->faces[mesh->face_count++] = face;
            }
            else if (face_token_count == 4) {
                srz_face_t a, b;
                for (int i = 0; i < 3; ++i) {
                    a.position_indices[i] = position_indices[i];
                    a.tex_coord_indices[i] = tex_coord_indices[i];
                    a.normal_indices[i] = normal_indices[i];

                    b.position_indices[i] = position_indices[(i + 2) % 4];
                    b.tex_coord_indices[i] = tex_coord_indices[(i + 2) % 4];
                    b.normal_indices[i] = normal_indices[(i + 2) % 4];
                }

                if (mesh->face_count >= face_max - 1) {
                    mesh->faces = realloc(mesh->faces, (face_max *= 2) * sizeof(srz_face_t));
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

void srz_mesh_free(srz_mesh_t* mesh) {
    free(mesh->positions);
    free(mesh->tex_coords);
    free(mesh->normals);
    free(mesh->faces);
}
#endif
#endif
