/*
offline-software-renderer
https://github.com/bu420/offline-software-renderer

Define 'OSR_IMPLEMENTATION' once.
Define 'OSR_SAVE_AND_LOAD' to enable convenience functions for saving to and loading from disk.

Credits to 'el nora' for sine and cosine implementations.
*/

#pragma once

#define OSR_PI 3.1415926f

typedef unsigned char osr_byte_t;

typedef union {
    struct {
        osr_byte_t x, y, z;
    };
    struct {
        osr_byte_t r, g, b;
    };
    osr_byte_t values[3];
} osr_byte3_t;

typedef union {
    struct {
        int x, y;
    };
    struct {
        int r, g;
    };
    int values[2];
} osr_int2_t;

typedef union {
    struct {
        float x, y;
    };
    struct {
        float r, g;
    };
    float values[2];
} osr_float2_t;

typedef union {
    struct {
        float x, y, z;
    };
    struct {
        float r, g, b;
    };
    float values[3];
} osr_float3_t;

typedef union {
    struct {
        float x, y, z, w;
    };
    struct {
        float r, g, b, a;
    };
    float values[4];
} osr_float4_t;

typedef struct {
    float elements[16];
} osr_matrix_t;

typedef struct {
    int w, h;
    osr_byte3_t* pixels;
} osr_image_t;

int osr_max(int a, int b);
int osr_min(int a, int b);
int osr_max3(int a, int b, int c);
int osr_min3(int a, int b, int c);
float osr_maxf(float a, float b);
float osr_minf(float a, float b);
float osr_absf(float n);
float osr_floorf(float n);
float osr_ceilf(float n);
float osr_roundf(float n);
float osr_sqrtf(float n);
float osr_fmodf(float a, float b);
float osr_sinf(float x);
float osr_cosf(float x);
float osr_tanf(float n);
osr_float3_t osr_normalize(osr_float3_t f);
osr_float3_t osr_cross(osr_float3_t a, osr_float3_t b);
float osr_dot(osr_float3_t a, osr_float3_t b);

void osr_matrix_init_zero(osr_matrix_t* matrix);
void osr_matrix_init_identity(osr_matrix_t* matrix);
float* osr_matrix_at(osr_matrix_t* matrix, int x, int y);
osr_matrix_t osr_matrix_mul(osr_matrix_t a, osr_matrix_t b);
osr_matrix_t osr_matrix_translate(osr_matrix_t matrix, osr_float3_t f);
// @param f radians.
osr_matrix_t osr_matrix_rotate(osr_matrix_t matrix, osr_float3_t f);
osr_matrix_t osr_matrix_scale(osr_matrix_t matrix, osr_float3_t f);
osr_float3_t osr_matrix_mul_float3(osr_matrix_t matrix, osr_float3_t f);
osr_float4_t osr_matrix_mul_float4(osr_matrix_t matrix, osr_float4_t f);
osr_float3_t osr_float3_mul_matrix(osr_float3_t f, osr_matrix_t matrix);
osr_float4_t osr_float4_mul_matrix(osr_float4_t f, osr_matrix_t matrix);
osr_matrix_t osr_create_view_matrix(osr_float3_t pos, osr_float3_t dir, osr_float3_t up);
// @param fov degrees.
osr_matrix_t osr_create_projection_matrix(float aspect, float fov, float near, float far);

// @param memory pointer to memory that the user must allocate and free themself.
void osr_image_init(osr_image_t* image, int w, int h, void* memory);
osr_byte3_t* osr_image_at(osr_image_t* image, int x, int y);
void osr_image_raster_line(osr_image_t* image, osr_int2_t a, osr_int2_t b, osr_byte3_t color);
void osr_image_raster_tri(osr_image_t* image, osr_int2_t a, osr_int2_t b, osr_int2_t c, osr_byte3_t color);
#ifdef OSR_SAVE_AND_LOAD
void osr_image_save_bmp(char const* filename, osr_image_t image);
#endif

#ifdef OSR_IMPLEMENTATION

#ifdef OSR_SAVE_AND_LOAD
#include <stdio.h>
#endif

int osr_max(int a, int b) {
    return a > b ? a : b;
}

int osr_min(int a, int b) {
    return a < b ? a : b;
}

int osr_max3(int a, int b, int c) {
    return osr_max(osr_max(a, b), c);
}

int osr_min3(int a, int b, int c) {
    return osr_min(osr_min(a, b), c);
}

float osr_maxf(float a, float b) {
    return a > b ? a : b;
}

float osr_minf(float a, float b) {
    return a < b ? a : b;
}

float osr_absf(float n) {
    return n >= 0 ? n : -n;
}

float osr_floorf(float n) {
    return n >= 0 ? (int)n : (int)n - 1;
}

float osr_ceilf(float n) {
    return -osr_floorf(-n);
}

float osr_roundf(float n) {
    return (n >= 0) ? osr_floorf(n + .5f) : osr_ceilf(n - .5f);
}

float osr_sqrtf(float n) {
    float lo = osr_minf(1, n);
    float hi = osr_maxf(1, n);
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

float osr_fmodf(float a, float b) {
    return a - (osr_roundf(a / b) * b);
}

float _osr_sf(float x) {
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

float _osr_cf(float x) {
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

float osr_sinf(float x) {
    int i = x / OSR_PI / 2;
    float a = x - 2 * OSR_PI * i;
    if (0 <= a && a < OSR_PI / 4) {
        return _osr_sf(a);
    }
    if (OSR_PI / 4 <= a && a < OSR_PI / 2) {
        return _osr_cf(OSR_PI / 2 - a);
    }
    if (OSR_PI / 2 <= a && a < 3 * OSR_PI / 4) {
        return _osr_cf(a - OSR_PI / 2);
    }
    if (3 * OSR_PI / 4 <= a && a < OSR_PI) {
        return _osr_sf(OSR_PI - a);
    }
    return -osr_sinf(a - OSR_PI);
}

float osr_cosf(float x) {
    return osr_sinf(x + OSR_PI / 2);
}

float osr_tanf(float n) {
    return osr_sinf(n) / osr_cosf(n);
}

osr_float3_t osr_normalize(osr_float3_t f) {
    float len = osr_sqrtf(f.x * f.x + f.y * f.y + f.z * f.z);
    return (osr_float3_t){f.x / len, f.y / len, f.z / len};
}

osr_float3_t osr_cross(osr_float3_t a, osr_float3_t b) {
    return (osr_float3_t){
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

float osr_dot(osr_float3_t a, osr_float3_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

void osr_matrix_init_zero(osr_matrix_t* matrix) {
    for (int i = 0; i < 16; ++i) {
        matrix->elements[i] = 0;
    }
}

void osr_matrix_init_identity(osr_matrix_t* matrix) {
    osr_matrix_init_zero(matrix);

    for (int i = 0; i < 4; ++i) {
        *osr_matrix_at(matrix, i, i) = 1;
    }
}

float* osr_matrix_at(osr_matrix_t* matrix, int x, int y) {
    return &matrix->elements[y * 4 + x];
}

osr_matrix_t osr_matrix_mul(osr_matrix_t a, osr_matrix_t b) {
    osr_matrix_t result;
    osr_matrix_init_zero(&result);

    for (int x = 0; x < 4; ++x) {
        for (int y = 0; y < 4; ++y) {
            for (int i = 0; i < 4; ++i) {
                *osr_matrix_at(&result, x, y) += *osr_matrix_at(&a, x, i) * *osr_matrix_at(&b, i, y);
            }
        }
    }

    return result;
}

osr_matrix_t osr_matrix_translate(osr_matrix_t matrix, osr_float3_t f) {
    osr_matrix_t m;
    osr_matrix_init_identity(&m);

    *osr_matrix_at(&m, 0, 3) = f.x;
    *osr_matrix_at(&m, 1, 3) = f.y;
    *osr_matrix_at(&m, 2, 3) = f.z;

    return osr_matrix_mul(matrix, m);
}

osr_matrix_t osr_matrix_rotate(osr_matrix_t matrix, osr_float3_t f) {
    if (f.x != 0) {
        float s = osr_sinf(f.x);
        float c = osr_cosf(f.x);

        osr_matrix_t x;
        osr_matrix_init_identity(&x);

        *osr_matrix_at(&x, 1, 1) = c;
        *osr_matrix_at(&x, 1, 2) = -s;
        *osr_matrix_at(&x, 2, 1) = s;
        *osr_matrix_at(&x, 2, 2) = c;

        matrix = osr_matrix_mul(matrix, x);
    }

    if (f.y != 0) {
        float s = osr_sinf(f.y);
        float c = osr_cosf(f.y);

        osr_matrix_t y;
        osr_matrix_init_identity(&y);

        *osr_matrix_at(&y, 0, 0) = c;
        *osr_matrix_at(&y, 0, 2) = s;
        *osr_matrix_at(&y, 2, 0) = -s;
        *osr_matrix_at(&y, 2, 2) = c;

        matrix = osr_matrix_mul(matrix, y);
    }

    if (f.z != 0) {
        float s = osr_sinf(f.z);
        float c = osr_cosf(f.z);

        osr_matrix_t z;
        osr_matrix_init_identity(&z);

        *osr_matrix_at(&z, 0, 0) = c;
        *osr_matrix_at(&z, 0, 1) = -s;
        *osr_matrix_at(&z, 1, 0) = s;
        *osr_matrix_at(&z, 1, 1) = c;

        matrix = osr_matrix_mul(matrix, z);
    }

    return matrix;
}

osr_matrix_t osr_matrix_scale(osr_matrix_t matrix, osr_float3_t f) {
    osr_matrix_t m;
    osr_matrix_init_identity(&m);

    *osr_matrix_at(&m, 0, 0) = f.x;
    *osr_matrix_at(&m, 1, 1) = f.y;
    *osr_matrix_at(&m, 2, 2) = f.z;

    return osr_matrix_mul(matrix, m);
}

osr_float3_t osr_matrix_mul_float3(osr_matrix_t matrix, osr_float3_t f) {
    osr_float4_t result;
    for (int x = 0; x < 4; ++x) {
        float i = 0;
        for (int y = 0; y < 3; ++y) {
            i += *osr_matrix_at(&matrix, x, y) * f.values[y];
        }
        result.values[x] = i + *osr_matrix_at(&matrix, x, 3);
    }
    
    if (result.w != 0) {
        result.x /= result.w;
        result.y /= result.w;
        result.z /= result.w;
    }

    return (osr_float3_t){result.x, result.y, result.z};
}

osr_float4_t osr_matrix_mul_float4(osr_matrix_t matrix, osr_float4_t f) {
    osr_float4_t result;
    for (int x = 0; x < 4; ++x) {
        float i = 0;
        for (int y = 0; y < 4; ++y) {
            i += *osr_matrix_at(&matrix, x, y) * f.values[y];
        }
        result.values[x] = i;
    }
    return result;
}

osr_float3_t osr_float3_mul_matrix(osr_float3_t f, osr_matrix_t matrix) {
    osr_float4_t result;
    for (int y = 0; y < 4; ++y) {
        float i = 0;
        for (int x = 0; x < 3; ++x) {
            i += *osr_matrix_at(&matrix, x, y) * f.values[x];
        }
        result.values[y] = i + *osr_matrix_at(&matrix, 3, y);
    }
    
    if (result.w != 0) {
        result.x /= result.w;
        result.y /= result.w;
        result.z /= result.w;
    }

    return (osr_float3_t){result.x, result.y, result.z};
}

osr_float4_t osr_float4_mul_matrix(osr_float4_t f, osr_matrix_t matrix) {
    osr_float4_t result;
    for (int y = 0; y < 4; ++y) {
        float i = 0;
        for (int x = 0; x < 4; ++x) {
            i += *osr_matrix_at(&matrix, x, y) * f.values[x];
        }
        result.values[y] = i;
    }
    return result;
}

osr_matrix_t osr_create_view_matrix(osr_float3_t pos, osr_float3_t dir, osr_float3_t up) {
    dir = osr_normalize(dir);
    osr_float3_t right = osr_normalize(osr_cross(up, dir));
    osr_float3_t local_up = osr_normalize(osr_cross(dir, right));

    osr_matrix_t m;
    osr_matrix_init_identity(&m);

    *osr_matrix_at(&m, 0, 0) = right.x;
    *osr_matrix_at(&m, 0, 1) = right.y;
    *osr_matrix_at(&m, 0, 2) = right.z;
    *osr_matrix_at(&m, 1, 0) = local_up.x;
    *osr_matrix_at(&m, 1, 1) = local_up.y;
    *osr_matrix_at(&m, 1, 2) = local_up.z;
    *osr_matrix_at(&m, 2, 0) = dir.x;
    *osr_matrix_at(&m, 2, 1) = dir.y;
    *osr_matrix_at(&m, 2, 2) = dir.z;
    *osr_matrix_at(&m, 0, 3) = -osr_dot(right, pos);
    *osr_matrix_at(&m, 1, 3) = -osr_dot(local_up, pos);
    *osr_matrix_at(&m, 2, 3) = -osr_dot(dir, pos);
    return m;
}

osr_matrix_t osr_create_projection_matrix(float aspect, float fov, float near, float far) {
    osr_matrix_t m;
    osr_matrix_init_zero(&m);

    float scale = 1.f / osr_tanf(fov * .5f / 180.f * OSR_PI);
    *osr_matrix_at(&m, 0, 0) = scale * aspect;
    *osr_matrix_at(&m, 1, 1) = scale;
    *osr_matrix_at(&m, 2, 2) = far / (far - near);
    *osr_matrix_at(&m, 3, 2) = (-far * near) / (far - near);
    *osr_matrix_at(&m, 2, 3) = 1;
    *osr_matrix_at(&m, 3, 3) = 0;
    return m;
}

void osr_image_init(osr_image_t* image, int w, int h, void* memory) {
    image->w = w;
    image->h = h;
    image->pixels = memory;
}

osr_byte3_t* osr_image_at(osr_image_t* image, int x, int y) {
    return &image->pixels[y * image->w + x];
}

void osr_image_raster_line(osr_image_t* image, osr_int2_t a, osr_int2_t b, osr_byte3_t color) {
    osr_int2_t diff = {b.x - a.x, b.y - a.y};
    osr_int2_t diff_abs = {osr_absf(diff.x), osr_absf(diff.y)};
    osr_int2_t p = {2 * diff_abs.y - diff_abs.x, 2 * diff_abs.x - diff_abs.y};
    osr_int2_t pos;
    osr_int2_t e;

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

        *osr_image_at(image, pos.x, pos.y) = color;

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
            *osr_image_at(image, pos.x, pos.y) = color;
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

        *osr_image_at(image, pos.x, pos.y) = color;

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
            *osr_image_at(image, pos.x, pos.y) = color;
        }
    }
}

int _osr_t(osr_int2_t a, osr_int2_t b, osr_int2_t c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

void osr_image_raster_tri(osr_image_t* image, osr_int2_t a, osr_int2_t b, osr_int2_t c, osr_byte3_t color) {
    // Bounding box.
    int min_x = osr_min3(a.x, b.x, c.x);
    int min_y = osr_min3(a.y, b.y, c.y);
    int max_x = osr_max3(a.x, b.x, c.x);
    int max_y = osr_max3(a.y, b.y, c.y);

    // Clip.
    min_x = osr_max(min_x, 0);
    min_y = osr_max(min_y, 0);
    max_x = osr_min(max_x, image->w - 1);
    max_y = osr_min(max_y, image->h - 1);

    // Rasterize.
    osr_int2_t p;
    for (p.x = min_x; p.x <= max_x; ++p.x) {
        for (p.y = min_y; p.y <= max_y; ++p.y) {
            int t0 = _osr_t(b, c, p);
            int t1 = _osr_t(c, a, p);
            int t2 = _osr_t(a, b, p);

            if (t0 >= 0 && t1 >= 0 && t2 >= 0) {
                *osr_image_at(image, p.x, p.y) = color;
            }
        }
    }
}

#ifdef OSR_SAVE_AND_LOAD
void osr_image_save_bmp(char const* filename, osr_image_t image) {
    osr_byte_t const padding[] = {0, 0, 0};
    int const padding_size = (4 - image.w * 3 % 4) % 4;
    int const stride = image.w * 3 + padding_size;
    int const file_size = 54 + stride * image.h;

    osr_byte_t const header[54] = {
        'B', 'M',
        file_size, file_size >> 8, file_size >> 16, file_size >> 24,
        0, 0, 0, 0,
        54, 0, 0, 0,
        40, 0, 0, 0,
        image.w, image.w >> 8, image.w >> 16, image.w >> 24,
        image.h, image.h >> 8, image.h >> 16, image.h >> 24,
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

    for (int y = 0; y < image.h; ++y) {
        for (int x = 0; x < image.w; ++x) {
            int i = (image.h - y) * image.h + x;
            osr_byte3_t pixel = image.pixels[i];

            // RBG to BGR.
            osr_byte3_t color = {.r = pixel.b, .g = pixel.g, .b = pixel.r};

            fwrite(&color, 3, 1, file);
        }

        fwrite(padding, 1, padding_size, file);
    }

    fclose(file);
}
#endif
#endif
