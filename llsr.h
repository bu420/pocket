#ifndef LLSR_H
#define LLSR_H

#define LLSR_PI 3.1415926f

typedef unsigned char llsr_byte_t;

typedef union {
    struct {
        llsr_byte_t x, y, z;
    };
    struct {
        llsr_byte_t r, g, b;
    };
    llsr_byte_t values[3];
} llsr_byte3_t;

typedef union {
    struct {
        llsr_byte_t x, y, z, w;
    };
    struct {
        llsr_byte_t r, g, b, a;
    };
    llsr_byte_t values[4];
} llsr_byte4_t;

typedef union {
    struct {
        int x, y;
    };
    struct {
        int r, g;
    };
    int values[2];
} llsr_int2_t;

typedef union {
    struct {
        float x, y;
    };
    struct {
        float r, g;
    };
    struct {
        float u, v;
    };
    float values[2];
} llsr_float2_t;

typedef union {
    struct {
        float x, y, z;
    };
    struct {
        float r, g, b;
    };
    struct {
        float u, v, w;
    };
    float values[3];
} llsr_float3_t;

typedef union {
    struct {
        float x, y, z, w;
    };
    struct {
        float r, g, b, a;
    };
    float values[4];
} llsr_float4_t;

typedef struct {
    float m[4][4];
} llsr_matrix_t;

typedef struct {
    int w, h;
    llsr_byte3_t* data;
} llsr_color_buffer_t;

typedef struct {
    int w, h;
    float* data;
} llsr_depth_buffer_t;

typedef struct {
    int w, h;
    llsr_byte4_t* data;
} llsr_texture_t;

typedef struct {
    llsr_int2_t start;
    llsr_int2_t end;
    llsr_int2_t current;
    llsr_int2_t delta;
    llsr_int2_t dir;
    int deviation;
} llsr_bresenham_line_t;

int llsr_max(int a, int b);
int llsr_min(int a, int b);
int llsr_max3(int a, int b, int c);
int llsr_min3(int a, int b, int c);
int llsr_abs(int n);
float llsr_maxf(float a, float b);
float llsr_minf(float a, float b);
float llsr_absf(float n);
int llsr_floorf(float n);
int llsr_ceilf(float n);
int llsr_roundf(float n);
float llsr_sqrtf(float n);
float llsr_fmodf(float a, float b);
float llsr_sinf(float x);
float llsr_cosf(float x);
float llsr_tanf(float n);
llsr_float3_t llsr_normalize(llsr_float3_t f);
llsr_float3_t llsr_cross(llsr_float3_t a, llsr_float3_t b);
float llsr_dot(llsr_float3_t a, llsr_float3_t b);
void llsr_swap(int* a, int* b);
void llsr_swapf(float* a, float* b);
void llsr_int2_swap(llsr_int2_t* a, llsr_int2_t* b);
void llsr_float2_swap(llsr_float2_t* a, llsr_float2_t* b);
void llsr_float3_swap(llsr_float3_t* a, llsr_float3_t* b);
llsr_byte3_t llsr_byte3_lerp(llsr_byte3_t a, llsr_byte3_t b, float amount);

void llsr_matrix_init_zero(llsr_matrix_t* matrix);
void llsr_matrix_init_identity(llsr_matrix_t* matrix);
llsr_matrix_t llsr_matrix_mul(llsr_matrix_t a, llsr_matrix_t b);
llsr_matrix_t llsr_matrix_translate(llsr_matrix_t matrix, llsr_float3_t f);
llsr_matrix_t llsr_matrix_rotate_x(llsr_matrix_t matrix, float a);
llsr_matrix_t llsr_matrix_rotate_y(llsr_matrix_t matrix, float a);
llsr_matrix_t llsr_matrix_rotate_z(llsr_matrix_t matrix, float a);
llsr_matrix_t llsr_matrix_scale(llsr_matrix_t matrix, llsr_float3_t f);
llsr_float3_t llsr_matrix_mul_float3(llsr_matrix_t matrix, llsr_float3_t f);
llsr_float4_t llsr_matrix_mul_float4(llsr_matrix_t matrix, llsr_float4_t f);
llsr_float3_t llsr_float3_mul_matrix(llsr_float3_t f, llsr_matrix_t matrix);
llsr_float4_t llsr_float4_mul_matrix(llsr_float4_t f, llsr_matrix_t matrix);

llsr_matrix_t llsr_look_at(llsr_float3_t pos, llsr_float3_t target, llsr_float3_t up);
llsr_matrix_t llsr_perspective(float aspect, float fov, float near, float far);

llsr_byte3_t* llsr_color_buffer_at(llsr_color_buffer_t* color_buffer, int x, int y);
void llsr_color_buffer_clear(llsr_color_buffer_t* color_buffer, llsr_byte3_t color);
float* llsr_depth_buffer_at(llsr_depth_buffer_t* depth_buffer, int x, int y);
void llsr_depth_buffer_clear(llsr_depth_buffer_t* depth_buffer);
llsr_byte4_t* llsr_texture_at(llsr_texture_t* texture, int x, int y);
llsr_byte4_t llsr_texture_sample(llsr_texture_t texture, float u, float v);

llsr_bresenham_line_t llsr_bresenham_line_create(llsr_int2_t start, llsr_int2_t end);
// @return true while it has not yet reached it's end position.
int llsr_bresenham_line_step(llsr_bresenham_line_t* line);
void llsr_raster_line(llsr_color_buffer_t* color_buffer, llsr_int2_t start, llsr_int2_t end, llsr_byte3_t start_color, llsr_byte3_t end_color);
void llsr_raster_triangle_2d(llsr_color_buffer_t* color_buffer, llsr_int2_t pos0, llsr_int2_t pos1, llsr_int2_t pos2, llsr_byte3_t color);
void llsr_raster_triangle_3d(llsr_color_buffer_t* color_buffer, llsr_depth_buffer_t* depth_buffer, llsr_float3_t pos0, llsr_float3_t pos1, llsr_float3_t pos2, llsr_byte3_t color);

#endif
