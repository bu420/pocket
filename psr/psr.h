#ifndef PSR_H
#define PSR_H

#define PSR_ARITHMETIC(out, a, b, nr_elements, op)          \
    for (int _i = 0; _i < nr_elements; _i++)                \
        out.values[_i] = a.values[_i] op b.values[_i];

#define PSR_ADD(out, a, b, nr_elements) PSR_ARITHMETIC(out, a, b, nr_elements, +)
#define PSR_SUB(out, a, b, nr_elements) PSR_ARITHMETIC(out, a, b, nr_elements, -)
#define PSR_MUL(out, a, b, nr_elements) PSR_ARITHMETIC(out, a, b, nr_elements, *)
#define PSR_DIV(out, a, b, nr_elements) PSR_ARITHMETIC(out, a, b, nr_elements, /)

typedef unsigned char psr_byte_t;

typedef union {
    struct {
        psr_byte_t x, y, z;
    };
    struct {
        psr_byte_t r, g, b;
    };
    psr_byte_t values[3];
} psr_byte3_t;

typedef union {
    struct {
        psr_byte_t x, y, z, w;
    };
    struct {
        psr_byte_t r, g, b, a;
    };
    psr_byte_t values[4];
} psr_byte4_t;

typedef union {
    struct {
        int x, y;
    };
    struct {
        int r, g;
    };
    int values[2];
} psr_int2_t;

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
} psr_float2_t;

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
} psr_float3_t;

typedef union {
    struct {
        float x, y, z, w;
    };
    struct {
        float r, g, b, a;
    };
    float values[4];
} psr_float4_t;

typedef struct {
    float m[3][3];
    struct {
        float m00, m10, m20, m01, m11, m21, m02, m12, m22;
    };
} psr_mat3_t;

typedef union {
    float m[4][4];
    struct {
        float m00, m10, m20, m30, m01, m11, m21, m31, m02, m12, m22, m32, m03, m13, m23, m33;
    };
} psr_mat4_t;

typedef struct {
    int x, y, w, h;
} psr_rect_t;

typedef struct {
    int w, h;
    psr_byte3_t* data;
} psr_color_buffer_t;

typedef struct {
    int w, h;
    float* data;
} psr_depth_buffer_t;

typedef enum {
    PSR_R5G6B5      = 1 << 0,
    PSR_R5G5B5A1    = 1 << 1, 
    PSR_R8G8B8      = 1 << 2,
    PSR_R8G8B8A8    = 1 << 3
} psr_color_depth_t;

typedef struct {
    int w, h;
    psr_color_depth_t color_depth;
    psr_byte_t* data;
} psr_image_t;

typedef struct {
    int position_indices[3];
    int tex_coord_indices[3];
    int normal_indices[3];
} psr_face_t;

typedef struct {
    psr_float3_t* positions;
    int position_count;
    psr_float2_t* tex_coords;
    int tex_coord_count;
    psr_float3_t* normals;
    int normal_count;
    psr_face_t* faces;
    int face_count;
} psr_mesh_t;

typedef struct psr_font_t psr_font_t;

psr_float3_t psr_normalize(psr_float3_t f);
psr_float3_t psr_cross(psr_float3_t a, psr_float3_t b);
float psr_dot(psr_float3_t a, psr_float3_t b);
psr_byte3_t psr_byte3_lerp(psr_byte3_t a, psr_byte3_t b, float amount);

// Matrix.

void psr_mat4_init_zero(psr_mat4_t* m);
void psr_mat4_init_identity(psr_mat4_t* m);
psr_mat4_t psr_mat4_mul(psr_mat4_t a, psr_mat4_t b);
psr_mat4_t psr_mat4_transpose(psr_mat4_t m);
psr_mat4_t psr_mat4_inverse(psr_mat4_t m);
psr_mat4_t psr_mat4_translate(psr_mat4_t mat4, psr_float3_t f);
psr_mat4_t psr_mat4_rotate_x(psr_mat4_t mat4, float a);
psr_mat4_t psr_mat4_rotate_y(psr_mat4_t mat4, float a);
psr_mat4_t psr_mat4_rotate_z(psr_mat4_t mat4, float a);
psr_mat4_t psr_mat4_scale(psr_mat4_t mat4, psr_float3_t f);
psr_float3_t psr_mat3_mul_float3(psr_mat3_t m, psr_float3_t f);
psr_float3_t psr_float3_mul_mat3(psr_float3_t f, psr_mat3_t m);
psr_float3_t psr_mat4_mul_float3(psr_mat4_t m, psr_float3_t f);
psr_float4_t psr_mat4_mul_float4(psr_mat4_t m, psr_float4_t f);
psr_float3_t psr_float3_mul_mat4(psr_float3_t f, psr_mat4_t m);
psr_float4_t psr_float4_mul_mat4(psr_float4_t f, psr_mat4_t m);
psr_mat3_t psr_mat4_to_mat3(psr_mat4_t m);

psr_mat4_t psr_look_at(psr_float3_t pos, psr_float3_t target, psr_float3_t up);
psr_mat4_t psr_perspective(float aspect, float fov, float near, float far);
psr_mat4_t psr_ortho();

// Buffer.

psr_color_buffer_t* psr_color_buffer_create(int w, int h);
void psr_color_buffer_destroy(psr_color_buffer_t* color_buffer);
psr_byte3_t* psr_color_buffer_at(psr_color_buffer_t* color_buffer, int x, int y);
void psr_color_buffer_clear(psr_color_buffer_t* color_buffer, psr_byte3_t color);

psr_depth_buffer_t* psr_depth_buffer_create(int w, int h);
void psr_depth_buffer_destroy(psr_depth_buffer_t* depth_buffer);
float* psr_depth_buffer_at(psr_depth_buffer_t* depth_buffer, int x, int y);
void psr_depth_buffer_clear(psr_depth_buffer_t* depth_buffer);

// Image.

psr_image_t* psr_image_create(psr_color_depth_t color_depth, int w, int h);
void psr_image_destroy(psr_image_t* image);
// Retuns the address of the first byte of the pixel.
psr_byte_t* psr_image_at(psr_image_t* image, int x, int y);

// Raster.

void psr_raster_line(psr_color_buffer_t* color_buffer, psr_int2_t start, psr_int2_t end, psr_byte3_t start_color, psr_byte3_t end_color);

void psr_raster_triangle_2d_color(psr_color_buffer_t* color_buffer, psr_int2_t pos0, psr_int2_t pos1, psr_int2_t pos2, psr_byte3_t color);
void psr_raster_triangle_2d_image(psr_color_buffer_t* color_buffer, psr_int2_t pos0, psr_int2_t pos1, psr_int2_t pos2);
void psr_raster_triangle_2d_callback(psr_color_buffer_t* color_buffer, psr_int2_t pos0, psr_int2_t pos1, psr_int2_t pos2, void (*callback)(psr_int2_t pixel_pos), void* user_data);

void psr_raster_triangle_3d(psr_color_buffer_t* color_buffer, psr_depth_buffer_t* depth_buffer, psr_float3_t pos0, psr_float3_t pos1, psr_float3_t pos2, psr_byte3_t color);

void psr_raster_image(psr_color_buffer_t* color_buffer, psr_image_t* image, psr_rect_t src, psr_rect_t dst);

void psr_raster_text(psr_color_buffer_t* color_buffer, char* text, psr_int2_t pos, psr_font_t* font, int scale);

// Asset I/O.

psr_image_t* psr_image_load_bmp(char* path, psr_color_depth_t color_depth);

void psr_save_bmp(char* path, psr_color_buffer_t* color_buffer);

psr_mesh_t* psr_mesh_load_obj(char* path);
void psr_mesh_free(psr_mesh_t* mesh);

// Parse AngelCode's bitmap font information.
psr_font_t* psr_font_load(char* image_path, char* info_path);
void psr_font_destroy(psr_font_t* font);

#endif
