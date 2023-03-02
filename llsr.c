#include "llsr.h"

int llsr_max(int a, int b) {
    return a > b ? a : b;
}

int llsr_min(int a, int b) {
    return a < b ? a : b;
}

int llsr_max3(int a, int b, int c) {
    return llsr_max(llsr_max(a, b), c);
}

int llsr_min3(int a, int b, int c) {
    return llsr_min(llsr_min(a, b), c);
}

int llsr_abs(int n) {
    return n >= 0 ? n : -n;
}

float llsr_maxf(float a, float b) {
    return a > b ? a : b;
}

float llsr_minf(float a, float b) {
    return a < b ? a : b;
}

float llsr_absf(float n) {
    return n >= 0 ? n : -n;
}

int llsr_floorf(float n) {
    return n >= 0 ? (int)n : (int)n - 1;
}

int llsr_ceilf(float n) {
    return -llsr_floorf(-n);
}

int llsr_roundf(float n) {
    return (n >= 0) ? llsr_floorf(n + .5f) : llsr_ceilf(n - .5f);
}

float llsr_sqrtf(float n) {
    float lo = llsr_minf(1, n);
    float hi = llsr_maxf(1, n);
    float mid;

    while (100 * lo * lo < n) {
        lo *= 10;
    }
    while (0.01 * hi * hi > n) {
        hi *= 0.1;
    }

    for (int i = 0; i < 100; i++) {
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

float llsr_fmodf(float a, float b) {
    return a - (llsr_roundf(a / b) * b);
}

float _llsr_sf(float x) {
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

float _llsr_cf(float x) {
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

float llsr_sinf(float x) {
    int i = x / LLSR_PI / 2;
    float a = x - 2 * LLSR_PI * i;
    if (0 <= a && a < LLSR_PI / 4) {
        return _llsr_sf(a);
    }
    if (LLSR_PI / 4 <= a && a < LLSR_PI / 2) {
        return _llsr_cf(LLSR_PI / 2 - a);
    }
    if (LLSR_PI / 2 <= a && a < 3 * LLSR_PI / 4) {
        return _llsr_cf(a - LLSR_PI / 2);
    }
    if (3 * LLSR_PI / 4 <= a && a < LLSR_PI) {
        return _llsr_sf(LLSR_PI - a);
    }
    return -llsr_sinf(a - LLSR_PI);
}

float llsr_cosf(float x) {
    return llsr_sinf(x + LLSR_PI / 2);
}

float llsr_tanf(float n) {
    return llsr_sinf(n) / llsr_cosf(n);
}

llsr_float3_t llsr_normalize(llsr_float3_t f) {
    float len = llsr_sqrtf(f.x * f.x + f.y * f.y + f.z * f.z);
    llsr_float3_t result = {f.x / len, f.y / len, f.z / len};
    return result;
}

llsr_float3_t llsr_cross(llsr_float3_t a, llsr_float3_t b) {
    llsr_float3_t result = {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
    return result;
}

float llsr_dot(llsr_float3_t a, llsr_float3_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

void llsr_swap(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

void llsr_swapf(float* a, float* b) {
    float temp = *a;
    *a = *b;
    *b = temp;
}

void llsr_int2_swap(llsr_int2_t* a, llsr_int2_t* b) {
    llsr_swap(&a->x, &b->x);
    llsr_swap(&a->y, &b->y);
}

void llsr_float2_swap(llsr_float2_t* a, llsr_float2_t* b) {
    llsr_swapf(&a->x, &b->x);
    llsr_swapf(&a->y, &b->y);
}

void llsr_float3_swap(llsr_float3_t* a, llsr_float3_t* b) {
    llsr_swapf(&a->x, &b->x);
    llsr_swapf(&a->y, &b->y);
    llsr_swapf(&a->z, &b->z);
}

llsr_byte3_t llsr_byte3_lerp(llsr_byte3_t a, llsr_byte3_t b, float amount) {
    llsr_byte3_t result;
    for (int i = 0; i < 3; i++) {
        result.values[i] = a.values[i] * (1 - amount) + b.values[i] * amount;
    }
    return result;
}

void llsr_matrix_init_zero(llsr_matrix_t* matrix) {
    for (int x = 0; x < 4; x++) {
        for (int y = 0; y < 4; y++) {
            matrix->m[x][y] = 0;
        }
    }
}

void llsr_matrix_init_identity(llsr_matrix_t* matrix) {
    llsr_matrix_init_zero(matrix);
    for (int i = 0; i < 4; i++) {
        matrix->m[i][i] = 1;
    }
}

llsr_matrix_t llsr_matrix_mul(llsr_matrix_t a, llsr_matrix_t b) {
    llsr_matrix_t result;
    llsr_matrix_init_zero(&result);

    for (int x = 0; x < 4; x++) {
        for (int y = 0; y < 4; y++) {
            for (int i = 0; i < 4; i++) {
                result.m[x][y] += a.m[x][i] * b.m[i][y];
            }
        }
    }
    return result;
}

llsr_matrix_t llsr_matrix_translate(llsr_matrix_t matrix, llsr_float3_t f) {
    llsr_matrix_t m;
    llsr_matrix_init_identity(&m);

    m.m[3][0] = f.x;
    m.m[3][1] = f.y;
    m.m[3][2] = f.z;

    return llsr_matrix_mul(matrix, m);
}

llsr_matrix_t llsr_matrix_rotate_x(llsr_matrix_t matrix, float a) {
    float s = llsr_sinf(a);
    float c = llsr_cosf(a);

    llsr_matrix_t x;
    llsr_matrix_init_identity(&x);

    x.m[1][1] = c;
    x.m[1][2] = -s;
    x.m[2][1] = s;
    x.m[2][2] = c;

    return llsr_matrix_mul(matrix, x);
}

llsr_matrix_t llsr_matrix_rotate_y(llsr_matrix_t matrix, float a) {
    float s = llsr_sinf(a);
    float c = llsr_cosf(a);

    llsr_matrix_t y;
    llsr_matrix_init_identity(&y);

    y.m[0][0] = c;
    y.m[0][2] = s;
    y.m[2][0] = -s;
    y.m[2][2] = c;

    return llsr_matrix_mul(matrix, y);
}

llsr_matrix_t llsr_matrix_rotate_z(llsr_matrix_t matrix, float a) {
    float s = llsr_sinf(a);
    float c = llsr_cosf(a);

    llsr_matrix_t z;
    llsr_matrix_init_identity(&z);

    z.m[0][0] = c;
    z.m[0][1] = -s;
    z.m[1][0] = s;
    z.m[1][1] = c;

    return llsr_matrix_mul(matrix, z);
}

llsr_matrix_t llsr_matrix_scale(llsr_matrix_t matrix, llsr_float3_t f) {
    llsr_matrix_t m;
    llsr_matrix_init_identity(&m);

    m.m[0][0] = f.x;
    m.m[1][1] = f.y;
    m.m[2][2] = f.z;

    return llsr_matrix_mul(matrix, m);
}

llsr_float3_t llsr_matrix_mul_float3(llsr_matrix_t matrix, llsr_float3_t f) {
    llsr_float4_t f4;
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

    llsr_float3_t result = {f4.x, f4.y, f4.z};
    return result;
}

llsr_float4_t llsr_matrix_mul_float4(llsr_matrix_t matrix, llsr_float4_t f) {
    llsr_float4_t result;
    for (int x = 0; x < 4; x++) {
        float i = 0;
        for (int y = 0; y < 4; y++) {
            i += matrix.m[x][y] * f.values[y];
        }
        result.values[x] = i;
    }
    return result;
}

llsr_float3_t llsr_float3_mul_matrix(llsr_float3_t f, llsr_matrix_t matrix) {
    llsr_float4_t f4;
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

    llsr_float3_t result = {f4.x, f4.y, f4.z};
    return result;
}

llsr_float4_t llsr_float4_mul_matrix(llsr_float4_t f, llsr_matrix_t matrix) {
    llsr_float4_t result;
    for (int y = 0; y < 4; y++) {
        float i = 0;
        for (int x = 0; x < 4; x++) {
            i += matrix.m[x][y] * f.values[x];
        }
        result.values[y] = i;
    }
    return result;
}

llsr_matrix_t llsr_look_at(llsr_float3_t pos, llsr_float3_t target, llsr_float3_t up) {
    llsr_float3_t diff = {target.x - pos.x, target.y - pos.y, target.z - pos.z};
    llsr_float3_t forward = llsr_normalize(diff);
    llsr_float3_t right = llsr_normalize(llsr_cross(forward, up));
    llsr_float3_t local_up = llsr_normalize(llsr_cross(right, forward));

    llsr_matrix_t m;
    llsr_matrix_init_identity(&m);
    m.m[0][0] = right.x;
    m.m[1][0] = right.y;
    m.m[2][0] = right.z;
    m.m[0][1] = local_up.x;
    m.m[1][1] = local_up.y;
    m.m[2][1] = local_up.z;
    m.m[0][2] = -forward.x;
    m.m[1][2] = -forward.y;
    m.m[2][2] = -forward.z;
    m.m[3][0] = -llsr_dot(right, pos);
    m.m[3][1] = -llsr_dot(local_up, pos);
    m.m[3][2] = llsr_dot(forward, pos);
    return m;
}

llsr_matrix_t llsr_perspective(float aspect, float fov, float near, float far) {    
    llsr_matrix_t m;
    llsr_matrix_init_zero(&m);

    float half_tan = llsr_tanf(fov / 2);

    m.m[0][0] = 1 / (half_tan * aspect);
    m.m[1][1] = 1 / half_tan;
    m.m[2][2] = -(far + near) / (far - near);
    m.m[2][3] = -1;
    m.m[3][2] = -(2 * far * near) / (far - near);
    return m;
}

llsr_byte3_t* llsr_color_buffer_at(llsr_color_buffer_t* color_buffer, int x, int y) {
    return &color_buffer->data[y * color_buffer->w + x];
}

void llsr_color_buffer_clear(llsr_color_buffer_t* color_buffer, llsr_byte3_t color) {
    for (int i = 0; i < color_buffer->w * color_buffer->h; i++) {
        color_buffer->data[i] = color;
    }
}

float* llsr_depth_buffer_at(llsr_depth_buffer_t* depth_buffer, int x, int y) {
    return &depth_buffer->data[y * depth_buffer->w + x];
}

void llsr_depth_buffer_clear(llsr_depth_buffer_t* depth_buffer) {
    for (int i = 0; i < depth_buffer->w * depth_buffer->h; i++) {
        depth_buffer->data[i] = 0;
    }
}

llsr_byte4_t* llsr_texture_at(llsr_texture_t* texture, int x, int y) {
    return &texture->data[y * texture->w + x];
}

llsr_byte4_t llsr_texture_sample(llsr_texture_t texture, float u, float v) {
    return *llsr_texture_at(&texture, (int)(u * texture.w), (int)(v * texture.h - 1));
}

llsr_bresenham_line_t llsr_bresenham_line_create(llsr_int2_t start, llsr_int2_t end) {
    llsr_bresenham_line_t line;
    line.start = start;
    line.end = end;
    line.current = start;
    line.delta = (llsr_int2_t){llsr_abs(end.x - start.x), -llsr_abs(end.y - start.y)};
    line.dir = (llsr_int2_t){start.x < end.x ? 1 : -1, start.y < end.y ? 1 : -1};
    line.deviation = line.delta.x + line.delta.y;
    return line;
}

int llsr_bresenham_line_step(llsr_bresenham_line_t* line) {
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
float _llsr_bresenham_line_inverse_lerp(llsr_int2_t start, llsr_int2_t end, llsr_int2_t current) {
    if (end.x - start.x > end.y - start.y) {
        return (current.x - start.x) / (float)(end.x - start.x);
    }
    return (current.y - start.y) / (float)(end.y - start.y);
}

void llsr_raster_line(llsr_color_buffer_t* color_buffer, llsr_int2_t start, llsr_int2_t end, llsr_byte3_t start_color, llsr_byte3_t end_color) {
    llsr_bresenham_line_t line = llsr_bresenham_line_create(start, end);
    do {
        llsr_byte3_t current_color = llsr_byte3_lerp(start_color, end_color, _llsr_bresenham_line_inverse_lerp(start, end, line.current));
        *llsr_color_buffer_at(color_buffer, line.current.x, line.current.y) = current_color;
    } 
    while (llsr_bresenham_line_step(&line));
}

void _llsr_sort_triangle_vertices_by_height(llsr_int2_t* pos0, llsr_int2_t* pos1, llsr_int2_t* pos2) {
    if (pos0->y > pos1->y) {
        llsr_int2_swap(pos0, pos1);
    }
    if (pos0->y > pos2->y) {
        llsr_int2_swap(pos0, pos2);
    }
    if (pos1->y > pos2->y) {
        llsr_int2_swap(pos1, pos2);
    }
}

int _llsr_bresenham_line_step_until_vertical_difference(llsr_bresenham_line_t* line) {
    int y = line->current.y;

    while (llsr_bresenham_line_step(line)) {
        if (line->current.y != y) {
            return 1;
        }
        y = line->current.y;
    }
    return 0;
}

// Either the top or bottom of the triangle must be flat (the lines must share the same start and end y position).
void _llsr_raster_triangle_2d_flat(llsr_color_buffer_t* color_buffer, llsr_bresenham_line_t line_a, llsr_bresenham_line_t line_b, llsr_byte3_t color) {
    for (int y = line_a.start.y; y <= line_a.end.y; y++) {
        int start = line_a.current.x;
        int end = line_b.current.x;

        if (start > end) {
            llsr_swap(&start, &end);
        }

        for (int x = start; x <= end; x++) {
            // TODO: find way to check bounds.
            if (x >= 0 && x < color_buffer->w && y >= 0 && y < color_buffer->h) {
                *llsr_color_buffer_at(color_buffer, x, y) = color;
            }
        }

        _llsr_bresenham_line_step_until_vertical_difference(&line_a);
        _llsr_bresenham_line_step_until_vertical_difference(&line_b);
    }
}

void llsr_raster_triangle_2d(llsr_color_buffer_t* color_buffer, llsr_int2_t pos0, llsr_int2_t pos1, llsr_int2_t pos2, llsr_byte3_t color) {
    _llsr_sort_triangle_vertices_by_height(&pos0, &pos1, &pos2);

    // Check if the top of the triangle is flat.
    if (pos0.y == pos1.y) {
        _llsr_raster_triangle_2d_flat(color_buffer, llsr_bresenham_line_create(pos0, pos2), llsr_bresenham_line_create(pos1, pos2), color);
    }
    // Check if the bottom is flat.
    else if (pos1.y == pos2.y) {
        _llsr_raster_triangle_2d_flat(color_buffer, llsr_bresenham_line_create(pos0, pos1), llsr_bresenham_line_create(pos0, pos2), color);
    }
    // Else plit triangle into two smaller triangles.
    else {
        llsr_int2_t pos3 = {(int)(pos0.x + ((float)(pos1.y - pos0.y) / (float)(pos2.y - pos0.y)) * (float)(pos2.x - pos0.x)), pos1.y};
        _llsr_raster_triangle_2d_flat(color_buffer, llsr_bresenham_line_create(pos0, pos1), llsr_bresenham_line_create(pos0, pos3), color);
        _llsr_raster_triangle_2d_flat(color_buffer, llsr_bresenham_line_create(pos3, pos2), llsr_bresenham_line_create(pos1, pos2), color);
    }
}

void llsr_raster_triangle_3d(llsr_color_buffer_t* color_buffer, llsr_depth_buffer_t* depth_buffer, llsr_float3_t pos0, llsr_float3_t pos1, llsr_float3_t pos2, llsr_byte3_t color) {
    
}
