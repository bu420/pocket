#ifndef POK_MATH_H
#define POK_MATH_H

#include <pok_core.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define POK_ABS(x) ((x) < 0 ? -(x) : (x))

typedef unsigned char Pok_Byte;

typedef union {
    struct {
        Pok_Byte x, y, z;
    };
    struct {
        Pok_Byte r, g, b;
    };
    Pok_Byte at[3];
} Pok_Byte3;

typedef union {
    struct {
        Pok_Byte x, y, z, w;
    };
    struct {
        Pok_Byte r, g, b, a;
    };
    Pok_Byte at[4];
} Pok_Byte4;

typedef union {
    struct {
        int x, y;
    };
    struct {
        int r, g;
    };
    int at[2];
} Pok_Int2;

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
    float at[2];
} Pok_Float2;

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
    float at[3];
} Pok_Float3;

typedef union {
    struct {
        float x, y, z, w;
    };
    struct {
        float r, g, b, a;
    };
    float at[4];
} Pok_Float4;

typedef union {
    float m[3][3];
    struct {
        float m00, m10, m20, m01, m11, m21, m02, m12, m22;
    };
} Pok_Mat3;

typedef union {
    float m[4][4];
    struct {
        float m00, m10, m20, m30, m01, m11, m21, m31, m02, m12, m22, m32, m03, m13, m23, m33;
    };
} Pok_Mat4;

typedef struct {
    int x, y, w, h;
} Pok_Rect;

Pok_Float3 Pok_Normalize(Pok_Float3 f);
Pok_Float3 Pok_Cross(Pok_Float3 a, Pok_Float3 b);
float Pok_Dot(Pok_Float3 a, Pok_Float3 b);
float Pok_Lerp(float a, float b, float amount);
float Pok_Clamp(float f, float min, float max);

void Pok_Mat4InitZero(Pok_Mat4* m);
void Pok_Mat4InitIdentity(Pok_Mat4* m);
Pok_Mat4 Pok_Mat4Mul(Pok_Mat4 a, Pok_Mat4 b);
Pok_Mat4 Pok_Mat4Transpose(Pok_Mat4 m);
Pok_Mat4 Pok_Mat4Inverse(Pok_Mat4 m);
Pok_Mat4 Pok_Mat4Translate(Pok_Mat4 m, Pok_Float3 f);
Pok_Mat4 Pok_Mat4RotateX(Pok_Mat4 m, float a);
Pok_Mat4 Pok_Mat4RotateY(Pok_Mat4 m, float a);
Pok_Mat4 Pok_Mat4RotateZ(Pok_Mat4 m, float a);
Pok_Mat4 Pok_Mat4Scale(Pok_Mat4 m, Pok_Float3 f);
Pok_Float3 Pok_Mat3MulFloat3(Pok_Mat3 m, Pok_Float3 f);
Pok_Float3 Pok_Float3MulMat3(Pok_Float3 f, Pok_Mat3 m);
Pok_Float3 Pok_Mat4MulFloat3(Pok_Mat4 m, Pok_Float3 f);
Pok_Float4 Pok_Mat4MulFloat4(Pok_Mat4 m, Pok_Float4 f);
Pok_Float3 Pok_Float3MulMat4(Pok_Float3 f, Pok_Mat4 m);
Pok_Float4 Pok_Float4MulMat4(Pok_Float4 f, Pok_Mat4 m);
Pok_Mat3 Pok_Mat4ToMat3(Pok_Mat4 m);

Pok_Mat4 Pok_LookAt(Pok_Float3 pos, Pok_Float3 target, Pok_Float3 up);
Pok_Mat4 Pok_Perspective(float aspect, float fov, float near, float far);
// TODO: implement.
Pok_Mat4 Pok_Ortho();

#endif
