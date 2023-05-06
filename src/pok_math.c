#include "pok_math.h"

float Pok_Radians(float degrees) {
    return degrees * (M_PI / 180.f);
}

float Pok_Degrees(float radians) {
    return radians * (180.f / M_PI);
}

Pok_Float3 Pok_Normalize(Pok_Float3 f) {
    float len = sqrtf(f.x * f.x + f.y * f.y + f.z * f.z);
    Pok_Float3 result = {f.x / len, f.y / len, f.z / len};
    return result;
}

Pok_Float3 Pok_Cross(Pok_Float3 a, Pok_Float3 b) {
    Pok_Float3 result = {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
    return result;
}

float Pok_Dot(Pok_Float3 a, Pok_Float3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

float Pok_Lerp(float a, float b, float amount) {
    return a * (1 - amount) + (b * amount);
}

float Pok_Clamp(float f, float min, float max) {
    float t = f < min ? min : f;
    return t > max ? max : t;
}

void Pok_Mat4InitZero(Pok_Mat4* m) {
    for (int x = 0; x < 4; x++) {
        for (int y = 0; y < 4; y++) {
            m->m[x][y] = 0;
        }
    }
}

void Pok_Mat4InitIdentity(Pok_Mat4* m) {
    Pok_Mat4InitZero(m);
    for (int i = 0; i < 4; i++) {
        m->m[i][i] = 1;
    }
}

Pok_Mat4 Pok_Mat4Mul(Pok_Mat4 a, Pok_Mat4 b) {
    Pok_Mat4 result;
    Pok_Mat4InitZero(&result);

    for (int x = 0; x < 4; x++) {
        for (int y = 0; y < 4; y++) {
            for (int i = 0; i < 4; i++) {
                result.m[x][y] += a.m[x][i] * b.m[i][y];
            }
        }
    }
    return result;
}

Pok_Mat4 Pok_Mat4Transpose(Pok_Mat4 m) {
    Pok_Mat4 r;
    for (int x = 0; x < 4; x++) {
        for (int y = 0; y < 4; y++) {
            r.m[y][x] = m.m[x][y];
        }
    }
    return r;
}

Pok_Mat4 Pok_Mat4Inverse(Pok_Mat4 m) {
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
    
    POK_ASSERT(det > 0, "Matrix is not invertible.");
    det = 1 / det;

    Pok_Mat4 result = {
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

Pok_Mat4 Pok_Mat4Translate(Pok_Mat4 mat4, Pok_Float3 f) {
    Pok_Mat4 m;
    Pok_Mat4InitIdentity(&m);

    m.m[3][0] = f.x;
    m.m[3][1] = f.y;
    m.m[3][2] = f.z;

    return Pok_Mat4Mul(mat4, m);
}

Pok_Mat4 Pok_Mat4RotateX(Pok_Mat4 m, float a) {
    float s = sinf(a);
    float c = cosf(a);

    Pok_Mat4 x;
    Pok_Mat4InitIdentity(&x);

    x.m[1][1] = c;
    x.m[1][2] = -s;
    x.m[2][1] = s;
    x.m[2][2] = c;

    return Pok_Mat4Mul(m, x);
}

Pok_Mat4 Pok_Mat4RotateY(Pok_Mat4 m, float a) {
    float s = sinf(a);
    float c = cosf(a);

    Pok_Mat4 y;
    Pok_Mat4InitIdentity(&y);

    y.m[0][0] = c;
    y.m[0][2] = s;
    y.m[2][0] = -s;
    y.m[2][2] = c;

    return Pok_Mat4Mul(m, y);
}

Pok_Mat4 Pok_Mat4RotateZ(Pok_Mat4 m, float a) {
    float s = sinf(a);
    float c = cosf(a);

    Pok_Mat4 z;
    Pok_Mat4InitIdentity(&z);

    z.m[0][0] = c;
    z.m[0][1] = -s;
    z.m[1][0] = s;
    z.m[1][1] = c;

    return Pok_Mat4Mul(m, z);
}

Pok_Mat4 Pok_Mat4Scale(Pok_Mat4 mat4, Pok_Float3 f) {
    Pok_Mat4 m;
    Pok_Mat4InitIdentity(&m);

    m.m[0][0] = f.x;
    m.m[1][1] = f.y;
    m.m[2][2] = f.z;

    return Pok_Mat4Mul(mat4, m);
}

Pok_Float3 Pok_Mat3MulFloat3(Pok_Mat3 m, Pok_Float3 f) {
    Pok_Float3 result;
    for (int x = 0; x < 3; x++) {
        float i = 0;
        for (int y = 0; y < 3; y++) {
            i += m.m[x][y] * f.at[y];
        }
        result.at[x] = i;
    }
    return result;
}

Pok_Float3 Pok_Float3MulMat3(Pok_Float3 f, Pok_Mat3 m) {
    Pok_Float3 result;
    for (int y = 0; y < 3; y++) {
        float i = 0;
        for (int x = 0; x < 3; x++) {
            i += m.m[x][y] * f.at[x];
        }
        result.at[y] = i;
    }
    return result;
}

Pok_Float3 Pok_Mat4MulFloat3(Pok_Mat4 m, Pok_Float3 f) {
    Pok_Float4 f4;
    for (int x = 0; x < 4; x++) {
        float i = 0;
        for (int y = 0; y < 3; y++) {
            i += m.m[x][y] * f.at[y];
        }
        f4.at[x] = i + m.m[x][3];
    }
    
    if (f4.w != 0) {
        f4.x /= f4.w;
        f4.y /= f4.w;
        f4.z /= f4.w;
    }

    Pok_Float3 result = {f4.x, f4.y, f4.z};
    return result;
}

Pok_Float4 Pok_Mat4MulFloat4(Pok_Mat4 m, Pok_Float4 f) {
    Pok_Float4 result;
    for (int x = 0; x < 4; x++) {
        float i = 0;
        for (int y = 0; y < 4; y++) {
            i += m.m[x][y] * f.at[y];
        }
        result.at[x] = i;
    }
    return result;
}

Pok_Float3 Pok_Float3MulMat4(Pok_Float3 f, Pok_Mat4 m) {
    Pok_Float4 f4;
    for (int y = 0; y < 4; y++) {
        float i = 0;
        for (int x = 0; x < 3; x++) {
            i += m.m[x][y] * f.at[x];
        }
        f4.at[y] = i + m.m[3][y];
    }
    
    if (f4.w != 0) {
        f4.x /= f4.w;
        f4.y /= f4.w;
        f4.z /= f4.w;
    }

    Pok_Float3 result = {f4.x, f4.y, f4.z};
    return result;
}

Pok_Float4 Pok_Float4MulMat4(Pok_Float4 f, Pok_Mat4 m) {
    Pok_Float4 result;
    for (int y = 0; y < 4; y++) {
        float i = 0;
        for (int x = 0; x < 4; x++) {
            i += m.m[x][y] * f.at[x];
        }
        result.at[y] = i;
    }
    return result;
}

Pok_Mat3 Pok_Mat4ToMat3(Pok_Mat4 m) {
    Pok_Mat3 r;
    for (int x = 0; x < 3; x++) {
        for (int y = 0; y < 3; y++) {
            r.m[x][y] = m.m[x][y];
        }
    }
    return r;
}

Pok_Mat4 Pok_LookAt(Pok_Float3 pos, Pok_Float3 target, Pok_Float3 up) {
    Pok_Float3 diff = {target.x - pos.x, target.y - pos.y, target.z - pos.z};
    Pok_Float3 forward = Pok_Normalize(diff);
    Pok_Float3 right = Pok_Normalize(Pok_Cross(forward, up));
    Pok_Float3 local_up = Pok_Normalize(Pok_Cross(right, forward));

    Pok_Mat4 m;
    Pok_Mat4InitIdentity(&m);
    m.m[0][0] = right.x;
    m.m[1][0] = right.y;
    m.m[2][0] = right.z;
    m.m[0][1] = local_up.x;
    m.m[1][1] = local_up.y;
    m.m[2][1] = local_up.z;
    m.m[0][2] = -forward.x;
    m.m[1][2] = -forward.y;
    m.m[2][2] = -forward.z;
    m.m[3][0] = -Pok_Dot(right, pos);
    m.m[3][1] = -Pok_Dot(local_up, pos);
    m.m[3][2] = Pok_Dot(forward, pos);
    return m;
}

Pok_Mat4 Pok_Perspective(float aspect, float fov, float near, float far) {    
    Pok_Mat4 m;
    Pok_Mat4InitZero(&m);

    float half_tan = tanf(fov / 2);

    m.m[0][0] = 1 / (half_tan * aspect);
    m.m[1][1] = 1 / half_tan;
    m.m[2][2] = -(far + near) / (far - near);
    m.m[2][3] = -1;
    m.m[3][2] = -(2 * far * near) / (far - near);
    return m;
}
