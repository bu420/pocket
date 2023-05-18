#ifndef POK_GFX_H
#define POK_GFX_H

#include <pok_core.h>
#include <pok_math.h>

#ifndef POK_MAX_ATTRIBUTES
#define POK_MAX_ATTRIBUTES 8
#endif

#define POK_ATTRIB_ARRAY(...) &(Pok_AttributeArray){.attributes = {__VA_ARGS__}}
#define POK_ATTRIB_1(a) (Pok_Attribute){.data = {.x = a.x}, .floatCount = 1}
#define POK_ATTRIB_2(a) (Pok_Attribute){.data = {.x = a.x, .y = a.y}, .floatCount = 2}
#define POK_ATTRIB_3(a) (Pok_Attribute){.data = {.x = a.x, .y = a.y, .z = a.z}, .floatCount = 3}
#define POK_ATTRIB_4(a) (Pok_Attribute){.data = {.x = a.x, .y = a.y, .z = a.z, .w = a.w}, .floatCount = 4}

#define POK_ATTRIB_TO_FLOAT1(a) (a.data.x)
#define POK_ATTRIB_TO_FLOAT2(a) (Pok_Float2){a.data.x, a.data.y}
#define POK_ATTRIB_TO_FLOAT3(a) (Pok_Float3){a.data.x, a.data.y, a.data.z}
#define POK_ATTRIB_TO_FLOAT4(a) (Pok_Float4){a.data.x, a.data.y, a.data.z, a.data.w}

typedef struct {
    int w, h;
    Pok_Byte3* data;
} Pok_ColorBuffer;

typedef struct {
    int w, h;
    float* data;
} Pok_DepthBuffer;

typedef enum {
    POK_R5G6B5      = 1 << 0,
    POK_R5G5B5A1    = 1 << 1, 
    POK_R8G8B8      = 1 << 2,
    POK_R8G8B8A8    = 1 << 3
} Pok_ColorDepth;

typedef struct {
    int w, h;
    Pok_ColorDepth colorDepth;
    Pok_Byte* data;
} Pok_Image;

typedef struct {
    Pok_Float4 data;
    int floatCount;
} Pok_Attribute;

typedef struct {
    Pok_Attribute attributes[POK_MAX_ATTRIBUTES];
} Pok_AttributeArray;

typedef struct {
    int positionIndices[3];
    int texCoordIndices[3];
    int normalIndices[3];
} Pok_Face;

typedef struct {
    Pok_Float3* positions;
    int positionCount;
    Pok_Float2* tex_coords;
    int texCoordCount;
    Pok_Float3* normals;
    int normalCount;
    Pok_Face* faces;
    int faceCount;
} Pok_Mesh;

typedef struct {
    int id;
    Pok_Rect src;
    Pok_Int2 offset;
    int xAdvance;
} Pok_CharInfo;

typedef struct {
    Pok_Image* image;
    Pok_CharInfo* charInfos;
    int charInfoCount;
    int size;
} Pok_Font;

typedef Pok_Byte3 (*Pok_PixelShaderCallback)(Pok_Int2 pixelPos, const Pok_AttributeArray* attributes, void* userData);

Pok_ColorBuffer* Pok_ColorBufferCreate(int w, int h);
void Pok_ColorBufferDestroy(Pok_ColorBuffer* colorBuffer);
void Pok_ColorBufferResize(Pok_ColorBuffer* colorBuffer, int w, int h);
Pok_Byte3* Pok_ColorBufferAt(Pok_ColorBuffer* colorBuffer, int x, int y);
void Pok_ColorBufferClear(Pok_ColorBuffer* colorBuffer, Pok_Byte3 color);

Pok_DepthBuffer* Pok_DepthBufferCreate(int w, int h);
void Pok_DepthBufferDestroy(Pok_DepthBuffer* depthBuffer);
void Pok_DepthBufferResize(Pok_DepthBuffer* depthBuffer, int w, int h);
float* Pok_DepthBufferAt(Pok_DepthBuffer* depthBuffer, int x, int y);
void Pok_DepthBufferClear(Pok_DepthBuffer* depthBuffer);

int Pok_GetPixelSize(Pok_ColorDepth colorDepth);
int Pok_GetChannels(Pok_ColorDepth colorDepth);

Pok_Image* Pok_ImageCreate(Pok_ColorDepth colorDepth, int w, int h);
void Pok_ImageDestroy(Pok_Image* image);
/**
 * @brief Retuns the address of the first byte of the pixel.
 */
Pok_Byte* Pok_ImageAt(Pok_Image* image, int x, int y);
Pok_Byte* Pok_ImageSample(Pok_Image* image, float u, float v);

/**
 * @param colorBuffer Optional.
 * @param depthBuffer Optional.
 */
void Pok_RenderTriangle3D(Pok_ColorBuffer* colorBuffer, 
                          Pok_DepthBuffer* depthBuffer, 
                          Pok_Float4 pos0, 
                          Pok_Float4 pos1, 
                          Pok_Float4 pos2, 
                          const Pok_AttributeArray* attributes0,
                          const Pok_AttributeArray* attributes1,
                          const Pok_AttributeArray* attributes2,
                          int attributeCount,
                          Pok_PixelShaderCallback pixelShader, 
                          void* userData);

void Pok_RenderImage(Pok_ColorBuffer* colorBuffer, Pok_Image* image, Pok_Rect src, Pok_Rect dst);

void Pok_RenderText(Pok_ColorBuffer* colorBuffer, const char* text, Pok_Int2 pos, Pok_Font* font, int scale);

Pok_Image* Pok_ImageLoadBMP(const char* path, Pok_ColorDepth colorDepth);

void Pok_SaveBMP(const char* path, const Pok_ColorBuffer* colorBuffer);

Pok_Mesh* Pok_MeshLoadOBJ(const char* path);
void Pok_MeshDestroy(Pok_Mesh* mesh);

/**
 * @brief Loads a corresponding info file to the bitmap font image.
 */
Pok_Font* Pok_FontLoad(Pok_Image* image, const char* infoPath);
void Pok_FontDestroy(Pok_Font* font);

#endif
