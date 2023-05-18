#include "pok_gfx.h"

#include <string.h>
#ifndef _NDEBUG
#include <stdio.h>
#endif

typedef enum {
    POK_LINE_INCREMENT_ANY,
    POK_LINE_INCREMENT_VERTICAL,
    POK_LINE_INCREMENT_HORIZONTAL
} Pok_LineIncrementType;

typedef struct {
    Pok_Int2 start;
    Pok_Int2 end;
    Pok_Int2 current;
    
    Pok_AttributeArray attributes;
    Pok_Float4 increments[POK_MAX_ATTRIBUTES];
    int attributeCount;

    int steps;
    Pok_LineIncrementType incrementType;

    // Data used with the bresenham algorithm.
    Pok_Int2 delta;
    Pok_Int2 dir;
    int error;
} Pok_Line;

typedef struct {
    Pok_Float3 actualPos;
    Pok_Int2 pixelPos;
} Pok_RenderPos;

typedef struct {
    Pok_Line line;

    Pok_Float3 actualStart;
    Pok_Float3 actualEnd;
    Pok_Float3 actualCurrent;

    Pok_Float3 posIncrement;
} Pok_Line3D;

Pok_ColorBuffer* Pok_ColorBufferCreate(int w, int h) {
    Pok_ColorBuffer* buf = Pok_AllocHeap(sizeof(Pok_ColorBuffer));
    buf->w = w;
    buf->h = h;
    buf->data = Pok_AllocHeap(w * h * sizeof(Pok_Byte3));
    return buf;
}

void Pok_ColorBufferDestroy(Pok_ColorBuffer* colorBuffer) {
    Pok_Free(colorBuffer->data);
    Pok_Free(colorBuffer);
}

void Pok_ColorBufferResize(Pok_ColorBuffer* colorBuffer, int w, int h) {
    colorBuffer->w = w;
    colorBuffer->h = h;
    Pok_Free(colorBuffer->data);
    colorBuffer->data = Pok_AllocHeap(w * h * sizeof(Pok_Byte3));
}

Pok_Byte3* Pok_ColorBufferAt(Pok_ColorBuffer* colorBuffer, int x, int y) {
    return &colorBuffer->data[y * colorBuffer->w + x];
}

void Pok_ColorBufferClear(Pok_ColorBuffer* colorBuffer, Pok_Byte3 color) {
    for (int i = 0; i < colorBuffer->w * colorBuffer->h; i++) {
        colorBuffer->data[i] = color;
    }
}

Pok_DepthBuffer* Pok_DepthBufferCreate(int w, int h) {
    Pok_DepthBuffer* buf = Pok_AllocHeap(sizeof(Pok_DepthBuffer));
    buf->w = w;
    buf->h = h;
    buf->data = Pok_AllocHeap(w * h * sizeof(float));
    return buf;
}

void Pok_DepthBufferDestroy(Pok_DepthBuffer* depthBuffer) {
    Pok_Free(depthBuffer->data);
    Pok_Free(depthBuffer);
}

void Pok_DepthBufferResize(Pok_DepthBuffer* depthBuffer, int w, int h) {
    depthBuffer->w = w;
    depthBuffer->h = h;
    Pok_Free(depthBuffer->data);
    depthBuffer->data = Pok_AllocHeap(w * h * sizeof(float));
}

float* Pok_DepthBufferAt(Pok_DepthBuffer* depthBuffer, int x, int y) {
    return &depthBuffer->data[y * depthBuffer->w + x];
}

void Pok_DepthBufferClear(Pok_DepthBuffer* depthBuffer) {
    for (int i = 0; i < depthBuffer->w * depthBuffer->h; i++) {
        depthBuffer->data[i] = 1;
    }
}

int Pok_GetPixelSize(Pok_ColorDepth colorDepth) {
    switch (colorDepth) {
    case POK_R5G6B5:
    case POK_R5G5B5A1:
        return 2;
    case POK_R8G8B8:
        return 3;
    case POK_R8G8B8A8:
        return 4;
    default:
        POK_ASSERT(false, "Unhandled case.");
        return -1;
    }
}

int Pok_GetChannels(Pok_ColorDepth colorDepth) {
    switch (colorDepth) {
    case POK_R5G6B5:
    case POK_R8G8B8:
        return 3;
    case POK_R5G5B5A1:
    case POK_R8G8B8A8:
        return 4;
    default:
        POK_ASSERT(false, "Unhandled case.");
        return -1;
    }
}

Pok_Image* Pok_ImageCreate(Pok_ColorDepth colorDepth, int w, int h) {
    Pok_Image* img = Pok_AllocHeap(sizeof(Pok_Image));
    img->w = w;
    img->h = h;
    img->colorDepth = colorDepth;
    img->data = Pok_AllocHeap(w * h * Pok_GetPixelSize(colorDepth));
    return img;
}

void Pok_ImageDestroy(Pok_Image* image) {
    Pok_Free(image->data);
    Pok_Free(image);
}

Pok_Byte* Pok_ImageAt(Pok_Image* image, int x, int y) {
    POK_ASSERT(x >= 0 && x < image->w && y >= 0 && y < image->h, "Outside bounds.");
    return &image->data[(y * image->w + x) * Pok_GetPixelSize(image->colorDepth)];
}

Pok_Byte* Pok_ImageSample(Pok_Image* image, float u, float v) {
    return Pok_ImageAt(
        image, 
        Pok_Clamp(u, 0, 1) * (image->w - 1), 
        Pok_Clamp(v, 0, 1) * (image->h - 1));
}

void Pok_LineInit(Pok_Line* line, 
                  Pok_Int2 start,
                  Pok_Int2 end,
                  const Pok_AttributeArray* startAttributes, 
                  const Pok_AttributeArray* endAttributes, 
                  int attributeCount,
                  Pok_LineIncrementType incrementType) {
    line->start = start;
    line->end = end;
    line->current = start;
    
    line->attributeCount = attributeCount;
    line->incrementType = incrementType;

    if (attributeCount > 0 && startAttributes) {
        memcpy(&line->attributes, startAttributes, sizeof(Pok_AttributeArray));
    }

    // Calculate number of steps.

    Pok_Int2 deltaAbs = {POK_ABS(end.x - start.x), POK_ABS(end.y - start.y)};

    switch (incrementType) {
        case POK_LINE_INCREMENT_ANY:        line->steps = (deltaAbs.x >= deltaAbs.y) ? deltaAbs.x : deltaAbs.y; break;
        case POK_LINE_INCREMENT_VERTICAL:   line->steps = deltaAbs.y; break;
        case POK_LINE_INCREMENT_HORIZONTAL: line->steps = deltaAbs.x; break;
        default:                            assert(!"Unhandled enum.");
    }

    // Setup increments.

    if (line->steps > 0) {
        for (int i = 0; i < attributeCount; i++) {
            for (int j = 0; j < line->attributes.attributes[i].floatCount; j++) {
                float delta = endAttributes->attributes[i].data.at[j] - startAttributes->attributes[i].data.at[j];
                line->increments[i].at[j] = delta / line->steps;
            }
        }
    }

    // Setup bresenham members.

    line->delta = (Pok_Int2){deltaAbs.x, -deltaAbs.y};
    line->dir = (Pok_Int2){start.x < end.x ? 1 : -1, start.y < end.y ? 1 : -1};
    line->error = line->delta.x + line->delta.y;
}

void Pok_Line3DInit(Pok_Line3D* line,
                    Pok_RenderPos start,
                    Pok_RenderPos end,
                    const Pok_AttributeArray* startAttributes,
                    const Pok_AttributeArray* endAttributes,
                    int attributeCount,
                    Pok_LineIncrementType incrementType) {
    Pok_LineInit(&line->line, start.pixelPos, end.pixelPos, startAttributes, endAttributes, attributeCount, incrementType);

    line->actualStart = start.actualPos;
    line->actualEnd = end.actualPos;
    line->actualCurrent = start.actualPos;

    // Setup position increment.

    for (int i = 0; i < 3; i++) {
        line->posIncrement.at[i] = (end.actualPos.at[i] - start.actualPos.at[i]) / line->line.steps;
    }
}

void Pok_AttributeArrayLerp(Pok_AttributeArray* out, const Pok_AttributeArray* a, const Pok_AttributeArray* b, int attributeCount, float amount) {
    POK_ASSERT(attributeCount <= POK_MAX_ATTRIBUTES, "Too many attributes.");
    
    for (int i = 0; i < attributeCount; i++) {
        POK_ASSERT(a->attributes[i].floatCount == b->attributes[i].floatCount, "Attribute arrays not matching.");

        for (int j = 0; j < a->attributes[i].floatCount; j++) {
            out->attributes[i].data.at[j] = Pok_Lerp(a->attributes[i].data.at[j], b->attributes[i].data.at[j], amount);
            out->attributes[i].floatCount = a->attributes[i].floatCount;
        }
    }
}

bool Pok_LineStep(Pok_Line* line) {
    Pok_Int2 prev = line->current;

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

        if (line->incrementType == POK_LINE_INCREMENT_HORIZONTAL) {
            if (prev.x != line->current.x) {
                break;
            }
        }
        else if (line->incrementType == POK_LINE_INCREMENT_VERTICAL) {
            if (prev.y != line->current.y) {
                break;
            }
        }
        else {
            break;
        }
    }

    // Increment attributes.
    for (int i = 0; i < line->attributeCount; i++) {
        for (int j = 0; j < line->attributes.attributes[i].floatCount; j++) {
            line->attributes.attributes[i].data.at[j] += line->increments[i].at[j];
        }
    }

    return true;
}

bool Pok_Line3DStep(Pok_Line3D* line) {
    if (Pok_LineStep(&line->line)) {
        for (int i = 0; i < 3; i++) {
            line->actualCurrent.at[i] += line->posIncrement.at[i];
        }
        return true;
    }
    return false;
}

void Pok_RenderTriangle3DBetweenTwoVerticalLines(Pok_ColorBuffer* colorBuffer, 
                                                 Pok_DepthBuffer* depthBuffer,
                                                 Pok_Line3D lineA,
                                                 Pok_Line3D lineB,
                                                 Pok_PixelShaderCallback pixelShader, 
                                                 void* userData) {    
    POK_ASSERT(lineA.line.incrementType == POK_LINE_INCREMENT_VERTICAL && lineB.line.incrementType == POK_LINE_INCREMENT_VERTICAL, 
        "Lines must be set to vertical increment.");
    
    if (lineA.line.start.x > lineB.line.start.x || lineA.line.end.x > lineB.line.end.x) {
        POK_SWAP(Pok_Line3D, lineA, lineB);
    }
    
    do {
        POK_ASSERT(lineA.line.current.y == lineB.line.current.y, "Something went wrong.");

        // Create a horizontal line between line a and line b.
        Pok_Line3D lineX;
        Pok_Line3DInit(&lineX, 
                       (Pok_RenderPos){.pixelPos = lineA.line.current, .actualPos = lineA.actualCurrent}, 
                       (Pok_RenderPos){.pixelPos = lineB.line.current, .actualPos = lineB.actualCurrent}, 
                       &lineA.line.attributes, 
                       &lineB.line.attributes, 
                       lineA.line.attributeCount, 
                       POK_LINE_INCREMENT_HORIZONTAL);

        do {
            int x = lineX.line.current.x;
            int y = lineX.line.current.y;
            
            // Update depth buffer.
            if (depthBuffer) {
                POK_ASSERT(x >= 0 && x < depthBuffer->w && y >= 0 && y < depthBuffer->h, "Outside bounds.");

                float z = lineX.actualCurrent.z;

                if (z < *Pok_DepthBufferAt(depthBuffer, x, y)) {
                    *Pok_DepthBufferAt(depthBuffer, x, y) = z;
                }
                // If pixel is invisible, skip color buffer update.
                else {
                    continue;
                }
            }

            // Update color buffer.
            if (colorBuffer) {
                POK_ASSERT(x >= 0 && x < colorBuffer->w && y >= 0 && y < colorBuffer->h, "Outside bounds.");

                // Invoke pixel shader.
                Pok_Byte3 color = 
                    pixelShader(lineX.line.current, &lineX.line.attributes, userData);

                *Pok_ColorBufferAt(colorBuffer, x, y) = color;
            }
        }
        while (Pok_Line3DStep(&lineX));
    }
    while (Pok_Line3DStep(&lineA) && Pok_Line3DStep(&lineB));
}

void Pok_RenderTriangle3DFlatTop(Pok_ColorBuffer* colorBuffer, 
                                 Pok_DepthBuffer* depthBuffer,
                                 Pok_RenderPos topLeftPos,
                                 Pok_RenderPos topRightPos,
                                 Pok_RenderPos bottomPos,
                                 const Pok_AttributeArray* topLeftAttribs,
                                 const Pok_AttributeArray* topRightAttribs,
                                 const Pok_AttributeArray* bottomAttribs,
                                 int attributeCount,
                                 Pok_PixelShaderCallback pixelShader, 
                                 void* userData) {
    POK_ASSERT(topLeftPos.pixelPos.y == topRightPos.pixelPos.y, "Y of top coordinates must be equal.");
    
    Pok_Line3D lineA;
    Pok_Line3DInit(&lineA, topLeftPos, bottomPos, topLeftAttribs, bottomAttribs, 
        attributeCount, POK_LINE_INCREMENT_VERTICAL);
    
    Pok_Line3D lineB;
    Pok_Line3DInit(&lineB, topRightPos, bottomPos, topRightAttribs, bottomAttribs, 
        attributeCount, POK_LINE_INCREMENT_VERTICAL);

    Pok_RenderTriangle3DBetweenTwoVerticalLines(colorBuffer, depthBuffer, 
        lineA, lineB, pixelShader, userData);
}

void Pok_RenderTriangle3DFlatBottom(Pok_ColorBuffer* colorBuffer, 
                                         Pok_DepthBuffer* depthBuffer,
                                         Pok_RenderPos topPos,
                                         Pok_RenderPos bottomLeftPos,
                                         Pok_RenderPos bottomRightPos,
                                         const Pok_AttributeArray* topAttribs,
                                         const Pok_AttributeArray* bottomLeftAttribs,
                                         const Pok_AttributeArray* bottomRightAttribs,
                                         int attributeCount,
                                         Pok_PixelShaderCallback pixelShader, 
                                         void* userData) {
    POK_ASSERT(bottomLeftPos.pixelPos.y == bottomRightPos.pixelPos.y, "Y of bottom coordinates must be equal.");
    
    Pok_Line3D lineA;
    Pok_Line3DInit(&lineA, topPos, bottomLeftPos, topAttribs, bottomLeftAttribs, 
        attributeCount, POK_LINE_INCREMENT_VERTICAL);
    
    Pok_Line3D lineB;
    Pok_Line3DInit(&lineB, topPos, bottomRightPos, topAttribs, bottomRightAttribs, 
        attributeCount, POK_LINE_INCREMENT_VERTICAL);

    Pok_RenderTriangle3DBetweenTwoVerticalLines(colorBuffer, depthBuffer, 
        lineA, lineB, pixelShader, userData);
}

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
                          void* userData) {
    // Cull against all 6 planes in homogeneous clip space.
    // If all vertices are outside view, discard triangle.
    if ((pos0.x > pos0.w && pos1.x > pos1.w && pos2.x > pos2.w) ||
        (pos0.x < -pos0.w && pos1.x < -pos1.w && pos2.x < -pos2.w) ||
        
        (pos0.y > pos0.w && pos1.y > pos1.w && pos2.y > pos2.w) ||
        (pos0.y < -pos0.w && pos1.y < -pos1.w && pos2.y < -pos2.w) ||
        
        (pos0.z > pos0.w && pos1.z > pos1.w && pos2.z > pos2.w) ||
        (pos0.z < 0 && pos1.z < 0 && pos2.z < 0)) {
        return;
    }

    // TODO: Clip.
    
    // W division (homogeneous clip space -> NDC space).
    for (int i = 0; i < 3; i++) {
        pos0.at[i] /= pos0.w;
        pos1.at[i] /= pos1.w;
        pos2.at[i] /= pos2.w;
    }

    // Viewport transformation. 
    // Scale from [-1, 1] to color buffer size.
    pos0.x = (pos0.x + 1) / 2.f * colorBuffer->w;
    pos0.y = (pos0.y + 1) / 2.f * colorBuffer->h;
    pos1.x = (pos1.x + 1) / 2.f * colorBuffer->w;
    pos1.y = (pos1.y + 1) / 2.f * colorBuffer->h;
    pos2.x = (pos2.x + 1) / 2.f * colorBuffer->w;
    pos2.y = (pos2.y + 1) / 2.f * colorBuffer->h;
    
    // Floor X and Y, otherwise there's missing pixel artifacts.
    for (int i = 0; i < 2; i++) {
        pos0.at[i] = floorf(pos0.at[i]);
        pos1.at[i] = floorf(pos1.at[i]);
        pos2.at[i] = floorf(pos2.at[i]);
    }
    
    Pok_RenderPos p0 = {.actualPos = {pos0.x, pos0.y, pos0.z}, .pixelPos = {pos0.x, pos0.y}};
    Pok_RenderPos p1 = {.actualPos = {pos1.x, pos1.y, pos1.z}, .pixelPos = {pos1.x, pos1.y}};
    Pok_RenderPos p2 = {.actualPos = {pos2.x, pos2.y, pos2.z}, .pixelPos = {pos2.x, pos2.y}};

    // HACK: discard triangle if outside view.
    /*if (p0.pixelPos.x < 0 || p0.pixelPos.x >= colorBuffer->w || p0.pixelPos.y < 0 || p0.pixelPos.y >= colorBuffer->h ||
        p1.pixelPos.x < 0 || p1.pixelPos.x >= colorBuffer->w || p1.pixelPos.y < 0 || p1.pixelPos.y >= colorBuffer->h ||
        p2.pixelPos.x < 0 || p2.pixelPos.x >= colorBuffer->w || p2.pixelPos.y < 0 || p2.pixelPos.y >= colorBuffer->h) {
        return;
    }*/
    
    // Sort vertices by height.
    if (p0.pixelPos.y > p1.pixelPos.y) {
        POK_SWAP(Pok_RenderPos, p0, p1);
        POK_SWAP(const Pok_AttributeArray*, attributes0, attributes1);
    }
    if (p0.pixelPos.y > p2.pixelPos.y) {
        POK_SWAP(Pok_RenderPos, p0, p2);
        POK_SWAP(const Pok_AttributeArray*, attributes0, attributes2);
    }
    if (p1.pixelPos.y > p2.pixelPos.y) {
        POK_SWAP(Pok_RenderPos, p1, p2);
        POK_SWAP(const Pok_AttributeArray*, attributes1, attributes2);
    }

    // Check if the top of the triangle is flat.
    if (p0.pixelPos.y == p1.pixelPos.y) {
        Pok_RenderTriangle3DFlatTop(colorBuffer, depthBuffer, 
            p0, p1, p2, attributes0, attributes1, attributes2, attributeCount,
            pixelShader, userData);
    }
    // Check if the bottom is flat.
    else if (p1.pixelPos.y == p2.pixelPos.y) {
        Pok_RenderTriangle3DFlatBottom(colorBuffer, depthBuffer,
            p0, p1, p2, attributes0, attributes1, attributes2, attributeCount, 
            pixelShader, userData);
    }
    // Else split into two smaller triangles.
    else {
        float alpha_split = 
            (p1.actualPos.y - p0.actualPos.y) / (p2.actualPos.y - p0.actualPos.y);

        Pok_RenderPos p3;
        p3.actualPos.y = p1.actualPos.y;
        p3.actualPos.x = Pok_Lerp(p0.actualPos.x, p2.actualPos.x, alpha_split);
        p3.actualPos.z = Pok_Lerp(p0.actualPos.z, p2.actualPos.z, alpha_split);
        p3.pixelPos = (Pok_Int2){(int)roundf(p3.actualPos.x), (int)roundf(p3.actualPos.y)};

        Pok_AttributeArray attributes3;
        Pok_AttributeArrayLerp(&attributes3, attributes0, attributes2, attributeCount, alpha_split);

        // Bottom (flat top).
        Pok_RenderTriangle3DFlatTop(colorBuffer, depthBuffer,
            p1, p3, p2, attributes1, &attributes3, attributes2, attributeCount, 
            pixelShader, userData);
        
        // Top (flat bottom).
        Pok_RenderTriangle3DFlatBottom(colorBuffer, depthBuffer,
            p0, p1, p3, attributes0, attributes1, &attributes3, attributeCount, 
            pixelShader, userData);
    }
}

void Pok_RenderImage(Pok_ColorBuffer* colorBuffer, Pok_Image* image, Pok_Rect src, Pok_Rect dst) {    
    Pok_Float2 scaleFactor = {src.w / (float)dst.w, src.h / (float)dst.h};
    
    for (int x = 0; x < dst.w; x++) {
        for (int y = 0; y < dst.h; y++) {
            Pok_Int2 sample = {src.x + (int)floorf(x * scaleFactor.x), src.y + (int)floorf(y * scaleFactor.y)};

            Pok_Byte* srcPixel = Pok_ImageAt(image, sample.x, sample.y);
            Pok_Byte3* dstPixel = Pok_ColorBufferAt(colorBuffer, dst.x + x, dst.y + y);

            switch (image->colorDepth) {
            case POK_R8G8B8:
                for (int i = 0; i < 3; i++) {
                    dstPixel->at[i] = *(srcPixel + i);
                }
                break;
            case POK_R8G8B8A8:
                // HACK: skip pixel if not fully opaque.
                if (*(srcPixel + 3) < 255) {
                    break;
                }
                for (int i = 0; i < 3; i++) {
                    dstPixel->at[i] = *(srcPixel + i);
                }
                break;
            default:
                POK_ASSERT(false, "Unhandled case.");
            }
        }
    }
}

void Pok_RenderText(Pok_ColorBuffer* colorBuffer, const char* text, Pok_Int2 pos, Pok_Font* font, int scale) {       
    for (size_t i = 0; i < strlen(text); i++) {
        Pok_CharInfo info = font->charInfos[(int)text[i]];

        Pok_Rect dst = {pos.x + info.offset.x * scale, pos.y + info.offset.y * scale, info.src.w * scale, info.src.h * scale};
        Pok_RenderImage(colorBuffer, font->image, info.src, dst);

        pos.x += info.xAdvance * scale;
    }
}
