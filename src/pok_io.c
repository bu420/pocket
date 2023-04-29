#include "pok_gfx.h"

#include <pok_util.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

Pok_Image* Pok_ImageLoadBMP(const char* path, Pok_ColorDepth colorDepth) {
    FILE* file = fopen(path, "rb");

    if (!file) {
        return NULL;
    }

    Pok_Byte header[54];
    fread(header, 1, 54, file);

    int fileSize = (header[2]) | (header[3] << 8) | (header[4] << 16) | (header[5] << 24);
    int headerSize = (header[10]) | (header[11] << 8) | (header[12] << 16) | (header[13] << 24);
    int width = (header[18]) | (header[19] << 8) | (header[20] << 16) | (header[21] << 24);
    int height = (header[22]) | (header[23] << 8) | (header[24] << 16) | (header[25] << 24);
    int pixelSize = header[28] / 8;
    int paddingSize = (4 - width * pixelSize % 4) % 4;

    POK_ASSERT(pixelSize == Pok_GetPixelSize(colorDepth), "Incompatible color depths.");

    Pok_Image* image = Pok_ImageCreate(colorDepth, width, height);

    Pok_Byte* pixelData = Pok_AllocHeap(fileSize - headerSize);
    fseek(file, headerSize, SEEK_SET);
    fread(pixelData, 1, fileSize - headerSize, file);
    fclose(file);

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int i = ((height - 1 - y) * (width + paddingSize) + x) * pixelSize;

            switch (colorDepth) {
            case POK_R8G8B8:
                for (int j = 0; j < 3; j++) {
                    *(Pok_ImageAt(image, x, y) + j) = *(pixelData + i + (2 - j));
                }
                break;
            case POK_R8G8B8A8:
                // Alpha.
                *(Pok_ImageAt(image, x, y) + 3) = *(pixelData + i + 3);
                // RGB.
                for (int j = 0; j < 3; j++) {
                    *(Pok_ImageAt(image, x, y) + j) = *(pixelData + i + (2 - j));
                }
                break;
            default:
                POK_ASSERT(false, "Unhandled case.");
            }
        }
    }
    Pok_Free(pixelData);

    return image;
}

void Pok_SaveBMP(const char* path, const Pok_ColorBuffer* colorBuffer) {
    const Pok_ColorBuffer* buf = colorBuffer;
    
    Pok_Byte padding[] = {0, 0, 0};
    int paddingSize = (4 - buf->w * 3 % 4) % 4;
    int stride = buf->w * 3 + paddingSize;
    int fileSize = 54 + stride * buf->h;

    Pok_Byte header[54] = {
        'B', 'M',
        fileSize, (fileSize >> 8), (fileSize >> 16), (fileSize >> 24),
        0, 0, 0, 0,
        54, 0, 0, 0,
        40, 0, 0, 0,
        buf->w, (buf->w >> 8), (buf->w >> 16), (buf->w >> 24),
        buf->h, (buf->h >> 8), (buf->h >> 16), (buf->h >> 24),
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
            Pok_Byte3 pixel = buf->data[i];
            Pok_Byte3 color = {.r = pixel.b, .g = pixel.g, .b = pixel.r};
            fwrite(&color, 3, 1, file);
        }

        fwrite(padding, 1, paddingSize, file);
    }

    fclose(file);
}

Pok_Mesh* Pok_MeshLoadOBJ(const char* path) {
    Pok_Mesh* mesh = Pok_AllocHeap(sizeof(Pok_Mesh));

    mesh->faceCount     = 0;
    mesh->positionCount = 0;
    mesh->texCoordCount = 0;
    mesh->normalCount   = 0;

    int faceMax     = 128;
    int positionMax = 128;
    int texCoordMax = 128;
    int normalMax   = 128;

    mesh->faces      = Pok_AllocHeap(faceMax     * sizeof(Pok_Face));
    mesh->positions  = Pok_AllocHeap(positionMax * sizeof(Pok_Float3));
    mesh->tex_coords = Pok_AllocHeap(texCoordMax * sizeof(Pok_Float2));
    mesh->normals    = Pok_AllocHeap(normalMax   * sizeof(Pok_Float3));

    // Load OBJ file.
    char* fileContent = Pok_LoadFileBinaryBlob(path);
    if (!fileContent) {
        return NULL;
    }

    // Split into lines.
    int lineCount;
    char** lines = Pok_StringSplit(fileContent, "\n", false, &lineCount);

    for (int i = 0; i < lineCount; i++) {
        // Split each line into tokens.
        int tokenCount;
        char** tokens = Pok_StringSplit(lines[i], " ", false, &tokenCount);

        if (tokenCount == 0) {
            continue;
        }

        // Check the first token of each line.
        // The first token of relevant lines contains either 'v' for 
        // vertex position, 'vt' for vertex texture coordinate, 
        // 'vn' for vertex normal or 'f' for face.

#define _PUSH_BACK(type, elementCount, count, max, array)   \
    assert(tokenCount >= 1 + elementCount && "Bad model."); \
    if (count == max)                                       \
        array = realloc(array, (max *= 2) * sizeof(type));  \
    for (int _i = 0; _i < elementCount; _i++)               \
        array[count].at[_i] = atof(tokens[_i + 1]);         \
    count++;

        if (strcmp(tokens[0], "v") == 0) {
            _PUSH_BACK(Pok_Float3, 3, mesh->positionCount, positionMax, mesh->positions);
        }
        else if (strcmp(tokens[0], "vt") == 0) {
            _PUSH_BACK(Pok_Float2, 2, mesh->texCoordCount, texCoordMax, mesh->tex_coords);
        }
        else if (strcmp(tokens[0], "vn") == 0) {
            _PUSH_BACK(Pok_Float3, 3, mesh->normalCount, normalMax, mesh->normals);
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

            int faceTokenCount = tokenCount - 1 > 4 ? 4 : tokenCount - 1;
            assert(faceTokenCount == 3 || faceTokenCount == 4);

            // Split each token into indices.

            int positionIndices[4];
            int texCoordIndices[4];
            int normalIndices[4];
            
            for (int j = 0; j < faceTokenCount; j++) {
                int indexCount;
                char** indices = Pok_StringSplit(tokens[j + 1], "/", true, &indexCount);
                POK_ASSERT(indexCount == 3, "Bad model.");

                // OBJ indices start at 1 and we want them to start at 0.
                positionIndices[j] = atoi(indices[0]) - 1;
                texCoordIndices[j] = atoi(indices[1]) - 1;
                normalIndices[j] = atoi(indices[2]) - 1;

                Pok_StringSplitFree(indices);
            }

            // Single face.
            if (faceTokenCount == 3) {
                Pok_Face face;

                for (int j = 0; j < 3; j++) {
                    face.positionIndices[j] = positionIndices[j];
                    face.texCoordIndices[j] = texCoordIndices[j];
                    face.normalIndices[j] = normalIndices[j];
                }

                if (mesh->faceCount == faceMax) {
                    mesh->faces = realloc(mesh->faces, (faceMax *= 2) * sizeof(Pok_Face));
                }
                mesh->faces[mesh->faceCount++] = face;
            }
            // 2 faces.
            // A face could theoretically contain 4 vertices but we keep
            // it simple and break it up into 2 faces with 3 vertices each.
            else if (faceTokenCount == 4) {
                Pok_Face faceA;
                Pok_Face faceB;

                for (int j = 0; j < 3; j++) {
                    faceA.positionIndices[j] = positionIndices[j];
                    faceA.texCoordIndices[j] = texCoordIndices[j];
                    faceA.normalIndices[j] = normalIndices[j];

                    faceB.positionIndices[j] = positionIndices[(j + 2) % 4];
                    faceB.texCoordIndices[j] = texCoordIndices[(j + 2) % 4];
                    faceB.normalIndices[j] = normalIndices[(j + 2) % 4];
                }

                if (mesh->faceCount >= faceMax - 1) {
                    mesh->faces = realloc(mesh->faces, (faceMax *= 2) * sizeof(Pok_Face));
                }
                mesh->faces[mesh->faceCount++] = faceA;
                mesh->faces[mesh->faceCount++] = faceB;
            }
        }

        Pok_StringSplitFree(tokens);
    }

    Pok_StringSplitFree(lines);
    Pok_Free(fileContent);

    return mesh;
}

void Pok_MeshDestroy(Pok_Mesh* mesh) {
    Pok_Free(mesh->positions);
    Pok_Free(mesh->tex_coords);
    Pok_Free(mesh->normals);
    Pok_Free(mesh->faces);
    Pok_Free(mesh);
}

int Pok_ParseInt(const char* str) {
    int count;
    char** tokens = Pok_StringSplit(str, "=", false, &count);
    assert(count == 2);
    int parsed = atoi(tokens[1]);
    Pok_StringSplitFree(tokens);
    return parsed;
}

Pok_Font* Pok_FontLoad(Pok_Image* image, const char* infoPath) { 
    Pok_Font* font = Pok_AllocHeap(sizeof(Pok_Font));
    
    font->image = image;
    font->charInfoCount = 256;
    font->charInfos = Pok_AllocHeap(256 * sizeof(Pok_CharInfo));
    font->size = -1;

    for (int i = 0; i < 256; i++) {
        font->charInfos[i].id = 0;
    }

    char* content = Pok_LoadFileBinaryBlob(infoPath);

    // Split into lines.
    int lineCount;
    char** lines = Pok_StringSplit(content, "\n", false, &lineCount);

    for (int i = 0; i < lineCount; i++) {
        // Split each line into words.
        int wordCount;
        char** words = Pok_StringSplit(lines[i], " ", false, &wordCount);

        if (wordCount == 0) {
            continue;
        }

        // Search for lines that start with either 'info' or 'char'.

        if (strcmp(words[0], "info") == 0) {
            font->size = Pok_ParseInt(words[2]);
        }
        else if (strcmp(words[0], "char") == 0) {
            assert(wordCount >= 9);

            int id = Pok_ParseInt(words[1]);

            if (id < 0 || id >= 256) {
                continue;
            }

            font->charInfos[id].id = id;
            font->charInfos[id].src.x = Pok_ParseInt(words[2]);
            font->charInfos[id].src.y = Pok_ParseInt(words[3]);
            font->charInfos[id].src.w = Pok_ParseInt(words[4]);
            font->charInfos[id].src.h = Pok_ParseInt(words[5]);
            font->charInfos[id].offset.x = Pok_ParseInt(words[6]);
            font->charInfos[id].offset.y = Pok_ParseInt(words[7]);
            font->charInfos[id].xAdvance = Pok_ParseInt(words[8]);
        }

        Pok_StringSplitFree(words);
    }

    Pok_StringSplitFree(lines);
    Pok_Free(content);
 
    return font;
}

void Pok_FontDestroy(Pok_Font* font) {
    Pok_Free(font->charInfos);
    Pok_Free(font);
}
