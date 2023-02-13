#ifndef LLSR_IO_H
#define LLSR_IO_H

#include "llsr.h"

typedef struct {
    int position_indices[3];
    int tex_coord_indices[3];
    int normal_indices[3];
} llsr_face_t;

typedef struct {
    llsr_float3_t* positions;
    llsr_float2_t* tex_coords;
    llsr_float3_t* normals;
    llsr_face_t* faces;
    int position_count;
    int tex_coord_count;
    int normal_count;
    int face_count;
} llsr_mesh_t;

void llsr_save_bmp(char const* filename, llsr_color_buffer_t color_buffer);
// Very basic OBJ loader.
// @return 0 on failure to open file and -1 on bad model.
int llsr_load_obj(char const* filename, llsr_mesh_t* mesh);
void llsr_mesh_free(llsr_mesh_t* mesh);

#endif
