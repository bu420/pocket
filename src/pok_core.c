#include "pok_core.h"

#include <stdlib.h>
#include <string.h>

void* Pok_AllocHeap(size_t size) {
    void* ptr = malloc(size);
	POK_ASSERT(ptr, "Memory allocation failed.");
    return ptr;
}

void Pok_Free(void* ptr) {
    free(ptr);
}
