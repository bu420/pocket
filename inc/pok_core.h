#ifndef POK_CORE_H
#define POK_CORE_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>

#define POK_ASSERT(expr, msg) \
    assert(expr && "Pocket: FATAL ERROR: " msg)

#define POK_SWAP(type, a, b) { type _temp = a; a = b; b = _temp; }

void* Pok_AllocHeap(size_t size);
void Pok_Free(void* ptr);

#endif
