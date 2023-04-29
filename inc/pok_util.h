#ifndef POK_UTIL_H
#define POK_UTIL_H

#include <pok_core.h>

/**
 * @brief Loads an entire file into a string. Remember to free return value.
 */
char* Pok_LoadFileBinaryBlob(const char* path);

/**
 * @brief Remember to call Pok_StringSplitFree() for return value.
 */
char** Pok_StringSplit(const char* str, 
                       const char* delim, 
                       bool count_consecutive_delimiters, 
                       int* out_count);

void Pok_StringSplitFree(char** tokens);

#endif
