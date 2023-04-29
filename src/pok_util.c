#include "pok_util.h"

#include <stdio.h>
#include <string.h>

char* Pok_LoadFileBinaryBlob(const char* path) {
    FILE* file = fopen(path, "rb");
    if (!file) {
        return NULL;
    }
    fseek(file, 0, SEEK_END);
    int size = ftell(file);
    fseek(file, 0, SEEK_SET);
    char* content = Pok_AllocHeap(size + 1);
    fread(content, 1, size, file);
    content[size] = '\0';
    fclose(file);
    return content;
}

char** Pok_StringSplit(const char* str, const char* delim, bool count_consecutive_delimiters, int* out_count) {
    const int str_len = strlen(str);
	const int delim_len = strlen(delim);

	char** tokens = Pok_AllocHeap((str_len / delim_len + 1) * sizeof(char*));
	int count = 0;

	const char* start = str;

	while (true) {
		// Find start of next delimiter (end of current token).
		const char* end = strstr(start, delim);

		// If there are more tokens.
		if (end) {
			int len = end - start;

            // Only copy token into array if delimiter is not consecutive
            // of if it is consecutive but flag is set to OK. 
            if (len > 0 || count_consecutive_delimiters) {
                char* token = Pok_AllocHeap(len + 1);
                memcpy(token, start, len);
                token[len] = '\0';

                tokens[count++] = token;
            }
		}
		// If this is the last token.
		else if (start < str + str_len) {
			tokens[count++] = _strdup(start);
		}

		// Break if this was the last token.
		if (!end) {
			break;
		}

		// Advance starting point.
		start = end + delim_len;
	}
	tokens[count] = NULL;

	if (out_count) {
		*out_count = count;
	}

	return tokens;
}

void Pok_StringSplitFree(char** tokens) {
	for (char** ptr = tokens; *ptr; ptr++) {
		Pok_Free(*ptr);
	}
	Pok_Free(tokens);
	tokens = NULL;
}
