#include <psr.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#define WIDTH 640
#define HEIGHT 640

psr_character_info_t font_info[256];

char** split(char* str, char* delims) {
    int count = 2;
    char** tokens = malloc(count * sizeof(char*));
    char* copy = strdup(str);

    char* token = strtok(copy, delims);
    int i = 0;
    while (token != NULL) {
        tokens[i] = strdup(token);
        token = strtok(NULL, delims);

        if (++i == count) {
            tokens = realloc(tokens, (count *= 2) * sizeof(char*));
        }
    }
    tokens[i] = NULL;
    return tokens;
}

int parse_int(char* str) {
    char** tokens = split(str, "=");
    return atoi(tokens[1]);
}

void load_font_info() {
    char* file_content;

    // Read file content into string.
    FILE* file = fopen("font_info.txt", "rb");
    fseek(file, 0, SEEK_END);
    int size = ftell(file);
    fseek(file, 0, SEEK_SET);
    file_content = malloc(size + 1);
    fread(file_content, 1, size, file);
    file_content[size] = '\0';
    fclose(file);

    // Split into lines.
    char** lines = split(file_content, "\n");
    for (char** line = lines; *line; line++) {
        // Split each line into words.
        char** words = split(*line, " ");

        // We only care about lines that start with "char".
        if (strcmp(words[0], "char") == 0) {
            int id = parse_int(words[1]);

            if (id >= 256) {
                continue;
            }

            font_info[id].src.x = parse_int(words[2]);
            font_info[id].src.y = parse_int(words[3]);
            font_info[id].src.w = parse_int(words[4]);
            font_info[id].src.h = parse_int(words[5]);
            font_info[id].offset.x = parse_int(words[6]);
            font_info[id].offset.y = parse_int(words[7]);
            font_info[id].x_advance = parse_int(words[8]);

            /*printf("%d: %d %d | %d %d | %d %d | %d\n", id, font_info[id].src.x, font_info[id].src.y,
                font_info[id].src.w, font_info[id].src.h, font_info[id].offset.x, font_info[id].offset.y,
                font_info[id].x_advance);*/
        }
    }
}

psr_character_info_t on_char_draw(int c, void* user_data) {
    assert(c < 256);
    return font_info[c];
}

int main() {
    load_font_info();

    psr_image_t background = psr_load_bmp("background.bmp", PSR_R8G8B8);
    psr_image_t font = psr_load_bmp("font.bmp", PSR_A8R8G8B8);

    psr_color_buffer_t color_buffer;
    psr_color_buffer_init(&color_buffer, background.w, background.h);
    psr_color_buffer_clear(&color_buffer, (psr_byte3_t){255, 255, 255});

    psr_raster_image(&color_buffer, background, (psr_rect_t){0, 0, background.w, background.h}, (psr_rect_t){80, 30, background.w / 2, background.h / 2});

    psr_raster_text(&color_buffer, "Hello World", (psr_int2_t){20, 40}, font, 72, 72, on_char_draw, NULL);
    psr_raster_text(&color_buffer, "Hello World 123 {}[]!?.:;", (psr_int2_t){200, 400}, font, 72, 24, on_char_draw, NULL);

    psr_save_bmp("output.bmp", color_buffer);
    printf("Done.");
}
