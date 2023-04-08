#include <psr/psr.h>
#include <stdio.h>
#include <assert.h>

#define WIDTH 640
#define HEIGHT 640

int main() {
    psr_image_t* background = psr_image_load_bmp("../assets/background.bmp", PSR_R8G8B8);
    assert(background);
    psr_image_t* font_image = psr_image_load_bmp("../assets/font_dos_vga.bmp", PSR_R8G8B8A8);
    assert(font_image);
    psr_font_t* font = psr_font_load(font_image, "../assets/font_dos_vga.txt");
    assert(font);

    psr_color_buffer_t* color_buffer = psr_color_buffer_create(WIDTH, HEIGHT);
    psr_color_buffer_clear(color_buffer, (psr_byte3_t){255, 255, 255});

    psr_raster_image(color_buffer, background, (psr_rect_t){0, 0, WIDTH, HEIGHT}, (psr_rect_t){80, 30, WIDTH / 2, HEIGHT / 2});

    psr_raster_text(color_buffer, "Hello World", (psr_int2_t){20, 40}, font, 5);
    psr_raster_text(color_buffer, "Hello World 123 {}[]!?.:;", (psr_int2_t){200, 400}, font, 1);

    psr_save_bmp("output.bmp", color_buffer);
    printf("Done.");
}
