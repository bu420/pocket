#include <pok.h>
#include <stdio.h>
#include <assert.h>

#define WIDTH 640
#define HEIGHT 640

int main() {
    Pok_Image* painting = Pok_ImageLoadBMP("assets/painting.bmp", POK_R8G8B8);
    assert(painting);
    Pok_Image* fontImage = Pok_ImageLoadBMP("assets/font_dos_vga.bmp", POK_R8G8B8A8);
    assert(fontImage);
    Pok_Font* font = Pok_FontLoad(fontImage, "assets/font_dos_vga.txt");
    assert(font);

    Pok_ColorBuffer* colorBuffer = Pok_ColorBufferCreate(WIDTH, HEIGHT);
    Pok_ColorBufferClear(colorBuffer, (Pok_Byte3){255, 255, 255});

    Pok_RenderImage(colorBuffer, painting, (Pok_Rect){0, 0, WIDTH, HEIGHT}, (Pok_Rect){80, 30, WIDTH / 2, HEIGHT / 2});

    Pok_RenderText(colorBuffer, "Hello World", (Pok_Int2){20, 40}, font, 5);
    Pok_RenderText(colorBuffer, "Hello World 123 {}[]!?.:;", (Pok_Int2){200, 400}, font, 1);

    Pok_SaveBMP("output.bmp", colorBuffer);
    printf("Done.");
}
