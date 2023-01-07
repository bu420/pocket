#define SRZ_IMPLEMENTATION
#define BOOTLEG_MINECRAFT_IMPLEMENTATION
#include <bootleg_minecraft.h>
#include <stdlib.h>
#include <SDL2/SDL.h>

#define WIDTH 640
#define HEIGHT 480

int main(int argc, char** argv) {
    bootleg_minecraft_t minecraft;
    bootleg_minecraft_init(&minecraft, WIDTH, HEIGHT, 
        malloc(WIDTH * HEIGHT * sizeof(srz_byte3_t)), 
        malloc(WIDTH * HEIGHT * sizeof(float)),
        malloc(4096 * sizeof(bootleg_minecraft_block_t)));

    SDL_Init(SDL_INIT_EVERYTHING);
    SDL_Window* window = SDL_CreateWindow("Bootleg Minecraft", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, 0);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, 0);

    while (1) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            switch (e.type) {
            case SDL_QUIT:
                return 0;
            }
        }

        bootleg_minecraft_update(&minecraft, SDL_GetTicks64(), 0);

        SDL_RenderClear(renderer);

        for (int x = 0; x < WIDTH; ++x) {
            for (int y = 0; y < HEIGHT; ++y) {
                srz_byte3_t* pixel = srz_color_buffer_at(&minecraft.color_buffer, x, y);
                SDL_SetRenderDrawColor(renderer, pixel->r, pixel->g, pixel->b, 255);
                SDL_RenderDrawPoint(renderer, x, y);
            }
        }

        SDL_RenderPresent(renderer);
    }
}
