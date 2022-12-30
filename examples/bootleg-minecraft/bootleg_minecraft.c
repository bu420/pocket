#define SRZ_IMPLEMENTATION
#include <srz.h>
#include <SDL2/SDL.h>

#define WIDTH 640
#define HEIGHT 480

int main(int argc, char** argv) {
    srz_byte_t cube_vertex_positions[36 * 3] = {
        0,0,0,0,1,0,1,1,0,0,0,0,1,1,0,1,0,0,
        1,0,0,1,1,0,1,1,1,1,0,0,1,1,1,1,0,1,
        1,0,1,1,1,1,0,1,1,1,0,1,0,1,1,0,0,1,
        0,0,1,0,1,1,0,1,0,0,0,1,0,1,0,0,0,0,
        0,1,0,0,1,1,1,1,1,0,1,0,1,1,1,1,1,0,
        1,0,1,0,0,1,0,0,0,1,0,1,0,0,0,1,0,0
    };

    char cube_vertex_normals[6 * 3] = {
        0, 0, 1,
        -1, 0, 0,
        0, 0, -1,
        1, 0, 0,
        0, -1, 0,
        0, 1, 0
    };

    srz_color_buffer_t color_buffer;
    color_buffer.w = WIDTH;
    color_buffer.h = HEIGHT;
    color_buffer.data = malloc(WIDTH * HEIGHT * sizeof(srz_byte3_t));

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

        SDL_RenderClear(renderer);

        uint64_t elapsed_time = SDL_GetTicks64();
        srz_byte3_t color = {(srz_sinf(elapsed_time / 1000.f) + 1) * 128, (srz_cosf(elapsed_time / 1000.f) + 1) * 128, 128};
        srz_color_buffer_clear(&color_buffer, color);

        for (int x = 0; x < WIDTH; ++x) {
            for (int y = 0; y < HEIGHT; ++y) {
                srz_byte3_t* pixel = srz_color_buffer_at(&color_buffer, x, y);

                SDL_SetRenderDrawColor(renderer, pixel->r, pixel->g, pixel->b, 255);
                SDL_RenderDrawPoint(renderer, x, y);
            }
        }

        SDL_RenderPresent(renderer);
    }
}
