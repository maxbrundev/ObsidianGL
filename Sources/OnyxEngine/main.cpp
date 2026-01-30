#define SDL_MAIN_USE_CALLBACKS 1
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

#include <ObsidianGL/ObsidianGL.h>

SDL_Window* window     = nullptr;
SDL_Renderer* renderer = nullptr;
SDL_Texture* texture   = nullptr;

unsigned int* pixels = nullptr;
ObsidianGL::FrameBuffer framebuffer;

SDL_AppResult SDL_AppInit(void**, int, char**)
{
	SDL_SetAppMetadata("OnyxEngine", "0.1", "com.obsidian.onyx");

	if (!SDL_Init(SDL_INIT_VIDEO))
	{
		SDL_Log("SDL_Init failed: %s", SDL_GetError());
		return SDL_APP_FAILURE;
	}

	const uint16_t width = 800;
	const uint16_t height = 600;

	if (!SDL_CreateWindowAndRenderer("OnyxEngine", width, height, SDL_WINDOW_RESIZABLE, &window, &renderer))
	{
		SDL_Log("Window/Renderer creation failed: %s", SDL_GetError());

		return SDL_APP_FAILURE;
	}

	texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, width, height);

	pixels = new uint32_t[width * height];

	framebuffer.Width = width;
	framebuffer.Height = height;
	framebuffer.Pixels = pixels;

	return SDL_APP_CONTINUE;
}

SDL_AppResult SDL_AppEvent(void*, SDL_Event* event)
{
	if (event->type == SDL_EVENT_QUIT)
		return SDL_APP_SUCCESS;

	return SDL_APP_CONTINUE;
}

SDL_AppResult SDL_AppIterate(void*)
{
	ObsidianGL::Test(&framebuffer);

	SDL_UpdateTexture(texture, nullptr, pixels, framebuffer.Width * sizeof(unsigned int));

	SDL_RenderClear(renderer);
	SDL_RenderTexture(renderer, texture, nullptr, nullptr);
	SDL_RenderPresent(renderer);

	return SDL_APP_CONTINUE;
}

void SDL_AppQuit(void*, SDL_AppResult)
{
	delete[] pixels;

	SDL_DestroyTexture(texture);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);

	SDL_Quit();
}
