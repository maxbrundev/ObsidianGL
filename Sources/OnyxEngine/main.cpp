#define SDL_MAIN_USE_CALLBACKS 1
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image/stb_image.h>

#include <ObsidianGL/ObsidianGL.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>

// Dirty main for quick tests

SDL_Window* window = nullptr;
SDL_Renderer* renderer = nullptr;
SDL_Texture* texture = nullptr;

constexpr uint16_t width = 800;
constexpr uint16_t height = 600;

glm::vec3 cameraPosition = glm::vec3(0.0f, 0.0f, -5.0f);
glm::vec3 cameraEuler = glm::vec3(0.0f, 0.0f, 0.0f);
glm::quat cameraRotation = glm::quat(glm::radians(cameraEuler));

float cameraMoveSpeed = 5.0f;
float cameraSensitivity = 0.2f;
bool rightMouseHeld = false;
bool firstMouse = true;
float cubeRotation = 0.0f;

uint64_t lastTime = 0;

uint32_t cubeTextureId = 0;

const char* TEXTURE_PATH = "dirt.bmp";

uint32_t LoadTexture(const char* p_path)
{
	int w, h, channels;
	unsigned char* pixels = stbi_load(p_path, &w, &h, &channels, 0);

	if (!pixels)
	{
		SDL_Log("Failed to load texture: %s (stbi: %s)", p_path, stbi_failure_reason());
		return 0;
	}

	uint32_t texId = ObsidianGL::GenTexture();
	ObsidianGL::BindTexture(texId);

	uint16_t format = (channels == 4) ? ObsidianGL::OGL_RGBA : ObsidianGL::OGL_RGB;
	ObsidianGL::TexImage2D(format, static_cast<uint32_t>(w), static_cast<uint32_t>(h), format, ObsidianGL::OGL_UNSIGNED_BYTE, pixels);

	ObsidianGL::TexParameteri(ObsidianGL::OGL_TEXTURE_MIN_FILTER, ObsidianGL::OGL_NEAREST);
	ObsidianGL::TexParameteri(ObsidianGL::OGL_TEXTURE_MAG_FILTER, ObsidianGL::OGL_NEAREST);
	ObsidianGL::TexParameteri(ObsidianGL::OGL_TEXTURE_WRAP_S, ObsidianGL::OGL_REPEAT);
	ObsidianGL::TexParameteri(ObsidianGL::OGL_TEXTURE_WRAP_T, ObsidianGL::OGL_REPEAT);

	ObsidianGL::BindTexture(0);

	stbi_image_free(pixels);

	SDL_Log("Loaded texture: %s (%dx%d, %d channels)", p_path, w, h, channels);

	return texId;
}

void DrawCube()
{
	ObsidianGL::Begin(ObsidianGL::OGL_QUADS);

	ObsidianGL::Color3f(1.0f, 0.0f, 0.0f);
	ObsidianGL::Normal3f(0.0f, 0.0f, 1.0f);
	ObsidianGL::Vertex3f(-0.5f, -0.5f, 0.5f);
	ObsidianGL::Vertex3f(-0.5f, 0.5f, 0.5f);
	ObsidianGL::Vertex3f(0.5f, 0.5f, 0.5f);
	ObsidianGL::Vertex3f(0.5f, -0.5f, 0.5f);

	ObsidianGL::Color3f(0.0f, 1.0f, 0.0f);
	ObsidianGL::Normal3f(0.0f, 0.0f, -1.0f);
	ObsidianGL::Vertex3f(0.5f, -0.5f, -0.5f);
	ObsidianGL::Vertex3f(0.5f, 0.5f, -0.5f);
	ObsidianGL::Vertex3f(-0.5f, 0.5f, -0.5f);
	ObsidianGL::Vertex3f(-0.5f, -0.5f, -0.5f);

	ObsidianGL::Color3f(0.0f, 0.0f, 1.0f);
	ObsidianGL::Normal3f(0.0f, 1.0f, 0.0f);
	ObsidianGL::Vertex3f(-0.5f, 0.5f, 0.5f);
	ObsidianGL::Vertex3f(-0.5f, 0.5f, -0.5f);
	ObsidianGL::Vertex3f(0.5f, 0.5f, -0.5f);
	ObsidianGL::Vertex3f(0.5f, 0.5f, 0.5f);

	ObsidianGL::Color3f(1.0f, 1.0f, 0.0f);
	ObsidianGL::Normal3f(0.0f, -1.0f, 0.0f);
	ObsidianGL::Vertex3f(-0.5f, -0.5f, -0.5f);
	ObsidianGL::Vertex3f(-0.5f, -0.5f, 0.5f);
	ObsidianGL::Vertex3f(0.5f, -0.5f, 0.5f);
	ObsidianGL::Vertex3f(0.5f, -0.5f, -0.5f);

	ObsidianGL::Color3f(1.0f, 0.0f, 1.0f);
	ObsidianGL::Normal3f(1.0f, 0.0f, 0.0f);
	ObsidianGL::Vertex3f(0.5f, -0.5f, 0.5f);
	ObsidianGL::Vertex3f(0.5f, 0.5f, 0.5f);
	ObsidianGL::Vertex3f(0.5f, 0.5f, -0.5f);
	ObsidianGL::Vertex3f(0.5f, -0.5f, -0.5f);

	ObsidianGL::Color3f(0.0f, 1.0f, 1.0f);
	ObsidianGL::Normal3f(-1.0f, 0.0f, 0.0f);
	ObsidianGL::Vertex3f(-0.5f, -0.5f, -0.5f);
	ObsidianGL::Vertex3f(-0.5f, 0.5f, -0.5f);
	ObsidianGL::Vertex3f(-0.5f, 0.5f, 0.5f);
	ObsidianGL::Vertex3f(-0.5f, -0.5f, 0.5f);

	ObsidianGL::End();
}

void DrawTexturedCube()
{
	ObsidianGL::Begin(ObsidianGL::OGL_QUADS);

	ObsidianGL::Color3f(1.0f, 1.0f, 1.0f);

	ObsidianGL::Normal3f(0.0f, 0.0f, 1.0f);
	ObsidianGL::TexCoord2f(0.0f, 1.0f); ObsidianGL::Vertex3f(-0.5f, -0.5f, 0.5f);
	ObsidianGL::TexCoord2f(0.0f, 0.0f); ObsidianGL::Vertex3f(-0.5f, 0.5f, 0.5f);
	ObsidianGL::TexCoord2f(1.0f, 0.0f); ObsidianGL::Vertex3f(0.5f, 0.5f, 0.5f);
	ObsidianGL::TexCoord2f(1.0f, 1.0f); ObsidianGL::Vertex3f(0.5f, -0.5f, 0.5f);

	ObsidianGL::Normal3f(0.0f, 0.0f, -1.0f);
	ObsidianGL::TexCoord2f(0.0f, 1.0f); ObsidianGL::Vertex3f(0.5f, -0.5f, -0.5f);
	ObsidianGL::TexCoord2f(0.0f, 0.0f); ObsidianGL::Vertex3f(0.5f, 0.5f, -0.5f);
	ObsidianGL::TexCoord2f(1.0f, 0.0f); ObsidianGL::Vertex3f(-0.5f, 0.5f, -0.5f);
	ObsidianGL::TexCoord2f(1.0f, 1.0f); ObsidianGL::Vertex3f(-0.5f, -0.5f, -0.5f);

	ObsidianGL::Normal3f(0.0f, 1.0f, 0.0f);
	ObsidianGL::TexCoord2f(0.0f, 1.0f); ObsidianGL::Vertex3f(-0.5f, 0.5f, 0.5f);
	ObsidianGL::TexCoord2f(0.0f, 0.0f); ObsidianGL::Vertex3f(-0.5f, 0.5f, -0.5f);
	ObsidianGL::TexCoord2f(1.0f, 0.0f); ObsidianGL::Vertex3f(0.5f, 0.5f, -0.5f);
	ObsidianGL::TexCoord2f(1.0f, 1.0f); ObsidianGL::Vertex3f(0.5f, 0.5f, 0.5f);

	ObsidianGL::Normal3f(0.0f, -1.0f, 0.0f);
	ObsidianGL::TexCoord2f(0.0f, 1.0f); ObsidianGL::Vertex3f(-0.5f, -0.5f, -0.5f);
	ObsidianGL::TexCoord2f(0.0f, 0.0f); ObsidianGL::Vertex3f(-0.5f, -0.5f, 0.5f);
	ObsidianGL::TexCoord2f(1.0f, 0.0f); ObsidianGL::Vertex3f(0.5f, -0.5f, 0.5f);
	ObsidianGL::TexCoord2f(1.0f, 1.0f); ObsidianGL::Vertex3f(0.5f, -0.5f, -0.5f);

	ObsidianGL::Normal3f(1.0f, 0.0f, 0.0f);
	ObsidianGL::TexCoord2f(0.0f, 1.0f); ObsidianGL::Vertex3f(0.5f, -0.5f, 0.5f);
	ObsidianGL::TexCoord2f(0.0f, 0.0f); ObsidianGL::Vertex3f(0.5f, 0.5f, 0.5f);
	ObsidianGL::TexCoord2f(1.0f, 0.0f); ObsidianGL::Vertex3f(0.5f, 0.5f, -0.5f);
	ObsidianGL::TexCoord2f(1.0f, 1.0f); ObsidianGL::Vertex3f(0.5f, -0.5f, -0.5f);

	ObsidianGL::Normal3f(-1.0f, 0.0f, 0.0f);
	ObsidianGL::TexCoord2f(0.0f, 1.0f); ObsidianGL::Vertex3f(-0.5f, -0.5f, -0.5f);
	ObsidianGL::TexCoord2f(0.0f, 0.0f); ObsidianGL::Vertex3f(-0.5f, 0.5f, -0.5f);
	ObsidianGL::TexCoord2f(1.0f, 0.0f); ObsidianGL::Vertex3f(-0.5f, 0.5f, 0.5f);
	ObsidianGL::TexCoord2f(1.0f, 1.0f); ObsidianGL::Vertex3f(-0.5f, -0.5f, 0.5f);

	ObsidianGL::End();
}

SDL_AppResult SDL_AppInit(void**, int, char**)
{
	SDL_SetAppMetadata("OnyxEngine", "0.1", "com.obsidian.onyx");

	if (!SDL_Init(SDL_INIT_VIDEO))
	{
		SDL_Log("SDL_Init failed: %s", SDL_GetError());
		return SDL_APP_FAILURE;
	}

	if (!SDL_CreateWindowAndRenderer("OnyxEngine", width, height, SDL_WINDOW_RESIZABLE, &window, &renderer))
	{
		SDL_Log("Window/Renderer creation failed: %s", SDL_GetError());
		return SDL_APP_FAILURE;
	}

	texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, width, height);

	ObsidianGL::Initialize(width, height);
	ObsidianGL::SetClearColor(0.1f, 0.1f, 0.1f, 1.0f);
	ObsidianGL::Enable(ObsidianGL::OGL_DEPTH_TEST);
	ObsidianGL::Enable(ObsidianGL::OGL_CULL_FACE);
	ObsidianGL::CullFace(ObsidianGL::OGL_BACK);

	cubeTextureId = LoadTexture(TEXTURE_PATH);

	lastTime = SDL_GetPerformanceCounter();

	return SDL_APP_CONTINUE;
}

SDL_AppResult SDL_AppEvent(void*, SDL_Event* event)
{
	if (event->type == SDL_EVENT_QUIT)
		return SDL_APP_SUCCESS;

	if (event->type == SDL_EVENT_MOUSE_BUTTON_DOWN && event->button.button == SDL_BUTTON_RIGHT)
	{
		rightMouseHeld = true;
		firstMouse = true;
		SDL_SetWindowRelativeMouseMode(window, true);
	}

	if (event->type == SDL_EVENT_MOUSE_BUTTON_UP && event->button.button == SDL_BUTTON_RIGHT)
	{
		rightMouseHeld = false;
		SDL_SetWindowRelativeMouseMode(window, false);
	}

	if (event->type == SDL_EVENT_MOUSE_MOTION && rightMouseHeld)
	{
		float xOffset = static_cast<float>(event->motion.xrel) * cameraSensitivity;
		float yOffset = static_cast<float>(event->motion.yrel) * cameraSensitivity;

		cameraEuler.y -= xOffset;
		cameraEuler.x += yOffset;

		cameraEuler.x = glm::clamp(cameraEuler.x, -89.0f, 89.0f);

		cameraRotation = glm::quat(glm::radians(cameraEuler));
	}

	if (event->type == SDL_EVENT_MOUSE_WHEEL)
	{
		cameraPosition += cameraRotation * glm::vec3(0.0f, 0.0f, 1.0f) * event->wheel.y;
	}

	return SDL_APP_CONTINUE;
}

SDL_AppResult SDL_AppIterate(void*)
{
	uint64_t currentTime = SDL_GetPerformanceCounter();
	float deltaTime = static_cast<float>(currentTime - lastTime) / static_cast<float>(SDL_GetPerformanceFrequency());
	lastTime = currentTime;

	const bool* keys = SDL_GetKeyboardState(nullptr);

	if (rightMouseHeld)
	{
		float velocity = cameraMoveSpeed * deltaTime;

		if (keys[SDL_SCANCODE_LSHIFT])
		{
			velocity *= 2.0f;
		}

		glm::vec3 forward = cameraRotation * glm::vec3(0.0f, 0.0f, 1.0f);
		glm::vec3 right = cameraRotation * glm::vec3(1.0f, 0.0f, 0.0f);
		glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);

		if (keys[SDL_SCANCODE_W]) cameraPosition += forward * velocity;
		if (keys[SDL_SCANCODE_S]) cameraPosition -= forward * velocity;
		if (keys[SDL_SCANCODE_A]) cameraPosition += right * velocity;
		if (keys[SDL_SCANCODE_D]) cameraPosition -= right * velocity;
		if (keys[SDL_SCANCODE_E]) cameraPosition += up * velocity;
		if (keys[SDL_SCANCODE_Q]) cameraPosition -= up * velocity;
	}

	cubeRotation += deltaTime * 45.0f;

	ObsidianGL::Clear(ObsidianGL::OGL_COLOR_BUFFER_BIT | ObsidianGL::OGL_DEPTH_BUFFER_BIT);

	ObsidianGL::MatrixMode(ObsidianGL::OGL_PROJECTION);
	ObsidianGL::LoadIdentity();
	ObsidianGL::Perspective(60.0f, static_cast<float>(width) / static_cast<float>(height), 0.1f, 1000.0f);

	ObsidianGL::MatrixMode(ObsidianGL::OGL_MODELVIEW);
	ObsidianGL::LoadIdentity();

	glm::vec3 forward = cameraRotation * glm::vec3(0.0f, 0.0f, 1.0f);
	glm::vec3 up = cameraRotation * glm::vec3(0.0f, 1.0f, 0.0f);
	glm::mat4 viewMatrix = glm::lookAt(cameraPosition, cameraPosition + forward, up);
	ObsidianGL::LoadMatrixf(&viewMatrix[0][0]);

	ObsidianGL::PushMatrix();
	ObsidianGL::Translatef(-2.0f, 0.0f, 0.0f);
	ObsidianGL::Rotatef(cubeRotation, 0.0f, 1.0f, 0.0f);
	ObsidianGL::Rotatef(cubeRotation * 0.7f, 1.0f, 0.0f, 0.0f);
	DrawCube();
	ObsidianGL::PopMatrix();

	ObsidianGL::EnableTexture2D();
	ObsidianGL::BindTexture(cubeTextureId);

	ObsidianGL::PushMatrix();
	ObsidianGL::Translatef(2.0f, 0.0f, 0.0f);
	DrawTexturedCube();
	ObsidianGL::PopMatrix();

	ObsidianGL::PolygonMode(ObsidianGL::OGL_FRONT, ObsidianGL::OGL_LINE);
	ObsidianGL::PushMatrix();
	ObsidianGL::Translatef(4.0f, 0.0f, 0.0f);
	DrawTexturedCube();
	ObsidianGL::PopMatrix();
	ObsidianGL::PolygonMode(ObsidianGL::OGL_FRONT, ObsidianGL::OGL_FILL);

	ObsidianGL::DisableTexture2D();
	ObsidianGL::BindTexture(0);

	ObsidianGL::Begin(ObsidianGL::OGL_TRIANGLES);
	ObsidianGL::Color3f(1.0f, 0.0f, 0.0f);
	ObsidianGL::Vertex3f(-1.0f, -1.0f, 0.0f);
	ObsidianGL::Color3f(0.0f, 1.0f, 0.0f);
	ObsidianGL::Vertex3f(1.0f, -1.0f, 0.0f);
	ObsidianGL::Color3f(0.0f, 0.0f, 1.0f);
	ObsidianGL::Vertex3f(0.0f, 1.0f, 0.0f);
	ObsidianGL::End();

	//ObsidianGL::Test();

	SDL_UpdateTexture(texture, nullptr, ObsidianGL::GetFrameBufferData(), width * sizeof(uint32_t));
	SDL_RenderClear(renderer);
	SDL_RenderTexture(renderer, texture, nullptr, nullptr);
	SDL_RenderPresent(renderer);

	return SDL_APP_CONTINUE;
}

void SDL_AppQuit(void*, SDL_AppResult)
{
	ObsidianGL::Shutdown();

	SDL_DestroyTexture(texture);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);

	SDL_Quit();
}