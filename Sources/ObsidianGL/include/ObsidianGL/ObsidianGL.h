#pragma once

#include <cstdint>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "API/Export.h"

namespace ObsidianGL
{
	constexpr uint16_t OGL_TRIANGLES = 0x0001;
	constexpr uint16_t OGL_QUADS = 0x0002;
	constexpr uint16_t OGL_TRIANGLE_STRIP = 0x0003;
	constexpr uint16_t OGL_TRIANGLE_FAN = 0x0004;

	constexpr uint16_t OGL_MODELVIEW = 0x0010;
	constexpr uint16_t OGL_PROJECTION = 0x0011;

	constexpr uint8_t OGL_DEPTH_TEST = 0x01;
	constexpr uint8_t OGL_CULL_FACE = 0x02;
	constexpr uint8_t OGL_DEPTH_WRITE = 0x04;

	constexpr uint16_t OGL_BACK = 0x0001;
	constexpr uint16_t OGL_FRONT = 0x0002;

	constexpr uint8_t OGL_COLOR_BUFFER_BIT = 0x01;
	constexpr uint8_t OGL_DEPTH_BUFFER_BIT = 0x02;

	OBSIDIANGL_API void Test();

	OBSIDIANGL_API void Initialize(uint16_t p_width, uint16_t p_height);
	OBSIDIANGL_API void Shutdown();
	OBSIDIANGL_API void SetClearColor(float p_r, float p_g, float p_b, float p_a);
	OBSIDIANGL_API void Clear(uint8_t p_flags);

	OBSIDIANGL_API void Enable(uint8_t p_state);
	OBSIDIANGL_API void Disable(uint8_t p_state);
	OBSIDIANGL_API void CullFace(uint16_t p_face);

	OBSIDIANGL_API void Begin(uint16_t p_mode);
	OBSIDIANGL_API void End();
	OBSIDIANGL_API void Vertex3f(float p_x, float p_y, float p_z);
	OBSIDIANGL_API void Color3f(float p_r, float p_g, float p_b);
	OBSIDIANGL_API void Color4f(float p_r, float p_g, float p_b, float p_a);
	OBSIDIANGL_API void Normal3f(float p_x, float p_y, float p_z);

	OBSIDIANGL_API void MatrixMode(uint16_t p_mode);
	OBSIDIANGL_API void LoadIdentity();
	OBSIDIANGL_API void PushMatrix();
	OBSIDIANGL_API void PopMatrix();
	OBSIDIANGL_API void LoadMatrixf(const float* p_matrix);
	OBSIDIANGL_API void MultMatrixf(const float* p_matrix);
	OBSIDIANGL_API void Translatef(float p_x, float p_y, float p_z);
	OBSIDIANGL_API void Rotatef(float p_angle, float p_x, float p_y, float p_z);
	OBSIDIANGL_API void Scalef(float p_x, float p_y, float p_z);
	OBSIDIANGL_API void Perspective(float p_fovY, float p_aspect, float p_near, float p_far);

	OBSIDIANGL_API uint32_t* GetFrameBufferData();
	OBSIDIANGL_API uint32_t GetFrameBufferPitch();
	OBSIDIANGL_API uint32_t GetFrameBufferWidth();
	OBSIDIANGL_API uint32_t GetFrameBufferHeight();
}