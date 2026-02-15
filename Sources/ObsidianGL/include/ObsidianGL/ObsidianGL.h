#pragma once

#include <cstdint>

#include <glm/gtc/matrix_transform.hpp>

#include "API/Export.h"

namespace ObsidianGL
{
	// Quick Immediate Mode OpenGL 1.1 API style to easily bridge Minecraft classic 0.30 LWJGL calls later.
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

	constexpr uint16_t OGL_FRONT_AND_BACK = 0x0003;

	constexpr uint16_t OGL_CW = 0x0900;
	constexpr uint16_t OGL_CCW = 0x0901;

	constexpr uint16_t OGL_FILL = 0x1B00;
	constexpr uint16_t OGL_LINE = 0x1B01;
	constexpr uint16_t OGL_POINT = 0x1B02;

	constexpr uint8_t OGL_COLOR_BUFFER_BIT = 0x01;
	constexpr uint8_t OGL_DEPTH_BUFFER_BIT = 0x02;

	constexpr uint16_t OGL_TEXTURE_2D = 0x0DE1;
	constexpr uint16_t OGL_TEXTURE_MAG_FILTER = 0x2800;
	constexpr uint16_t OGL_TEXTURE_MIN_FILTER = 0x2801;
	constexpr uint16_t OGL_TEXTURE_WRAP_S = 0x2802;
	constexpr uint16_t OGL_TEXTURE_WRAP_T = 0x2803;
	constexpr uint16_t OGL_NEAREST = 0x2600;
	constexpr uint16_t OGL_LINEAR = 0x2601;
	constexpr uint16_t OGL_REPEAT = 0x2901;
	constexpr uint16_t OGL_CLAMP_TO_EDGE = 0x812F;
	constexpr uint16_t OGL_RGB = 0x1907;
	constexpr uint16_t OGL_RGBA = 0x1908;
	constexpr uint16_t OGL_UNSIGNED_BYTE = 0x1401;

	OBSIDIANGL_API void Test();

	OBSIDIANGL_API void Initialize(uint16_t p_width, uint16_t p_height);
	OBSIDIANGL_API void Shutdown();
	OBSIDIANGL_API void SetClearColor(float p_r, float p_g, float p_b, float p_a);
	OBSIDIANGL_API void Clear(uint8_t p_flags);

	OBSIDIANGL_API void Enable(uint8_t p_state);
	OBSIDIANGL_API void Disable(uint8_t p_state);
	OBSIDIANGL_API void EnableTexture2D();
	OBSIDIANGL_API void DisableTexture2D();
	OBSIDIANGL_API void CullFace(uint16_t p_face);
	OBSIDIANGL_API void FrontFace(uint16_t p_mode);
	OBSIDIANGL_API void PolygonMode(uint16_t p_face, uint16_t p_mode);
	OBSIDIANGL_API void PointSize(float p_size);
	OBSIDIANGL_API void LineWidth(float p_width);

	OBSIDIANGL_API void Begin(uint16_t p_mode);
	OBSIDIANGL_API void End();
	OBSIDIANGL_API void Vertex3f(float p_x, float p_y, float p_z);
	OBSIDIANGL_API void Color3f(float p_r, float p_g, float p_b);
	OBSIDIANGL_API void Color4f(float p_r, float p_g, float p_b, float p_a);
	OBSIDIANGL_API void Normal3f(float p_x, float p_y, float p_z);
	OBSIDIANGL_API void TexCoord2f(float p_s, float p_t);

	OBSIDIANGL_API uint32_t GenTexture();
	OBSIDIANGL_API void DeleteTexture(uint32_t p_textureId);
	OBSIDIANGL_API void BindTexture(uint32_t p_textureId);
	OBSIDIANGL_API void TexImage2D(uint16_t p_internalFormat, uint32_t p_width, uint32_t p_height, uint16_t p_format, uint16_t p_type, const void* p_data);
	OBSIDIANGL_API void TexParameteri(uint16_t p_param, uint16_t p_value);

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
	OBSIDIANGL_API void Ortho(double p_left, double p_right, double p_bottom, double p_top, double p_near, double p_far);

	OBSIDIANGL_API uint32_t* GetFrameBufferData();
	OBSIDIANGL_API uint32_t GetFrameBufferPitch();
	OBSIDIANGL_API uint32_t GetFrameBufferWidth();
	OBSIDIANGL_API uint32_t GetFrameBufferHeight();
}