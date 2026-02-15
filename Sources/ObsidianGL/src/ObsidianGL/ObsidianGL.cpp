#include "ObsidianGL/ObsidianGL.h"

#include <algorithm>
#include <vector>

#include <glm/gtc/type_ptr.hpp>

#include "ObsidianGL/API/SIMDConfig.h"

#include "ObsidianGL/Buffers/FrameBuffer.h"

#include "ObsidianGL/Geometry/Clipper.h"
#include "ObsidianGL/Geometry/Triangle.h"

#include "ObsidianGL/Graphics/Color.h"
#include "ObsidianGL/Graphics/Texture.h"

#include "ObsidianGL/Managers/TextureManager.h"

// The SIMD parts of the pipeline is so complex that extending them with new features and avoiding avx2 cache pressure is overly challenging
// It's a nice exercise, but it requires heavily self-documented code to continue, I wil address that asap.

// TODO:
// - Reduce AVX2 cache pressure
// - Backbuffer / Frontbuffer
// - MSAA
// - Program -> Better architecture than AmberGL, need to store attributes, interpolate varyings...
// - Tile-Based
// - Stencil Buffer
// - Modern API

namespace
{
	constexpr uint8_t MATRIX_STACK_MAX_DEPTH = 32;

	struct MatrixStack
	{
		alignas(32) glm::mat4 Data[MATRIX_STACK_MAX_DEPTH];

		uint8_t Index = 0;

		MatrixStack()
		{
			Data[0] = glm::mat4(1.0f);
		}

		glm::mat4& Top()
		{
			return Data[Index];
		}

		const glm::mat4& Top() const
		{
			return Data[Index];
		}

		void Push()
		{
			if (Index < MATRIX_STACK_MAX_DEPTH - 1)
			{
				Data[Index + 1] = Data[Index];

				++Index;
			}
		}

		void Pop()
		{
			if (Index > 0)
			{
				--Index;
			}
		}
	};

	struct RenderState
	{
		ObsidianGL::Buffers::ColorBuffer* ColorBuffer = nullptr;
		ObsidianGL::Buffers::DepthBuffer* DepthBuffer = nullptr;

		uint8_t State = ObsidianGL::OGL_DEPTH_TEST | ObsidianGL::OGL_DEPTH_WRITE;
		uint16_t CullFace = ObsidianGL::OGL_BACK;
		uint16_t FrontFace = ObsidianGL::OGL_CCW;
		uint16_t PolygonModeFront = ObsidianGL::OGL_FILL;
		uint16_t PolygonModeBack = ObsidianGL::OGL_FILL;
		float PointSize = 1.0f;
		float LineWidth = 1.0f;

		MatrixStack ModelViewStack;
		MatrixStack ProjectionStack;
		uint16_t CurrentMatrixMode = ObsidianGL::OGL_MODELVIEW;

		bool InBeginEnd = false;
		uint16_t CurrentPrimitiveMode = 0;
		glm::vec4 CurrentColor = glm::vec4(1.0f);
		glm::vec3 CurrentNormal = glm::vec3(0.0f, 0.0f, 1.0f);
		glm::vec2 CurrentTexCoord = glm::vec2(0.0f);

		bool TextureEnabled = false;
		uint32_t BoundTextureId = 0;
	};

	struct ImmediateVertex
	{
		glm::vec3 Position;
		glm::vec4 Color;
		glm::vec3 Normal;
		glm::vec2 TexCoord;
	};

	RenderState CurrentRenderState;
	std::vector<ImmediateVertex> VertexBuffer;
	ObsidianGL::Managers::TextureManager TextureManager;

	MatrixStack& GetCurrentStack()
	{
		return (CurrentRenderState.CurrentMatrixMode == ObsidianGL::OGL_PROJECTION) ? CurrentRenderState.ProjectionStack : CurrentRenderState.ModelViewStack;
	}

	const ObsidianGL::Graphics::Texture* GetActiveTexture()
	{
		if (!CurrentRenderState.TextureEnabled || CurrentRenderState.BoundTextureId == 0)
			return nullptr;

		const ObsidianGL::Graphics::Texture* tex = TextureManager.Get(CurrentRenderState.BoundTextureId);
		return (tex && tex->Data) ? tex : nullptr;
	}

	void RasterizeTriangle3D_Scalar(const glm::vec4& p_vertex0, const glm::vec4& p_vertex1, const glm::vec4& p_vertex2, const glm::vec4& p_color0, const glm::vec4& p_color1, const glm::vec4& p_color2, const glm::vec2& p_texCoord0, const glm::vec2& p_texCoord1, const glm::vec2& p_texCoord2)
	{
		if (!CurrentRenderState.ColorBuffer)
			return;

		const uint32_t width = CurrentRenderState.ColorBuffer->Width;
		const uint32_t height = CurrentRenderState.ColorBuffer->Height;

		const float invW0 = 1.0f / p_vertex0.w;
		float invW1 = 1.0f / p_vertex1.w;
		float invW2 = 1.0f / p_vertex2.w;

		const glm::vec3 ndc0 = glm::vec3(p_vertex0) * invW0;
		glm::vec3 ndc1 = glm::vec3(p_vertex1) * invW1;
		glm::vec3 ndc2 = glm::vec3(p_vertex2) * invW2;

		const float halfW = static_cast<float>(width) * 0.5f;
		const float halfH = static_cast<float>(height) * 0.5f;

		const glm::vec2 screen0((ndc0.x + 1.0f) * halfW, (1.0f - ndc0.y) * halfH);
		glm::vec2 screen1((ndc1.x + 1.0f) * halfW, (1.0f - ndc1.y) * halfH);
		glm::vec2 screen2((ndc2.x + 1.0f) * halfW, (1.0f - ndc2.y) * halfH);

		const float area = (screen1.x - screen0.x) * (screen2.y - screen0.y) - (screen1.y - screen0.y) * (screen2.x - screen0.x);

		if (std::abs(area) < 0.5f)
			return;

		if (CurrentRenderState.State & ObsidianGL::OGL_CULL_FACE)
		{
			const bool isFrontFacing = (CurrentRenderState.FrontFace == ObsidianGL::OGL_CCW) ? (area > 0.0f) : (area < 0.0f);

			if ((CurrentRenderState.CullFace == ObsidianGL::OGL_BACK && !isFrontFacing) || (CurrentRenderState.CullFace == ObsidianGL::OGL_FRONT && isFrontFacing))
				return;
		}

		const bool swapWinding = (area < 0.0f);

		if (swapWinding)
		{
			std::swap(screen1, screen2);
			std::swap(ndc1, ndc2);
			std::swap(invW1, invW2);
		}

		const auto& color1 = swapWinding ? p_color2 : p_color1;
		const auto& color2 = swapWinding ? p_color1 : p_color2;
		const auto& texCoord1 = swapWinding ? p_texCoord2 : p_texCoord1;
		const auto& texCoord2 = swapWinding ? p_texCoord1 : p_texCoord2;

		ObsidianGL::Geometry::Triangle triangle(screen0, screen1, screen2);

		if (triangle.InvFixedDoubledArea == 0.0f)
			return;

		const int32_t minX = std::max(triangle.BoundingBox.MinX, 0);
		const int32_t minY = std::max(triangle.BoundingBox.MinY, 0);
		const int32_t maxX = std::min(triangle.BoundingBox.MaxX, static_cast<int32_t>(width - 1));
		const int32_t maxY = std::min(triangle.BoundingBox.MaxY, static_cast<int32_t>(height - 1));

		if (maxX < minX || maxY < minY)
			return;

		const float depth0 = ndc0.z * 0.5f + 0.5f;
		const float depth1 = ndc1.z * 0.5f + 0.5f;
		const float depth2 = ndc2.z * 0.5f + 0.5f;

		const glm::vec4 c0w = p_color0 * invW0;
		const glm::vec4 c1w = color1 * invW1;
		const glm::vec4 c2w = color2 * invW2;

		const glm::vec2 t0w = p_texCoord0 * invW0;
		const glm::vec2 t1w = texCoord1 * invW1;
		const glm::vec2 t2w = texCoord2 * invW2;

		const bool depthTest = (CurrentRenderState.State & ObsidianGL::OGL_DEPTH_TEST) && CurrentRenderState.DepthBuffer;
		const bool depthWrite = (CurrentRenderState.State & ObsidianGL::OGL_DEPTH_WRITE) && CurrentRenderState.DepthBuffer;

		const ObsidianGL::Graphics::Texture* texture = GetActiveTexture();

		for (int32_t y = minY; y <= maxY; ++y)
		{
			RGBA8* colorRow = CurrentRenderState.ColorBuffer->GetRow(static_cast<uint32_t>(y));
			Depth* depthRow = (depthTest || depthWrite) ? CurrentRenderState.DepthBuffer->GetRow(static_cast<uint32_t>(y)) : nullptr;

			for (int32_t x = minX; x <= maxX; ++x)
			{
				float w;
				float u;
				float v;

				if (!triangle.ComputeBarycentricCoordinates(x, y, w, u, v))
					continue;

				const float oneOverW = w * invW0 + u * invW1 + v * invW2;
				const float correctedW = 1.0f / oneOverW;

				const float depth = w * depth0 + u * depth1 + v * depth2;

				if (depth < 0.0f || depth > 1.0f)
					continue;

				if (depthTest && depthRow)
				{
					if (depth >= depthRow[x])
						continue;
				}

				const glm::vec4 colorInterpolated = w * c0w + u * c1w + v * c2w;
				glm::vec4 color = colorInterpolated * correctedW;

				if (texture)
				{
					const glm::vec2 texInterpolated = w * t0w + u * t1w + v * t2w;
					const glm::vec2 texCoord = texInterpolated * correctedW;

					uint32_t texel = texture->Sample(texCoord.x, texCoord.y);

					float r;
					float g;
					float b;
					float a;

					ObsidianGL::Graphics::UnpackColorToFloat(texel, r, g, b, a);

					color.r *= r;
					color.g *= g;
					color.b *= b;
					color.a *= a;
				}

				colorRow[x] = ObsidianGL::Graphics::PackColor(color);

				if (depthWrite && depthRow)
				{
					depthRow[x] = depth;
				}
			}
		}
	}

#ifdef OBSIDIAN_USE_AVX2
	void RasterizeTriangle3D_AVX2(const glm::vec4& p_vertex0, const glm::vec4& p_vertex1, const glm::vec4& p_vertex2, const glm::vec4& p_color0, const glm::vec4& p_color1, const glm::vec4& p_color2, const glm::vec2& p_texCoord0, const glm::vec2& p_texCoord1, const glm::vec2& p_texCoord2)
	{
		if (!CurrentRenderState.ColorBuffer)
			return;

		const uint32_t width = CurrentRenderState.ColorBuffer->Width;
		const uint32_t height = CurrentRenderState.ColorBuffer->Height;

		const float invW0 = 1.0f / p_vertex0.w;
		float invW1 = 1.0f / p_vertex1.w;
		float invW2 = 1.0f / p_vertex2.w;

		const glm::vec3 ndc0 = glm::vec3(p_vertex0) * invW0;
		glm::vec3 ndc1 = glm::vec3(p_vertex1) * invW1;
		glm::vec3 ndc2 = glm::vec3(p_vertex2) * invW2;

		const float halfW = static_cast<float>(width) * 0.5f;
		const float halfH = static_cast<float>(height) * 0.5f;

		const glm::vec2 screen0((ndc0.x + 1.0f) * halfW, (1.0f - ndc0.y) * halfH);
		glm::vec2 screen1((ndc1.x + 1.0f) * halfW, (1.0f - ndc1.y) * halfH);
		glm::vec2 screen2((ndc2.x + 1.0f) * halfW, (1.0f - ndc2.y) * halfH);

		const float area = (screen1.x - screen0.x) * (screen2.y - screen0.y) - (screen1.y - screen0.y) * (screen2.x - screen0.x);

		if (std::abs(area) < 0.5f)
			return;

		if (CurrentRenderState.State & ObsidianGL::OGL_CULL_FACE)
		{
			const bool isFrontFacing = (CurrentRenderState.FrontFace == ObsidianGL::OGL_CCW) ? (area > 0.0f) : (area < 0.0f);

			if ((CurrentRenderState.CullFace == ObsidianGL::OGL_BACK && !isFrontFacing) || (CurrentRenderState.CullFace == ObsidianGL::OGL_FRONT && isFrontFacing))
				return;
		}

		const bool swapWinding = (area < 0.0f);

		if (swapWinding)
		{
			std::swap(screen1, screen2);
			std::swap(ndc1, ndc2);
			std::swap(invW1, invW2);
		}

		const auto& color1 = swapWinding ? p_color2 : p_color1;
		const auto& color2 = swapWinding ? p_color1 : p_color2;

		const auto& texCoord1 = swapWinding ? p_texCoord2 : p_texCoord1;
		const auto& texCoord2 = swapWinding ? p_texCoord1 : p_texCoord2;

		ObsidianGL::Geometry::Triangle triangle(screen0, screen1, screen2);

		if (triangle.InvFixedDoubledArea == 0.0f)
			return;

		const int32_t minX = std::max(triangle.BoundingBox.MinX, 0);
		const int32_t minY = std::max(triangle.BoundingBox.MinY, 0);
		const int32_t maxX = std::min(triangle.BoundingBox.MaxX, static_cast<int32_t>(width - 1));
		const int32_t maxY = std::min(triangle.BoundingBox.MaxY, static_cast<int32_t>(height - 1));

		if (maxX < minX || maxY < minY)
			return;

		const float depth0 = ndc0.z * 0.5f + 0.5f;
		const float depth1 = ndc1.z * 0.5f + 0.5f;
		const float depth2 = ndc2.z * 0.5f + 0.5f;

		const glm::vec4 c0w = p_color0 * invW0;
		const glm::vec4 c1w = color1 * invW1;
		const glm::vec4 c2w = color2 * invW2;

		const glm::vec2 t0w = p_texCoord0 * invW0;
		const glm::vec2 t1w = texCoord1 * invW1;
		const glm::vec2 t2w = texCoord2 * invW2;

		const bool depthTest = (CurrentRenderState.State & ObsidianGL::OGL_DEPTH_TEST) && CurrentRenderState.DepthBuffer;
		const bool depthWrite = (CurrentRenderState.State & ObsidianGL::OGL_DEPTH_WRITE) && CurrentRenderState.DepthBuffer;

		const ObsidianGL::Graphics::Texture* texture = GetActiveTexture();

		const __m256 vInvW0 = _mm256_set1_ps(invW0);
		const __m256 vInvW1 = _mm256_set1_ps(invW1);
		const __m256 vInvW2 = _mm256_set1_ps(invW2);
		const __m256 vZero = _mm256_setzero_ps();
		const __m256 vOne = _mm256_set1_ps(1.0f);

		const auto baryConst = triangle.PrepareSIMDBarycentricSetup8();

		for (int32_t y = minY; y <= maxY; ++y)
		{
			RGBA8* colorRow = CurrentRenderState.ColorBuffer->GetRow(static_cast<uint32_t>(y));
			Depth* depthRow = (depthTest || depthWrite) ? CurrentRenderState.DepthBuffer->GetRow(static_cast<uint32_t>(y)) : nullptr;

			int32_t x = minX;

			for (; x + 7 <= maxX; x += 8)
			{
				__m256 vBaryU;
				__m256 vBaryV;
				__m256 vBaryW;
				__m256 insideMask;

				triangle.ComputeBarycentricCoordinates8(baryConst, x, y, vBaryW, vBaryU, vBaryV, insideMask);

				int insideBits = _mm256_movemask_ps(insideMask);

				if (insideBits == 0)
					continue;

				__m256 vOneOverW = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(vBaryW, vInvW0), _mm256_mul_ps(vBaryU, vInvW1)), _mm256_mul_ps(vBaryV, vInvW2));
				__m256 vCorrectedW = _mm256_div_ps(vOne, vOneOverW);

				__m256 vDepth = _mm256_add_ps(_mm256_add_ps(
					_mm256_mul_ps(vBaryW, _mm256_set1_ps(depth0)),
					_mm256_mul_ps(vBaryU, _mm256_set1_ps(depth1))),
					_mm256_mul_ps(vBaryV, _mm256_set1_ps(depth2)));

				__m256 depthValidMask = _mm256_and_ps(_mm256_cmp_ps(vDepth, vZero, _CMP_GE_OQ), _mm256_cmp_ps(vDepth, vOne, _CMP_LE_OQ));

				insideMask = _mm256_and_ps(insideMask, depthValidMask);

				if (depthTest && depthRow)
				{
					__m256 vBufferDepth = _mm256_loadu_ps(&depthRow[x]);

					__m256 depthPassMask = _mm256_cmp_ps(vDepth, vBufferDepth, _CMP_LT_OQ);

					insideMask = _mm256_and_ps(insideMask, depthPassMask);
				}

				insideBits = _mm256_movemask_ps(insideMask);

				if (insideBits == 0)
					continue;

				__m256i vTexels = _mm256_setzero_si256();

				if (texture)
				{
					__m256 vTexU = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(
						_mm256_mul_ps(vBaryW, _mm256_set1_ps(t0w.x)),
						_mm256_mul_ps(vBaryU, _mm256_set1_ps(t1w.x))),
						_mm256_mul_ps(vBaryV, _mm256_set1_ps(t2w.x))), vCorrectedW);

					__m256 vTexV = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(
						_mm256_mul_ps(vBaryW, _mm256_set1_ps(t0w.y)),
						_mm256_mul_ps(vBaryU, _mm256_set1_ps(t1w.y))),
						_mm256_mul_ps(vBaryV, _mm256_set1_ps(t2w.y))), vCorrectedW);

					alignas(32) float texU[8];
					alignas(32) float texV[8];

					_mm256_store_ps(texU, vTexU);
					_mm256_store_ps(texV, vTexV);

					alignas(32) uint32_t texels[8];

					texture->SampleNearest8_AVX2(texU, texV, texels);
					vTexels = _mm256_load_si256(reinterpret_cast<const __m256i*>(texels));
				}

				__m256 vChan = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(
					_mm256_mul_ps(vBaryW, _mm256_set1_ps(c0w.r)),
					_mm256_mul_ps(vBaryU, _mm256_set1_ps(c1w.r))),
					_mm256_mul_ps(vBaryV, _mm256_set1_ps(c2w.r))), vCorrectedW);

				if (texture)
					vChan = _mm256_mul_ps(vChan, _mm256_mul_ps(_mm256_cvtepi32_ps(_mm256_and_si256(_mm256_srli_epi32(vTexels, 24), _mm256_set1_epi32(0xFF))), _mm256_set1_ps(1.0f / 255.0f)));

				__m256i vColor = _mm256_slli_epi32(_mm256_cvtps_epi32(
					_mm256_mul_ps(_mm256_min_ps(_mm256_max_ps(vChan, vZero), vOne), _mm256_set1_ps(255.0f))), 24);

				vChan = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(
					_mm256_mul_ps(vBaryW, _mm256_set1_ps(c0w.g)),
					_mm256_mul_ps(vBaryU, _mm256_set1_ps(c1w.g))),
					_mm256_mul_ps(vBaryV, _mm256_set1_ps(c2w.g))), vCorrectedW);

				if (texture)
					vChan = _mm256_mul_ps(vChan, _mm256_mul_ps(_mm256_cvtepi32_ps(_mm256_and_si256(_mm256_srli_epi32(vTexels, 16), _mm256_set1_epi32(0xFF))), _mm256_set1_ps(1.0f / 255.0f)));

				vColor = _mm256_or_si256(vColor, _mm256_slli_epi32(_mm256_cvtps_epi32(
					_mm256_mul_ps(_mm256_min_ps(_mm256_max_ps(vChan, vZero), vOne), _mm256_set1_ps(255.0f))), 16));

				vChan = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(
					_mm256_mul_ps(vBaryW, _mm256_set1_ps(c0w.b)),
					_mm256_mul_ps(vBaryU, _mm256_set1_ps(c1w.b))),
					_mm256_mul_ps(vBaryV, _mm256_set1_ps(c2w.b))), vCorrectedW);

				if (texture)
					vChan = _mm256_mul_ps(vChan, _mm256_mul_ps(_mm256_cvtepi32_ps(_mm256_and_si256(_mm256_srli_epi32(vTexels, 8), _mm256_set1_epi32(0xFF))), _mm256_set1_ps(1.0f / 255.0f)));

				vColor = _mm256_or_si256(vColor, _mm256_slli_epi32(_mm256_cvtps_epi32(
					_mm256_mul_ps(_mm256_min_ps(_mm256_max_ps(vChan, vZero), vOne), _mm256_set1_ps(255.0f))), 8));

				vChan = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(
					_mm256_mul_ps(vBaryW, _mm256_set1_ps(c0w.a)),
					_mm256_mul_ps(vBaryU, _mm256_set1_ps(c1w.a))),
					_mm256_mul_ps(vBaryV, _mm256_set1_ps(c2w.a))), vCorrectedW);

				if (texture)
					vChan = _mm256_mul_ps(vChan, _mm256_mul_ps(_mm256_cvtepi32_ps(_mm256_and_si256(vTexels, _mm256_set1_epi32(0xFF))), _mm256_set1_ps(1.0f / 255.0f)));

				vColor = _mm256_or_si256(vColor, _mm256_cvtps_epi32(
					_mm256_mul_ps(_mm256_min_ps(_mm256_max_ps(vChan, vZero), vOne), _mm256_set1_ps(255.0f))));


				if (insideBits == 0xFF)
				{
					_mm256_storeu_si256(reinterpret_cast<__m256i*>(&colorRow[x]), vColor);

					if (depthWrite && depthRow)
					{
						_mm256_storeu_ps(&depthRow[x], vDepth);
					}
				}
				else
				{
					__m256 oldColors = _mm256_castsi256_ps(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(&colorRow[x])));
					__m256 newColors = _mm256_castsi256_ps(vColor);
					__m256 blended = _mm256_blendv_ps(oldColors, newColors, insideMask);

					_mm256_storeu_si256(reinterpret_cast<__m256i*>(&colorRow[x]), _mm256_castps_si256(blended));

					if (depthWrite && depthRow)
					{
						__m256 oldDepth = _mm256_loadu_ps(&depthRow[x]);
						__m256 newDepth = _mm256_blendv_ps(oldDepth, vDepth, insideMask);
						_mm256_storeu_ps(&depthRow[x], newDepth);
					}
				}
			}

			// Scalar fallback for remaining pixels
			for (; x <= maxX; ++x)
			{
				float w;
				float u;
				float v;

				if (!triangle.ComputeBarycentricCoordinates(x, y, w, u, v))
					continue;

				const float oneOverW = w * invW0 + u * invW1 + v * invW2;
				const float correctedW = 1.0f / oneOverW;
				const float depth = w * depth0 + u * depth1 + v * depth2;

				if (depth < 0.0f || depth > 1.0f)
					continue;

				if (depthTest && depthRow && depth >= depthRow[x])
					continue;

				const glm::vec4 colorInterp = w * c0w + u * c1w + v * c2w;
				glm::vec4 color = colorInterp * correctedW;

				if (texture)
				{
					const glm::vec2 texInterp = w * t0w + u * t1w + v * t2w;
					const glm::vec2 texCoord = texInterp * correctedW;

					uint32_t texel = texture->Sample(texCoord.x, texCoord.y);

					float r;
					float g;
					float b;
					float a;

					ObsidianGL::Graphics::UnpackColorToFloat(texel, r, g, b, a);

					color.r *= r;
					color.g *= g;
					color.b *= b;
					color.a *= a;
				}

				colorRow[x] = ObsidianGL::Graphics::PackColor(color);

				if (depthWrite && depthRow)
				{
					depthRow[x] = depth;
				}
			}
		}
	}
#endif

#if defined(OBSIDIAN_USE_SSE) || defined(OBSIDIAN_USE_SSE41)
	void RasterizeTriangle3D_SSE(const glm::vec4& p_vertex0, const glm::vec4& p_vertex1, const glm::vec4& p_vertex2, const glm::vec4& p_color0, const glm::vec4& p_color1, const glm::vec4& p_color2, const glm::vec2& p_texCoord0, const glm::vec2& p_texCoord1, const glm::vec2& p_texCoord2)
	{
		if (!CurrentRenderState.ColorBuffer)
			return;

		const uint32_t width = CurrentRenderState.ColorBuffer->Width;
		const uint32_t height = CurrentRenderState.ColorBuffer->Height;

		const float invW0 = 1.0f / p_vertex0.w;
		float invW1 = 1.0f / p_vertex1.w;
		float invW2 = 1.0f / p_vertex2.w;

		const glm::vec3 ndc0 = glm::vec3(p_vertex0) * invW0;
		glm::vec3 ndc1 = glm::vec3(p_vertex1) * invW1;
		glm::vec3 ndc2 = glm::vec3(p_vertex2) * invW2;

		const float halfW = static_cast<float>(width) * 0.5f;
		const float halfH = static_cast<float>(height) * 0.5f;

		const glm::vec2 screen0((ndc0.x + 1.0f) * halfW, (1.0f - ndc0.y) * halfH);

		glm::vec2 screen1((ndc1.x + 1.0f) * halfW, (1.0f - ndc1.y) * halfH);
		glm::vec2 screen2((ndc2.x + 1.0f) * halfW, (1.0f - ndc2.y) * halfH);

		const float area = (screen1.x - screen0.x) * (screen2.y - screen0.y) - (screen1.y - screen0.y) * (screen2.x - screen0.x);

		if (std::abs(area) < 0.5f)
			return;

		if (CurrentRenderState.State & ObsidianGL::OGL_CULL_FACE)
		{
			const bool isFrontFacing = (CurrentRenderState.FrontFace == ObsidianGL::OGL_CCW) ? (area > 0.0f) : (area < 0.0f);

			if ((CurrentRenderState.CullFace == ObsidianGL::OGL_BACK && !isFrontFacing) || (CurrentRenderState.CullFace == ObsidianGL::OGL_FRONT && isFrontFacing))
				return;
		}

		const bool swapWinding = (area < 0.0f);

		if (swapWinding)
		{
			std::swap(screen1, screen2);
			std::swap(ndc1, ndc2);
			std::swap(invW1, invW2);
		}

		const auto& color1 = swapWinding ? p_color2 : p_color1;
		const auto& color2 = swapWinding ? p_color1 : p_color2;
		const auto& texCoord1 = swapWinding ? p_texCoord2 : p_texCoord1;
		const auto& texCoord2 = swapWinding ? p_texCoord1 : p_texCoord2;

		ObsidianGL::Geometry::Triangle triangle(screen0, screen1, screen2);

		if (triangle.FixedDoubledArea == 0.0f)
			return;

		const int32_t minX = std::max(triangle.BoundingBox.MinX, 0);
		const int32_t minY = std::max(triangle.BoundingBox.MinY, 0);
		const int32_t maxX = std::min(triangle.BoundingBox.MaxX, static_cast<int32_t>(width - 1));
		const int32_t maxY = std::min(triangle.BoundingBox.MaxY, static_cast<int32_t>(height - 1));

		if (maxX < minX || maxY < minY)
			return;

		const float depth0 = ndc0.z * 0.5f + 0.5f;
		const float depth1 = ndc1.z * 0.5f + 0.5f;
		const float depth2 = ndc2.z * 0.5f + 0.5f;

		const glm::vec4 c0w = p_color0 * invW0;
		const glm::vec4 c1w = color1 * invW1;
		const glm::vec4 c2w = color2 * invW2;

		const glm::vec2 t0w = p_texCoord0 * invW0;
		const glm::vec2 t1w = texCoord1 * invW1;
		const glm::vec2 t2w = texCoord2 * invW2;

		const bool depthTest = (CurrentRenderState.State & ObsidianGL::OGL_DEPTH_TEST) && CurrentRenderState.DepthBuffer;
		const bool depthWrite = (CurrentRenderState.State & ObsidianGL::OGL_DEPTH_WRITE) && CurrentRenderState.DepthBuffer;

		const ObsidianGL::Graphics::Texture* texture = GetActiveTexture();

		const __m128 vInvW0 = _mm_set1_ps(invW0);
		const __m128 vInvW1 = _mm_set1_ps(invW1);
		const __m128 vInvW2 = _mm_set1_ps(invW2);
		const __m128 vZero = _mm_setzero_ps();
		const __m128 vOne = _mm_set1_ps(1.0f);

		const auto baryConst = triangle.PrepareSIMDBarycentricSetup4();

		for (int32_t y = minY; y <= maxY; ++y)
		{
			RGBA8* colorRow = CurrentRenderState.ColorBuffer->GetRow(static_cast<uint32_t>(y));
			Depth* depthRow = (depthTest || depthWrite) ? CurrentRenderState.DepthBuffer->GetRow(static_cast<uint32_t>(y)) : nullptr;

			int32_t x = minX;

			for (; x + 3 <= maxX; x += 4)
			{
				__m128 vBaryU;
				__m128 vBaryV;
				__m128 vBaryW;
				__m128 insideMask;

				triangle.ComputeBarycentricCoordinates4(baryConst, x, y, vBaryW, vBaryU, vBaryV, insideMask);

				int insideBits = _mm_movemask_ps(insideMask);

				if (insideBits == 0)
					continue;

				__m128 vOneOverW = _mm_add_ps(_mm_add_ps(_mm_mul_ps(vBaryW, vInvW0), _mm_mul_ps(vBaryU, vInvW1)), _mm_mul_ps(vBaryV, vInvW2));
				__m128 vCorrectedW = _mm_div_ps(vOne, vOneOverW);

				__m128 vDepth = _mm_add_ps(_mm_add_ps(
					_mm_mul_ps(vBaryW, _mm_set1_ps(depth0)),
					_mm_mul_ps(vBaryU, _mm_set1_ps(depth1))),
					_mm_mul_ps(vBaryV, _mm_set1_ps(depth2)));

				__m128 depthValidMask = _mm_and_ps(_mm_cmpge_ps(vDepth, vZero), _mm_cmple_ps(vDepth, vOne));

				insideMask = _mm_and_ps(insideMask, depthValidMask);

				if (depthTest && depthRow)
				{
					__m128 vBufferDepth = _mm_loadu_ps(&depthRow[x]);
					__m128 depthPassMask = _mm_cmplt_ps(vDepth, vBufferDepth);

					insideMask = _mm_and_ps(insideMask, depthPassMask);
				}

				insideBits = _mm_movemask_ps(insideMask);

				if (insideBits == 0)
					continue;

				__m128i vTexels = _mm_setzero_si128();

				if (texture)
				{
					__m128 vTexU = _mm_mul_ps(_mm_add_ps(_mm_add_ps(
						_mm_mul_ps(vBaryW, _mm_set1_ps(t0w.x)),
						_mm_mul_ps(vBaryU, _mm_set1_ps(t1w.x))),
						_mm_mul_ps(vBaryV, _mm_set1_ps(t2w.x))), vCorrectedW);
					__m128 vTexV = _mm_mul_ps(_mm_add_ps(_mm_add_ps(
						_mm_mul_ps(vBaryW, _mm_set1_ps(t0w.y)),
						_mm_mul_ps(vBaryU, _mm_set1_ps(t1w.y))),
						_mm_mul_ps(vBaryV, _mm_set1_ps(t2w.y))), vCorrectedW);

					alignas(16) float texU[4];
					alignas(16) float texV[4];

					_mm_store_ps(texU, vTexU);
					_mm_store_ps(texV, vTexV);

					alignas(16) uint32_t texels[4];

					texture->SampleNearest4_SSE(texU, texV, texels);
					vTexels = _mm_load_si128(reinterpret_cast<const __m128i*>(texels));
				}

				__m128 vChan = _mm_mul_ps(_mm_add_ps(_mm_add_ps(
					_mm_mul_ps(vBaryW, _mm_set1_ps(c0w.r)),
					_mm_mul_ps(vBaryU, _mm_set1_ps(c1w.r))),
					_mm_mul_ps(vBaryV, _mm_set1_ps(c2w.r))), vCorrectedW);

				if (texture)
					vChan = _mm_mul_ps(vChan, _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(_mm_srli_epi32(vTexels, 24), _mm_set1_epi32(0xFF))), _mm_set1_ps(1.0f / 255.0f)));

				__m128i vColor = _mm_slli_epi32(_mm_cvtps_epi32(
					_mm_mul_ps(_mm_min_ps(_mm_max_ps(vChan, vZero), vOne), _mm_set1_ps(255.0f))), 24);

				vChan = _mm_mul_ps(_mm_add_ps(_mm_add_ps(
					_mm_mul_ps(vBaryW, _mm_set1_ps(c0w.g)),
					_mm_mul_ps(vBaryU, _mm_set1_ps(c1w.g))),
					_mm_mul_ps(vBaryV, _mm_set1_ps(c2w.g))), vCorrectedW);

				if (texture)
					vChan = _mm_mul_ps(vChan, _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(_mm_srli_epi32(vTexels, 16), _mm_set1_epi32(0xFF))), _mm_set1_ps(1.0f / 255.0f)));

				vColor = _mm_or_si128(vColor, _mm_slli_epi32(_mm_cvtps_epi32(
					_mm_mul_ps(_mm_min_ps(_mm_max_ps(vChan, vZero), vOne), _mm_set1_ps(255.0f))), 16));

				vChan = _mm_mul_ps(_mm_add_ps(_mm_add_ps(
					_mm_mul_ps(vBaryW, _mm_set1_ps(c0w.b)),
					_mm_mul_ps(vBaryU, _mm_set1_ps(c1w.b))),
					_mm_mul_ps(vBaryV, _mm_set1_ps(c2w.b))), vCorrectedW);

				if (texture)
					vChan = _mm_mul_ps(vChan, _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(_mm_srli_epi32(vTexels, 8), _mm_set1_epi32(0xFF))), _mm_set1_ps(1.0f / 255.0f)));

				vColor = _mm_or_si128(vColor, _mm_slli_epi32(_mm_cvtps_epi32(
					_mm_mul_ps(_mm_min_ps(_mm_max_ps(vChan, vZero), vOne), _mm_set1_ps(255.0f))), 8));

				vChan = _mm_mul_ps(_mm_add_ps(_mm_add_ps(
					_mm_mul_ps(vBaryW, _mm_set1_ps(c0w.a)),
					_mm_mul_ps(vBaryU, _mm_set1_ps(c1w.a))),
					_mm_mul_ps(vBaryV, _mm_set1_ps(c2w.a))), vCorrectedW);

				if (texture)
					vChan = _mm_mul_ps(vChan, _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(vTexels, _mm_set1_epi32(0xFF))), _mm_set1_ps(1.0f / 255.0f)));

				vColor = _mm_or_si128(vColor, _mm_cvtps_epi32(
					_mm_mul_ps(_mm_min_ps(_mm_max_ps(vChan, vZero), vOne), _mm_set1_ps(255.0f))));

#ifdef OBSIDIAN_USE_SSE41
				if (insideBits == 0x0F)
				{
					_mm_storeu_si128(reinterpret_cast<__m128i*>(&colorRow[x]), vColor);

					if (depthWrite && depthRow)
					{
						_mm_storeu_ps(&depthRow[x], vDepth);
					}
				}
				else
				{
					__m128 oldColors = _mm_castsi128_ps(_mm_loadu_si128(reinterpret_cast<const __m128i*>(&colorRow[x])));
					__m128 newColors = _mm_castsi128_ps(vColor);
					__m128 blended = _mm_blendv_ps(oldColors, newColors, insideMask);

					_mm_storeu_si128(reinterpret_cast<__m128i*>(&colorRow[x]), _mm_castps_si128(blended));

					if (depthWrite && depthRow)
					{
						__m128 oldDepth = _mm_loadu_ps(&depthRow[x]);
						__m128 newDepth = _mm_blendv_ps(oldDepth, vDepth, insideMask);
						_mm_storeu_ps(&depthRow[x], newDepth);
					}
				}
#else
				alignas(16) uint32_t colors[4];
				alignas(16) float depths[4];

				_mm_store_si128(reinterpret_cast<__m128i*>(colors), vColor);
				_mm_store_ps(depths, vDepth);

				for (int i = 0; i < 4; ++i)
				{
					if (insideBits & (1 << i))
					{
						colorRow[x + i] = colors[i];

						if (depthWrite && depthRow)
						{
							depthRow[x + i] = depths[i];
						}
					}
				}
#endif
			}

			// Scalar fallback for remaining pixels
			for (; x <= maxX; ++x)
			{
				float w;
				float u;
				float v;

				if (!triangle.ComputeBarycentricCoordinates(x, y, w, u, v))
					continue;

				const float oneOverW = w * invW0 + u * invW1 + v * invW2;
				const float correctedW = 1.0f / oneOverW;
				const float depth = w * depth0 + u * depth1 + v * depth2;

				if (depth < 0.0f || depth > 1.0f)
					continue;

				if (depthTest && depthRow && depth >= depthRow[x])
					continue;

				const glm::vec4 colorInterpolated = w * c0w + u * c1w + v * c2w;
				glm::vec4 color = colorInterpolated * correctedW;

				if (texture)
				{
					const glm::vec2 texInterpolated = w * t0w + u * t1w + v * t2w;
					const glm::vec2 texCoord = texInterpolated * correctedW;

					uint32_t texel = texture->Sample(texCoord.x, texCoord.y);

					float r;
					float g;
					float b;
					float a;

					ObsidianGL::Graphics::UnpackColorToFloat(texel, r, g, b, a);

					color.r *= r;
					color.g *= g;
					color.b *= b;
					color.a *= a;
				}

				colorRow[x] = ObsidianGL::Graphics::PackColor(color);

				if (depthWrite && depthRow)
				{
					depthRow[x] = depth;
				}
			}
		}
	}
#endif

	void RasterizeTriangle3D(const glm::vec4& p_vertex0, const glm::vec4& p_vertex1, const glm::vec4& p_vertex2, const glm::vec4& p_color0, const glm::vec4& p_color1, const glm::vec4& p_color2, const glm::vec2& p_texCoord0, const glm::vec2& p_texCoord1, const glm::vec2& p_texCoord2)
	{
#ifdef OBSIDIAN_USE_AVX2
		RasterizeTriangle3D_AVX2(p_vertex0, p_vertex1, p_vertex2, p_color0, p_color1, p_color2, p_texCoord0, p_texCoord1, p_texCoord2);
#elif defined(OBSIDIAN_USE_SSE) || defined(OBSIDIAN_USE_SSE41)
		RasterizeTriangle3D_SSE(p_vertex0, p_vertex1, p_vertex2, p_color0, p_color1, p_color2, p_texCoord0, p_texCoord1, p_texCoord2);
#else
		RasterizeTriangle3D_Scalar(p_vertex0, p_vertex1, p_vertex2, p_color0, p_color1, p_color2, p_texCoord0, p_texCoord1, p_texCoord2);
#endif
	}

	void RasterizePoint(const glm::vec4& p_clipPosition, const glm::vec4& p_color, const glm::vec2& p_texCoord0, const glm::vec2& p_texCoord1)
	{
		if (!CurrentRenderState.ColorBuffer)
			return;

		const float w = p_clipPosition.w;

		if (w <= 0.0f)
			return;

		const float absX = std::abs(p_clipPosition.x);
		const float absY = std::abs(p_clipPosition.y);
		const float absZ = std::abs(p_clipPosition.z);

		if (absX > w || absY > w || absZ > w)
			return;

		const uint32_t width = CurrentRenderState.ColorBuffer->Width;
		const uint32_t height = CurrentRenderState.ColorBuffer->Height;

		const float invW = 1.0f / w;
		const float ndcX = p_clipPosition.x * invW;
		const float ndcY = p_clipPosition.y * invW;
		const float ndcZ = p_clipPosition.z * invW;

		const float halfW = static_cast<float>(width) * 0.5f;
		const float halfH = static_cast<float>(height) * 0.5f;

		const float screenX = (ndcX + 1.0f) * halfW;
		const float screenY = (1.0f - ndcY) * halfH;

		const float depth = ndcZ * 0.5f + 0.5f;

		if (depth < 0.0f || depth > 1.0f)
			return;

		const bool depthTest = (CurrentRenderState.State & ObsidianGL::OGL_DEPTH_TEST) && CurrentRenderState.DepthBuffer;
		const bool depthWrite = (CurrentRenderState.State & ObsidianGL::OGL_DEPTH_WRITE) && CurrentRenderState.DepthBuffer;

		const ObsidianGL::Graphics::Texture* texture = GetActiveTexture();

		glm::vec4 finalColor = p_color;

		if (texture != nullptr)
		{
			uint32_t texel0 = texture->Sample(p_texCoord0.x, p_texCoord0.y);

			float r0;
			float g0;
			float b0;
			float a0;

			ObsidianGL::Graphics::UnpackColorToFloat(texel0, r0, g0, b0, a0);
			finalColor *= glm::vec4(r0, g0, b0, a0);

			uint32_t texel1 = texture->Sample(p_texCoord1.x, p_texCoord1.y);

			float r1;
			float g1;
			float b1;
			float a1;

			ObsidianGL::Graphics::UnpackColorToFloat(texel1, r1, g1, b1, a1);
			finalColor *= glm::vec4(r1, g1, b1, a1);
		}

		const uint32_t packedColor = ObsidianGL::Graphics::PackColor(finalColor);

		const float pointSize = CurrentRenderState.PointSize;
		const float halfSize = pointSize * 0.5f;

		const int32_t startX = std::max(0, static_cast<int32_t>(std::floor(screenX - halfSize)));
		const int32_t endX = std::min(static_cast<int32_t>(width - 1), static_cast<int32_t>(std::ceil(screenX + halfSize)));
		const int32_t startY = std::max(0, static_cast<int32_t>(std::floor(screenY - halfSize)));
		const int32_t endY = std::min(static_cast<int32_t>(height - 1), static_cast<int32_t>(std::ceil(screenY + halfSize)));

		if (startX > endX || startY > endY)
			return;

		for (int32_t y = startY; y <= endY; ++y)
		{
			RGBA8* colorRow = CurrentRenderState.ColorBuffer->GetRow(static_cast<uint32_t>(y));
			Depth* depthRow = (depthTest || depthWrite) ? CurrentRenderState.DepthBuffer->GetRow(static_cast<uint32_t>(y)) : nullptr;

			int32_t x = startX;

#ifdef OBSIDIAN_USE_AVX2
			const __m256 vDepth = _mm256_set1_ps(depth);
			const __m256i vPackedColor = _mm256_set1_epi32(static_cast<int32_t>(packedColor));

			for (; x + 7 <= endX; x += 8)
			{
				if (depthTest && depthRow)
				{
					__m256 vBufferDepth = _mm256_loadu_ps(&depthRow[x]);
					__m256 depthPassMask = _mm256_cmp_ps(vDepth, vBufferDepth, _CMP_LT_OQ);

					int passBits = _mm256_movemask_ps(depthPassMask);

					if (passBits == 0)
						continue;

					if (passBits == 0xFF)
					{
						_mm256_storeu_si256(reinterpret_cast<__m256i*>(&colorRow[x]), vPackedColor);

						if (depthWrite)
						{
							_mm256_storeu_ps(&depthRow[x], vDepth);
						}
					}
					else
					{
						__m256 oldColors = _mm256_castsi256_ps(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(&colorRow[x])));
						__m256 newColors = _mm256_castsi256_ps(vPackedColor);
						__m256 blended = _mm256_blendv_ps(oldColors, newColors, depthPassMask);

						_mm256_storeu_si256(reinterpret_cast<__m256i*>(&colorRow[x]), _mm256_castps_si256(blended));

						if (depthWrite)
						{
							__m256 oldDepth = _mm256_loadu_ps(&depthRow[x]);
							__m256 newDepth = _mm256_blendv_ps(oldDepth, vDepth, depthPassMask);
							_mm256_storeu_ps(&depthRow[x], newDepth);
						}
					}
				}
				else
				{
					_mm256_storeu_si256(reinterpret_cast<__m256i*>(&colorRow[x]), vPackedColor);

					if (depthWrite && depthRow)
					{
						_mm256_storeu_ps(&depthRow[x], vDepth);
					}
				}
			}
#elif defined(OBSIDIAN_USE_SSE) || defined(OBSIDIAN_USE_SSE41)
			const __m128 vDepth = _mm_set1_ps(depth);
			const __m128i vPackedColor = _mm_set1_epi32(static_cast<int32_t>(packedColor));

			for (; x + 3 <= endX; x += 4)
			{
				if (depthTest && depthRow)
				{
					__m128 vBufferDepth = _mm_loadu_ps(&depthRow[x]);
					__m128 depthPassMask = _mm_cmplt_ps(vDepth, vBufferDepth);

					int passBits = _mm_movemask_ps(depthPassMask);

					if (passBits == 0)
						continue;

#ifdef OBSIDIAN_USE_SSE41
					if (passBits == 0x0F)
					{
						_mm_storeu_si128(reinterpret_cast<__m128i*>(&colorRow[x]), vPackedColor);

						if (depthWrite)
						{
							_mm_storeu_ps(&depthRow[x], vDepth);
						}
					}
					else
					{
						__m128 oldColors = _mm_castsi128_ps(_mm_loadu_si128(reinterpret_cast<const __m128i*>(&colorRow[x])));
						__m128 newColors = _mm_castsi128_ps(vPackedColor);
						__m128 blended = _mm_blendv_ps(oldColors, newColors, depthPassMask);

						_mm_storeu_si128(reinterpret_cast<__m128i*>(&colorRow[x]), _mm_castps_si128(blended));

						if (depthWrite)
						{
							__m128 oldDepth = _mm_loadu_ps(&depthRow[x]);
							__m128 newDepth = _mm_blendv_ps(oldDepth, vDepth, depthPassMask);
							_mm_storeu_ps(&depthRow[x], newDepth);
						}
					}
#else
					for (int i = 0; i < 4; ++i)
					{
						if (passBits & (1 << i))
						{
							colorRow[x + i] = packedColor;

							if (depthWrite)
							{
								depthRow[x + i] = depth;
							}
						}
					}
#endif
				}
				else
				{
					_mm_storeu_si128(reinterpret_cast<__m128i*>(&colorRow[x]), vPackedColor);

					if (depthWrite && depthRow)
					{
						_mm_storeu_ps(&depthRow[x], vDepth);
					}
				}
			}
#endif

			for (; x <= endX; ++x)
			{
				if (depthTest && depthRow && depth >= depthRow[x])
					continue;

				colorRow[x] = packedColor;

				if (depthWrite && depthRow)
				{
					depthRow[x] = depth;
				}
			}
		}
	}

	void RasterizeLine(const glm::vec4& p_clipPos0, const glm::vec4& p_clipPos1, const glm::vec4& p_color0, const glm::vec4& p_color1, const glm::vec2& p_texCoord0, const glm::vec2& p_texCoord1)
	{
		if (!CurrentRenderState.ColorBuffer)
			return;

		const float w0 = p_clipPos0.w;
		const float w1 = p_clipPos1.w;

		if (w0 <= 0.0f && w1 <= 0.0f)
			return;

		if ((p_clipPos0.x < -w0 && p_clipPos1.x < -w1)
			|| (p_clipPos0.x > w0 && p_clipPos1.x > w1)
			|| (p_clipPos0.y < -w0 && p_clipPos1.y < -w1)
			|| (p_clipPos0.y > w0 && p_clipPos1.y > w1)
			|| (p_clipPos0.z < -w0 && p_clipPos1.z < -w1)
			|| (p_clipPos0.z > w0 && p_clipPos1.z > w1))
			return;

		const uint32_t width = CurrentRenderState.ColorBuffer->Width;
		const uint32_t height = CurrentRenderState.ColorBuffer->Height;

		auto ProjectToScreen = [&](const glm::vec4& clipPos, glm::vec2& screen, float& depth, float& invW) -> bool
			{
				if (clipPos.w <= 0.0f)
					return false;

				invW = 1.0f / clipPos.w;
				const float ndcX = clipPos.x * invW;
				const float ndcY = clipPos.y * invW;
				const float ndcZ = clipPos.z * invW;

				screen.x = (ndcX + 1.0f) * static_cast<float>(width) * 0.5f;
				screen.y = (1.0f - ndcY) * static_cast<float>(height) * 0.5f;
				depth = ndcZ * 0.5f + 0.5f;

				return true;
			};

		glm::vec2 screen0;
		glm::vec2 screen1;
		float depth0;
		float depth1;
		float invW0;
		float invW1;

		if (!ProjectToScreen(p_clipPos0, screen0, depth0, invW0))
			return;
		if (!ProjectToScreen(p_clipPos1, screen1, depth1, invW1))
			return;

		const bool depthTest = (CurrentRenderState.State & ObsidianGL::OGL_DEPTH_TEST) && CurrentRenderState.DepthBuffer;
		const bool depthWrite = (CurrentRenderState.State & ObsidianGL::OGL_DEPTH_WRITE) && CurrentRenderState.DepthBuffer;

		const ObsidianGL::Graphics::Texture* texture = GetActiveTexture();

		glm::vec4 finalColor0 = p_color0;
		glm::vec4 finalColor1 = p_color1;

		if (texture != nullptr)
		{
			uint32_t texel0 = texture->Sample(p_texCoord0.x, p_texCoord0.y);

			float r0;
			float g0;
			float b0;
			float a0;

			ObsidianGL::Graphics::UnpackColorToFloat(texel0, r0, g0, b0, a0);
			finalColor0 *= glm::vec4(r0, g0, b0, a0);

			uint32_t texel1 = texture->Sample(p_texCoord1.x, p_texCoord1.y);

			float r1;
			float g1;
			float b1;
			float a1;

			ObsidianGL::Graphics::UnpackColorToFloat(texel1, r1, g1, b1, a1);
			finalColor1 *= glm::vec4(r1, g1, b1, a1);
		}

		const float lineWidth = CurrentRenderState.LineWidth;
		const float halfWidth = lineWidth * 0.5f;
		const float halfWidthSq = halfWidth * halfWidth;

		const float deltaX = screen1.x - screen0.x;
		const float deltaY = screen1.y - screen0.y;
		const float lineLengthSq = deltaX * deltaX + deltaY * deltaY;

		if (lineLengthSq < 0.0001f)
		{
			const float midDepth = (depth0 + depth1) * 0.5f;
			const glm::vec4 midColor = (finalColor0 + finalColor1) * 0.5f;
			const uint32_t packedColor = ObsidianGL::Graphics::PackColor(midColor);

			const int32_t startX = std::max(0, static_cast<int32_t>(std::floor(screen0.x - halfWidth)));
			const int32_t endX = std::min(static_cast<int32_t>(width - 1), static_cast<int32_t>(std::ceil(screen0.x + halfWidth)));
			const int32_t startY = std::max(0, static_cast<int32_t>(std::floor(screen0.y - halfWidth)));
			const int32_t endY = std::min(static_cast<int32_t>(height - 1), static_cast<int32_t>(std::ceil(screen0.y + halfWidth)));

			for (int32_t y = startY; y <= endY; ++y)
			{
				RGBA8* colorRow = CurrentRenderState.ColorBuffer->GetRow(static_cast<uint32_t>(y));
				Depth* depthRow = (depthTest || depthWrite) ? CurrentRenderState.DepthBuffer->GetRow(static_cast<uint32_t>(y)) : nullptr;

				for (int32_t x = startX; x <= endX; ++x)
				{
					const float dx = static_cast<float>(x) + 0.5f - screen0.x;
					const float dy = static_cast<float>(y) + 0.5f - screen0.y;

					if (dx * dx + dy * dy > halfWidthSq)
						continue;

					if (depthTest && depthRow && midDepth >= depthRow[x])
						continue;

					colorRow[x] = packedColor;

					if (depthWrite && depthRow)
					{
						depthRow[x] = midDepth;
					}
				}
			}

			return;
		}

		const float invLineLengthSq = 1.0f / lineLengthSq;

		const int32_t screenWidthMax = static_cast<int32_t>(width - 1);
		const int32_t screenHeightMax = static_cast<int32_t>(height - 1);

		const int32_t startX = std::max(0, static_cast<int32_t>(std::floor(std::min(screen0.x, screen1.x) - halfWidth - 1.0f)));
		const int32_t endX = std::min(screenWidthMax, static_cast<int32_t>(std::ceil(std::max(screen0.x, screen1.x) + halfWidth + 1.0f)));
		const int32_t startY = std::max(0, static_cast<int32_t>(std::floor(std::min(screen0.y, screen1.y) - halfWidth - 1.0f)));
		const int32_t endY = std::min(screenHeightMax, static_cast<int32_t>(std::ceil(std::max(screen0.y, screen1.y) + halfWidth + 1.0f)));

		if (startX > endX || startY > endY)
			return;

		for (int32_t y = startY; y <= endY; ++y)
		{
			RGBA8* colorRow = CurrentRenderState.ColorBuffer->GetRow(static_cast<uint32_t>(y));
			Depth* depthRow = (depthTest || depthWrite) ? CurrentRenderState.DepthBuffer->GetRow(static_cast<uint32_t>(y)) : nullptr;

			const float pixelCenterY = static_cast<float>(y) + 0.5f;
			const float pyMinusS0Y = pixelCenterY - screen0.y;

			int32_t x = startX;

#ifdef OBSIDIAN_USE_AVX2
			const __m256 vDeltaX = _mm256_set1_ps(deltaX);
			const __m256 vDeltaY = _mm256_set1_ps(deltaY);

			const __m256 vInvLineLenSq = _mm256_set1_ps(invLineLengthSq);

			const __m256 vScreen0X = _mm256_set1_ps(screen0.x);
			const __m256 vScreen0Y = _mm256_set1_ps(screen0.y);

			const __m256 vHalfWidthSq = _mm256_set1_ps(halfWidthSq);

			const __m256 vDepth0 = _mm256_set1_ps(depth0);
			const __m256 vDepth1 = _mm256_set1_ps(depth1);

			const __m256 vPyMinusS0Y = _mm256_set1_ps(pyMinusS0Y);

			const __m256 vColor0R = _mm256_set1_ps(finalColor0.r);
			const __m256 vColor0G = _mm256_set1_ps(finalColor0.g);
			const __m256 vColor0B = _mm256_set1_ps(finalColor0.b);
			const __m256 vColor0A = _mm256_set1_ps(finalColor0.a);
			const __m256 vColor1R = _mm256_set1_ps(finalColor1.r);
			const __m256 vColor1G = _mm256_set1_ps(finalColor1.g);
			const __m256 vColor1B = _mm256_set1_ps(finalColor1.b);
			const __m256 vColor1A = _mm256_set1_ps(finalColor1.a);

			for (; x + 7 <= endX; x += 8)
			{
				__m256 vPixelX = _mm256_add_ps(_mm256_add_ps(_mm256_set1_ps(static_cast<float>(x)), _mm256_set_ps(7.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f)), _mm256_set1_ps(0.5f));

				__m256 vPxMinusS0X = _mm256_sub_ps(vPixelX, vScreen0X);

				__m256 vT = _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(vPxMinusS0X, vDeltaX), _mm256_mul_ps(vPyMinusS0Y, vDeltaY)), vInvLineLenSq);

				vT = _mm256_min_ps(_mm256_max_ps(vT, _mm256_setzero_ps()), _mm256_set1_ps(1.0f));

				__m256 vClosestX = _mm256_add_ps(vScreen0X, _mm256_mul_ps(vT, vDeltaX));
				__m256 vClosestY = _mm256_add_ps(vScreen0Y, _mm256_mul_ps(vT, vDeltaY));

				__m256 vDx = _mm256_sub_ps(vPixelX, vClosestX);
				__m256 vDy = _mm256_sub_ps(_mm256_set1_ps(pixelCenterY), vClosestY);
				__m256 vDistSq = _mm256_add_ps(_mm256_mul_ps(vDx, vDx), _mm256_mul_ps(vDy, vDy));

				__m256 insideMask = _mm256_cmp_ps(vDistSq, vHalfWidthSq, _CMP_LE_OQ);

				int insideBits = _mm256_movemask_ps(insideMask);

				if (insideBits == 0)
					continue;

				__m256 vDepth = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(_mm256_set1_ps(1.0f), vT), vDepth0), _mm256_mul_ps(vT, vDepth1));

				__m256 depthValidMask = _mm256_and_ps(_mm256_cmp_ps(vDepth, _mm256_setzero_ps(), _CMP_GE_OQ), _mm256_cmp_ps(vDepth, _mm256_set1_ps(1.0f), _CMP_LE_OQ));

				insideMask = _mm256_and_ps(insideMask, depthValidMask);

				if (depthTest && depthRow)
				{
					__m256 vBufferDepth = _mm256_loadu_ps(&depthRow[x]);
					__m256 depthPassMask = _mm256_cmp_ps(vDepth, vBufferDepth, _CMP_LT_OQ);

					insideMask = _mm256_and_ps(insideMask, depthPassMask);
				}

				insideBits = _mm256_movemask_ps(insideMask);

				if (insideBits == 0)
					continue;

				__m256 vOneMinusT = _mm256_sub_ps(_mm256_set1_ps(1.0f), vT);

				__m256 vR = _mm256_add_ps(_mm256_mul_ps(vOneMinusT, vColor0R), _mm256_mul_ps(vT, vColor1R));
				__m256 vG = _mm256_add_ps(_mm256_mul_ps(vOneMinusT, vColor0G), _mm256_mul_ps(vT, vColor1G));
				__m256 vB = _mm256_add_ps(_mm256_mul_ps(vOneMinusT, vColor0B), _mm256_mul_ps(vT, vColor1B));
				__m256 vA = _mm256_add_ps(_mm256_mul_ps(vOneMinusT, vColor0A), _mm256_mul_ps(vT, vColor1A));

				__m256i vColor = _mm256_slli_epi32(_mm256_cvtps_epi32(_mm256_mul_ps(_mm256_min_ps(_mm256_max_ps(vR, _mm256_setzero_ps()), _mm256_set1_ps(1.0f)), _mm256_set1_ps(255.0f))), 24);
				vColor = _mm256_or_si256(vColor, _mm256_slli_epi32(_mm256_cvtps_epi32(_mm256_mul_ps(_mm256_min_ps(_mm256_max_ps(vG, _mm256_setzero_ps()), _mm256_set1_ps(1.0f)), _mm256_set1_ps(255.0f))), 16));
				vColor = _mm256_or_si256(vColor, _mm256_slli_epi32(_mm256_cvtps_epi32(_mm256_mul_ps(_mm256_min_ps(_mm256_max_ps(vB, _mm256_setzero_ps()), _mm256_set1_ps(1.0f)), _mm256_set1_ps(255.0f))), 8));
				vColor = _mm256_or_si256(vColor, _mm256_cvtps_epi32(_mm256_mul_ps(_mm256_min_ps(_mm256_max_ps(vA, _mm256_setzero_ps()), _mm256_set1_ps(1.0f)), _mm256_set1_ps(255.0f))));

				if (insideBits == 0xFF)
				{
					_mm256_storeu_si256(reinterpret_cast<__m256i*>(&colorRow[x]), vColor);

					if (depthWrite && depthRow)
					{
						_mm256_storeu_ps(&depthRow[x], vDepth);
					}
				}
				else
				{
					__m256 oldColors = _mm256_castsi256_ps(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(&colorRow[x])));
					__m256 newColors = _mm256_castsi256_ps(vColor);
					__m256 blended = _mm256_blendv_ps(oldColors, newColors, insideMask);

					_mm256_storeu_si256(reinterpret_cast<__m256i*>(&colorRow[x]), _mm256_castps_si256(blended));

					if (depthWrite && depthRow)
					{
						__m256 oldDepth = _mm256_loadu_ps(&depthRow[x]);
						__m256 newDepth = _mm256_blendv_ps(oldDepth, vDepth, insideMask);
						_mm256_storeu_ps(&depthRow[x], newDepth);
					}
				}
			}
#elif defined(OBSIDIAN_USE_SSE) || defined(OBSIDIAN_USE_SSE41)
			const __m128 vDeltaX = _mm_set1_ps(deltaX);
			const __m128 vDeltaY = _mm_set1_ps(deltaY);

			const __m128 vInvLineLenSq = _mm_set1_ps(invLineLengthSq);

			const __m128 vScreen0X = _mm_set1_ps(screen0.x);
			const __m128 vScreen0Y = _mm_set1_ps(screen0.y);

			const __m128 vHalfWidthSq = _mm_set1_ps(halfWidthSq);

			const __m128 vDepth0 = _mm_set1_ps(depth0);
			const __m128 vDepth1 = _mm_set1_ps(depth1);

			const __m128 vZero = _mm_setzero_ps();

			const __m128 vPyMinusS0Y = _mm_set1_ps(pyMinusS0Y);

			const __m128 vColor0R = _mm_set1_ps(finalColor0.r);
			const __m128 vColor0G = _mm_set1_ps(finalColor0.g);
			const __m128 vColor0B = _mm_set1_ps(finalColor0.b);
			const __m128 vColor0A = _mm_set1_ps(finalColor0.a);
			const __m128 vColor1R = _mm_set1_ps(finalColor1.r);
			const __m128 vColor1G = _mm_set1_ps(finalColor1.g);
			const __m128 vColor1B = _mm_set1_ps(finalColor1.b);
			const __m128 vColor1A = _mm_set1_ps(finalColor1.a);

			for (; x + 3 <= endX; x += 4)
			{
				__m128 vPixelX = _mm_add_ps(_mm_add_ps(_mm_set1_ps(static_cast<float>(x)), _mm_set_ps(3.0f, 2.0f, 1.0f, 0.0f)), _mm_set1_ps(0.5f));

				__m128 vPxMinusS0X = _mm_sub_ps(vPixelX, vScreen0X);

				__m128 vT = _mm_mul_ps(_mm_add_ps(_mm_mul_ps(vPxMinusS0X, vDeltaX), _mm_mul_ps(vPyMinusS0Y, vDeltaY)), vInvLineLenSq);

				vT = _mm_min_ps(_mm_max_ps(vT, vZero), _mm_set1_ps(1.0f));

				__m128 vClosestX = _mm_add_ps(vScreen0X, _mm_mul_ps(vT, vDeltaX));
				__m128 vClosestY = _mm_add_ps(vScreen0Y, _mm_mul_ps(vT, vDeltaY));

				__m128 vDx = _mm_sub_ps(vPixelX, vClosestX);
				__m128 vDy = _mm_sub_ps(_mm_set1_ps(pixelCenterY), vClosestY);
				__m128 vDistSq = _mm_add_ps(_mm_mul_ps(vDx, vDx), _mm_mul_ps(vDy, vDy));

				__m128 insideMask = _mm_cmple_ps(vDistSq, vHalfWidthSq);

				int insideBits = _mm_movemask_ps(insideMask);

				if (insideBits == 0)
					continue;

				__m128 vDepth = _mm_add_ps(_mm_mul_ps(_mm_sub_ps(_mm_set1_ps(1.0f), vT), vDepth0), _mm_mul_ps(vT, vDepth1));

				__m128 depthValidMask = _mm_and_ps(_mm_cmpge_ps(vDepth, vZero), _mm_cmple_ps(vDepth, _mm_set1_ps(1.0f)));

				insideMask = _mm_and_ps(insideMask, depthValidMask);

				if (depthTest && depthRow)
				{
					__m128 vBufferDepth = _mm_loadu_ps(&depthRow[x]);
					__m128 depthPassMask = _mm_cmplt_ps(vDepth, vBufferDepth);

					insideMask = _mm_and_ps(insideMask, depthPassMask);
				}

				insideBits = _mm_movemask_ps(insideMask);

				if (insideBits == 0)
					continue;

				__m128 vOneMinusT = _mm_sub_ps(_mm_set1_ps(1.0f), vT);

				__m128 vR = _mm_add_ps(_mm_mul_ps(vOneMinusT, vColor0R), _mm_mul_ps(vT, vColor1R));
				__m128 vG = _mm_add_ps(_mm_mul_ps(vOneMinusT, vColor0G), _mm_mul_ps(vT, vColor1G));
				__m128 vB = _mm_add_ps(_mm_mul_ps(vOneMinusT, vColor0B), _mm_mul_ps(vT, vColor1B));
				__m128 vA = _mm_add_ps(_mm_mul_ps(vOneMinusT, vColor0A), _mm_mul_ps(vT, vColor1A));

				__m128i vColor = _mm_slli_epi32(_mm_cvtps_epi32(_mm_mul_ps(_mm_min_ps(_mm_max_ps(vR, vZero), _mm_set1_ps(1.0f)), _mm_set1_ps(255.0f))), 24);
				vColor = _mm_or_si128(vColor, _mm_slli_epi32(_mm_cvtps_epi32(_mm_mul_ps(_mm_min_ps(_mm_max_ps(vG, vZero), _mm_set1_ps(1.0f)), _mm_set1_ps(255.0f))), 16));
				vColor = _mm_or_si128(vColor, _mm_slli_epi32(_mm_cvtps_epi32(_mm_mul_ps(_mm_min_ps(_mm_max_ps(vB, vZero), _mm_set1_ps(1.0f)), _mm_set1_ps(255.0f))), 8));
				vColor = _mm_or_si128(vColor, _mm_cvtps_epi32(_mm_mul_ps(_mm_min_ps(_mm_max_ps(vA, vZero), _mm_set1_ps(1.0f)), _mm_set1_ps(255.0f))));

#ifdef OBSIDIAN_USE_SSE41
				if (insideBits == 0x0F)
				{
					_mm_storeu_si128(reinterpret_cast<__m128i*>(&colorRow[x]), vColor);

					if (depthWrite && depthRow)
					{
						_mm_storeu_ps(&depthRow[x], vDepth);
					}
				}
				else
				{
					__m128 oldColors = _mm_castsi128_ps(_mm_loadu_si128(reinterpret_cast<const __m128i*>(&colorRow[x])));
					__m128 newColors = _mm_castsi128_ps(vColor);
					__m128 blended = _mm_blendv_ps(oldColors, newColors, insideMask);

					_mm_storeu_si128(reinterpret_cast<__m128i*>(&colorRow[x]), _mm_castps_si128(blended));

					if (depthWrite && depthRow)
					{
						__m128 oldDepth = _mm_loadu_ps(&depthRow[x]);
						__m128 newDepth = _mm_blendv_ps(oldDepth, vDepth, insideMask);
						_mm_storeu_ps(&depthRow[x], newDepth);
					}
				}
#else
				alignas(16) uint32_t colors[4];
				alignas(16) float depths[4];

				_mm_store_si128(reinterpret_cast<__m128i*>(colors), vColor);
				_mm_store_ps(depths, vDepth);

				for (int i = 0; i < 4; ++i)
				{
					if (insideBits & (1 << i))
					{
						colorRow[x + i] = colors[i];

						if (depthWrite && depthRow)
						{
							depthRow[x + i] = depths[i];
						}
					}
				}
#endif
			}
#endif

			for (; x <= endX; ++x)
			{
				const float pixelCenterX = static_cast<float>(x) + 0.5f;

				const float pxMinusS0X = pixelCenterX - screen0.x;

				float t = (pxMinusS0X * deltaX + pyMinusS0Y * deltaY) * invLineLengthSq;
				t = std::max(0.0f, std::min(1.0f, t));

				const float closestX = screen0.x + t * deltaX;
				const float closestY = screen0.y + t * deltaY;

				const float dx = pixelCenterX - closestX;
				const float dy = pixelCenterY - closestY;
				const float distSq = dx * dx + dy * dy;

				if (distSq > halfWidthSq)
					continue;

				const float depth = depth0 * (1.0f - t) + depth1 * t;

				if (depth < 0.0f || depth > 1.0f)
					continue;

				if (depthTest && depthRow && depth >= depthRow[x])
					continue;

				const glm::vec4 color = finalColor0 * (1.0f - t) + finalColor1 * t;

				colorRow[x] = ObsidianGL::Graphics::PackColor(color);

				if (depthWrite && depthRow)
				{
					depthRow[x] = depth;
				}
			}
		}
	}

	void RasterizeTriangleWireframe(const glm::vec4& p_clipPos0, const glm::vec4& p_clipPos1, const glm::vec4& p_clipPos2, const glm::vec4& p_color0, const glm::vec4& p_color1, const glm::vec4& p_color2, const glm::vec2& p_texCoord0, const glm::vec2& p_texCoord1, const glm::vec2& p_texCoord2)
	{
		RasterizeLine(p_clipPos0, p_clipPos1, p_color0, p_color1, p_texCoord0, p_texCoord1);
		RasterizeLine(p_clipPos1, p_clipPos2, p_color1, p_color2, p_texCoord1, p_texCoord2);
		RasterizeLine(p_clipPos2, p_clipPos0, p_color2, p_color0, p_texCoord2, p_texCoord0);
	}

	void RasterizeTrianglePoints(const glm::vec4& p_clipPos0, const glm::vec4& p_clipPos1, const glm::vec4& p_clipPos2, const glm::vec4& p_color0, const glm::vec4& p_color1, const glm::vec4& p_color2, const glm::vec2& p_texCoord0, const glm::vec2& p_texCoord1, const glm::vec2& p_texCoord2)
	{
		RasterizePoint(p_clipPos0, p_color0, p_texCoord0, p_texCoord1);
		RasterizePoint(p_clipPos1, p_color1, p_texCoord1, p_texCoord2);
		RasterizePoint(p_clipPos2, p_color2, p_texCoord2, p_texCoord0);
	}

	uint16_t DeterminePolygonMode(const glm::vec2& p_screen0, const glm::vec2& p_screen1, const glm::vec2& p_screen2)
	{
		const float area = (p_screen1.x - p_screen0.x) * (p_screen2.y - p_screen0.y) - (p_screen1.y - p_screen0.y) * (p_screen2.x - p_screen0.x);

		const bool isFrontFacing = (CurrentRenderState.FrontFace == ObsidianGL::OGL_CCW) ? (area > 0.0f) : (area < 0.0f);

		return isFrontFacing ? CurrentRenderState.PolygonModeFront : CurrentRenderState.PolygonModeBack;
	}

	void DispatchTriangle(const glm::vec4& p_clipPos0, const glm::vec4& p_clipPos1, const glm::vec4& p_clipPos2, const glm::vec4& p_color0, const glm::vec4& p_color1, const glm::vec4& p_color2, const glm::vec2& p_texCoord0, const glm::vec2& p_texCoord1, const glm::vec2& p_texCoord2)
	{
		if (CurrentRenderState.PolygonModeFront == ObsidianGL::OGL_FILL && CurrentRenderState.PolygonModeBack == ObsidianGL::OGL_FILL)
		{
			RasterizeTriangle3D(p_clipPos0, p_clipPos1, p_clipPos2, p_color0, p_color1, p_color2, p_texCoord0, p_texCoord1, p_texCoord2);

			return;
		}

		const float w0 = p_clipPos0.w;
		const float w1 = p_clipPos1.w;
		const float w2 = p_clipPos2.w;

		if (w0 <= 0.0f || w1 <= 0.0f || w2 <= 0.0f)
		{
			switch (CurrentRenderState.PolygonModeFront)
			{
			case ObsidianGL::OGL_FILL:
				RasterizeTriangle3D(p_clipPos0, p_clipPos1, p_clipPos2, p_color0, p_color1, p_color2, p_texCoord0, p_texCoord1, p_texCoord2);
				break;

			case ObsidianGL::OGL_LINE:
				RasterizeTriangleWireframe(p_clipPos0, p_clipPos1, p_clipPos2, p_color0, p_color1, p_color2, p_texCoord0, p_texCoord1, p_texCoord2);
				break;

			case ObsidianGL::OGL_POINT:
				RasterizeTrianglePoints(p_clipPos0, p_clipPos1, p_clipPos2, p_color0, p_color1, p_color2, p_texCoord0, p_texCoord1, p_texCoord2);
				break;
			}

			return;
		}

		const uint32_t screenWidth = CurrentRenderState.ColorBuffer->Width;
		const uint32_t screenHeight = CurrentRenderState.ColorBuffer->Height;

		const float halfW = static_cast<float>(screenWidth) * 0.5f;
		const float halfH = static_cast<float>(screenHeight) * 0.5f;

		const glm::vec2 s0((p_clipPos0.x / w0 + 1.0f) * halfW, (1.0f - p_clipPos0.y / w0) * halfH);
		const glm::vec2 s1((p_clipPos1.x / w1 + 1.0f) * halfW, (1.0f - p_clipPos1.y / w1) * halfH);
		const glm::vec2 s2((p_clipPos2.x / w2 + 1.0f) * halfW, (1.0f - p_clipPos2.y / w2) * halfH);

		const uint16_t mode = DeterminePolygonMode(s0, s1, s2);

		if (mode == ObsidianGL::OGL_LINE)
		{
			RasterizeTriangleWireframe(p_clipPos0, p_clipPos1, p_clipPos2, p_color0, p_color1, p_color2, p_texCoord0, p_texCoord1, p_texCoord2);
		}
		else if (mode == ObsidianGL::OGL_POINT)
		{
			RasterizeTrianglePoints(p_clipPos0, p_clipPos1, p_clipPos2, p_color0, p_color1, p_color2, p_texCoord0, p_texCoord1, p_texCoord2);
		}
		else
		{
			RasterizeTriangle3D(p_clipPos0, p_clipPos1, p_clipPos2, p_color0, p_color1, p_color2, p_texCoord0, p_texCoord1, p_texCoord2);
		}
	}

	void ClipAndRasterizeTriangle(const glm::vec4& p_clipPosition0, const glm::vec4& p_clipPosition1, const glm::vec4& p_clipPosition2, const glm::vec4& p_color0, const glm::vec4& p_color1, const glm::vec4& p_color2, const glm::vec3& p_normal0, const glm::vec3& p_normal1, const glm::vec3& p_normal2, const glm::vec2& p_texCoord0, const glm::vec2& p_texCoord1, const glm::vec2& p_texCoord2)
	{
		uint8_t code0 = ObsidianGL::Geometry::Clipper::ComputeOutcode(p_clipPosition0);
		uint8_t code1 = ObsidianGL::Geometry::Clipper::ComputeOutcode(p_clipPosition1);
		uint8_t code2 = ObsidianGL::Geometry::Clipper::ComputeOutcode(p_clipPosition2);

		if ((code0 | code1 | code2) == ObsidianGL::Geometry::Clipper::OUTCODE_INSIDE)
		{
			DispatchTriangle(p_clipPosition0, p_clipPosition1, p_clipPosition2, p_color0, p_color1, p_color2, p_texCoord0, p_texCoord1, p_texCoord2);
			return;
		}

		if ((code0 & code1 & code2) != ObsidianGL::Geometry::Clipper::OUTCODE_INSIDE)
		{
			return;
		}

		ObsidianGL::Geometry::ClipVertex v0{ .Position = p_clipPosition0, .Color = p_color0, .Normal = p_normal0, .TexCoord = p_texCoord0 };
		ObsidianGL::Geometry::ClipVertex v1{ .Position = p_clipPosition1, .Color = p_color1, .Normal = p_normal1, .TexCoord = p_texCoord1 };
		ObsidianGL::Geometry::ClipVertex v2{ .Position = p_clipPosition2, .Color = p_color2, .Normal = p_normal2, .TexCoord = p_texCoord2 };

		ObsidianGL::Geometry::ClipPolygon clippedPoly;
		if (!ObsidianGL::Geometry::Clipper::ClipTriangle(v0, v1, v2, clippedPoly))
		{
			return;
		}

		for (uint8_t i = 1; i + 1 < clippedPoly.VertexCount; ++i)
		{
			DispatchTriangle(
				clippedPoly.Vertices[0].Position,
				clippedPoly.Vertices[i].Position,
				clippedPoly.Vertices[i + 1].Position,
				clippedPoly.Vertices[0].Color,
				clippedPoly.Vertices[i].Color,
				clippedPoly.Vertices[i + 1].Color,
				clippedPoly.Vertices[0].TexCoord,
				clippedPoly.Vertices[i].TexCoord,
				clippedPoly.Vertices[i + 1].TexCoord
			);
		}
	}

	static std::vector<glm::vec4> clipPositions;

	void FlushVertices()
	{
		if (VertexBuffer.empty())
			return;

		const glm::mat4 viewProjectionMatrix = CurrentRenderState.ProjectionStack.Top() * CurrentRenderState.ModelViewStack.Top();

		clipPositions.clear();
		clipPositions.reserve(VertexBuffer.size());

		for (const auto& immediateVertex : VertexBuffer)
		{
			clipPositions.push_back(viewProjectionMatrix * glm::vec4(immediateVertex.Position, 1.0f));
		}

		auto ProcessTriangle = [&](size_t index0, size_t index1, size_t index2)
			{
				ClipAndRasterizeTriangle(
					clipPositions[index0], clipPositions[index1], clipPositions[index2],
					VertexBuffer[index0].Color, VertexBuffer[index1].Color, VertexBuffer[index2].Color,
					VertexBuffer[index0].Normal, VertexBuffer[index1].Normal, VertexBuffer[index2].Normal,
					VertexBuffer[index0].TexCoord, VertexBuffer[index1].TexCoord, VertexBuffer[index2].TexCoord
				);
			};

		switch (CurrentRenderState.CurrentPrimitiveMode)
		{
		case ObsidianGL::OGL_TRIANGLES:
			for (size_t i = 0; i + 2 < VertexBuffer.size(); i += 3)
			{
				ProcessTriangle(i, i + 1, i + 2);
			}
			break;
		case ObsidianGL::OGL_QUADS:
			for (size_t i = 0; i + 3 < VertexBuffer.size(); i += 4)
			{
				ProcessTriangle(i, i + 1, i + 2);
				ProcessTriangle(i, i + 2, i + 3);
			}
			break;
		case ObsidianGL::OGL_TRIANGLE_STRIP:
			for (size_t i = 0; i + 2 < VertexBuffer.size(); ++i)
			{
				if (i % 2 == 0)
				{
					ProcessTriangle(i, i + 1, i + 2);
				}
				else
				{
					ProcessTriangle(i + 1, i, i + 2);
				}
			}
			break;
		case ObsidianGL::OGL_TRIANGLE_FAN:
			for (size_t i = 1; i + 1 < VertexBuffer.size(); ++i)
			{
				ProcessTriangle(0, i, i + 1);
			}
			break;
		default:
			for (size_t i = 0; i + 2 < VertexBuffer.size(); i += 3)
			{
				ProcessTriangle(i, i + 1, i + 2);
			}
			break;
		}

		VertexBuffer.clear();
	}
}

void ObsidianGL::Test()
{
	const uint32_t width = CurrentRenderState.ColorBuffer->Width;
	const uint32_t height = CurrentRenderState.ColorBuffer->Height;
	const float invHeight = 1.0f / static_cast<float>(height - 1);
	const float invWidth = 1.0f / static_cast<float>(width - 1);

	for (uint16_t y = 0; y < CurrentRenderState.ColorBuffer->Height; ++y)
	{
		const float v = static_cast<float>(y) * invHeight;

		uint8_t  g = static_cast<uint8_t>(v * 255.0f);

		RGBA8* row = CurrentRenderState.ColorBuffer->GetRow(y);

		for (uint16_t x = 0; x < CurrentRenderState.ColorBuffer->Width; ++x)
		{
			const float u = static_cast<float>(x) * invWidth;

			uint8_t  r = static_cast<uint8_t>(u * 255.0f);
			uint8_t  b = static_cast<uint8_t>((1.0f - u) * 255.0f);

			row[x] = (r << 24) | (g << 16) | (b << 8) | 255;
		}
	}
}

void ObsidianGL::Initialize(uint16_t p_width, uint16_t p_height)
{
	CurrentRenderState.ColorBuffer = new Buffers::ColorBuffer(p_width, p_height);
	CurrentRenderState.DepthBuffer = new Buffers::DepthBuffer(p_width, p_height);
	CurrentRenderState.ColorBuffer->Clear();
	CurrentRenderState.DepthBuffer->Clear();
}

void ObsidianGL::Shutdown()
{
	delete CurrentRenderState.ColorBuffer;
	delete CurrentRenderState.DepthBuffer;
	CurrentRenderState.ColorBuffer = nullptr;
	CurrentRenderState.DepthBuffer = nullptr;
}

void ObsidianGL::SetClearColor(float p_r, float p_g, float p_b, float p_a)
{
	if (CurrentRenderState.ColorBuffer)
	{
		CurrentRenderState.ColorBuffer->SetClearColor(p_r, p_g, p_b, p_a);
	}
}

void ObsidianGL::Clear(uint8_t p_flags)
{
	if ((p_flags & OGL_COLOR_BUFFER_BIT) && CurrentRenderState.ColorBuffer)
	{
		CurrentRenderState.ColorBuffer->Clear();
	}

	if ((p_flags & OGL_DEPTH_BUFFER_BIT) && CurrentRenderState.DepthBuffer)
	{
		CurrentRenderState.DepthBuffer->Clear();
	}
}

void ObsidianGL::Enable(uint8_t p_state)
{
	CurrentRenderState.State |= p_state;
}

void ObsidianGL::Disable(uint8_t p_state)
{
	CurrentRenderState.State &= ~p_state;
}

void ObsidianGL::EnableTexture2D()
{
	CurrentRenderState.TextureEnabled = true;
}

void ObsidianGL::DisableTexture2D()
{
	CurrentRenderState.TextureEnabled = false;
}

void ObsidianGL::CullFace(uint16_t p_face)
{
	CurrentRenderState.CullFace = p_face;
}

void ObsidianGL::FrontFace(uint16_t p_mode)
{
	CurrentRenderState.FrontFace = p_mode;
}

void ObsidianGL::PolygonMode(uint16_t p_face, uint16_t p_mode)
{
	if (p_face == OGL_FRONT || p_face == OGL_FRONT_AND_BACK)
	{
		CurrentRenderState.PolygonModeFront = p_mode;
	}

	if (p_face == OGL_BACK || p_face == OGL_FRONT_AND_BACK)
	{
		CurrentRenderState.PolygonModeBack = p_mode;
	}
}

void ObsidianGL::PointSize(float p_size)
{
	CurrentRenderState.PointSize = std::max(p_size, 1.0f);
}

void ObsidianGL::LineWidth(float p_width)
{
	CurrentRenderState.LineWidth = std::max(p_width, 1.0f);
}

void ObsidianGL::Begin(uint16_t p_mode)
{
	CurrentRenderState.InBeginEnd = true;
	CurrentRenderState.CurrentPrimitiveMode = p_mode;
	VertexBuffer.clear();
}

void ObsidianGL::End()
{
	CurrentRenderState.InBeginEnd = false;
	FlushVertices();
}

void ObsidianGL::Vertex3f(float p_x, float p_y, float p_z)
{
	if (!CurrentRenderState.InBeginEnd)
		return;

	ImmediateVertex immediateVertex;
	immediateVertex.Position = glm::vec3(p_x, p_y, p_z);
	immediateVertex.Color = CurrentRenderState.CurrentColor;
	immediateVertex.Normal = CurrentRenderState.CurrentNormal;
	immediateVertex.TexCoord = CurrentRenderState.CurrentTexCoord;
	VertexBuffer.push_back(immediateVertex);
}

void ObsidianGL::Color3f(float p_r, float p_g, float p_b)
{
	CurrentRenderState.CurrentColor = glm::vec4(p_r, p_g, p_b, 1.0f);
}

void ObsidianGL::Color4f(float p_r, float p_g, float p_b, float p_a)
{
	CurrentRenderState.CurrentColor = glm::vec4(p_r, p_g, p_b, p_a);
}

void ObsidianGL::Normal3f(float p_x, float p_y, float p_z)
{
	CurrentRenderState.CurrentNormal = glm::vec3(p_x, p_y, p_z);
}

void ObsidianGL::TexCoord2f(float p_s, float p_t)
{
	CurrentRenderState.CurrentTexCoord = glm::vec2(p_s, p_t);
}

void ObsidianGL::MatrixMode(uint16_t p_mode)
{
	CurrentRenderState.CurrentMatrixMode = p_mode;
}

void ObsidianGL::LoadIdentity()
{
	GetCurrentStack().Top() = glm::mat4(1.0f);
}

void ObsidianGL::PushMatrix()
{
	GetCurrentStack().Push();
}

void ObsidianGL::PopMatrix()
{
	GetCurrentStack().Pop();
}

void ObsidianGL::LoadMatrixf(const float* p_matrix)
{
	GetCurrentStack().Top() = glm::make_mat4(p_matrix);
}

void ObsidianGL::MultMatrixf(const float* p_matrix)
{
	auto& stack = GetCurrentStack();
	stack.Top() = stack.Top() * glm::make_mat4(p_matrix);
}

void ObsidianGL::Translatef(float p_x, float p_y, float p_z)
{
	auto& stack = GetCurrentStack();
	stack.Top() = glm::translate(stack.Top(), glm::vec3(p_x, p_y, p_z));
}

void ObsidianGL::Rotatef(float p_angle, float p_x, float p_y, float p_z)
{
	auto& stack = GetCurrentStack();
	stack.Top() = glm::rotate(stack.Top(), glm::radians(p_angle), glm::vec3(p_x, p_y, p_z));
}

void ObsidianGL::Scalef(float p_x, float p_y, float p_z)
{
	auto& stack = GetCurrentStack();
	stack.Top() = glm::scale(stack.Top(), glm::vec3(p_x, p_y, p_z));
}

void ObsidianGL::Perspective(float p_fovY, float p_aspect, float p_near, float p_far)
{
	auto& stack = GetCurrentStack();
	stack.Top() = stack.Top() * glm::perspective(glm::radians(p_fovY), p_aspect, p_near, p_far);
}

void ObsidianGL::Ortho(double p_left, double p_right, double p_bottom, double p_top, double p_near, double p_far)
{
	auto& stack = GetCurrentStack();
	stack.Top() = stack.Top() * glm::ortho(static_cast<float>(p_left), static_cast<float>(p_right), static_cast<float>(p_bottom), static_cast<float>(p_top), static_cast<float>(p_near), static_cast<float>(p_far));
}

uint32_t* ObsidianGL::GetFrameBufferData()
{
	return CurrentRenderState.ColorBuffer ? CurrentRenderState.ColorBuffer->Data : nullptr;
}

uint32_t ObsidianGL::GetFrameBufferPitch()
{
	return CurrentRenderState.ColorBuffer ? CurrentRenderState.ColorBuffer->Width * sizeof(uint32_t) : 0;
}

uint32_t ObsidianGL::GetFrameBufferWidth()
{
	return CurrentRenderState.ColorBuffer ? CurrentRenderState.ColorBuffer->Width : 0;
}

uint32_t ObsidianGL::GetFrameBufferHeight()
{
	return CurrentRenderState.ColorBuffer ? CurrentRenderState.ColorBuffer->Height : 0;
}

uint32_t ObsidianGL::GenTexture()
{
	return TextureManager.Generate();
}

void ObsidianGL::DeleteTexture(uint32_t p_textureId)
{
	if (CurrentRenderState.BoundTextureId == p_textureId)
	{
		CurrentRenderState.BoundTextureId = 0;
	}

	TextureManager.Delete(p_textureId);
}

void ObsidianGL::BindTexture(uint32_t p_textureId)
{
	CurrentRenderState.BoundTextureId = p_textureId;
}

void ObsidianGL::TexImage2D(uint16_t p_internalFormat, uint32_t p_width, uint32_t p_height, uint16_t p_format, uint16_t p_type, const void* p_data)
{
	(void)p_internalFormat;
	(void)p_type;

	ObsidianGL::Graphics::Texture* tex = TextureManager.Get(CurrentRenderState.BoundTextureId);

	if (!tex || !p_data)
		return;

	tex->Upload(p_width, p_height, p_format, static_cast<const uint8_t*>(p_data));
}

void ObsidianGL::TexParameteri(uint16_t p_param, uint16_t p_value)
{
	ObsidianGL::Graphics::Texture* tex = TextureManager.Get(CurrentRenderState.BoundTextureId);

	if (!tex)
		return;

	switch (p_param)
	{
	case OGL_TEXTURE_MAG_FILTER:
		tex->MagFilter = p_value;
		break;
	case OGL_TEXTURE_MIN_FILTER:
		tex->MinFilter = p_value;
		break;
	case OGL_TEXTURE_WRAP_S:
		tex->WrapS = p_value;
		break;
	case OGL_TEXTURE_WRAP_T:
		tex->WrapT = p_value;
		break;
	}
}