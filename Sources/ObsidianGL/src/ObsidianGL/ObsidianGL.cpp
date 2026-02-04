#include "ObsidianGL/ObsidianGL.h"

#include <algorithm>
#include <vector>

#include <glm/gtc/type_ptr.hpp>

#include "ObsidianGL/Geometry/Triangle.h"
#include "ObsidianGL/Geometry/Clipper.h"
#include "ObsidianGL/Graphics/Color.h"
#include "ObsidianGL/Buffers/FrameBuffer.h"

// ============================================================================
// #define OBSIDIAN_FORCE_SCALAR // Pure scalar - no SIMD
// #define OBSIDIAN_FORCE_SSE    // SSE only (4 pixels at a time)
// #define OBSIDIAN_FORCE_AVX2   // AVX2 (8 pixels at a time) - requires AVX2 CPU
// ============================================================================
// If none are defined, it auto-selects based on compiler flags:
// - __AVX2__ defined -> AVX2
// - Otherwise -> SSE
// ============================================================================

#if defined(OBSIDIAN_FORCE_SCALAR)
	// No SIMD includes needed for scalar
#elif defined(OBSIDIAN_FORCE_SSE)
#include <emmintrin.h>
#elif defined(OBSIDIAN_FORCE_AVX2)
#include <immintrin.h>
#elif defined(__AVX2__)
#include <immintrin.h>
#define OBSIDIAN_AUTO_AVX2
#else
#include <emmintrin.h>
#define OBSIDIAN_AUTO_SSE
#endif

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

		MatrixStack ModelViewStack;
		MatrixStack ProjectionStack;
		uint16_t CurrentMatrixMode = ObsidianGL::OGL_MODELVIEW;

		bool InBeginEnd = false;
		uint16_t CurrentPrimitiveMode = 0;
		glm::vec4 CurrentColor = glm::vec4(1.0f);
		glm::vec3 CurrentNormal = glm::vec3(0.0f, 0.0f, 1.0f);
	};

	struct ImmediateVertex
	{
		glm::vec3 Position;
		glm::vec4 Color;
		glm::vec3 Normal;
	};

	RenderState CurrentRenderState;
	std::vector<ImmediateVertex> VertexBuffer;

	MatrixStack& GetCurrentStack()
	{
		return (CurrentRenderState.CurrentMatrixMode == ObsidianGL::OGL_PROJECTION) ? CurrentRenderState.ProjectionStack : CurrentRenderState.ModelViewStack;
	}

	void RasterizeTriangle3D_Scalar(const glm::vec4& p_vertex0, const glm::vec4& p_vertex1, const glm::vec4& p_vertex2, const glm::vec4& p_color0, const glm::vec4& p_color1, const glm::vec4& p_color2)
	{
		if (!CurrentRenderState.ColorBuffer)
			return;

		const uint32_t width = CurrentRenderState.ColorBuffer->Width;
		const uint32_t height = CurrentRenderState.ColorBuffer->Height;

		const float invW0 = 1.0f / p_vertex0.w;
		const float invW1 = 1.0f / p_vertex1.w;
		const float invW2 = 1.0f / p_vertex2.w;

		const glm::vec3 ndc0 = glm::vec3(p_vertex0) * invW0;
		const glm::vec3 ndc1 = glm::vec3(p_vertex1) * invW1;
		const glm::vec3 ndc2 = glm::vec3(p_vertex2) * invW2;

		const float halfW = static_cast<float>(width) * 0.5f;
		const float halfH = static_cast<float>(height) * 0.5f;

		const glm::vec2 screen0((ndc0.x + 1.0f) * halfW, (1.0f - ndc0.y) * halfH);
		const glm::vec2 screen1((ndc1.x + 1.0f) * halfW, (1.0f - ndc1.y) * halfH);
		const glm::vec2 screen2((ndc2.x + 1.0f) * halfW, (1.0f - ndc2.y) * halfH);

		const float area = (screen1.x - screen0.x) * (screen2.y - screen0.y) - (screen1.y - screen0.y) * (screen2.x - screen0.x);

		if (std::abs(area) < 0.5f)
			return;

		if (CurrentRenderState.State & ObsidianGL::OGL_CULL_FACE)
		{
			if ((CurrentRenderState.CullFace == ObsidianGL::OGL_BACK && area < 0.0f)
				|| (CurrentRenderState.CullFace == ObsidianGL::OGL_FRONT && area > 0.0f))
				return;
		}

		ObsidianGL::Geometry::Triangle triangle(screen0, screen1, screen2);

		if (triangle.InvDenom == 0.0f)
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
		const glm::vec4 c1w = p_color1 * invW1;
		const glm::vec4 c2w = p_color2 * invW2;

		const bool depthTest = (CurrentRenderState.State & ObsidianGL::OGL_DEPTH_TEST) && CurrentRenderState.DepthBuffer;
		const bool depthWrite = (CurrentRenderState.State & ObsidianGL::OGL_DEPTH_WRITE) && CurrentRenderState.DepthBuffer;

		for (int32_t y = minY; y <= maxY; ++y)
		{
			const float py = static_cast<float>(y) + 0.5f;

			RGBA8* colorRow = CurrentRenderState.ColorBuffer->GetRow(static_cast<uint32_t>(y));
			Depth* depthRow = (depthTest || depthWrite) ? CurrentRenderState.DepthBuffer->GetRow(static_cast<uint32_t>(y)) : nullptr;

			for (int32_t x = minX; x <= maxX; ++x)
			{
				const float pxf = static_cast<float>(x) + 0.5f;

				const glm::vec3 bary = triangle.GetBarycentricCoordinates(pxf, py);

				if (bary.x >= 0.0f && bary.y >= 0.0f && bary.z >= 0.0f)
				{
					const float w = bary.x;
					const float u = bary.y;
					const float v = bary.z;

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
					const glm::vec4 color = colorInterpolated * correctedW;

					colorRow[x] = ObsidianGL::Graphics::PackColor(color);

					if (depthWrite && depthRow)
					{
						depthRow[x] = depth;
					}
				}
			}
		}
	}

#if defined(OBSIDIAN_FORCE_SSE) || defined(OBSIDIAN_AUTO_SSE)
	void RasterizeTriangle3D_SSE(const glm::vec4& p_vertex0, const glm::vec4& p_vertex1, const glm::vec4& p_vertex2, const glm::vec4& p_color0, const glm::vec4& p_color1, const glm::vec4& p_color2)
	{
		if (!CurrentRenderState.ColorBuffer)
			return;

		const uint32_t width = CurrentRenderState.ColorBuffer->Width;
		const uint32_t height = CurrentRenderState.ColorBuffer->Height;

		const float invW0 = 1.0f / p_vertex0.w;
		const float invW1 = 1.0f / p_vertex1.w;
		const float invW2 = 1.0f / p_vertex2.w;

		const glm::vec3 ndc0 = glm::vec3(p_vertex0) * invW0;
		const glm::vec3 ndc1 = glm::vec3(p_vertex1) * invW1;
		const glm::vec3 ndc2 = glm::vec3(p_vertex2) * invW2;

		const float halfW = static_cast<float>(width) * 0.5f;
		const float halfH = static_cast<float>(height) * 0.5f;

		const glm::vec2 screen0((ndc0.x + 1.0f) * halfW, (1.0f - ndc0.y) * halfH);
		const glm::vec2 screen1((ndc1.x + 1.0f) * halfW, (1.0f - ndc1.y) * halfH);
		const glm::vec2 screen2((ndc2.x + 1.0f) * halfW, (1.0f - ndc2.y) * halfH);

		const float area = (screen1.x - screen0.x) * (screen2.y - screen0.y) - (screen1.y - screen0.y) * (screen2.x - screen0.x);

		if (std::abs(area) < 0.5f)
			return;

		if (CurrentRenderState.State & ObsidianGL::OGL_CULL_FACE)
		{
			if ((CurrentRenderState.CullFace == ObsidianGL::OGL_BACK && area < 0.0f)
				|| (CurrentRenderState.CullFace == ObsidianGL::OGL_FRONT && area > 0.0f))
				return;
		}

		ObsidianGL::Geometry::Triangle triangle(screen0, screen1, screen2);

		if (triangle.InvDenom == 0.0f)
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
		const glm::vec4 c1w = p_color1 * invW1;
		const glm::vec4 c2w = p_color2 * invW2;

		const bool depthTest = (CurrentRenderState.State & ObsidianGL::OGL_DEPTH_TEST) && CurrentRenderState.DepthBuffer;
		const bool depthWrite = (CurrentRenderState.State & ObsidianGL::OGL_DEPTH_WRITE) && CurrentRenderState.DepthBuffer;

		const __m128 vInvW0 = _mm_set1_ps(invW0);
		const __m128 vInvW1 = _mm_set1_ps(invW1);
		const __m128 vInvW2 = _mm_set1_ps(invW2);
		const __m128 vDepth0 = _mm_set1_ps(depth0);
		const __m128 vDepth1 = _mm_set1_ps(depth1);
		const __m128 vDepth2 = _mm_set1_ps(depth2);
		const __m128 vZero = _mm_setzero_ps();
		const __m128 vOne = _mm_set1_ps(1.0f);
		const __m128 v255 = _mm_set1_ps(255.0f);

		const __m128 vC0R = _mm_set1_ps(c0w.r);
		const __m128 vC0G = _mm_set1_ps(c0w.g);
		const __m128 vC0B = _mm_set1_ps(c0w.b);
		const __m128 vC0A = _mm_set1_ps(c0w.a);
		const __m128 vC1R = _mm_set1_ps(c1w.r);
		const __m128 vC1G = _mm_set1_ps(c1w.g);
		const __m128 vC1B = _mm_set1_ps(c1w.b);
		const __m128 vC1A = _mm_set1_ps(c1w.a);
		const __m128 vC2R = _mm_set1_ps(c2w.r);
		const __m128 vC2G = _mm_set1_ps(c2w.g);
		const __m128 vC2B = _mm_set1_ps(c2w.b);
		const __m128 vC2A = _mm_set1_ps(c2w.a);

		alignas(16) float px[4];
		alignas(16) float baryU[4];
		alignas(16) float baryV[4];
		alignas(16) float baryW[4];

		for (int32_t y = minY; y <= maxY; ++y)
		{
			const float py = static_cast<float>(y) + 0.5f;

			RGBA8* colorRow = CurrentRenderState.ColorBuffer->GetRow(static_cast<uint32_t>(y));
			Depth* depthRow = (depthTest || depthWrite) ? CurrentRenderState.DepthBuffer->GetRow(static_cast<uint32_t>(y)) : nullptr;

			int32_t x = minX;

			for (; x + 3 <= maxX; x += 4)
			{
				px[0] = static_cast<float>(x + 0) + 0.5f;
				px[1] = static_cast<float>(x + 1) + 0.5f;
				px[2] = static_cast<float>(x + 2) + 0.5f;
				px[3] = static_cast<float>(x + 3) + 0.5f;

				triangle.GetBarycentricCoordinates4(px, py, baryU, baryV, baryW);

				__m128 vBaryW = _mm_load_ps(baryW);
				__m128 vBaryU = _mm_load_ps(baryU);
				__m128 vBaryV = _mm_load_ps(baryV);

				__m128 maskW = _mm_cmpge_ps(vBaryW, vZero);
				__m128 maskU = _mm_cmpge_ps(vBaryU, vZero);
				__m128 maskV = _mm_cmpge_ps(vBaryV, vZero);
				__m128 insideMask = _mm_and_ps(_mm_and_ps(maskW, maskU), maskV);

				int insideBits = _mm_movemask_ps(insideMask);

				if (insideBits == 0)
					continue;

				__m128 vOneOverW = _mm_add_ps(_mm_add_ps(_mm_mul_ps(vBaryW, vInvW0), _mm_mul_ps(vBaryU, vInvW1)), _mm_mul_ps(vBaryV, vInvW2));
				__m128 vCorrectedW = _mm_div_ps(vOne, vOneOverW);

				__m128 vDepth = _mm_add_ps(_mm_add_ps(_mm_mul_ps(vBaryW, vDepth0), _mm_mul_ps(vBaryU, vDepth1)), _mm_mul_ps(vBaryV, vDepth2));

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

				__m128 vR = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(vBaryW, vC0R), _mm_mul_ps(vBaryU, vC1R)), _mm_mul_ps(vBaryV, vC2R)), vCorrectedW);
				__m128 vG = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(vBaryW, vC0G), _mm_mul_ps(vBaryU, vC1G)), _mm_mul_ps(vBaryV, vC2G)), vCorrectedW);
				__m128 vB = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(vBaryW, vC0B), _mm_mul_ps(vBaryU, vC1B)), _mm_mul_ps(vBaryV, vC2B)), vCorrectedW);
				__m128 vA = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(vBaryW, vC0A), _mm_mul_ps(vBaryU, vC1A)), _mm_mul_ps(vBaryV, vC2A)), vCorrectedW);

				vR = _mm_mul_ps(_mm_min_ps(_mm_max_ps(vR, vZero), vOne), v255);
				vG = _mm_mul_ps(_mm_min_ps(_mm_max_ps(vG, vZero), vOne), v255);
				vB = _mm_mul_ps(_mm_min_ps(_mm_max_ps(vB, vZero), vOne), v255);
				vA = _mm_mul_ps(_mm_min_ps(_mm_max_ps(vA, vZero), vOne), v255);

				__m128i vRi = _mm_cvtps_epi32(vR);
				__m128i vGi = _mm_cvtps_epi32(vG);
				__m128i vBi = _mm_cvtps_epi32(vB);
				__m128i vAi = _mm_cvtps_epi32(vA);

				__m128i vColor = _mm_or_si128(_mm_or_si128(_mm_slli_epi32(vRi, 24), _mm_slli_epi32(vGi, 16)), _mm_or_si128(_mm_slli_epi32(vBi, 8), vAi));

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
			}

			// Scalar fallback
			for (; x <= maxX; ++x)
			{
				const float pxf = static_cast<float>(x) + 0.5f;
				const glm::vec3 bary = triangle.GetBarycentricCoordinates(pxf, py);

				if (bary.x >= 0.0f && bary.y >= 0.0f && bary.z >= 0.0f)
				{
					const float w = bary.x, u = bary.y, v = bary.z;
					const float oneOverW = w * invW0 + u * invW1 + v * invW2;
					const float correctedW = 1.0f / oneOverW;
					const float depth = w * depth0 + u * depth1 + v * depth2;

					if (depth < 0.0f || depth > 1.0f)
						continue;

					if (depthTest && depthRow && depth >= depthRow[x])
						continue;

					const glm::vec4 colorInterpolated = w * c0w + u * c1w + v * c2w;
					const glm::vec4 color = colorInterpolated * correctedW;

					colorRow[x] = ObsidianGL::Graphics::PackColor(color);

					if (depthWrite && depthRow)
					{
						depthRow[x] = depth;
					}
				}
			}
		}
	}
#endif

#if defined(OBSIDIAN_FORCE_AVX2) || defined(OBSIDIAN_AUTO_AVX2)
	void RasterizeTriangle3D_AVX2(const glm::vec4& p_vertex0, const glm::vec4& p_vertex1, const glm::vec4& p_vertex2, const glm::vec4& p_color0, const glm::vec4& p_color1, const glm::vec4& p_color2)
	{
		if (!CurrentRenderState.ColorBuffer)
			return;

		const uint32_t width = CurrentRenderState.ColorBuffer->Width;
		const uint32_t height = CurrentRenderState.ColorBuffer->Height;

		const float invW0 = 1.0f / p_vertex0.w;
		const float invW1 = 1.0f / p_vertex1.w;
		const float invW2 = 1.0f / p_vertex2.w;

		const glm::vec3 ndc0 = glm::vec3(p_vertex0) * invW0;
		const glm::vec3 ndc1 = glm::vec3(p_vertex1) * invW1;
		const glm::vec3 ndc2 = glm::vec3(p_vertex2) * invW2;

		const float halfW = static_cast<float>(width) * 0.5f;
		const float halfH = static_cast<float>(height) * 0.5f;

		const glm::vec2 screen0((ndc0.x + 1.0f) * halfW, (1.0f - ndc0.y) * halfH);
		const glm::vec2 screen1((ndc1.x + 1.0f) * halfW, (1.0f - ndc1.y) * halfH);
		const glm::vec2 screen2((ndc2.x + 1.0f) * halfW, (1.0f - ndc2.y) * halfH);

		const float area = (screen1.x - screen0.x) * (screen2.y - screen0.y) - (screen1.y - screen0.y) * (screen2.x - screen0.x);

		if (std::abs(area) < 0.5f)
			return;

		if (CurrentRenderState.State & ObsidianGL::OGL_CULL_FACE)
		{
			if ((CurrentRenderState.CullFace == ObsidianGL::OGL_BACK && area < 0.0f)
				|| (CurrentRenderState.CullFace == ObsidianGL::OGL_FRONT && area > 0.0f))
				return;
		}

		ObsidianGL::Geometry::Triangle triangle(screen0, screen1, screen2);

		if (triangle.InvDenom == 0.0f)
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
		const glm::vec4 c1w = p_color1 * invW1;
		const glm::vec4 c2w = p_color2 * invW2;

		const bool depthTest = (CurrentRenderState.State & ObsidianGL::OGL_DEPTH_TEST) && CurrentRenderState.DepthBuffer;
		const bool depthWrite = (CurrentRenderState.State & ObsidianGL::OGL_DEPTH_WRITE) && CurrentRenderState.DepthBuffer;

		const __m256 vInvW0 = _mm256_set1_ps(invW0);
		const __m256 vInvW1 = _mm256_set1_ps(invW1);
		const __m256 vInvW2 = _mm256_set1_ps(invW2);
		const __m256 vDepth0 = _mm256_set1_ps(depth0);
		const __m256 vDepth1 = _mm256_set1_ps(depth1);
		const __m256 vDepth2 = _mm256_set1_ps(depth2);
		const __m256 vZero = _mm256_setzero_ps();
		const __m256 vOne = _mm256_set1_ps(1.0f);
		const __m256 v255 = _mm256_set1_ps(255.0f);

		const __m256 vC0R = _mm256_set1_ps(c0w.r);
		const __m256 vC0G = _mm256_set1_ps(c0w.g);
		const __m256 vC0B = _mm256_set1_ps(c0w.b);
		const __m256 vC0A = _mm256_set1_ps(c0w.a);
		const __m256 vC1R = _mm256_set1_ps(c1w.r);
		const __m256 vC1G = _mm256_set1_ps(c1w.g);
		const __m256 vC1B = _mm256_set1_ps(c1w.b);
		const __m256 vC1A = _mm256_set1_ps(c1w.a);
		const __m256 vC2R = _mm256_set1_ps(c2w.r);
		const __m256 vC2G = _mm256_set1_ps(c2w.g);
		const __m256 vC2B = _mm256_set1_ps(c2w.b);
		const __m256 vC2A = _mm256_set1_ps(c2w.a);

		alignas(32) float px[8];
		alignas(32) float baryU[8];
		alignas(32) float baryV[8];
		alignas(32) float baryW[8];

		for (int32_t y = minY; y <= maxY; ++y)
		{
			const float py = static_cast<float>(y) + 0.5f;

			RGBA8* colorRow = CurrentRenderState.ColorBuffer->GetRow(static_cast<uint32_t>(y));

			Depth* depthRow = (depthTest || depthWrite) ? CurrentRenderState.DepthBuffer->GetRow(static_cast<uint32_t>(y)) : nullptr;

			int32_t x = minX;

			for (; x + 7 <= maxX; x += 8)
			{
				px[0] = static_cast<float>(x + 0) + 0.5f;
				px[1] = static_cast<float>(x + 1) + 0.5f;
				px[2] = static_cast<float>(x + 2) + 0.5f;
				px[3] = static_cast<float>(x + 3) + 0.5f;
				px[4] = static_cast<float>(x + 4) + 0.5f;
				px[5] = static_cast<float>(x + 5) + 0.5f;
				px[6] = static_cast<float>(x + 6) + 0.5f;
				px[7] = static_cast<float>(x + 7) + 0.5f;

				triangle.GetBarycentricCoordinates8(px, py, baryU, baryV, baryW);

				__m256 vBaryW = _mm256_load_ps(baryW);
				__m256 vBaryU = _mm256_load_ps(baryU);
				__m256 vBaryV = _mm256_load_ps(baryV);

				__m256 maskW = _mm256_cmp_ps(vBaryW, vZero, _CMP_GE_OQ);
				__m256 maskU = _mm256_cmp_ps(vBaryU, vZero, _CMP_GE_OQ);
				__m256 maskV = _mm256_cmp_ps(vBaryV, vZero, _CMP_GE_OQ);

				__m256 insideMask = _mm256_and_ps(_mm256_and_ps(maskW, maskU), maskV);

				int insideBits = _mm256_movemask_ps(insideMask);

				if (insideBits == 0)
					continue;

				__m256 vOneOverW = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(vBaryW, vInvW0), _mm256_mul_ps(vBaryU, vInvW1)), _mm256_mul_ps(vBaryV, vInvW2));
				__m256 vCorrectedW = _mm256_div_ps(vOne, vOneOverW);

				__m256 vDepth = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(vBaryW, vDepth0), _mm256_mul_ps(vBaryU, vDepth1)), _mm256_mul_ps(vBaryV, vDepth2));

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

				__m256 vR = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(vBaryW, vC0R), _mm256_mul_ps(vBaryU, vC1R)), _mm256_mul_ps(vBaryV, vC2R)), vCorrectedW);
				__m256 vG = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(vBaryW, vC0G), _mm256_mul_ps(vBaryU, vC1G)), _mm256_mul_ps(vBaryV, vC2G)), vCorrectedW);
				__m256 vB = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(vBaryW, vC0B), _mm256_mul_ps(vBaryU, vC1B)), _mm256_mul_ps(vBaryV, vC2B)), vCorrectedW);
				__m256 vA = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(vBaryW, vC0A), _mm256_mul_ps(vBaryU, vC1A)), _mm256_mul_ps(vBaryV, vC2A)), vCorrectedW);

				vR = _mm256_mul_ps(_mm256_min_ps(_mm256_max_ps(vR, vZero), vOne), v255);
				vG = _mm256_mul_ps(_mm256_min_ps(_mm256_max_ps(vG, vZero), vOne), v255);
				vB = _mm256_mul_ps(_mm256_min_ps(_mm256_max_ps(vB, vZero), vOne), v255);
				vA = _mm256_mul_ps(_mm256_min_ps(_mm256_max_ps(vA, vZero), vOne), v255);

				__m256i vRi = _mm256_cvtps_epi32(vR);
				__m256i vGi = _mm256_cvtps_epi32(vG);
				__m256i vBi = _mm256_cvtps_epi32(vB);
				__m256i vAi = _mm256_cvtps_epi32(vA);

				__m256i vColor = _mm256_or_si256(_mm256_or_si256(_mm256_slli_epi32(vRi, 24), _mm256_slli_epi32(vGi, 16)), _mm256_or_si256(_mm256_slli_epi32(vBi, 8), vAi));

				alignas(32) uint32_t colors[8];
				alignas(32) float depths[8];

				_mm256_store_si256(reinterpret_cast<__m256i*>(colors), vColor);

				_mm256_store_ps(depths, vDepth);

				for (int i = 0; i < 8; ++i)
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
			}

			// Scalar fallback
			for (; x <= maxX; ++x)
			{
				const float pxf = static_cast<float>(x) + 0.5f;

				const glm::vec3 bary = triangle.GetBarycentricCoordinates(pxf, py);

				if (bary.x >= 0.0f && bary.y >= 0.0f && bary.z >= 0.0f)
				{
					const float w = bary.x, u = bary.y, v = bary.z;
					const float oneOverW = w * invW0 + u * invW1 + v * invW2;
					const float correctedW = 1.0f / oneOverW;
					const float depth = w * depth0 + u * depth1 + v * depth2;

					if (depth < 0.0f || depth > 1.0f)
						continue;

					if (depthTest && depthRow && depth >= depthRow[x])
						continue;

					const glm::vec4 colorInterp = w * c0w + u * c1w + v * c2w;
					const glm::vec4 color = colorInterp * correctedW;

					colorRow[x] = ObsidianGL::Graphics::PackColor(color);

					if (depthWrite && depthRow)
					{
						depthRow[x] = depth;
					}
				}
			}
		}
	}
#endif

	inline void RasterizeTriangle3D(const glm::vec4& p_vertex0, const glm::vec4& p_vertex1, const glm::vec4& p_vertex2, const glm::vec4& p_color0, const glm::vec4& p_color1, const glm::vec4& p_color2)
	{
#if defined(OBSIDIAN_FORCE_SCALAR)
		RasterizeTriangle3D_Scalar(p_vertex0, p_vertex1, p_vertex2, p_color0, p_color1, p_color2);
#elif defined(OBSIDIAN_FORCE_SSE)
		RasterizeTriangle3D_SSE(p_vertex0, p_vertex1, p_vertex2, p_color0, p_color1, p_color2);
#elif defined(OBSIDIAN_FORCE_AVX2)
		RasterizeTriangle3D_AVX2(p_vertex0, p_vertex1, p_vertex2, p_color0, p_color1, p_color2);
#elif defined(OBSIDIAN_AUTO_AVX2)
		RasterizeTriangle3D_AVX2(p_vertex0, p_vertex1, p_vertex2, p_color0, p_color1, p_color2);
#elif defined(OBSIDIAN_AUTO_SSE)
		RasterizeTriangle3D_SSE(p_vertex0, p_vertex1, p_vertex2, p_color0, p_color1, p_color2);
#else
		RasterizeTriangle3D_Scalar(p_vertex0, p_vertex1, p_vertex2, p_color0, p_color1, p_color2);
#endif
	}

	void ClipAndRasterizeTriangle(const glm::vec4& p_clipPosition0, const glm::vec4& p_clipPosition1, const glm::vec4& p_clipPosition2, const glm::vec4& p_color0, const glm::vec4& p_color1, const glm::vec4& p_color2, const glm::vec3& p_normal0, const glm::vec3& p_normal1, const glm::vec3& p_normal2)
	{
		uint8_t code0 = ObsidianGL::Geometry::Clipper::ComputeOutcode(p_clipPosition0);
		uint8_t code1 = ObsidianGL::Geometry::Clipper::ComputeOutcode(p_clipPosition1);
		uint8_t code2 = ObsidianGL::Geometry::Clipper::ComputeOutcode(p_clipPosition2);

		if ((code0 | code1 | code2) == ObsidianGL::Geometry::Clipper::OUTCODE_INSIDE)
		{
			RasterizeTriangle3D(p_clipPosition0, p_clipPosition1, p_clipPosition2, p_color0, p_color1, p_color2);
			return;
		}

		if ((code0 & code1 & code2) != ObsidianGL::Geometry::Clipper::OUTCODE_INSIDE)
		{
			return;
		}

		ObsidianGL::Geometry::ClipVertex v0 { .Position = p_clipPosition0, .Color = p_color0, .Normal = p_normal0 };
		ObsidianGL::Geometry::ClipVertex v1 { .Position = p_clipPosition1, .Color = p_color1, .Normal = p_normal1 };
		ObsidianGL::Geometry::ClipVertex v2 { .Position = p_clipPosition2, .Color = p_color2, .Normal = p_normal2 };

		ObsidianGL::Geometry::ClipPolygon clippedPoly;
		if (!ObsidianGL::Geometry::Clipper::ClipTriangle(v0, v1, v2, clippedPoly))
		{
			return;
		}

		for (uint8_t i = 1; i + 1 < clippedPoly.VertexCount; ++i)
		{
			RasterizeTriangle3D(
				clippedPoly.Vertices[0].Position,
				clippedPoly.Vertices[i].Position,
				clippedPoly.Vertices[i + 1].Position,
				clippedPoly.Vertices[0].Color,
				clippedPoly.Vertices[i].Color,
				clippedPoly.Vertices[i + 1].Color
			);
		}
	}

	void FlushVertices()
	{
		if (VertexBuffer.empty()) 
			return;

		const glm::mat4 mvp = CurrentRenderState.ProjectionStack.Top() * CurrentRenderState.ModelViewStack.Top();

		std::vector<glm::vec4> clipPositions;

		clipPositions.reserve(VertexBuffer.size());

		for (const auto& immediateVertex : VertexBuffer)
		{
			clipPositions.push_back(mvp * glm::vec4(immediateVertex.Position, 1.0f));
		}

		auto ProcessTriangle = [&](size_t index0, size_t index1, size_t index2)
			{
				ClipAndRasterizeTriangle(
					clipPositions[index0], clipPositions[index1], clipPositions[index2],
					VertexBuffer[index0].Color, VertexBuffer[index1].Color, VertexBuffer[index2].Color,
					VertexBuffer[index0].Normal, VertexBuffer[index1].Normal, VertexBuffer[index2].Normal
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

void ObsidianGL::CullFace(uint16_t p_face)
{
	CurrentRenderState.CullFace = p_face;
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