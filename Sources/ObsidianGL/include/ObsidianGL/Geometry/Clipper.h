#pragma once

#include <cstdint>
#include <cmath>

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#ifdef __AVX2__
#include <immintrin.h>
#else
#include <emmintrin.h>
#endif

namespace ObsidianGL::Geometry
{
	constexpr uint8_t MAX_CLIP_VERTICES = 12;

	struct ClipVertex
	{
		alignas(16) glm::vec4 Position;
		alignas(16) glm::vec4 Color;
		alignas(16) glm::vec3 Normal;
	};

	struct ClipPolygon
	{
		alignas(32) ClipVertex Vertices[MAX_CLIP_VERTICES];
		uint8_t VertexCount = 0;
	};

	enum ClipPlane : uint8_t
	{
		CLIP_PLANE_NEAR   = 0,
		CLIP_PLANE_FAR    = 1,
		CLIP_PLANE_LEFT   = 2,
		CLIP_PLANE_RIGHT  = 3,
		CLIP_PLANE_BOTTOM = 4,
		CLIP_PLANE_TOP    = 5,
		CLIP_PLANE_COUNT  = 6
	};

	class Clipper
	{
	public:
		// Outcode bits for Cohen-Sutherland style trivial rejection
		static constexpr uint8_t OUTCODE_INSIDE = 0x00;
		static constexpr uint8_t OUTCODE_NEAR   = 0x01;
		static constexpr uint8_t OUTCODE_FAR    = 0x02;
		static constexpr uint8_t OUTCODE_LEFT   = 0x04;
		static constexpr uint8_t OUTCODE_RIGHT  = 0x08;
		static constexpr uint8_t OUTCODE_BOTTOM = 0x10;
		static constexpr uint8_t OUTCODE_TOP    = 0x20;

		static inline uint8_t ComputeOutcode(const glm::vec4& p_vertex)
		{
			uint8_t code = OUTCODE_INSIDE;

			if (p_vertex.z < -p_vertex.w) code |= OUTCODE_NEAR;
			if (p_vertex.z >  p_vertex.w) code |= OUTCODE_FAR;
			if (p_vertex.x < -p_vertex.w) code |= OUTCODE_LEFT;
			if (p_vertex.x >  p_vertex.w) code |= OUTCODE_RIGHT;
			if (p_vertex.y < -p_vertex.w) code |= OUTCODE_BOTTOM;
			if (p_vertex.y >  p_vertex.w) code |= OUTCODE_TOP;

			return code;
		}

#ifdef __AVX2__
		// AVX2: Compute outcodes for 2 vertices at once (8 floats per vertex = 2 vec4s)
		static inline void ComputeOutcodes2_AVX2(const glm::vec4& p_vertex0, const glm::vec4& p_vertex1, uint8_t& outcode0, uint8_t& outcode1)
		{
			__m256 vertices = _mm256_set_ps(p_vertex1.w, p_vertex1.z, p_vertex1.y, p_vertex1.x, p_vertex0.w, p_vertex0.z, p_vertex0.y, p_vertex0.x);

			__m256 x = _mm256_shuffle_ps(vertices, vertices, _MM_SHUFFLE(0, 0, 0, 0));
			__m256 y = _mm256_shuffle_ps(vertices, vertices, _MM_SHUFFLE(1, 1, 1, 1));
			__m256 z = _mm256_shuffle_ps(vertices, vertices, _MM_SHUFFLE(2, 2, 2, 2));
			__m256 w = _mm256_shuffle_ps(vertices, vertices, _MM_SHUFFLE(3, 3, 3, 3));
			__m256 negW = _mm256_sub_ps(_mm256_setzero_ps(), w);

			__m256 nearTest   = _mm256_cmp_ps(z, negW, _CMP_LT_OQ); // z < -w
			__m256 farTest    = _mm256_cmp_ps(z, w, _CMP_GT_OQ);    // z > w
			__m256 leftTest   = _mm256_cmp_ps(x, negW, _CMP_LT_OQ); // x < -w
			__m256 rightTest  = _mm256_cmp_ps(x, w, _CMP_GT_OQ);    // x > w
			__m256 bottomTest = _mm256_cmp_ps(y, negW, _CMP_LT_OQ); // y < -w
			__m256 topTest    = _mm256_cmp_ps(y, w, _CMP_GT_OQ);    // y > w

			int nearMask   = _mm256_movemask_ps(nearTest);
			int farMask    = _mm256_movemask_ps(farTest);
			int leftMask   = _mm256_movemask_ps(leftTest);
			int rightMask  = _mm256_movemask_ps(rightTest);
			int bottomMask = _mm256_movemask_ps(bottomTest);
			int topMask    = _mm256_movemask_ps(topTest);

			outcode0 = ((nearMask & 1) ? OUTCODE_NEAR : 0) |
			           ((farMask & 1) ? OUTCODE_FAR : 0) |
			           ((leftMask & 1) ? OUTCODE_LEFT : 0) |
			           ((rightMask & 1) ? OUTCODE_RIGHT : 0) |
			           ((bottomMask & 1) ? OUTCODE_BOTTOM : 0) |
			           ((topMask & 1) ? OUTCODE_TOP : 0);

			outcode1 = ((nearMask & 0x10) ? OUTCODE_NEAR : 0) |
			           ((farMask & 0x10) ? OUTCODE_FAR : 0) |
			           ((leftMask & 0x10) ? OUTCODE_LEFT : 0) |
			           ((rightMask & 0x10) ? OUTCODE_RIGHT : 0) |
			           ((bottomMask & 0x10) ? OUTCODE_BOTTOM : 0) |
			           ((topMask & 0x10) ? OUTCODE_TOP : 0);
		}
#endif

		static inline float ComputeClipDistance(const glm::vec4& p_vertex, ClipPlane p_plane)
		{
			switch (p_plane)
			{
				case CLIP_PLANE_NEAR:   return p_vertex.z + p_vertex.w; // z >= -w => z + w >= 0
				case CLIP_PLANE_FAR:    return p_vertex.w - p_vertex.z; // z <= w  => w - z >= 0
				case CLIP_PLANE_LEFT:   return p_vertex.x + p_vertex.w; // x >= -w => x + w >= 0
				case CLIP_PLANE_RIGHT:  return p_vertex.w - p_vertex.x; // x <= w  => w - x >= 0
				case CLIP_PLANE_BOTTOM: return p_vertex.y + p_vertex.w; // y >= -w => y + w >= 0
				case CLIP_PLANE_TOP:    return p_vertex.w - p_vertex.y; // y <= w  => w - y >= 0
				default: return 0.0f;
			}
		}

#ifdef __AVX2__
		// AVX2: Compute clip distances for 8 vertices against one plane
		static inline __m256 ComputeClipDistances8_AVX2(const float* p_positionX, const float* p_positionY, const float* p_positionZ, const float* p_positionW, ClipPlane p_plane)
		{
			__m256 x = _mm256_loadu_ps(p_positionX);
			__m256 y = _mm256_loadu_ps(p_positionY);
			__m256 z = _mm256_loadu_ps(p_positionZ);
			__m256 w = _mm256_loadu_ps(p_positionW);

			switch (p_plane)
			{
				case CLIP_PLANE_NEAR:   return _mm256_add_ps(z, w);
				case CLIP_PLANE_FAR:    return _mm256_sub_ps(w, z);
				case CLIP_PLANE_LEFT:   return _mm256_add_ps(x, w);
				case CLIP_PLANE_RIGHT:  return _mm256_sub_ps(w, x);
				case CLIP_PLANE_BOTTOM: return _mm256_add_ps(y, w);
				case CLIP_PLANE_TOP:    return _mm256_sub_ps(w, y);
				default: return _mm256_setzero_ps();
			}
		}
#else
		// SSE: Compute clip distances for 4 vertices against one plane
		static inline __m128 ComputeClipDistances4_SSE(const float* p_positionX, const float* p_positionY, const float* p_positionZ, const float* p_positionW, ClipPlane p_plane)
		{
			__m128 x = _mm_loadu_ps(p_positionX);
			__m128 y = _mm_loadu_ps(p_positionY);
			__m128 z = _mm_loadu_ps(p_positionZ);
			__m128 w = _mm_loadu_ps(p_positionW);

			switch (p_plane)
			{
				case CLIP_PLANE_NEAR:   return _mm_add_ps(z, w);
				case CLIP_PLANE_FAR:    return _mm_sub_ps(w, z);
				case CLIP_PLANE_LEFT:   return _mm_add_ps(x, w);
				case CLIP_PLANE_RIGHT:  return _mm_sub_ps(w, x);
				case CLIP_PLANE_BOTTOM: return _mm_add_ps(y, w);
				case CLIP_PLANE_TOP:    return _mm_sub_ps(w, y);
				default: return _mm_setzero_ps();
			}
		}
#endif

		static inline ClipVertex InterpolateVertex(const ClipVertex& p_clipVertex0, const ClipVertex& p_clipVertex1, float p_t)
		{
			ClipVertex result;

			result.Position.x = p_clipVertex0.Position.x + p_t * (p_clipVertex1.Position.x - p_clipVertex0.Position.x);
			result.Position.y = p_clipVertex0.Position.y + p_t * (p_clipVertex1.Position.y - p_clipVertex0.Position.y);
			result.Position.z = p_clipVertex0.Position.z + p_t * (p_clipVertex1.Position.z - p_clipVertex0.Position.z);
			result.Position.w = p_clipVertex0.Position.w + p_t * (p_clipVertex1.Position.w - p_clipVertex0.Position.w);

			result.Color.r = p_clipVertex0.Color.r + p_t * (p_clipVertex1.Color.r - p_clipVertex0.Color.r);
			result.Color.g = p_clipVertex0.Color.g + p_t * (p_clipVertex1.Color.g - p_clipVertex0.Color.g);
			result.Color.b = p_clipVertex0.Color.b + p_t * (p_clipVertex1.Color.b - p_clipVertex0.Color.b);
			result.Color.a = p_clipVertex0.Color.a + p_t * (p_clipVertex1.Color.a - p_clipVertex0.Color.a);

			result.Normal.x = p_clipVertex0.Normal.x + p_t * (p_clipVertex1.Normal.x - p_clipVertex0.Normal.x);
			result.Normal.y = p_clipVertex0.Normal.y + p_t * (p_clipVertex1.Normal.y - p_clipVertex0.Normal.y);
			result.Normal.z = p_clipVertex0.Normal.z + p_t * (p_clipVertex1.Normal.z - p_clipVertex0.Normal.z);

			return result;
		}

#ifdef __AVX2__
		// AVX2: Interpolate 4 floats at once
		static inline __m128 InterpolateVec4_AVX2(const __m128& p_a, const __m128& p_b, float p_t)
		{
			__m128 vt = _mm_set1_ps(p_t);
			__m128 delta = _mm_sub_ps(p_b, p_a);
			return _mm_add_ps(p_a, _mm_mul_ps(delta, vt));
		}

		static inline ClipVertex InterpolateVertex_AVX2(const ClipVertex& p_clipVertex0, const ClipVertex& p_clipVertex1, float t)
		{
			ClipVertex result;

			__m128 position0 = _mm_loadu_ps(&p_clipVertex0.Position.x);
			__m128 position1 = _mm_loadu_ps(&p_clipVertex1.Position.x);
			__m128 positionResult = InterpolateVec4_AVX2(position0, position1, t);
			_mm_storeu_ps(&result.Position.x, positionResult);

			__m128 color0 = _mm_loadu_ps(&p_clipVertex0.Color.r);
			__m128 color1 = _mm_loadu_ps(&p_clipVertex1.Color.r);
			__m128 colorResult = InterpolateVec4_AVX2(color0, color1, t);
			_mm_storeu_ps(&result.Color.r, colorResult);

			__m128 normal0 = _mm_set_ps(0.0f, p_clipVertex0.Normal.z, p_clipVertex0.Normal.y, p_clipVertex0.Normal.x);
			__m128 normal1 = _mm_set_ps(0.0f, p_clipVertex1.Normal.z, p_clipVertex1.Normal.y, p_clipVertex1.Normal.x);
			__m128 normalResult = InterpolateVec4_AVX2(normal0, normal1, t);

			alignas(16) float normalData[4];
			_mm_store_ps(normalData, normalResult);
			result.Normal.x = normalData[0];
			result.Normal.y = normalData[1];
			result.Normal.z = normalData[2];

			return result;
		}
#else
		// SSE: Interpolate vec4
		static inline __m128 InterpolateVec4_SSE(const __m128& p_a, const __m128& p_b, float p_t)
		{
			__m128 vt = _mm_set1_ps(p_t);
			__m128 delta = _mm_sub_ps(p_b, p_a);
			return _mm_add_ps(p_a, _mm_mul_ps(delta, vt));
		}

		static inline ClipVertex InterpolateVertex_SSE(const ClipVertex& p_clipVertex0, const ClipVertex& p_clipVertex1, float p_t)
		{
			ClipVertex result;

			__m128 position0 = _mm_loadu_ps(&p_clipVertex0.Position.x);
			__m128 position1 = _mm_loadu_ps(&p_clipVertex1.Position.x);
			__m128 positionResult = InterpolateVec4_SSE(position0, position1, p_t);
			_mm_storeu_ps(&result.Position.x, positionResult);

			__m128 color0 = _mm_loadu_ps(&p_clipVertex0.Color.r);
			__m128 color1 = _mm_loadu_ps(&p_clipVertex1.Color.r);
			__m128 colorResult = InterpolateVec4_SSE(color0, color1, p_t);
			_mm_storeu_ps(&result.Color.r, colorResult);

			__m128 normal0 = _mm_set_ps(0.0f, p_clipVertex0.Normal.z, p_clipVertex0.Normal.y, p_clipVertex0.Normal.x);
			__m128 normal1 = _mm_set_ps(0.0f, p_clipVertex1.Normal.z, p_clipVertex1.Normal.y, p_clipVertex1.Normal.x);
			__m128 normalResult = InterpolateVec4_SSE(normal0, normal1, p_t);

			alignas(16) float normalData[4];
			_mm_store_ps(normalData, normalResult);
			result.Normal.x = normalData[0];
			result.Normal.y = normalData[1];
			result.Normal.z = normalData[2];

			return result;
		}
#endif

		// Clip polygon against a single plane using Sutherland-Hodgman algorithm
		static inline void ClipPolygonAgainstPlane(ClipPolygon& p_clipPolygon, ClipPlane p_clipPlane)
		{
			if (p_clipPolygon.VertexCount == 0)
				return;

			ClipVertex outputVertices[MAX_CLIP_VERTICES];
			uint8_t outputCount = 0;

			const ClipVertex* previousClipVertex = &p_clipPolygon.Vertices[p_clipPolygon.VertexCount - 1];
			float previousDistance = ComputeClipDistance(previousClipVertex->Position, p_clipPlane);

			for (uint8_t i = 0; i < p_clipPolygon.VertexCount; ++i)
			{
				const ClipVertex* currentVertex = &p_clipPolygon.Vertices[i];

				float currentDistance = ComputeClipDistance(currentVertex->Position, p_clipPlane);

				if ((previousDistance >= 0.0f) != (currentDistance >= 0.0f))
				{
					float t = previousDistance / (previousDistance - currentDistance);

#ifdef __AVX2__
					outputVertices[outputCount++] = InterpolateVertex_AVX2(*previousClipVertex, *currentVertex, t);
#else
					outputVertices[outputCount++] = InterpolateVertex_SSE(*previousClipVertex, *currentVertex, t);
#endif
				}

				if (currentDistance >= 0.0f)
				{
					outputVertices[outputCount++] = *currentVertex;
				}

				previousClipVertex = currentVertex;
				previousDistance = currentDistance;
			}

			p_clipPolygon.VertexCount = outputCount;

			for (uint8_t i = 0; i < outputCount; ++i)
			{
				p_clipPolygon.Vertices[i] = outputVertices[i];
			}
		}

		static inline bool ClipTriangle(const ClipVertex& p_clipVertex0, const ClipVertex& p_clipVertex1, const ClipVertex& p_clipVertex2, ClipPolygon& p_outClipPolygon)
		{
			uint8_t code0 = ComputeOutcode(p_clipVertex0.Position);
			uint8_t code1 = ComputeOutcode(p_clipVertex1.Position);
			uint8_t code2 = ComputeOutcode(p_clipVertex2.Position);

			if ((code0 | code1 | code2) == OUTCODE_INSIDE)
			{
				p_outClipPolygon.Vertices[0] = p_clipVertex0;
				p_outClipPolygon.Vertices[1] = p_clipVertex1;
				p_outClipPolygon.Vertices[2] = p_clipVertex2;
				p_outClipPolygon.VertexCount = 3;
				return true;
			}

			if ((code0 & code1 & code2) != OUTCODE_INSIDE)
			{
				p_outClipPolygon.VertexCount = 0;
				return false;
			}

			p_outClipPolygon.Vertices[0] = p_clipVertex0;
			p_outClipPolygon.Vertices[1] = p_clipVertex1;
			p_outClipPolygon.Vertices[2] = p_clipVertex2;
			p_outClipPolygon.VertexCount = 3;

			uint8_t combinedCode = code0 | code1 | code2;

			if (combinedCode & OUTCODE_NEAR)
				ClipPolygonAgainstPlane(p_outClipPolygon, CLIP_PLANE_NEAR);

			if (p_outClipPolygon.VertexCount == 0) 
				return false;

			if (combinedCode & OUTCODE_FAR)
				ClipPolygonAgainstPlane(p_outClipPolygon, CLIP_PLANE_FAR);

			if (p_outClipPolygon.VertexCount == 0) 
				return false;

			if (combinedCode & OUTCODE_LEFT)
				ClipPolygonAgainstPlane(p_outClipPolygon, CLIP_PLANE_LEFT);

			if (p_outClipPolygon.VertexCount == 0) 
				return false;

			if (combinedCode & OUTCODE_RIGHT)
				ClipPolygonAgainstPlane(p_outClipPolygon, CLIP_PLANE_RIGHT);

			if (p_outClipPolygon.VertexCount == 0) 
				return false;

			if (combinedCode & OUTCODE_BOTTOM)
				ClipPolygonAgainstPlane(p_outClipPolygon, CLIP_PLANE_BOTTOM);

			if (p_outClipPolygon.VertexCount == 0) 
				return false;

			if (combinedCode & OUTCODE_TOP)
				ClipPolygonAgainstPlane(p_outClipPolygon, CLIP_PLANE_TOP);

			return p_outClipPolygon.VertexCount >= 3;
		}
	};
}
