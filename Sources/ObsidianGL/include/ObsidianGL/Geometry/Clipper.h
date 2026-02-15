#pragma once

#include <cstdint>
#include <cmath>

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include "ObsidianGL/API/SIMDConfig.h"

namespace ObsidianGL::Geometry
{
	constexpr uint8_t MAX_CLIP_VERTICES = 12;

	struct ClipVertex
	{
		alignas(32) glm::vec4 Position;
		alignas(32) glm::vec4 Color;
		alignas(32) glm::vec3 Normal;
		glm::vec2 TexCoord;
	};

	struct ClipPolygon
	{
		alignas(32) ClipVertex Vertices[MAX_CLIP_VERTICES];
		uint8_t VertexCount = 0;
	};

	enum ClipPlane : uint8_t
	{
		CLIP_PLANE_NEAR = 0,
		CLIP_PLANE_FAR = 1,
		CLIP_PLANE_LEFT = 2,
		CLIP_PLANE_RIGHT = 3,
		CLIP_PLANE_BOTTOM = 4,
		CLIP_PLANE_TOP = 5,
		CLIP_PLANE_COUNT = 6
	};

	class Clipper
	{
	public:
		// Outcode bits for Cohen-Sutherland style trivial rejection
		static constexpr uint8_t OUTCODE_INSIDE = 0x00;
		static constexpr uint8_t OUTCODE_NEAR = 0x01;
		static constexpr uint8_t OUTCODE_FAR = 0x02;
		static constexpr uint8_t OUTCODE_LEFT = 0x04;
		static constexpr uint8_t OUTCODE_RIGHT = 0x08;
		static constexpr uint8_t OUTCODE_BOTTOM = 0x10;
		static constexpr uint8_t OUTCODE_TOP = 0x20;

		static uint8_t ComputeOutcode(const glm::vec4& p_vertex)
		{
			uint8_t code = OUTCODE_INSIDE;

			if (p_vertex.z < -p_vertex.w) code |= OUTCODE_NEAR;
			if (p_vertex.z > p_vertex.w) code |= OUTCODE_FAR;
			if (p_vertex.x < -p_vertex.w) code |= OUTCODE_LEFT;
			if (p_vertex.x > p_vertex.w) code |= OUTCODE_RIGHT;
			if (p_vertex.y < -p_vertex.w) code |= OUTCODE_BOTTOM;
			if (p_vertex.y > p_vertex.w) code |= OUTCODE_TOP;

			return code;
		}

		static float ComputeClipDistance(const glm::vec4& p_vertex, ClipPlane p_plane)
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


#ifdef OBSIDIAN_USE_AVX2
		static __m128 InterpolateVec4_AVX2(const __m128& p_a, const __m128& p_b, float p_t)
		{
			__m128 vt = _mm_set1_ps(p_t);
			__m128 delta = _mm_sub_ps(p_b, p_a);
			return _mm_add_ps(p_a, _mm_mul_ps(delta, vt));
		}

		static ClipVertex InterpolateVertex_AVX2(const ClipVertex& p_clipVertex0, const ClipVertex& p_clipVertex1, float t)
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

			result.TexCoord.x = p_clipVertex0.TexCoord.x + t * (p_clipVertex1.TexCoord.x - p_clipVertex0.TexCoord.x);
			result.TexCoord.y = p_clipVertex0.TexCoord.y + t * (p_clipVertex1.TexCoord.y - p_clipVertex0.TexCoord.y);

			return result;
		}
#else
		static __m128 InterpolateVec4_SSE(const __m128& p_a, const __m128& p_b, float p_t)
		{
			__m128 vt = _mm_set1_ps(p_t);
			__m128 delta = _mm_sub_ps(p_b, p_a);
			return _mm_add_ps(p_a, _mm_mul_ps(delta, vt));
		}

		static ClipVertex InterpolateVertex_SSE(const ClipVertex& p_clipVertex0, const ClipVertex& p_clipVertex1, float p_t)
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

			result.TexCoord.x = p_clipVertex0.TexCoord.x + p_t * (p_clipVertex1.TexCoord.x - p_clipVertex0.TexCoord.x);
			result.TexCoord.y = p_clipVertex0.TexCoord.y + p_t * (p_clipVertex1.TexCoord.y - p_clipVertex0.TexCoord.y);

			return result;
		}
#endif

		// Clip polygon against a single plane using Sutherland-Hodgman algorithm
		static void ClipPolygonAgainstPlane(ClipPolygon& p_clipPolygon, ClipPlane p_clipPlane)
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

#ifdef OBSIDIAN_USE_AVX2
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

		static bool ClipTriangle(const ClipVertex& p_clipVertex0, const ClipVertex& p_clipVertex1, const ClipVertex& p_clipVertex2, ClipPolygon& p_outClipPolygon)
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