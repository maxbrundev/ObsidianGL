#pragma once

#include <algorithm>
#include <array>
#include <cmath>

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#ifdef __AVX2__
#include <immintrin.h>
#else
#include <emmintrin.h>
#endif

#include "ObsidianGL/Geometry/BoundingBox2D.h"

namespace ObsidianGL::Geometry
{
	struct Triangle
	{
		alignas(16) float VerticesX[4];
		alignas(16) float VerticesY[4];
		BoundingBox2D BoundingBox;

		float V0X;
		float V0Y;
		float V1X;
		float V1Y;
		float D00;
		float D01;
		float D11;
		float InvDenom;

		float Edge0DeltaX;
		float Edge0DeltaY;
		float Edge1DeltaX;
		float Edge1DeltaY;
		float Edge2DeltaX;
		float Edge2DeltaY;

		Triangle(const glm::vec2& p_vertex0, const glm::vec2& p_vertex1, const glm::vec2& p_vertex2)
		{
			VerticesX[0] = p_vertex0.x;
			VerticesY[0] = p_vertex0.y;
			VerticesX[1] = p_vertex1.x;
			VerticesY[1] = p_vertex1.y;
			VerticesX[2] = p_vertex2.x;
			VerticesY[2] = p_vertex2.y;
			// Duplicate vertex 2 instead of 0.0f to prevent min/max pollution
			VerticesX[3] = p_vertex2.x;
			VerticesY[3] = p_vertex2.y;

			ComputeBoundingBox();
			PrecomputeBarycentric();
			PrecomputeEdgeDeltas();
		}

		void ComputeBoundingBox()
		{
			// Scalar bounding box - more reliable than SIMD horizontal reductions
			float minX = std::min({ VerticesX[0], VerticesX[1], VerticesX[2] });
			float maxX = std::max({ VerticesX[0], VerticesX[1], VerticesX[2] });
			float minY = std::min({ VerticesY[0], VerticesY[1], VerticesY[2] });
			float maxY = std::max({ VerticesY[0], VerticesY[1], VerticesY[2] });

			BoundingBox.MinX = static_cast<int32_t>(std::floor(minX));
			BoundingBox.MinY = static_cast<int32_t>(std::floor(minY));
			BoundingBox.MaxX = static_cast<int32_t>(std::ceil(maxX));
			BoundingBox.MaxY = static_cast<int32_t>(std::ceil(maxY));
		}

		void PrecomputeBarycentric()
		{
			V0X = VerticesX[2] - VerticesX[0];
			V0Y = VerticesY[2] - VerticesY[0];
			V1X = VerticesX[1] - VerticesX[0];
			V1Y = VerticesY[1] - VerticesY[0];

			D00 = V0X * V0X + V0Y * V0Y;
			D01 = V0X * V1X + V0Y * V1Y;
			D11 = V1X * V1X + V1Y * V1Y;

			float denom = D00 * D11 - D01 * D01;
			InvDenom = (denom != 0.0f) ? 1.0f / denom : 0.0f;
		}

		void PrecomputeEdgeDeltas()
		{
			Edge0DeltaX = VerticesY[1] - VerticesY[0];
			Edge0DeltaY = VerticesX[0] - VerticesX[1];
			Edge1DeltaX = VerticesY[2] - VerticesY[1];
			Edge1DeltaY = VerticesX[1] - VerticesX[2];
			Edge2DeltaX = VerticesY[0] - VerticesY[2];
			Edge2DeltaY = VerticesX[2] - VerticesX[0];
		}

		glm::vec3 GetBarycentricCoordinates(float p_positionX, float p_positionY) const
		{
			float x = p_positionX - VerticesX[0];
			float y = p_positionY - VerticesY[0];

			float d02 = V0X * x + V0Y * y;
			float d12 = V1X * x + V1Y * y;

			float u = (D00 * d12 - D01 * d02) * InvDenom;
			float v = (D11 * d02 - D01 * d12) * InvDenom;

			return glm::vec3(1.0f - u - v, u, v);
		}

#ifdef __AVX2__
		void GetBarycentricCoordinates8(const float* p_positionX, float p_positionY, float* p_outU, float* p_outV, float* p_outW) const
		{
			__m256 vx0 = _mm256_set1_ps(VerticesX[0]);
			__m256 vy0 = _mm256_set1_ps(VerticesY[0]);
			__m256 v0x = _mm256_set1_ps(V0X);
			__m256 v0y = _mm256_set1_ps(V0Y);
			__m256 v1x = _mm256_set1_ps(V1X);
			__m256 v1y = _mm256_set1_ps(V1Y);
			__m256 d00 = _mm256_set1_ps(D00);
			__m256 d01 = _mm256_set1_ps(D01);
			__m256 d11 = _mm256_set1_ps(D11);
			__m256 invDenom = _mm256_set1_ps(InvDenom);
			__m256 one = _mm256_set1_ps(1.0f);

			__m256 pxVec = _mm256_loadu_ps(p_positionX);
			__m256 pyVec = _mm256_set1_ps(p_positionY);

			__m256 x = _mm256_sub_ps(pxVec, vx0);
			__m256 y = _mm256_sub_ps(pyVec, vy0);

			__m256 d02 = _mm256_add_ps(_mm256_mul_ps(v0x, x), _mm256_mul_ps(v0y, y));
			__m256 d12 = _mm256_add_ps(_mm256_mul_ps(v1x, x), _mm256_mul_ps(v1y, y));

			__m256 u = _mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(d00, d12), _mm256_mul_ps(d01, d02)), invDenom);
			__m256 v = _mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(d11, d02), _mm256_mul_ps(d01, d12)), invDenom);
			__m256 w = _mm256_sub_ps(_mm256_sub_ps(one, u), v);

			_mm256_storeu_ps(p_outW, w);
			_mm256_storeu_ps(p_outU, u);
			_mm256_storeu_ps(p_outV, v);
		}
#else
		void GetBarycentricCoordinates4(const float* p_positionX, float p_positionY, float* p_outU, float* p_outV, float* p_outW) const
		{
			__m128 vx0 = _mm_set1_ps(VerticesX[0]);
			__m128 vy0 = _mm_set1_ps(VerticesY[0]);
			__m128 v0x = _mm_set1_ps(V0X);
			__m128 v0y = _mm_set1_ps(V0Y);
			__m128 v1x = _mm_set1_ps(V1X);
			__m128 v1y = _mm_set1_ps(V1Y);
			__m128 d00 = _mm_set1_ps(D00);
			__m128 d01 = _mm_set1_ps(D01);
			__m128 d11 = _mm_set1_ps(D11);
			__m128 invDenom = _mm_set1_ps(InvDenom);
			__m128 one = _mm_set1_ps(1.0f);

			__m128 pxVec = _mm_loadu_ps(p_positionX);
			__m128 pyVec = _mm_set1_ps(p_positionY);

			__m128 x = _mm_sub_ps(pxVec, vx0);
			__m128 y = _mm_sub_ps(pyVec, vy0);

			__m128 d02 = _mm_add_ps(_mm_mul_ps(v0x, x), _mm_mul_ps(v0y, y));
			__m128 d12 = _mm_add_ps(_mm_mul_ps(v1x, x), _mm_mul_ps(v1y, y));

			__m128 u = _mm_mul_ps(_mm_sub_ps(_mm_mul_ps(d00, d12), _mm_mul_ps(d01, d02)), invDenom);
			__m128 v = _mm_mul_ps(_mm_sub_ps(_mm_mul_ps(d11, d02), _mm_mul_ps(d01, d12)), invDenom);
			__m128 w = _mm_sub_ps(_mm_sub_ps(one, u), v);

			_mm_storeu_ps(p_outW, w);
			_mm_storeu_ps(p_outU, u);
			_mm_storeu_ps(p_outV, v);
		}
#endif

		float ComputeEdge(const glm::vec2& p_vertex0, const glm::vec2& p_vertex1, const glm::vec2& p_vertex2) const
		{
			return (p_vertex1.x - p_vertex0.x) * (p_vertex2.y - p_vertex0.y) - (p_vertex1.y - p_vertex0.y) * (p_vertex2.x - p_vertex0.x);
		}

		float ComputeArea() const
		{
			return (VerticesX[1] - VerticesX[0]) * (VerticesY[2] - VerticesY[0]) - (VerticesY[1] - VerticesY[0]) * (VerticesX[2] - VerticesX[0]);
		}

		float ComputeMaxDepthSlope(const std::array<glm::vec4, 3>& p_vertices) const
		{
			float dx = p_vertices[1].x - p_vertices[0].x;
			float dy = p_vertices[2].y - p_vertices[0].y;

			float dzdx = (dx != 0.0f) ? std::abs((p_vertices[1].z - p_vertices[0].z) / dx) : 0.0f;
			float dzdy = (dy != 0.0f) ? std::abs((p_vertices[2].z - p_vertices[0].z) / dy) : 0.0f;

			return std::max(dzdx, dzdy);
		}
	};
}