#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>

#include <glm/vec2.hpp>

#include "ObsidianGL/API/SIMDConfig.h"

#include "ObsidianGL/Geometry/BoundingBox2D.h"

namespace ObsidianGL::Geometry
{
	struct Triangle
	{
		static constexpr int32_t SUBPIXEL_SHIFT = 8;
		static constexpr int32_t SUBPIXEL_SCALE = 1 << SUBPIXEL_SHIFT;
		static constexpr int32_t SUBPIXEL_SCALE_HALF = SUBPIXEL_SCALE >> 1;

		BoundingBox2D BoundingBox;

		// TODO: glm::ivec2 instead...?
		int32_t FixedPointVertex0X;
		int32_t FixedPointVertex0Y;
		int32_t FixedPointVertex1X;
		int32_t FixedPointVertex1Y;
		int32_t FixedPointVertex2X;
		int32_t FixedPointVertex2Y;

		int32_t EdgeABDeltaX;
		int32_t EdgeABDeltaY;
		int32_t EdgeBCDeltaX;
		int32_t EdgeBCDeltaY;
		int32_t EdgeCADeltaX;
		int32_t EdgeCADeltaY;

		int64_t FillBiasEdgeAB;
		int64_t FillBiasEdgeBC;
		int64_t FillBiasEdgeCA;

		float FixedDoubledArea;
		float InvFixedDoubledArea;

		inline float GetArea() const
		{
			return std::abs(FixedDoubledArea) / (2.0f * static_cast<float>(SUBPIXEL_SCALE));
		}

		static inline int32_t FloatToFixedPoint(float p_value)
		{
			return static_cast<int32_t>(std::floor(p_value * static_cast<float>(SUBPIXEL_SCALE) + 0.5f));
		}

		static inline int32_t FixedPointCeil(int32_t p_fixedValue)
		{
			return (p_fixedValue + SUBPIXEL_SCALE - 1) >> SUBPIXEL_SHIFT;
		}

		static inline int32_t FixedPointFloor(int32_t p_fixedValue)
		{
			return p_fixedValue >> SUBPIXEL_SHIFT;
		}

		static inline bool IsTopOrLeftEdge(int32_t p_edgeDeltaX, int32_t p_edgeDeltaY)
		{
			return (p_edgeDeltaY > 0) || (p_edgeDeltaY == 0 && p_edgeDeltaX < 0);
		}

		Triangle(const glm::vec2& p_screenVertexA, const glm::vec2& p_screenVertexB, const glm::vec2& p_screenVertexC)
		{
			FixedPointVertex0X = FloatToFixedPoint(p_screenVertexA.x);
			FixedPointVertex0Y = FloatToFixedPoint(p_screenVertexA.y);
			FixedPointVertex1X = FloatToFixedPoint(p_screenVertexB.x);
			FixedPointVertex1Y = FloatToFixedPoint(p_screenVertexB.y);
			FixedPointVertex2X = FloatToFixedPoint(p_screenVertexC.x);
			FixedPointVertex2Y = FloatToFixedPoint(p_screenVertexC.y);

			EdgeABDeltaX = FixedPointVertex1X - FixedPointVertex0X;
			EdgeABDeltaY = FixedPointVertex1Y - FixedPointVertex0Y;
			EdgeBCDeltaX = FixedPointVertex2X - FixedPointVertex1X;
			EdgeBCDeltaY = FixedPointVertex2Y - FixedPointVertex1Y;
			EdgeCADeltaX = FixedPointVertex0X - FixedPointVertex2X;
			EdgeCADeltaY = FixedPointVertex0Y - FixedPointVertex2Y;

			FillBiasEdgeAB = IsTopOrLeftEdge(EdgeABDeltaX, EdgeABDeltaY) ? 0 : -1;
			FillBiasEdgeBC = IsTopOrLeftEdge(EdgeBCDeltaX, EdgeBCDeltaY) ? 0 : -1;
			FillBiasEdgeCA = IsTopOrLeftEdge(EdgeCADeltaX, EdgeCADeltaY) ? 0 : -1;

			int64_t fixedCrossProduct = (static_cast<int64_t>(EdgeABDeltaX) * EdgeBCDeltaY - static_cast<int64_t>(EdgeABDeltaY) * EdgeBCDeltaX) >> SUBPIXEL_SHIFT;
			FixedDoubledArea = static_cast<float>(fixedCrossProduct);
			InvFixedDoubledArea = (fixedCrossProduct != 0) ? 1.0f / FixedDoubledArea : 0.0f;

			int32_t fixedMinX = std::min({ FixedPointVertex0X, FixedPointVertex1X, FixedPointVertex2X });
			int32_t fixedMaxX = std::max({ FixedPointVertex0X, FixedPointVertex1X, FixedPointVertex2X });
			int32_t fixedMinY = std::min({ FixedPointVertex0Y, FixedPointVertex1Y, FixedPointVertex2Y });
			int32_t fixedMaxY = std::max({ FixedPointVertex0Y, FixedPointVertex1Y, FixedPointVertex2Y });

			BoundingBox.MinX = FixedPointFloor(fixedMinX);
			BoundingBox.MinY = FixedPointFloor(fixedMinY);
			BoundingBox.MaxX = FixedPointCeil(fixedMaxX);
			BoundingBox.MaxY = FixedPointCeil(fixedMaxY);
		}

		int64_t ComputeEdgeFunction(int32_t p_edgeDeltaX, int32_t p_edgeDeltaY, int32_t p_fixedPointX, int32_t p_fixedPointY, int32_t p_edgeOriginX, int32_t p_edgeOriginY, int64_t p_fillBias)
		{
			int64_t value = static_cast<int64_t>(p_edgeDeltaX) * (p_fixedPointY - p_edgeOriginY) - static_cast<int64_t>(p_edgeDeltaY) * (p_fixedPointX - p_edgeOriginX);
			return (value + p_fillBias) >> SUBPIXEL_SHIFT;
		}

		bool ComputeBarycentricCoordinates(int32_t p_pixelX, int32_t p_pixelY, float& p_outWeightA, float& p_outWeightB, float& p_outWeightC)
		{
			int32_t fixedPixelX = (p_pixelX << SUBPIXEL_SHIFT) + SUBPIXEL_SCALE_HALF;
			int32_t fixedPixelY = (p_pixelY << SUBPIXEL_SHIFT) + SUBPIXEL_SCALE_HALF;

			int64_t edgeValueAB = ComputeEdgeFunction(EdgeABDeltaX, EdgeABDeltaY, fixedPixelX, fixedPixelY, FixedPointVertex0X, FixedPointVertex0Y, FillBiasEdgeAB);
			int64_t edgeValueBC = ComputeEdgeFunction(EdgeBCDeltaX, EdgeBCDeltaY, fixedPixelX, fixedPixelY, FixedPointVertex1X, FixedPointVertex1Y, FillBiasEdgeBC);
			int64_t edgeValueCA = ComputeEdgeFunction(EdgeCADeltaX, EdgeCADeltaY, fixedPixelX, fixedPixelY, FixedPointVertex2X, FixedPointVertex2Y, FillBiasEdgeCA);

			if ((edgeValueAB | edgeValueBC | edgeValueCA) < 0)
				return false;

			p_outWeightA = static_cast<float>(edgeValueBC) * InvFixedDoubledArea;
			p_outWeightB = static_cast<float>(edgeValueCA) * InvFixedDoubledArea;
			p_outWeightC = 1.0f - p_outWeightA - p_outWeightB;

			return true;
		}

#ifdef OBSIDIAN_USE_AVX2
		struct SIMDBarycentricSetup8
		{
			int32_t EdgeABRowOrigin;
			int32_t EdgeBCRowOrigin;
			int32_t EdgeCARowOrigin;

			int32_t EdgeABStepY;
			int32_t EdgeBCStepY;
			int32_t EdgeCAStepY;

			int32_t EdgeABStepX;
			int32_t EdgeBCStepX;
			int32_t EdgeCAStepX;

			__m256i SIMDEdgeABStepX;
			__m256i SIMDEdgeBCStepX;
			__m256i SIMDEdgeCAStepX;

			__m256 SIMDInvFixedDoubledArea;
		};

		SIMDBarycentricSetup8 PrepareSIMDBarycentricSetup8() const
		{
			SIMDBarycentricSetup8 setup;

			setup.EdgeABStepX = -EdgeABDeltaY;
			setup.EdgeBCStepX = -EdgeBCDeltaY;
			setup.EdgeCAStepX = -EdgeCADeltaY;

			setup.EdgeABStepY = EdgeABDeltaX;
			setup.EdgeBCStepY = EdgeBCDeltaX;
			setup.EdgeCAStepY = EdgeCADeltaX;

			setup.EdgeABRowOrigin = static_cast<int32_t>((static_cast<int64_t>(EdgeABDeltaX) * (SUBPIXEL_SCALE_HALF - FixedPointVertex0Y) - static_cast<int64_t>(EdgeABDeltaY) * (SUBPIXEL_SCALE_HALF - FixedPointVertex0X) + FillBiasEdgeAB) >> SUBPIXEL_SHIFT);
			setup.EdgeBCRowOrigin = static_cast<int32_t>((static_cast<int64_t>(EdgeBCDeltaX) * (SUBPIXEL_SCALE_HALF - FixedPointVertex1Y) - static_cast<int64_t>(EdgeBCDeltaY) * (SUBPIXEL_SCALE_HALF - FixedPointVertex1X) + FillBiasEdgeBC) >> SUBPIXEL_SHIFT);
			setup.EdgeCARowOrigin = static_cast<int32_t>((static_cast<int64_t>(EdgeCADeltaX) * (SUBPIXEL_SCALE_HALF - FixedPointVertex2Y) - static_cast<int64_t>(EdgeCADeltaY) * (SUBPIXEL_SCALE_HALF - FixedPointVertex2X) + FillBiasEdgeCA) >> SUBPIXEL_SHIFT);

			setup.SIMDEdgeABStepX = _mm256_set1_epi32(setup.EdgeABStepX);
			setup.SIMDEdgeBCStepX = _mm256_set1_epi32(setup.EdgeBCStepX);
			setup.SIMDEdgeCAStepX = _mm256_set1_epi32(setup.EdgeCAStepX);

			setup.SIMDInvFixedDoubledArea = _mm256_set1_ps(InvFixedDoubledArea);

			return setup;
		}

		void ComputeBarycentricCoordinates8(const SIMDBarycentricSetup8& p_setup, int32_t p_startPixelX, int32_t p_pixelY, __m256& p_outWeightA, __m256& p_outWeightB, __m256& p_outWeightC, __m256& p_outInsideMask) const
		{
			const __m256i PIXEL_OFFSETS = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);

			int32_t rowEdgeAB = static_cast<int32_t>(p_setup.EdgeABRowOrigin + static_cast<int64_t>(p_pixelY) * p_setup.EdgeABStepY);
			int32_t rowEdgeBC = static_cast<int32_t>(p_setup.EdgeBCRowOrigin + static_cast<int64_t>(p_pixelY) * p_setup.EdgeBCStepY);
			int32_t rowEdgeCA = static_cast<int32_t>(p_setup.EdgeCARowOrigin + static_cast<int64_t>(p_pixelY) * p_setup.EdgeCAStepY);

			int32_t edgeABAtX = rowEdgeAB + p_startPixelX * p_setup.EdgeABStepX;
			int32_t edgeBCAtX = rowEdgeBC + p_startPixelX * p_setup.EdgeBCStepX;
			int32_t edgeCAAtX = rowEdgeCA + p_startPixelX * p_setup.EdgeCAStepX;

			__m256i vEdgeAB = _mm256_add_epi32(_mm256_set1_epi32(edgeABAtX), _mm256_mullo_epi32(PIXEL_OFFSETS, p_setup.SIMDEdgeABStepX));
			__m256i vEdgeBC = _mm256_add_epi32(_mm256_set1_epi32(edgeBCAtX), _mm256_mullo_epi32(PIXEL_OFFSETS, p_setup.SIMDEdgeBCStepX));
			__m256i vEdgeCA = _mm256_add_epi32(_mm256_set1_epi32(edgeCAAtX), _mm256_mullo_epi32(PIXEL_OFFSETS, p_setup.SIMDEdgeCAStepX));

			__m256i orAll = _mm256_or_si256(_mm256_or_si256(vEdgeBC, vEdgeCA), vEdgeAB);
			__m256i insideInt = _mm256_cmpgt_epi32(orAll, _mm256_set1_epi32(-1));
			p_outInsideMask = _mm256_castsi256_ps(insideInt);

			__m256 weightA = _mm256_mul_ps(_mm256_cvtepi32_ps(vEdgeBC), p_setup.SIMDInvFixedDoubledArea);
			__m256 weightB = _mm256_mul_ps(_mm256_cvtepi32_ps(vEdgeCA), p_setup.SIMDInvFixedDoubledArea);
			__m256 weightC = _mm256_sub_ps(_mm256_set1_ps(1.0f), _mm256_add_ps(weightA, weightB));

			p_outWeightA = _mm256_blendv_ps(_mm256_set1_ps(-1.0f), weightA, p_outInsideMask);
			p_outWeightB = _mm256_blendv_ps(_mm256_set1_ps(-1.0f), weightB, p_outInsideMask);
			p_outWeightC = _mm256_blendv_ps(_mm256_set1_ps(-1.0f), weightC, p_outInsideMask);
		}
#endif

#if defined(OBSIDIAN_USE_SSE) || defined(OBSIDIAN_USE_SSE41)
		struct SIMDBarycentricSetup4
		{
			int32_t EdgeABRowOrigin;
			int32_t EdgeBCRowOrigin;
			int32_t EdgeCARowOrigin;

			int32_t EdgeABStepY;
			int32_t EdgeBCStepY;
			int32_t EdgeCAStepY;

			int32_t EdgeABStepX;
			int32_t EdgeBCStepX;
			int32_t EdgeCAStepX;

			__m128i SIMDEdgeABStepX;
			__m128i SIMDEdgeBCStepX;
			__m128i SIMDEdgeCAStepX;

			__m128 SIMDInvFixedDoubledArea;
		};

		SIMDBarycentricSetup4 PrepareSIMDBarycentricSetup4() const
		{
			SIMDBarycentricSetup4 setup;

			setup.EdgeABStepX = -EdgeABDeltaY;
			setup.EdgeBCStepX = -EdgeBCDeltaY;
			setup.EdgeCAStepX = -EdgeCADeltaY;

			setup.EdgeABStepY = EdgeABDeltaX;
			setup.EdgeBCStepY = EdgeBCDeltaX;
			setup.EdgeCAStepY = EdgeCADeltaX;

			setup.EdgeABRowOrigin = static_cast<int32_t>((static_cast<int64_t>(EdgeABDeltaX) * (SUBPIXEL_SCALE_HALF - FixedPointVertex0Y) - static_cast<int64_t>(EdgeABDeltaY) * (SUBPIXEL_SCALE_HALF - FixedPointVertex0X) + FillBiasEdgeAB) >> SUBPIXEL_SHIFT);
			setup.EdgeBCRowOrigin = static_cast<int32_t>((static_cast<int64_t>(EdgeBCDeltaX) * (SUBPIXEL_SCALE_HALF - FixedPointVertex1Y) - static_cast<int64_t>(EdgeBCDeltaY) * (SUBPIXEL_SCALE_HALF - FixedPointVertex1X) + FillBiasEdgeBC) >> SUBPIXEL_SHIFT);
			setup.EdgeCARowOrigin = static_cast<int32_t>((static_cast<int64_t>(EdgeCADeltaX) * (SUBPIXEL_SCALE_HALF - FixedPointVertex2Y) - static_cast<int64_t>(EdgeCADeltaY) * (SUBPIXEL_SCALE_HALF - FixedPointVertex2X) + FillBiasEdgeCA) >> SUBPIXEL_SHIFT);

			setup.SIMDEdgeABStepX = _mm_set1_epi32(setup.EdgeABStepX);
			setup.SIMDEdgeBCStepX = _mm_set1_epi32(setup.EdgeBCStepX);
			setup.SIMDEdgeCAStepX = _mm_set1_epi32(setup.EdgeCAStepX);

			setup.SIMDInvFixedDoubledArea = _mm_set1_ps(InvFixedDoubledArea);

			return setup;
		}

		void ComputeBarycentricCoordinates4(const SIMDBarycentricSetup4& p_setup, int32_t p_startPixelX, int32_t p_pixelY, __m128& p_outWeightA, __m128& p_outWeightB, __m128& p_outWeightC, __m128& p_outInsideMask) const
		{
			const __m128i PIXEL_OFFSETS = _mm_set_epi32(3, 2, 1, 0);

			int32_t rowEdgeBC = static_cast<int32_t>(p_setup.EdgeBCRowOrigin + static_cast<int64_t>(p_pixelY) * p_setup.EdgeBCStepY);
			int32_t rowEdgeCA = static_cast<int32_t>(p_setup.EdgeCARowOrigin + static_cast<int64_t>(p_pixelY) * p_setup.EdgeCAStepY);
			int32_t rowEdgeAB = static_cast<int32_t>(p_setup.EdgeABRowOrigin + static_cast<int64_t>(p_pixelY) * p_setup.EdgeABStepY);

			int32_t edgeBCAtX = rowEdgeBC + p_startPixelX * p_setup.EdgeBCStepX;
			int32_t edgeCAAtX = rowEdgeCA + p_startPixelX * p_setup.EdgeCAStepX;
			int32_t edgeABAtX = rowEdgeAB + p_startPixelX * p_setup.EdgeABStepX;

#ifdef OBSIDIAN_USE_SSE41
			__m128i vEdgeBC = _mm_add_epi32(_mm_set1_epi32(edgeBCAtX), _mm_mullo_epi32(PIXEL_OFFSETS, p_setup.SIMDEdgeBCStepX));
			__m128i vEdgeCA = _mm_add_epi32(_mm_set1_epi32(edgeCAAtX), _mm_mullo_epi32(PIXEL_OFFSETS, p_setup.SIMDEdgeCAStepX));
			__m128i vEdgeAB = _mm_add_epi32(_mm_set1_epi32(edgeABAtX), _mm_mullo_epi32(PIXEL_OFFSETS, p_setup.SIMDEdgeABStepX));
#else
		
			int32_t stepBC = p_setup.EdgeBCStepX;
			int32_t stepCA = p_setup.EdgeCAStepX;
			int32_t stepAB = p_setup.EdgeABStepX;

			__m128i vEdgeBC = _mm_set_epi32(edgeBCAtX + 3 * stepBC, edgeBCAtX + 2 * stepBC, edgeBCAtX + stepBC, edgeBCAtX);
			__m128i vEdgeCA = _mm_set_epi32(edgeCAAtX + 3 * stepCA, edgeCAAtX + 2 * stepCA, edgeCAAtX + stepCA, edgeCAAtX);
			__m128i vEdgeAB = _mm_set_epi32(edgeABAtX + 3 * stepAB, edgeABAtX + 2 * stepAB, edgeABAtX + stepAB, edgeABAtX);
#endif

			__m128i orAll = _mm_or_si128(_mm_or_si128(vEdgeBC, vEdgeCA), vEdgeAB);
			__m128i insideInt = _mm_cmpgt_epi32(orAll, _mm_set1_epi32(-1));

			p_outInsideMask = _mm_castsi128_ps(insideInt);

			__m128 weightA = _mm_mul_ps(_mm_cvtepi32_ps(vEdgeBC), p_setup.SIMDInvFixedDoubledArea);
			__m128 weightB = _mm_mul_ps(_mm_cvtepi32_ps(vEdgeCA), p_setup.SIMDInvFixedDoubledArea);
			__m128 weightC = _mm_sub_ps(_mm_set1_ps(1.0f), _mm_add_ps(weightA, weightB));

#ifdef OBSIDIAN_USE_SSE41
			p_outWeightA = _mm_blendv_ps(_mm_set1_ps(-1.0f), weightA, p_outInsideMask);
			p_outWeightB = _mm_blendv_ps(_mm_set1_ps(-1.0f), weightB, p_outInsideMask);
			p_outWeightC = _mm_blendv_ps(_mm_set1_ps(-1.0f), weightC, p_outInsideMask);
#else
			p_outWeightA = _mm_or_ps(_mm_and_ps(p_outInsideMask, weightA), _mm_andnot_ps(p_outInsideMask, _mm_set1_ps(-1.0f)));
			p_outWeightB = _mm_or_ps(_mm_and_ps(p_outInsideMask, weightB), _mm_andnot_ps(p_outInsideMask, _mm_set1_ps(-1.0f)));
			p_outWeightC = _mm_or_ps(_mm_and_ps(p_outInsideMask, weightC), _mm_andnot_ps(p_outInsideMask, _mm_set1_ps(-1.0f)));
#endif
		}
#endif
	};
}