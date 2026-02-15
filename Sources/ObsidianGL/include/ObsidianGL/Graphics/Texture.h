#pragma once

#include <cstdint>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>

#include "ObsidianGL/API/SIMDConfig.h"

#include "ObsidianGL/Utils/Memory/Common.h"

namespace ObsidianGL::Graphics
{
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

	struct Texture
	{
		uint32_t* Data = nullptr;

		uint32_t Width = 0;
		uint32_t Height = 0;

		uint16_t MagFilter = OGL_NEAREST;
		uint16_t MinFilter = OGL_NEAREST;
		uint16_t WrapS = OGL_REPEAT;
		uint16_t WrapT = OGL_REPEAT;

		int32_t WidthMask = 0;
		int32_t HeightMask = 0;

		bool IsWidthPow2 = false;
		bool IsHeightPow2 = false;

		~Texture()
		{
			if (Data != nullptr)
			{
				Utils::Memory::AlignedFree(Data);
				Data = nullptr;
			}
		}

		void Upload(uint32_t p_width, uint32_t p_height, uint16_t p_format, const uint8_t* p_pixels)
		{
			if (Data != nullptr)
			{
				Utils::Memory::AlignedFree(Data);
			}

			Width = p_width;
			Height = p_height;

			IsWidthPow2 = (p_width & (p_width - 1)) == 0 && p_width > 0;
			IsHeightPow2 = (p_height & (p_height - 1)) == 0 && p_height > 0;

			WidthMask = static_cast<int32_t>(p_width - 1);
			HeightMask = static_cast<int32_t>(p_height - 1);

			uint32_t pixelCount = p_width * p_height;
			Data = static_cast<uint32_t*>(Utils::Memory::AlignedAlloc(pixelCount * sizeof(uint32_t), 32));

			if (p_format == OGL_RGBA)
			{
				for (uint32_t i = 0; i < pixelCount; ++i)
				{
					uint8_t r = p_pixels[i * 4 + 0];
					uint8_t g = p_pixels[i * 4 + 1];
					uint8_t b = p_pixels[i * 4 + 2];
					uint8_t a = p_pixels[i * 4 + 3];
					Data[i] = (static_cast<uint32_t>(r) << 24) | (static_cast<uint32_t>(g) << 16) | (static_cast<uint32_t>(b) << 8) | static_cast<uint32_t>(a);
				}
			}
			else if (p_format == OGL_RGB)
			{
				for (uint32_t i = 0; i < pixelCount; ++i)
				{
					uint8_t r = p_pixels[i * 3 + 0];
					uint8_t g = p_pixels[i * 3 + 1];
					uint8_t b = p_pixels[i * 3 + 2];
					Data[i] = (static_cast<uint32_t>(r) << 24) | (static_cast<uint32_t>(g) << 16) | (static_cast<uint32_t>(b) << 8) | 0xFF;
				}
			}
		}

		int32_t WrapCoordX(int32_t p_coord) const
		{
			if (WrapS == OGL_REPEAT)
			{
				if (IsWidthPow2)
				{
					return p_coord & WidthMask;
				}

				p_coord = p_coord % static_cast<int32_t>(Width);

				return (p_coord < 0) ? p_coord + static_cast<int32_t>(Width) : p_coord;
			}

			if (p_coord < 0)
				return 0;

			if (p_coord >= static_cast<int32_t>(Width))
			{
				return static_cast<int32_t>(Width) - 1;
			}

			return p_coord;
		}

		int32_t WrapCoordY(int32_t p_coord) const
		{
			if (WrapT == OGL_REPEAT)
			{
				if (IsHeightPow2)
				{
					return p_coord & HeightMask;
				}

				p_coord = p_coord % static_cast<int32_t>(Height);

				return (p_coord < 0) ? p_coord + static_cast<int32_t>(Height) : p_coord;
			}

			if (p_coord < 0)
				return 0;

			if (p_coord >= static_cast<int32_t>(Height))
			{
				return static_cast<int32_t>(Height) - 1;
			}

			return p_coord;
		}

		uint32_t SampleNearest(float p_u, float p_v) const
		{
			int32_t texelX = static_cast<int32_t>(p_u * Width);
			int32_t texelY = static_cast<int32_t>(p_v * Height);

			texelX = WrapCoordX(texelX);
			texelY = WrapCoordY(texelY);

			return Data[texelY * Width + texelX];
		}

		uint32_t SampleBilinear(float p_u, float p_v) const
		{
			float texelCenterX = p_u * Width - 0.5f;
			float texelCenterY = p_v * Height - 0.5f;

			int32_t texelFloorX = static_cast<int32_t>(std::floor(texelCenterX));
			int32_t texelFloorY = static_cast<int32_t>(std::floor(texelCenterY));
			int32_t texelCeilX = texelFloorX + 1;
			int32_t texelCeilY = texelFloorY + 1;

			float blendFactorX = texelCenterX - static_cast<float>(texelFloorX);
			float blendFactorY = texelCenterY - static_cast<float>(texelFloorY);

			texelFloorX = WrapCoordX(texelFloorX);
			texelFloorY = WrapCoordY(texelFloorY);
			texelCeilX = WrapCoordX(texelCeilX);
			texelCeilY = WrapCoordY(texelCeilY);

			uint32_t topLeftTexel = Data[texelFloorY * Width + texelFloorX];
			uint32_t topRightTexel = Data[texelFloorY * Width + texelCeilX];
			uint32_t bottomLeftTexel = Data[texelCeilY * Width + texelFloorX];
			uint32_t bottomRightTexel = Data[texelCeilY * Width + texelCeilX];

			float topLeftWeight = (1.0f - blendFactorX) * (1.0f - blendFactorY);
			float topRightWeight = blendFactorX * (1.0f - blendFactorY);
			float bottomLeftWeight = (1.0f - blendFactorX) * blendFactorY;
			float bottomRightWeight = blendFactorX * blendFactorY;

			float blendedRed = static_cast<float>((topLeftTexel >> 24) & 0xFF) * topLeftWeight + static_cast<float>((topRightTexel >> 24) & 0xFF) * topRightWeight + static_cast<float>((bottomLeftTexel >> 24) & 0xFF) * bottomLeftWeight + static_cast<float>((bottomRightTexel >> 24) & 0xFF) * bottomRightWeight;
			float blendedGreen = static_cast<float>((topLeftTexel >> 16) & 0xFF) * topLeftWeight + static_cast<float>((topRightTexel >> 16) & 0xFF) * topRightWeight + static_cast<float>((bottomLeftTexel >> 16) & 0xFF) * bottomLeftWeight + static_cast<float>((bottomRightTexel >> 16) & 0xFF) * bottomRightWeight;
			float blendedBlue = static_cast<float>((topLeftTexel >> 8) & 0xFF) * topLeftWeight + static_cast<float>((topRightTexel >> 8) & 0xFF) * topRightWeight + static_cast<float>((bottomLeftTexel >> 8) & 0xFF) * bottomLeftWeight + static_cast<float>((bottomRightTexel >> 8) & 0xFF) * bottomRightWeight;
			float blendedAlpha = static_cast<float>(topLeftTexel & 0xFF) * topLeftWeight + static_cast<float>(topRightTexel & 0xFF) * topRightWeight + static_cast<float>(bottomLeftTexel & 0xFF) * bottomLeftWeight + static_cast<float>(bottomRightTexel & 0xFF) * bottomRightWeight;

			uint8_t finalRed = static_cast<uint8_t>(std::clamp(blendedRed, 0.0f, 255.0f));
			uint8_t finalGreen = static_cast<uint8_t>(std::clamp(blendedGreen, 0.0f, 255.0f));
			uint8_t finalBlue = static_cast<uint8_t>(std::clamp(blendedBlue, 0.0f, 255.0f));
			uint8_t finalAlpha = static_cast<uint8_t>(std::clamp(blendedAlpha, 0.0f, 255.0f));

			return (static_cast<uint32_t>(finalRed) << 24) | (static_cast<uint32_t>(finalGreen) << 16) | (static_cast<uint32_t>(finalBlue) << 8) | static_cast<uint32_t>(finalAlpha);
		}

		uint32_t Sample(float p_u, float p_v) const
		{
			if (MagFilter == OGL_LINEAR || MinFilter == OGL_LINEAR)
			{
				return SampleBilinear(p_u, p_v);
			}

			return SampleNearest(p_u, p_v);
		}

#ifdef OBSIDIAN_USE_AVX2
		void SampleNearest8_AVX2(const float* p_u, const float* p_v, uint32_t* p_outTexels) const
		{
			__m256i texelCoordsX = _mm256_cvttps_epi32(_mm256_mul_ps(_mm256_loadu_ps(p_u), _mm256_set1_ps(Width)));
			__m256i texelCoordsY = _mm256_cvttps_epi32(_mm256_mul_ps(_mm256_loadu_ps(p_v), _mm256_set1_ps(Height)));

			bool wrapXResolved = true;
			bool wrapYResolved = true;

			if (WrapS == OGL_REPEAT && IsWidthPow2)
			{
				texelCoordsX = _mm256_and_si256(texelCoordsX, _mm256_set1_epi32(WidthMask));
			}
			else if (WrapS == OGL_CLAMP_TO_EDGE)
			{
				texelCoordsX = _mm256_max_epi32(_mm256_setzero_si256(), _mm256_min_epi32(texelCoordsX, _mm256_set1_epi32(static_cast<int32_t>(Width - 1))));
			}
			else
			{
				wrapXResolved = false;
			}

			if (WrapT == OGL_REPEAT && IsHeightPow2)
			{
				texelCoordsY = _mm256_and_si256(texelCoordsY, _mm256_set1_epi32(HeightMask));
			}
			else if (WrapT == OGL_CLAMP_TO_EDGE)
			{
				texelCoordsY = _mm256_max_epi32(_mm256_setzero_si256(), _mm256_min_epi32(texelCoordsY, _mm256_set1_epi32(static_cast<int32_t>(Height - 1))));
			}
			else
			{
				wrapYResolved = false;
			}

			if (!wrapXResolved || !wrapYResolved)
			{
				alignas(32) int32_t scalarCoordsX[8];
				alignas(32) int32_t scalarCoordsY[8];

				_mm256_store_si256(reinterpret_cast<__m256i*>(scalarCoordsX), texelCoordsX);
				_mm256_store_si256(reinterpret_cast<__m256i*>(scalarCoordsY), texelCoordsY);

				for (int i = 0; i < 8; ++i)
				{
					if (!wrapXResolved)
					{
						scalarCoordsX[i] = scalarCoordsX[i] % Width;

						if (scalarCoordsX[i] < 0)
						{
							scalarCoordsX[i] += Width;
						}
					}

					if (!wrapYResolved)
					{
						scalarCoordsY[i] = scalarCoordsY[i] % Height;

						if (scalarCoordsY[i] < 0)
						{
							scalarCoordsY[i] += Height;
						}
					}
				}

				texelCoordsX = _mm256_load_si256(reinterpret_cast<const __m256i*>(scalarCoordsX));
				texelCoordsY = _mm256_load_si256(reinterpret_cast<const __m256i*>(scalarCoordsY));
			}

			__m256i texelIndices = _mm256_add_epi32(_mm256_mullo_epi32(texelCoordsY, _mm256_set1_epi32(static_cast<int32_t>(Width))), texelCoordsX);
			__m256i sampledTexels = _mm256_i32gather_epi32(reinterpret_cast<const int*>(Data), texelIndices, 4);

			_mm256_storeu_si256(reinterpret_cast<__m256i*>(p_outTexels), sampledTexels);
		}

#endif

#if defined(OBSIDIAN_USE_SSE) || defined(OBSIDIAN_USE_SSE41)
		void SampleNearest4_SSE(const float* p_u, const float* p_v, uint32_t* p_outTexels) const
		{
			__m128i texelCoordsX = _mm_cvttps_epi32(_mm_mul_ps(_mm_loadu_ps(p_u), _mm_set1_ps(Width)));
			__m128i texelCoordsY = _mm_cvttps_epi32(_mm_mul_ps(_mm_loadu_ps(p_v), _mm_set1_ps(Height)));

			bool wrapXResolved = (WrapS == OGL_REPEAT && IsWidthPow2);
			bool wrapYResolved = (WrapT == OGL_REPEAT && IsHeightPow2);

			if (wrapXResolved)
			{
				texelCoordsX = _mm_and_si128(texelCoordsX, _mm_set1_epi32(WidthMask));
			}

			if (wrapYResolved)
			{
				texelCoordsY = _mm_and_si128(texelCoordsY, _mm_set1_epi32(HeightMask));
			}

#ifdef OBSIDIAN_USE_SSE41
			if (!wrapXResolved && WrapS == OGL_CLAMP_TO_EDGE)
			{
				texelCoordsX = _mm_max_epi32(_mm_setzero_si128(), _mm_min_epi32(texelCoordsX, _mm_set1_epi32(static_cast<int32_t>(Width - 1))));

				wrapXResolved = true;
			}

			if (!wrapYResolved && WrapT == OGL_CLAMP_TO_EDGE)
			{
				texelCoordsY = _mm_max_epi32(_mm_setzero_si128(), _mm_min_epi32(texelCoordsY, _mm_set1_epi32(static_cast<int32_t>(Height - 1))));

				wrapYResolved = true;
			}

			if (wrapXResolved && wrapYResolved)
			{
				__m128i texelIndices = _mm_add_epi32(_mm_mullo_epi32(texelCoordsY, _mm_set1_epi32(static_cast<int32_t>(Width))), texelCoordsX);

				alignas(16) int32_t linearIndices[4];
				_mm_store_si128(reinterpret_cast<__m128i*>(linearIndices), texelIndices);

				p_outTexels[0] = Data[linearIndices[0]];
				p_outTexels[1] = Data[linearIndices[1]];
				p_outTexels[2] = Data[linearIndices[2]];
				p_outTexels[3] = Data[linearIndices[3]];
			}
			else
			{
				alignas(16) int32_t scalarCoordsX[4];
				alignas(16) int32_t scalarCoordsY[4];

				_mm_store_si128(reinterpret_cast<__m128i*>(scalarCoordsX), texelCoordsX);
				_mm_store_si128(reinterpret_cast<__m128i*>(scalarCoordsY), texelCoordsY);

				int32_t signedWidth = static_cast<int32_t>(Width);
				int32_t signedHeight = static_cast<int32_t>(Height);

				for (int i = 0; i < 4; ++i)
				{
					int32_t wrappedX = scalarCoordsX[i];
					int32_t wrappedY = scalarCoordsY[i];

					if (!wrapXResolved)
					{
						wrappedX = wrappedX % signedWidth;

						if (wrappedX < 0)
						{
							wrappedX += signedWidth;
						}
					}

					if (!wrapYResolved)
					{
						wrappedY = wrappedY % signedHeight;

						if (wrappedY < 0)
						{
							wrappedY += signedHeight;
						}
					}

					p_outTexels[i] = Data[wrappedY * Width + wrappedX];
				}
			}
#else
			alignas(16) int32_t scalarCoordsX[4];
			alignas(16) int32_t scalarCoordsY[4];

			_mm_store_si128(reinterpret_cast<__m128i*>(scalarCoordsX), texelCoordsX);
			_mm_store_si128(reinterpret_cast<__m128i*>(scalarCoordsY), texelCoordsY);

			if (wrapXResolved && wrapYResolved)
			{
				p_outTexels[0] = Data[scalarCoordsY[0] * Width + scalarCoordsX[0]];
				p_outTexels[1] = Data[scalarCoordsY[1] * Width + scalarCoordsX[1]];
				p_outTexels[2] = Data[scalarCoordsY[2] * Width + scalarCoordsX[2]];
				p_outTexels[3] = Data[scalarCoordsY[3] * Width + scalarCoordsX[3]];
			}
			else
			{
				for (int i = 0; i < 4; ++i)
				{
					int32_t wrappedX = scalarCoordsX[i];
					int32_t wrappedY = scalarCoordsY[i];

					if (!wrapXResolved)
					{
						if (WrapS == OGL_REPEAT)
						{
							wrappedX = wrappedX % Width;

							if (wrappedX < 0)
							{
								wrappedX += Width;
							}
						}
						else
						{
							// OGL_CLAMP_TO_EDGE
							if (wrappedX < 0)
							{
								wrappedX = 0;
							}

							else if (wrappedX >= Width)
							{
								wrappedX = Width - 1;
							}
						}
					}

					if (!wrapYResolved)
					{
						if (WrapT == OGL_REPEAT)
						{
							wrappedY = wrappedY % Height;

							if (wrappedY < 0)
							{
								wrappedY += Height;
							}
						}
						else
						{
							// OGL_CLAMP_TO_EDGE
							if (wrappedY < 0)
							{
								wrappedY = 0;
							}
							else if (wrappedY >= Height)
							{
								wrappedY = Height - 1;
							}
						}
					}

					p_outTexels[i] = Data[wrappedY * Width + wrappedX];
				}
			}
#endif
		}
#endif
	};
}