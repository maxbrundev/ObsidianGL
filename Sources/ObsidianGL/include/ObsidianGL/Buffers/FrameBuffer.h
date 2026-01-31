#pragma once

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <limits>

#ifdef _WIN32
#include <malloc.h>
#endif

#ifdef __AVX2__
#include <immintrin.h>
#else
#include <emmintrin.h>
#endif

typedef uint32_t RGBA8;
typedef float Depth;

namespace ObsidianGL::Buffers
{
	inline void* AlignedAlloc(size_t p_size, size_t p_alignment)
	{
		size_t alignedSize = (p_size + p_alignment - 1) & ~(p_alignment - 1);

#ifdef _WIN32
		return _aligned_malloc(alignedSize, p_alignment);
#else
		return std::aligned_alloc(p_alignment, alignedSize);
#endif
	}

	inline void AlignedFree(void* p_ptr)
	{
#ifdef _WIN32
		_aligned_free(p_ptr);
#else
		std::free(p_ptr);
#endif
	}

	inline uint32_t FloatToUint32(float p_value)
	{
		uint32_t bits;
		std::memcpy(&bits, &p_value, sizeof(float));
		return bits;
	}

	inline float Uint32ToFloat(uint32_t p_bits)
	{
		float value;
		std::memcpy(&value, &p_bits, sizeof(float));
		return value;
	}

	inline void SIMDFill(uint32_t* p_data, size_t p_size, uint32_t p_value)
	{
#ifdef __AVX2__
		__m256i fillValue = _mm256_set1_epi32(static_cast<int>(p_value));

		size_t blockCount = p_size / 16;

		for (size_t i = 0; i < blockCount; i++)
		{
			_mm256_store_si256(reinterpret_cast<__m256i*>(p_data + i * 16), fillValue);
			_mm256_store_si256(reinterpret_cast<__m256i*>(p_data + i * 16 + 8), fillValue);
		}

		for (size_t i = blockCount * 16; i < p_size; i++)
		{
			p_data[i] = p_value;
		}
#else
		__m128i fillValue = _mm_set1_epi32(static_cast<int>(p_value));

		size_t blockCount = p_size / 16;

		for (size_t i = 0; i < blockCount; i++)
		{
			_mm_store_si128(reinterpret_cast<__m128i*>(p_data + i * 16), fillValue);
			_mm_store_si128(reinterpret_cast<__m128i*>(p_data + i * 16 + 4), fillValue);
			_mm_store_si128(reinterpret_cast<__m128i*>(p_data + i * 16 + 8), fillValue);
			_mm_store_si128(reinterpret_cast<__m128i*>(p_data + i * 16 + 12), fillValue);
		}

		for (size_t i = blockCount * 16; i < p_size; i++)
		{
			p_data[i] = p_value;
		}
#endif
	}

	template<typename T>
	struct FrameBufferBase
	{
		static_assert(std::is_same_v<T, RGBA8> || std::is_same_v<T, Depth>, "FrameBuffer supports RGBA8 (uint32_t), Depth (float).");

		static constexpr size_t Alignment = 32;

		uint32_t Width;
		uint32_t Height;
		uint32_t Size;
		T* Data;

		FrameBufferBase(uint32_t p_width, uint32_t p_height) :
		Width(p_width), 
		Height(p_height), 
		Size(Width* Height), 
		Data(static_cast<T*>(AlignedAlloc(Size * sizeof(T), Alignment)))
		{
		}

		FrameBufferBase(const FrameBufferBase&) = delete;
		FrameBufferBase& operator=(const FrameBufferBase&) = delete;

		~FrameBufferBase()
		{
			AlignedFree(Data);
		}

		FrameBufferBase& operator=(FrameBufferBase&& p_other) noexcept
		{
			if (this != &p_other)
			{
				AlignedFree(Data);

				Width = p_other.Width;
				Height = p_other.Height;
				Size = p_other.Size;
				Data = p_other.Data;

				p_other.Data = nullptr;
				p_other.Width = 0;
				p_other.Height = 0;
				p_other.Size = 0;
			}

			return *this;
		}

		void SetPixel(uint32_t p_x, uint32_t p_y, T p_value)
		{
			if constexpr (std::is_same_v<T, RGBA8>)
			{
				Data[p_y * Width + p_x] = p_value;
			}
			else if constexpr (std::is_same_v<T, Depth>)
			{
				reinterpret_cast<uint32_t*>(Data)[p_y * Width + p_x] = FloatToUint32(p_value);
			}
		}

		T GetPixel(uint32_t p_x, uint32_t p_y) const
		{
			if constexpr (std::is_same_v<T, RGBA8>)
			{
				return Data[p_y * Width + p_x];
			}
			
			// T is Depth.
			return Uint32ToFloat(reinterpret_cast<uint32_t*>(Data)[p_y * Width + p_x]);
		}

		T* GetRow(uint32_t p_y)
		{
			return &Data[p_y * Width];
		}
	};

	template<typename T>
	struct FrameBuffer;

	template<>
	struct FrameBuffer<RGBA8> : FrameBufferBase<RGBA8>
	{
		uint32_t ClearColor = 0x000000FF;

		FrameBuffer(uint32_t p_width, uint32_t p_height) :
			FrameBufferBase<RGBA8>(p_width, p_height)
		{
		}

		void SetClearColor(float p_r, float p_g, float p_b, float p_a)
		{
			uint8_t r = static_cast<uint8_t>(std::clamp(p_r, 0.0f, 1.0f) * 255.0f);
			uint8_t g = static_cast<uint8_t>(std::clamp(p_g, 0.0f, 1.0f) * 255.0f);
			uint8_t b = static_cast<uint8_t>(std::clamp(p_b, 0.0f, 1.0f) * 255.0f);
			uint8_t a = static_cast<uint8_t>(std::clamp(p_a, 0.0f, 1.0f) * 255.0f);

			ClearColor = (r << 24) | (g << 16) | (b << 8) | a;
		}

		void Clear()
		{
			SIMDFill(Data, Size, ClearColor);
		}

		void Resize(uint32_t p_width, uint32_t p_height)
		{
			if (Width == p_width && Height == p_height)
				return;

			AlignedFree(Data);

			Width = p_width;
			Height = p_height;
			Size = Width * Height;

			Data = static_cast<RGBA8*>(AlignedAlloc(Size * sizeof(RGBA8), Alignment));
			Clear();
		}
	};

	template<>
	struct FrameBuffer<Depth> : FrameBufferBase<Depth>
	{
		FrameBuffer(uint32_t p_width, uint32_t p_height) :
			FrameBufferBase<Depth>(p_width, p_height)
		{
		}

		void Clear()
		{
			uint32_t depthMax = FloatToUint32(std::numeric_limits<float>::max());
			SIMDFill(reinterpret_cast<uint32_t*>(Data), Size, depthMax);
		}

		void Resize(uint32_t p_width, uint32_t p_height)
		{
			if (Width == p_width && Height == p_height)
				return;

			AlignedFree(Data);

			Width = p_width;
			Height = p_height;
			Size = Width * Height;

			Data = static_cast<Depth*>(AlignedAlloc(Size * sizeof(Depth), Alignment));
			Clear();
		}
	};

	using ColorBuffer = FrameBuffer<RGBA8>;
	using DepthBuffer = FrameBuffer<Depth>;
}
