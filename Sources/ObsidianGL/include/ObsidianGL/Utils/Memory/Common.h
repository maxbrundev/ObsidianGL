#pragma once

#ifdef _WIN32
#include <malloc.h>
#endif

namespace ObsidianGL::Utils::Memory
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
}
