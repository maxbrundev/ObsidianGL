#pragma once

// #define OBSIDIAN_FORCE_SCALAR
// #define OBSIDIAN_FORCE_SSE
// #define OBSIDIAN_FORCE_SSE41
// #define OBSIDIAN_FORCE_AVX2

#if defined(OBSIDIAN_FORCE_AVX2)
#define OBSIDIAN_USE_AVX2
#elif defined(OBSIDIAN_FORCE_SSE41)
#define OBSIDIAN_USE_SSE41
#elif defined(OBSIDIAN_FORCE_SSE)
#define OBSIDIAN_USE_SSE
#elif defined(OBSIDIAN_FORCE_SCALAR)
#define OBSIDIAN_USE_SCALAR

#elif defined(__AVX2__)
#define OBSIDIAN_USE_AVX2
#elif defined(__SSE4_1__)
#define OBSIDIAN_USE_SSE41
#else
#define OBSIDIAN_USE_SSE
#endif

#if defined(OBSIDIAN_USE_AVX2)
#include <immintrin.h>
#elif defined(OBSIDIAN_USE_SSE41)
#include <smmintrin.h>
#elif defined(OBSIDIAN_USE_SSE)
#include <emmintrin.h>
#endif
