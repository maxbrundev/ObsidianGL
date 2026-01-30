#pragma once

#if defined(_WIN32)
	#pragma warning(disable : 4251)
	#pragma warning(disable : 4275)
	#if defined(OBSIDIANGL_BUILD)
		#define OBSIDIANGL_API __declspec(dllexport)
	#else
		#define OBSIDIANGL_API __declspec(dllimport)
	#endif
#elif defined(__GNUC__) && __GNUC__ >= 4
	#define OBSIDIANGL_API __attribute__((visibility("default")))
#else
	#define OBSIDIANGL_API
#endif