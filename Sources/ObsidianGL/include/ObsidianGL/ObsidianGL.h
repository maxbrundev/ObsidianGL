#pragma once

#include <cstdint>

#include "API/Export.h"

namespace ObsidianGL
{
	struct FrameBuffer
	{
		uint16_t Width;
		uint16_t Height;
		uint32_t* Pixels;
	};

	OBSIDIANGL_API void Test(const FrameBuffer* p_outBuffer);
}
