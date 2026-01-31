#include "ObsidianGL/ObsidianGL.h"

void ObsidianGL::Test(Buffers::FrameBuffer<RGBA8>& p_outBuffer)
{
	const uint32_t width = p_outBuffer.Width;
	const uint32_t height = p_outBuffer.Height;
	const float invHeight = 1.0f / static_cast<float>(height - 1);
	const float invWidth = 1.0f / static_cast<float>(width - 1);

	for (uint16_t y = 0; y < p_outBuffer.Height; ++y)
	{
		const float v = static_cast<float>(y) * invHeight;

		uint8_t  g = static_cast<uint8_t>(v * 255.0f);

		RGBA8* row = p_outBuffer.GetRow(y);

		for (uint16_t x = 0; x < p_outBuffer.Width; ++x)
		{
			const float u = static_cast<float>(x) * invWidth;

			uint8_t  r = static_cast<uint8_t>(u * 255.0f);
			uint8_t  b = static_cast<uint8_t>((1.0f - u) * 255.0f);

			row[x] = (r << 24) | (g << 16) | (b << 8) | 255;
		}
	}
}
