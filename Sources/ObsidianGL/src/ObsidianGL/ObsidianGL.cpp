#include "ObsidianGL/ObsidianGL.h"

void ObsidianGL::Test(const FrameBuffer* p_outBuffer)
{
	const uint16_t width = p_outBuffer->Width;
	const uint16_t height = p_outBuffer->Height;

	for (uint16_t y = 0; y < height; ++y)
	{
		const float v = static_cast<float>(y) / static_cast<float>(height - 1);

		for (uint16_t x = 0; x < width; ++x)
		{
			const float u = static_cast<float>(x) / static_cast<float>(width - 1);

			uint8_t  r = static_cast<uint8_t>(u * 255.0f);
			uint8_t  g = static_cast<uint8_t>(v * 255.0f);
			uint8_t  b = static_cast<uint8_t>((1.0f - u) * 255.0f);
			uint8_t  a = 255;

			p_outBuffer->Pixels[y * width + x] = (r << 24) | (g << 16) | (b << 8) | a;
		}
	}
}
