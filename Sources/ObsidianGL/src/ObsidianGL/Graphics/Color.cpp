#include "ObsidianGL/Graphics/Color.h"

#include <glm/common.hpp>

#include "ObsidianGL/Utils/Bitwise/Common.h"

std::uint32_t ObsidianGL::Graphics::PackColor(const glm::vec4& p_color)
{
	return Utils::Bitwise::PackBytes(
		static_cast<uint8_t>(glm::clamp(p_color.r, 0.0f, 1.0f) * 255.0f),
		static_cast<uint8_t>(glm::clamp(p_color.g, 0.0f, 1.0f) * 255.0f),
		static_cast<uint8_t>(glm::clamp(p_color.b, 0.0f, 1.0f) * 255.0f),
		static_cast<uint8_t>(glm::clamp(p_color.a, 0.0f, 1.0f) * 255.0f)
	);
}

void ObsidianGL::Graphics::UnpackColorToFloat(uint32_t p_packed, float& p_r, float& p_g, float& p_b, float& p_a)
{
	constexpr float inv255 = 1.0f / 255.0f;

	p_r = static_cast<float>((p_packed >> 24) & 0xFF) * inv255;
	p_g = static_cast<float>((p_packed >> 16) & 0xFF) * inv255;
	p_b = static_cast<float>((p_packed >> 8) & 0xFF) * inv255;
	p_a = static_cast<float>(p_packed & 0xFF) * inv255;
}
