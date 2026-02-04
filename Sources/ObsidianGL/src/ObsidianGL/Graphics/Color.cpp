#include "ObsidianGL/Graphics/Color.h"

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
