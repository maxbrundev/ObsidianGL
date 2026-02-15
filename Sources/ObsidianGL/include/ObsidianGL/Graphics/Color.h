#pragma once

#include <cstdint>

#include <glm/vec4.hpp>

namespace ObsidianGL::Graphics
{
	std::uint32_t PackColor(const glm::vec4& p_color);

	void UnpackColorToFloat(uint32_t p_packed, float& p_r, float& p_g, float& p_b, float& p_a);
}