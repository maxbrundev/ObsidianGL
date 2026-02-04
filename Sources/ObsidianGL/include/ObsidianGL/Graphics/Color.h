#pragma once

#include <cstdint>

#include <glm/common.hpp>
#include <glm/vec4.hpp>

namespace ObsidianGL::Graphics
{
	std::uint32_t PackColor(const glm::vec4& p_color);
}