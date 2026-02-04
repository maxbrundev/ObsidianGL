#pragma once

#include <cstdint>

namespace ObsidianGL::Utils::Bitwise
{
	uint32_t PackBytes(const uint8_t p_r, const uint8_t p_g, const uint8_t p_b, const uint8_t p_a);

	void UnpackBytes(const uint32_t p_value, uint8_t& p_b3, uint8_t& p_b2, uint8_t& p_b1, uint8_t& p_b0);

	uint32_t FloatToUint32(const float p_value);

	float Uint32ToFloat(const uint32_t p_bits);
}
