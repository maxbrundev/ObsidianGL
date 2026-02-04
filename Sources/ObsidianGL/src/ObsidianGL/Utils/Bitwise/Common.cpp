#include "ObsidianGL/Utils/Bitwise/Common.h"

#include <memory>

uint32_t ObsidianGL::Utils::Bitwise::PackBytes(const uint8_t p_b3, const uint8_t p_b2, const uint8_t p_b1, const uint8_t p_b0)
{
	return (static_cast<uint32_t>(p_b3) << 24) | (static_cast<uint32_t>(p_b2) << 16) | (static_cast<uint32_t>(p_b1) << 8) | static_cast<uint32_t>(p_b0);
}

void ObsidianGL::Utils::Bitwise::UnpackBytes(const uint32_t p_value, uint8_t& p_b3, uint8_t& p_b2, uint8_t& p_b1, uint8_t& p_b0)
{
	p_b3 = static_cast<uint8_t>((p_value >> 24) & 0xFF);
	p_b2 = static_cast<uint8_t>((p_value >> 16) & 0xFF);
	p_b1 = static_cast<uint8_t>((p_value >> 8) & 0xFF);
	p_b0 = static_cast<uint8_t>(p_value & 0xFF);
}

uint32_t ObsidianGL::Utils::Bitwise::FloatToUint32(const float p_value)
{
	uint32_t bits;
	std::memcpy(&bits, &p_value, sizeof(float));
	return bits;
}

float ObsidianGL::Utils::Bitwise::Uint32ToFloat(const uint32_t p_bits)
{
	float value;
	std::memcpy(&value, &p_bits, sizeof(float));
	return value;
}
