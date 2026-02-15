#pragma once

#include <cstdint>
#include <vector>

#include "ObsidianGL/Graphics/Texture.h"

namespace ObsidianGL::Managers
{
	class TextureManager
	{
	public:
		TextureManager() = default;
		~TextureManager();

		uint32_t Generate();

		void Delete(uint32_t p_id);

		Graphics::Texture* Get(uint32_t p_id) const;

	private:
		std::vector<Graphics::Texture*> m_textures;
		uint32_t m_nextID = 1;
	};
}
