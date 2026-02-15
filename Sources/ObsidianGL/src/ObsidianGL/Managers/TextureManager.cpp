#include "ObsidianGL/Managers/TextureManager.h"

ObsidianGL::Managers::TextureManager::~TextureManager()
{
	for (auto* texture : m_textures)
	{
		delete texture;
	}
}

uint32_t ObsidianGL::Managers::TextureManager::Generate()
{
	uint32_t id = m_nextID++;

	if (id >= m_textures.size())
	{
		m_textures.resize(id + 1, nullptr);
	}

	m_textures[id] = new Graphics::Texture();

	return id;
}

void ObsidianGL::Managers::TextureManager::Delete(uint32_t p_id)
{
	if (p_id > 0 && p_id < m_textures.size() && m_textures[p_id])
	{
		delete m_textures[p_id];
		m_textures[p_id] = nullptr;
	}
}

ObsidianGL::Graphics::Texture* ObsidianGL::Managers::TextureManager::Get(uint32_t p_id) const
{
	if (p_id > 0 && p_id < m_textures.size())
	{
		return m_textures[p_id];
	}

	return nullptr;
}
