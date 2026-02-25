# ObsidianGL

# Description
Real-time SIMD-optimized CPU software renderer with a modern API
A complete rewrite of [AmberGL](https://github.com/maxbrundev/Rasterizer) from scratch, designed to also provide an immediate mode capable of rendering Minecraft Classic 0.30 (OpenGL 1.1), for science.

# Features
- OpenGL 1.1 immediate-mode API (Begin/End, matrix stacks, texture binding...)
- Triangle rasterization (barycentric, subpixel-accurate with top-left fill rule, no coplanar triangle issues)
- Line rasterization (configurable width)
- Point rasterization (configurable size)
- Three SIMD paths: AVX2 (8px), SSE/SSE4.1 (4px), scalar
- SIMD-accelerated buffer clears
- Texture mapping (nearest + bilinear filtering, repeat + clamp-to-edge wrapping)
- SIMD batch texture sampling
- Perspective-correct interpolation (colors, texture coordinates)
- Sutherland-Hodgman frustum clipping (6 planes, Cohen-Sutherland trivial reject)
- 32-bit float depth buffer (depth test / depth write)
- Face culling (front/back, CW/CCW)
- Fill, wireframe, and point polygon modes
- Perspective / orthographic projection

Coming soon:
- Tile-based multithreaded rasterization
- Tile debug overlay (heat map, grid, per-tile visibility control)

# Video and Screenshots

Minecraft 0.30 running on CPU using ObsidianGL and a Java bridge. I removed all GL11 function calls and now just blit the ObsidianGL framebuffer to OpenGL when Display is called. A lot of rendering issues remain, but at least it runs!

https://github.com/user-attachments/assets/f22f4969-07eb-465c-918e-3db4355673d9

[![perf showcase](https://img.youtube.com/vi/gS9VbkYYQCY/maxresdefault.jpg?cache=2026)](https://www.youtube.com/watch?v=gS9VbkYYQCY)

<img width="802" height="632" alt="image" src="https://github.com/user-attachments/assets/ab2821b5-c3cb-417e-a119-fc1fbd323f59" />
<img width="802" height="632" alt="Screenshot 2026-02-16 205801" src="https://github.com/user-attachments/assets/40bdd92c-c0f3-4c04-9f1c-7b68e9c839f9" />
<img width="802" height="632" alt="Screenshot 2026-02-16 205822" src="https://github.com/user-attachments/assets/72b73707-0f49-41f4-9878-00ccb6b11f03" />

*Amber preserves, Onyx structures, Obsidian cuts.*
