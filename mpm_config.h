#pragma once

constexpr int quality = 16;
constexpr int n_grid = 8 * quality;
constexpr int n_particles = 32 * quality * quality;
constexpr float dt = 2e-3 / quality;
constexpr float dx = 1.0f / n_grid;
constexpr float inv_dx = n_grid;
