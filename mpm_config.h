#pragma once

constexpr int kQuality = 8;
constexpr int kGridSize = 8 * kQuality;
constexpr int kGridRowSize = kGridSize + 1;
constexpr int kParticleNumber = 32 * kQuality * kQuality;
constexpr float kDeltaT = 2e-3f / kQuality;
constexpr float kDeltaX = 1.0f / kGridSize;
constexpr float kInvDeltaX = kGridSize;
