#pragma once

using Real = double;

constexpr int kQuality = 8;
constexpr int kGridSize = 8 * kQuality;
constexpr int kGridRowSize = kGridSize + 1;
constexpr int kParticleNumber = 32 * kQuality * kQuality;
constexpr Real kDeltaT = 2e-3 / kQuality;
constexpr Real kDeltaX = 1.0 / kGridSize;
constexpr Real kInvDeltaX = kGridSize;
