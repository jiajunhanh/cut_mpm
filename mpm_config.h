#pragma once

using Real = double;
// using Real = float;

constexpr int kQuality = 8;
constexpr int kGridSize = 8 * kQuality;
constexpr int kRowSize = kGridSize + 1;
constexpr int kParticleNumber = 32 * kQuality * kQuality;
constexpr Real kDeltaT = 2e-3 / kQuality;
constexpr Real kDeltaX = 1.0 / kGridSize;
constexpr Real kKernelRange = 1.5;
constexpr Real kMargin = kDeltaX / 64;
constexpr int kValidBit = 30;
// constexpr Real kInvDeltaX = kGridSize;

// using Vec2i = Eigen::Vector2i;
using Mat23 = Eigen::Matrix<Real, 2, 3>;
using Vec2 = std::conditional_t<std::is_same_v<Real, float>, Eigen::Vector2f,
                                Eigen::Vector2d>;
using Vec3 = std::conditional_t<std::is_same_v<Real, float>, Eigen::Vector3f,
                                Eigen::Vector3d>;
using Mat2 = std::conditional_t<std::is_same_v<Real, float>, Eigen::Matrix2f,
                                Eigen::Matrix2d>;
using Mat3 = std::conditional_t<std::is_same_v<Real, float>, Eigen::Matrix3f,
                                Eigen::Matrix3d>;
