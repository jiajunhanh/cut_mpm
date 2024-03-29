#pragma once

using Real = float;

// constexpr int kQuality = 8;
// constexpr int kGridSize = 8 * kQuality;
// constexpr int kRowSize = kGridSize + 1;
// constexpr int kParticleNumber = 32 * kQuality * kQuality;
// constexpr Real kDeltaT = 2e-3 / kQuality;
// constexpr Real kDeltaX = 1.0 / kGridSize;
constexpr Real kKernelRange = 1.5;
constexpr Real kInf = std::numeric_limits<Real>::infinity();
// constexpr Real kMargin = kDeltaX / 32;
//  constexpr Real kInvDeltaX = kGridSize;

using Vec2i = Eigen::Vector2i;
using Mat23 = Eigen::Matrix<Real, 2, 3>;
using Vec2 = std::conditional_t<std::is_same_v<Real, float>, Eigen::Vector2f,
                                Eigen::Vector2d>;
using Vec3 = std::conditional_t<std::is_same_v<Real, float>, Eigen::Vector3f,
                                Eigen::Vector3d>;
using Mat2 = std::conditional_t<std::is_same_v<Real, float>, Eigen::Matrix2f,
                                Eigen::Matrix2d>;
using Mat3 = std::conditional_t<std::is_same_v<Real, float>, Eigen::Matrix3f,
                                Eigen::Matrix3d>;

inline Real interpolate(Real x) {
    x = std::abs(x);
    if (x < 0.5) return Real{0.75} - x * x;
    if (x < 1.5) return Real{0.5} * (Real{1.5} - x) * (Real{1.5} - x);
    return 0.0;
    // return std::max(Real{0}, Real{1.5} - x);
}

inline Real interpolate(Vec2 x) {
    return interpolate(x.x()) * interpolate(x.y());
}

inline Vec2 project(const Vec2& normal, const Vec2& v) {
    return v - std::min(Real{0}, normal.dot(v)) * normal;
}