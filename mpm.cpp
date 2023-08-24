#include "mpm.h"

#include <Eigen/SVD>
#include <random>

constexpr float kParticleVolume = kDeltaX * kDeltaX * 0.25f;
constexpr float kParticleDensity = 1.0f;
constexpr float kParticleMass = kParticleVolume * kParticleDensity;
constexpr float kYoungModulus = 1e3f;
constexpr float kPoissonRatio = 0.2f;
constexpr float kMu0 = kYoungModulus / (2.0f * (1.0f + kPoissonRatio));
constexpr float kLambda0 =
    kYoungModulus * kPoissonRatio /
    ((1.0f + kPoissonRatio) * (1.0f - 2.0f * kPoissonRatio));
const MPM::Vec2 kGravity(0.0f, 1.0f);

static float random_float() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<float> dis(0.0f, 1.0f);
    return dis(gen);
}

static auto svd(const MPM::Mat2 &m) {
    Eigen::JacobiSVD<MPM::Mat2, Eigen::ComputeFullU | Eigen::ComputeFullV> svd(
        m);
    return std::tuple{svd.matrixU(), svd.singularValues(), svd.matrixV()};
}

MPM::MPM(const std::shared_ptr<HalfEdgeMesh> &cut_mesh_)
    : cut_mesh_(cut_mesh_) {
    int n_vertices = static_cast<int>(cut_mesh_->vertices().size());
    grid_nodes_.resize(n_vertices);
    for (auto [i, v] = std::tuple(0, begin(cut_mesh_->vertices()));
         i < n_vertices; ++i, ++v) {
        grid_nodes_[i].vertex = v;
    }
}

void MPM::initialize() {
    particles_.clear();
    particles_.resize(kParticleNumber);
    for (auto &p : particles_) {
        p.x.x() = random_float() * 0.35f + 0.2f;
        p.x.y() = random_float() * 0.35f + 0.2f;
    }
}

void MPM::update() {
    for (auto &g : grid_nodes_) {
        g.m = 0.0f;
        g.v.setZero();
    }
    // P2G
    for (auto &p : particles_) {
        Vec2i base = ((p.x * kInvDeltaX).array() - 0.5f).cast<int>();
        Vec2 fx = p.x * kInvDeltaX - base.cast<float>();
        Vec2 weights[3] = {0.5f * (1.5f - fx.array()).square(),
                           0.75f - (fx.array() - 1.0f).square(),
                           0.5f * (fx.array() - 0.5f).square()};
        p.F = (Mat2::Identity() + kDeltaT * p.C) * p.F;
        auto hardening_coefficient = 0.5f;
        auto mu = kMu0 * hardening_coefficient;
        auto lambda = kLambda0 * hardening_coefficient;
        auto [U, sig, V] = svd(p.F);
        float J = p.F.determinant();
        Mat2 PF = 2.0f * mu * (p.F - U * V.transpose()) * p.F.transpose() +
                  Mat2::Identity() * lambda * J * (J - 1.0f);
        Mat3 M_inv = p.M_inv;
        Mat23 stress = Mat23::Zero();
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j)
                stress.row(i) += M_inv.row(j + 1) * PF(i, j);
        }
        stress *= -kDeltaT * kParticleVolume;
        Mat2 affine = kParticleMass * p.C;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                auto offset = Vec2i(i, j);
                Vec2 distance = (offset.cast<float>() - fx);
                auto weight = weights[i][0] * weights[j][1];
                distance *= kDeltaX;
                int x = base.x() + i;
                int y = base.y() + j;
                auto &g = grid_nodes_[y * (kGridSize + 1) + x];
                g.v += weight * (kParticleMass * p.v + affine * distance);
                g.v += weight * stress * Vec3(1.0f, distance.x(), distance.y());
                g.m += weight * kParticleMass;
            }
        }
    }
    // Grid update
    for (auto &g : grid_nodes_) {
        auto id = g.vertex->id;
        auto x = id % (kGridSize + 1);
        auto y = id / (kGridSize + 1);
        if (g.m <= 0.0f) {
            continue;
        }
        g.v /= g.m;
        g.v += kDeltaT * kGravity * 30.0f;
        if (x < 3 && g.v[0] < 0.0f) {
            g.v[0] = 0.0f;
        }
        if (x > kGridSize - 3 && g.v[0] > 0.0f) {
            g.v[0] = 0.0f;
        }
        if (y < 3 && g.v[1] < 0.0f) {
            g.v[1] = 0.0f;
        }
        if (y > kGridSize - 3 && g.v[1] > 0.0f) {
            g.v[1] = 0.0f;
        }
    }
    // G2P
    for (auto &p : particles_) {
        Vec2i base = ((p.x * kInvDeltaX).array() - 0.5f).cast<int>();
        Vec2 fx = p.x * kInvDeltaX - base.cast<float>();
        Vec2 weights[3] = {0.5f * (1.5f - fx.array()).square(),
                           0.75f - (fx.array() - 1.0f).square(),
                           0.5f * (fx.array() - 0.5f).square()};
        Vec2 new_v = Vec2::Zero();
        Mat2 new_C = Mat2::Zero();
        Mat3 new_M = Mat3::Zero();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                auto offset = Vec2i(i, j);
                Vec2 distance = (offset.cast<float>() - fx) * kDeltaX;
                int x = base.x() + i;
                int y = base.y() + j;
                const auto &g = grid_nodes_[y * (kGridSize + 1) + x];
                auto weight = weights[i][0] * weights[j][1];
                new_v += weight * g.v;
                new_C += weight * (g.v * distance.transpose());
                Vec3 P(1.0f, distance.x(), distance.y());
                new_M += weight * P * P.transpose();
            }
        }
        p.v = new_v;
        p.M_inv = new_M.inverse();
        p.C = new_C * new_M.block<2, 2>(1, 1).inverse();
        p.x += kDeltaT * p.v;
    }
}