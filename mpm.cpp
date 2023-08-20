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
const MPM::Vector kGravity(0.0f, 1.0f);

static float random_float() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<float> dis(0.0f, 1.0f);
    return dis(gen);
}

MPM::MPM(const std::shared_ptr<HalfEdgeMesh> &cut_mesh_)
    : cut_mesh_(cut_mesh_) {
    grid_nodes_.resize(cut_mesh_->vertices().size());
    for (auto [i, v] =
             std::tuple{static_cast<size_t>(0), cut_mesh_->vertices().begin()};
         i < grid_nodes_.size() && v != cut_mesh_->vertices().end(); ++i, ++v) {
        grid_nodes_[i].vertex = v;
    }
}

void MPM::initialize() {
    particles_.clear();
    particles_.resize(kNumberOfParticles);
    for (auto &p : particles_) {
        p.position.x() = random_float() * 0.4f + 0.2f;
        p.position.y() = random_float() * 0.4f + 0.2f;
    }
}

void MPM::update() {
    for (auto &g : grid_nodes_) {
        g.mass = 0.0f;
        g.velocity.setZero();
    }
    // P2G
    for (auto &p : particles_) {
        Vectori base = ((p.position * kInvDeltaX).array() - 0.5f).cast<int>();
        Vector fx = p.position * kInvDeltaX - base.cast<float>();
        Vector weights[3] = {0.5f * (1.5f - fx.array()).square(),
                             0.75f - (fx.array() - 1.0f).square(),
                             0.5f * (fx.array() - 0.5f).square()};
        p.deformation_gradient =
            (Matrix::Identity() + kDeltaT * p.affine_matrix) *
            p.deformation_gradient;
        auto hardening_coefficient =
            std::exp(10.0f * (1.0f - p.deformation_jacobian));
        hardening_coefficient = 0.3f;
        auto mu = kMu0 * hardening_coefficient;
        auto lambda = kLambda0 * hardening_coefficient;
        // mu = 0.0f;
        Eigen::JacobiSVD<Matrix, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(
            p.deformation_gradient);
        Vector sig = svd.singularValues();
        float plastic_deformation = 1.0f;
        for (int d = 0; d < 2; ++d) {
            auto new_sig = sig[d];
            p.deformation_jacobian *= sig[d] / new_sig;
            sig[d] = new_sig;
            plastic_deformation *= new_sig;
        }
        // p.deformation_gradient =
        //     Matrix::Identity() * std::sqrt(plastic_deformation);
        Matrix stress = 2.0f * mu *
                            (p.deformation_gradient -
                             svd.matrixU() * svd.matrixV().transpose()) *
                            p.deformation_gradient.transpose() +
                        Matrix::Identity() * lambda * plastic_deformation *
                            (plastic_deformation - 1.0f);
        stress = (-kDeltaT * kParticleVolume * 4.0f * kInvDeltaX * kInvDeltaX) *
                 stress;
        Matrix affine = stress + kParticleMass * p.affine_matrix;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                auto offset = Vectori(i, j);
                Vector distance = (offset.cast<float>() - fx) * kDeltaX;
                auto weight = weights[i][0] * weights[j][1];
                Vectori grid_idx = base + offset;
                auto &g = grid_nodes_[grid_idx[0] * kGridSize + grid_idx[1]];
                g.velocity +=
                    weight * (kParticleMass * p.velocity + affine * distance);
                g.mass += weight * kParticleMass;
            }
        }
    }
    // Grid update
    for (auto &g : grid_nodes_) {
        auto id = g.vertex->id;
        auto x = id / kGridSize;
        auto y = id % kGridSize;
        if (g.mass <= 0.0f) {
            continue;
        }
        g.velocity /= g.mass;
        g.velocity += kDeltaT * kGravity * 30.0f;
        if (x < 3 && g.velocity[0] < 0.0f) {
            g.velocity[0] = 0.0f;
        }
        if (x > kGridSize - 3 && g.velocity[0] > 0.0f) {
            g.velocity[0] = 0.0f;
        }
        if (y < 3 && g.velocity[1] < 0.0f) {
            g.velocity[1] = 0.0f;
        }
        if (y > kGridSize - 3 && g.velocity[1] > 0.0f) {
            g.velocity[1] = 0.0f;
        }
    }
    // G2P
    for (auto &p : particles_) {
        Vectori base = ((p.position * kInvDeltaX).array() - 0.5f).cast<int>();
        Vector fx = p.position * kInvDeltaX - base.cast<float>();
        Vector weights[3] = {0.5f * (1.5f - fx.array()).square(),
                             0.75f - (fx.array() - 1.0f).square(),
                             0.5f * (fx.array() - 0.5f).square()};
        Vector new_velocity = Vector::Zero();
        Matrix new_affine_velocity = Matrix::Zero();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                auto offset = Vectori(i, j);
                Vector distance = offset.cast<float>() - fx;
                Vectori grid_idx = base + offset;
                const auto &g =
                    grid_nodes_[grid_idx[0] * kGridSize + grid_idx[1]];
                auto weight = weights[i][0] * weights[j][1];
                new_velocity += weight * g.velocity;
                new_affine_velocity += 4.0f * kInvDeltaX * weight * g.velocity *
                                       distance.transpose();
            }
            p.velocity = new_velocity;
            p.affine_matrix = new_affine_velocity;
            p.position += kDeltaT * p.velocity;
        }
    }
}