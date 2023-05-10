#include "mpm.h"
#include <Eigen/SVD>
#include <random>

constexpr float particle_volume = dx * dx * 0.25f;
constexpr float particle_density = 1.0f;
constexpr float particle_mass = particle_volume * particle_density;
constexpr float young_modulus = 1e3f;
constexpr float poisson_ratio = 0.2f;
constexpr float mu_0 = young_modulus / (2.0f * (1.0f + poisson_ratio));
constexpr float lambda_0 =
    young_modulus * poisson_ratio /
    ((1.0f + poisson_ratio) * (1.0f - 2.0f * poisson_ratio));

static float random_float() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<float> dis(0.0f, 1.0f);
    return dis(gen);
}

MPM::MPM(const std::shared_ptr<HalfEdgeMesh> &cut_mesh_) : cut_mesh(cut_mesh_) {
    grid.resize(cut_mesh->vertices.size());
    for (auto [i, v] =
             std::tuple{static_cast<size_t>(0), cut_mesh->vertices.begin()};
         i < grid.size() && v != cut_mesh->vertices.end(); ++i, ++v) {
        grid[i].vertex = v;
    }
}

void MPM::initialize() {
    particles.clear();
    particles.resize(n_particles);
    for (auto &p : particles) {
        p.position.x() = random_float() * 0.4f + 0.2f;
        p.position.y() = random_float() * 0.4f + 0.2f;
    }
}

void MPM::update() {
    for (auto &g : grid) {
        g.mass = 0.0f;
        g.velocity.setZero();
    }
    for (auto &p : particles) {
        Vectori base = ((p.position * inv_dx).array() - 0.5f).cast<int>();
        Vector fx = p.position * inv_dx - base.cast<float>();
        Vector weights[3] = {0.5f * (1.5f - fx.array()).square(),
                             0.75f - (fx.array() - 1.0f).square(),
                             0.5f * (0.5f - fx.array()).square()};
        p.deformation_gradient =
            (Matrix::Identity().array() + dt * p.affine_velocity.array())
                .matrix() *
            p.deformation_gradient;
        auto hardening_coefficient =
            std::exp(10.0f * (1.0f - p.plastic_deformation));
        // auto mu = mu_0 * hardening_coefficient;
        auto mu = 0.0f;
        auto lambda = lambda_0 * hardening_coefficient;
        Eigen::JacobiSVD<Matrix, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(
            p.deformation_gradient);
        Vector sig = svd.singularValues();
        float plastic_deformation = 1.0f;
        for (int d = 0; d < 2; ++d) {
            auto new_sig = sig[d];
            p.plastic_deformation *= sig[d] / new_sig;
            sig[d] = new_sig;
            plastic_deformation *= new_sig;
        }
        p.deformation_gradient =
            Matrix::Identity() * std::sqrt(plastic_deformation);
        Matrix stress = 2.0f * mu *
                            (p.deformation_gradient -
                             svd.matrixU() * svd.matrixV().transpose()) *
                            p.deformation_gradient.transpose() +
                        Matrix::Identity() * lambda * plastic_deformation *
                            (plastic_deformation - 1.0f);
        stress = (-dt * particle_volume * 4.0f * inv_dx * inv_dx) * stress;
        Matrix affine = stress + particle_mass * p.affine_velocity;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                auto offset = Vectori(i, j);
                Vector dpos = (offset.cast<float>() - fx) * dx;
                auto weight = weights[i][0] * weights[j][1];
                Vectori grid_idx = base + offset;
                auto &g = grid[grid_idx[0] * n_grid + grid_idx[1]];
                g.velocity +=
                    weight * (particle_mass * p.velocity + affine * dpos);
                g.mass += weight * particle_mass;
            }
        }
    }
}