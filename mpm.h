#pragma once

#include <Eigen/Dense>
#include <queue>
#include <vector>

#include "cut_mesh.h"
#include "mpm_config.h"

class MPM {
   public:
    struct Particle {
        Vec2 x = Vec2::Zero();
        Vec2 v = Vec2::Zero();
        Mat2 F = Mat2::Identity();
        Mat23 C = Mat23::Zero();
        Mat3 M_inv = Mat3::Zero();
    };

    struct GridNode {
        Real m = 0.0;
        Vec2 v = Vec2::Zero();
        CutMesh::VertexRef vertex;
        std::vector<CutMesh::FaceRef> faces;
    };

    MPM(const std::shared_ptr<CutMesh>& cut_mesh, int quality, int material);
    void initialize();
    void update();
    [[nodiscard]] const std::vector<Particle>& particles() const {
        return particles_;
    }

   private:
    int quality_ = 8;
    int grid_size_ = 8 * quality_;
    int row_size_ = grid_size_ + 1;
    int n_particles_ = 32 * quality_ * quality_;
    int material_ = 0;
    Real delta_t_ = Real{2e-3} / static_cast<Real>(quality_);
    Real delta_x_ = Real{1.0} / static_cast<Real>(grid_size_);
    Real margin_ = delta_x_ / 32;
    Real particle_volume_ = delta_x_ * delta_x_ * Real{0.25};
    Real particle_density_ = 1.0;
    Real k_particle_mass_ = particle_volume_ * particle_density_;

    std::shared_ptr<CutMesh> cut_mesh_;
    std::vector<Particle> particles_;
    std::vector<GridNode> nodes_;
};
