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

        Real weight_sum = 0;
        std::vector<Real> weights;
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

    int quality_;
    int grid_size_;
    int row_size_;
    int n_particles_;
    int material_;
    Real delta_t_;
    Real delta_x_;
    Real margin_;
    Real particle_volume_;
    Real particle_density_;
    Real particle_mass_;
    Real inv_delta_x_;

   private:
    std::shared_ptr<CutMesh> cut_mesh_;
    std::vector<Particle> particles_;
    std::vector<GridNode> nodes_;

    GridNode& grid(int x, int y) { return nodes_[y * row_size_ + x]; }
    [[nodiscard]] const GridNode& grid(int x, int y) const {
        return nodes_[y * row_size_ + x];
    }
};
