#pragma once

#include <Eigen/Dense>
#include <memory>
#include <vector>

#include "cut_mesh.h"
#include "mpm_config.h"

class MPM {
   public:
    using Vector = Eigen::Vector2f;
    using Vectori = Eigen::Vector2i;
    using Matrix = Eigen::Matrix2f;

    struct Particle {
        Vector position{};
        Vector velocity{};
        Matrix deformation_gradient = Matrix::Identity();
        Matrix affine_velocity = Matrix::Identity();
        float plastic_deformation = 1.0f;
    };

    struct Grid {
        float mass = 0.0f;
        Vector velocity{};
        HalfEdgeMesh::VertexRef vertex{};
    };

    std::shared_ptr<HalfEdgeMesh> cut_mesh;
    std::vector<Particle> particles;
    std::vector<Grid> grids;

    explicit MPM(const std::shared_ptr<HalfEdgeMesh> &cut_mesh_);
    void initialize();
    void update();
};