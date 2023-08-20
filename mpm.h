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
        Vector position;
        Vector velocity;
        Matrix deformation_gradient = Matrix::Identity();
        float deformation_jacobian = 1.0f;
        Matrix affine_matrix = Matrix::Identity();
    };

    struct GridNode {
        float mass = 0.0f;
        Vector velocity;
        HalfEdgeMesh::VertexRef vertex;
    };

    explicit MPM(const std::shared_ptr<HalfEdgeMesh>& cut_mesh_);
    void initialize();
    void update();
    [[nodiscard]] const std::vector<Particle>& particles() const {
        return particles_;
    }

   private:
    std::shared_ptr<HalfEdgeMesh> cut_mesh_;
    std::vector<Particle> particles_;
    std::vector<GridNode> grid_nodes_;
};