#pragma once

#include <Eigen/Dense>
#include <memory>
#include <vector>

#include "cut_mesh.h"
#include "mpm_config.h"

class MPM {
   public:
    using Vec2 = Eigen::Vector2f;
    using Vec2i = Eigen::Vector2i;
    using Vec3 = Eigen::Vector3f;
    using Mat2 = Eigen::Matrix2f;
    using Mat3 = Eigen::Matrix3f;
    using Mat23 = Eigen::Matrix<float, 2, 3>;

    struct Particle {
        Vec2 x = Vec2::Zero();
        Vec2 v = Vec2::Zero();
        Mat2 F = Mat2::Identity();
        Mat23 C = Mat23::Zero();
        Mat3 M_inv = Mat3::Zero();
    };

    struct GridNode {
        float m = 0.0f;
        Vec2 v = Vec2::Zero();
        HalfEdgeMesh::VertexRef vertex{};
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