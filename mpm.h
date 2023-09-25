#pragma once

#include <Eigen/Dense>
#include <memory>
#include <vector>

#include "cut_mesh.h"
#include "mpm_config.h"

class MPM {
   public:
    using Vec2i = Eigen::Vector2i;
    using Mat23 = Eigen::Matrix<Real, 2, 3>;
    using Vec2 = std::conditional_t<std::is_same_v<Real, float>,
                                    Eigen::Vector2f, Eigen::Vector2d>;
    using Vec3 = std::conditional_t<std::is_same_v<Real, float>,
                                    Eigen::Vector3f, Eigen::Vector3d>;
    using Mat2 = std::conditional_t<std::is_same_v<Real, float>,
                                    Eigen::Matrix2f, Eigen::Matrix2d>;
    using Mat3 = std::conditional_t<std::is_same_v<Real, float>,
                                    Eigen::Matrix3f, Eigen::Matrix3d>;

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
        CutMesh::VertexRef vertex{};
        std::vector<CutMesh::FaceRef> faces{};
    };

    explicit MPM(const std::shared_ptr<CutMesh>& cut_mesh_);
    void initialize();
    void update();
    [[nodiscard]] const std::vector<Particle>& particles() const {
        return particles_;
    }

   private:
    std::shared_ptr<CutMesh> cut_mesh_;
    std::vector<Particle> particles_;
    std::vector<GridNode> nodes_;
};
