#pragma once

#include <Eigen/Dense>
#include <memory>
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
    std::unordered_map<int, int> node_of_vertex_;

    [[nodiscard]] std::vector<int> get_neighbor_nodes(const Vec2& x) const;
};
