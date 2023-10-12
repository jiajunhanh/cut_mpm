#pragma once

#include <Eigen/Dense>
#include <memory>
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
        Mat2 C = Mat2::Zero();
    };

    struct GridNode {
        Real m = 0.0;
        Vec2 v = Vec2::Zero();
        CutMesh::VertexRef vertex;
        std::vector<CutMesh::FaceRef> faces;
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

    std::vector<bool> face_visited_;
    std::vector<bool> is_neighbor_node_;
    std::vector<CutMesh::FaceRef> neighbor_faces_;
    std::queue<CutMesh::FaceRef> face_queue_;

    std::vector<int> get_neighbor_nodes(const Vec2& x);
    void get_neighbor_nodes_in_place(const Vec2& x,
                                     std::vector<int>& neighbor_nodes);
};
