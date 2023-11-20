#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "mpm_config.h"

class CutMesh {
   public:
    struct HalfEdge;
    struct Vertex;
    struct Edge;
    struct Face;

    using HalfEdgeRef = std::vector<HalfEdge>::iterator;
    using VertexRef = std::vector<Vertex>::iterator;
    using EdgeRef = std::vector<Edge>::iterator;
    using FaceRef = std::vector<Face>::iterator;

    struct HalfEdge {
        HalfEdgeRef twin;
        HalfEdgeRef next;
        EdgeRef edge;
        VertexRef vertex;
        FaceRef face;

        int id{};
        bool is_boundary{};

        void set_tnvef(const HalfEdgeRef& twin_, const HalfEdgeRef& next_,
                       const VertexRef& vertex_, const EdgeRef& edge_,
                       const FaceRef& face_) {
            twin = twin_;
            next = next_;
            vertex = vertex_;
            edge = edge_;
            face = face_;
        }
    };

    struct Edge {
        HalfEdgeRef half_edge;
        int id{};
    };

    struct Vertex {
        HalfEdgeRef half_edge;
        int id{};
        int grid_id{};
        Vec2 position{};
        std::array<bool, 2> on_edge{};
    };

    struct Face {
        HalfEdgeRef half_edge;
        int id{};
        Vec2 center = Vec2::Zero();
        std::vector<int> neighbor_nodes;
        std::vector<int> neighbor_node_sides;
        void calculate_center();
        [[nodiscard]] bool enclose(const Vec2& x) const;
    };

    struct Grid {
        VertexRef vertex;
        std::vector<FaceRef> faces;
    };

    auto& half_edges() { return half_edges_; }
    auto& vertices() { return vertices_; }
    auto& edges() { return edges_; }
    auto& faces() { return faces_; }
    auto& grids() { return grids_; }
    Grid& grid(int r, int c) { return grids_[r * kGridSize + c]; }
    [[nodiscard]] const auto& half_edges() const { return half_edges_; }
    [[nodiscard]] const auto& vertices() const { return vertices_; }
    [[nodiscard]] const auto& edges() const { return edges_; }
    [[nodiscard]] const auto& faces() const { return faces_; }
    [[nodiscard]] const auto& grids() const { return grids_; }
    [[nodiscard]] const Grid& grid(int r, int c) const {
        return grids_[r * kGridSize + c];
    }

    HalfEdgeRef emplace_half_edge();
    VertexRef emplace_vertex();
    EdgeRef emplace_edge();
    FaceRef emplace_face();

    // void erase_half_edge(HalfEdgeRef&& h);
    // void erase_vertex(VertexRef&& v);
    // void erase_edge(EdgeRef&& e);
    // void erase_face(FaceRef&& f);

    [[nodiscard]] FaceRef get_enclosing_face(Vec2 center) const;
    void calculate_neighbor_nodes_of_faces();

   private:
    std::vector<HalfEdge> half_edges_;
    std::vector<Vertex> vertices_;
    std::vector<Edge> edges_;
    std::vector<Face> faces_;
    std::vector<Grid> grids_;
    // std::stack<HalfEdgeRef> recycled_half_edges_;
    // std::stack<VertexRef> recycled_vertices_;
    // std::stack<EdgeRef> recycled_edges_;
    // std::stack<FaceRef> recycled_faces_;
};

CutMesh construct_cut_mesh(const std::vector<std::array<Real, 2>>& vertices,
                           const std::vector<std::array<int, 2>>& edges);
