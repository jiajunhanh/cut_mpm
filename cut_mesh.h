#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include <unordered_set>
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

        int id = 0;
        bool is_boundary = false;
        Vec2 normal = Vec2::Zero();

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
        int id = 0;
    };

    struct Vertex {
        HalfEdgeRef half_edge;
        int id = 0;
        // int grid_id = 0;
        Vec2 position{};

        bool on_boundary = false;
        Vec2 normal = Vec2::Zero();
    };

    struct Face {
        HalfEdgeRef half_edge;
        int id = 0;
        Vec2 center = Vec2::Zero();
        std::vector<int> neighbor_nodes;
        std::vector<bool> neighbor_node_sides;
        std::vector<int> neighbor_boundaries;
        void calculate_center();
        [[nodiscard]] bool enclose(const Vec2& x) const;
    };

    struct Grid {
        VertexRef vertex;
        std::vector<FaceRef> faces;
        bool near_boundary = false;
    };

    CutMesh() = default;
    CutMesh(std::vector<std::array<Real, 2>> vertices, int quality);

    auto& half_edges() { return half_edges_; }
    auto& vertices() { return vertices_; }
    auto& edges() { return edges_; }
    auto& faces() { return faces_; }
    auto& grids() { return grids_; }
    Grid& grid(int r, int c) { return grids_[r * grid_size_ + c]; }
    [[nodiscard]] const auto& half_edges() const { return half_edges_; }
    [[nodiscard]] const auto& vertices() const { return vertices_; }
    [[nodiscard]] const auto& edges() const { return edges_; }
    [[nodiscard]] const auto& faces() const { return faces_; }
    [[nodiscard]] const auto& grids() const { return grids_; }
    [[nodiscard]] const Grid& grid(int r, int c) const {
        return grids_[r * grid_size_ + c];
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
    void calculate_neighbor_nodes_and_boundaries_of_faces();
    void calculate_node_normals();

   private:
    struct CutVertex {
        int c[2] = {};
        Real q[2] = {};
        bool b[2] = {};

        [[nodiscard]] Real coord(int d) const {
            return static_cast<Real>(c[d]) + q[d];
        }
    };

    struct CutEdge {
       public:
        int i = 0;
        int j = 0;
        bool is_boundary = false;
        CutEdge(int i_, int j_, bool is_boundary_)
            : i(i_), j(j_), is_boundary(is_boundary_) {}
    };

    int grid_size_;
    int row_size_;
    int n_grid_nodes_;
    Real delta_x_;
    Real margin_;

    std::vector<HalfEdge> half_edges_;
    std::vector<Vertex> vertices_;
    std::vector<Edge> edges_;
    std::vector<Face> faces_;
    std::vector<Grid> grids_;
    // std::stack<HalfEdgeRef> recycled_half_edges_;
    // std::stack<VertexRef> recycled_vertices_;
    // std::stack<EdgeRef> recycled_edges_;
    // std::stack<FaceRef> recycled_faces_;

    void lerp(CutVertex& v, const CutVertex& vi, const CutVertex& vj, Real t,
              int d) const;
    void add_grid_edges(std::vector<CutVertex>& cut_vertices,
                        std::vector<CutEdge>& cut_edges,
                        const std::unordered_set<uint64_t>& edge_hash) const;
    std::pair<std::vector<CutVertex>, std::vector<CutEdge>>
    compute_cut_vertices_and_edges(
        const std::vector<std::array<Real, 2>>& vertices);
    [[nodiscard]] bool is_neighbor_face(CutMesh::FaceRef face0,
                                        CutMesh::FaceRef face1) const;
};
