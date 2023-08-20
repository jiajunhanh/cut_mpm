#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <list>
#include <stack>
#include <utility>
#include <vector>

class HalfEdgeMesh {
   public:
    struct HalfEdge;
    struct Vertex;
    struct Edge;
    struct Face;

    using HalfEdgeRef = std::list<HalfEdge>::iterator;
    using VertexRef = std::list<Vertex>::iterator;
    using EdgeRef = std::list<Edge>::iterator;
    using FaceRef = std::list<Face>::iterator;

    struct HalfEdge {
        HalfEdgeRef twin;
        HalfEdgeRef next;
        EdgeRef edge;
        VertexRef vertex;
        FaceRef face;

        int id{};

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
        Eigen::Vector2f position{};
    };

    struct Face {
        HalfEdgeRef half_edge;
        int id{};
    };

    auto& half_edges() { return half_edges_; }
    auto& vertices() { return vertices_; }
    auto& edges() { return edges_; }
    auto& faces() { return faces_; }
    [[nodiscard]] const auto& half_edges() const { return half_edges_; }
    [[nodiscard]] const auto& vertices() const { return vertices_; }
    [[nodiscard]] const auto& edges() const { return edges_; }
    [[nodiscard]] const auto& faces() const { return faces_; }

    HalfEdgeRef emplace_half_edge();
    VertexRef emplace_vertex();
    EdgeRef emplace_edge();
    FaceRef emplace_face();

    void erase_half_edge(HalfEdgeRef&& h);
    void erase_vertex(VertexRef&& v);
    void erase_edge(EdgeRef&& e);
    void erase_face(FaceRef&& f);

   private:
    int next_id{};
    std::list<HalfEdge> half_edges_;
    std::list<Vertex> vertices_;
    std::list<Edge> edges_;
    std::list<Face> faces_;
    std::stack<HalfEdgeRef> recycled_half_edges_;
    std::stack<VertexRef> recycled_vertices_;
    std::stack<EdgeRef> recycled_edges_;
    std::stack<FaceRef> recycled_faces_;
};