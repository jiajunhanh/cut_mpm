#pragma once

#include <iostream>
#include <list>

class HalfEdgeMesh {
    size_t next_id{};

  public:
    class HalfEdge;
    class Vertex;
    class Edge;
    class Face;

    using HalfEdgeRef = std::list<HalfEdge>::iterator;
    using VertexRef = std::list<Vertex>::iterator;
    using EdgeRef = std::list<Edge>::iterator;
    using FaceRef = std::list<Face>::iterator;

    class HalfEdge {
      public:
        HalfEdgeRef twin;
        HalfEdgeRef next;
        EdgeRef edge;
        VertexRef vertex;
        FaceRef face;

        size_t id{};

        void set_tnvef(HalfEdgeRef twin_, HalfEdgeRef next_, VertexRef vertex_,
                       EdgeRef edge_, FaceRef face_) {
            twin = twin_;
            next = next_;
            vertex = vertex_;
            edge = edge_;
            face = face_;
        }
    };

    class Edge {
      public:
        HalfEdgeRef half_edge;
        size_t id{};
    };

    class Vertex {
      public:
        HalfEdgeRef half_edge;
        size_t id{};
    };

    class Face {
      public:
        HalfEdgeRef half_edge;
        size_t id{};
    };

    std::list<HalfEdge> half_edges;
    std::list<Vertex> vertices;
    std::list<Edge> edges;
    std::list<Face> faces;

    HalfEdgeRef emplace_half_edge();
    VertexRef emplace_vertex();
    EdgeRef emplace_edge();
    FaceRef emplace_face();

    void free_half_edge(HalfEdgeRef h);
    void free_vertex(VertexRef v);
    void free_edge(EdgeRef e);
    void free_face(FaceRef f);
};