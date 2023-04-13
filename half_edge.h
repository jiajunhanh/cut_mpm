#pragma once

#include <memory>

namespace hfeg {

class Edge;
class Vertex;
class Face;

class HalfEdge {
    std::weak_ptr<HalfEdge> twin;
    std::weak_ptr<HalfEdge> next;
    std::weak_ptr<Edge> edge;
    std::weak_ptr<Vertex> vertex;
    std::weak_ptr<Face> face;
};

class Edge {
    std::weak_ptr<HalfEdge> half_edge;
};

class Vertex {
    std::weak_ptr<HalfEdge> half_edge;
};

class Face {
    std::weak_ptr<HalfEdge> half_edge;
};

} // namespace hfeg