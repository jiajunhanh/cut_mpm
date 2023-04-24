#include "half_edge.h"
HalfEdgeMesh::HalfEdgeRef HalfEdgeMesh::emplace_half_edge() {
    half_edges.emplace_back();
    auto h = std::prev(half_edges.end());
    h->twin = half_edges.end();
    h->next = half_edges.end();
    h->edge = edges.end();
    h->vertex = vertices.end();
    h->face = faces.end();
    h->id = next_id++;
    return h;
}

HalfEdgeMesh::VertexRef HalfEdgeMesh::emplace_vertex() {
    vertices.emplace_back();
    auto v = std::prev(vertices.end());
    v->half_edge = half_edges.end();
    v->id = next_id++;
    return v;
}

HalfEdgeMesh::EdgeRef HalfEdgeMesh::emplace_edge() {
    edges.emplace_back();
    auto e = std::prev(edges.end());
    e->half_edge = half_edges.end();
    e->id = next_id++;
    return e;
}

HalfEdgeMesh::FaceRef HalfEdgeMesh::emplace_face() {
    faces.emplace_back();
    auto f = std::prev(faces.end());
    f->half_edge = half_edges.end();
    f->id = next_id++;
    return f;
}

void HalfEdgeMesh::HalfEdgeMesh::free_half_edge(HalfEdgeRef h) {
    half_edges.erase(h);
}

void HalfEdgeMesh::HalfEdgeMesh::free_vertex(VertexRef v) { vertices.erase(v); }

void HalfEdgeMesh::HalfEdgeMesh::free_edge(EdgeRef e) { edges.erase(e); }

void HalfEdgeMesh::HalfEdgeMesh::free_face(FaceRef f) { faces.erase(f); }
