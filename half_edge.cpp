#include "half_edge.h"

HalfEdgeMesh::HalfEdgeRef HalfEdgeMesh::emplace_half_edge() {
    HalfEdgeRef h;
    if (recycled_half_edges_.empty()) {
        h = half_edges_.emplace(half_edges_.end());
    } else {
        h = recycled_half_edges_.top();
        recycled_half_edges_.pop();
    }
    h->twin = half_edges_.end();
    h->next = half_edges_.end();
    h->edge = edges_.end();
    h->vertex = vertices_.end();
    h->face = faces_.end();
    h->id = next_id++;
    return h;
}

HalfEdgeMesh::VertexRef HalfEdgeMesh::emplace_vertex() {
    VertexRef v;
    if (recycled_vertices_.empty()) {
        v = vertices_.emplace(vertices_.end());
    } else {
        v = recycled_vertices_.top();
        recycled_vertices_.pop();
    }
    v->half_edge = half_edges_.end();
    v->id = next_id++;
    return v;
}

HalfEdgeMesh::EdgeRef HalfEdgeMesh::emplace_edge() {
    EdgeRef e;
    if (recycled_edges_.empty()) {
        e = edges_.emplace(edges_.end());
    } else {
        e = recycled_edges_.top();
        recycled_edges_.pop();
    }
    e->half_edge = half_edges_.end();
    e->id = next_id++;
    return e;
}

HalfEdgeMesh::FaceRef HalfEdgeMesh::emplace_face() {
    FaceRef f;
    if (recycled_faces_.empty()) {
        f = faces_.emplace(faces_.end());
    } else {
        f = recycled_faces_.top();
        recycled_faces_.pop();
    }
    f->half_edge = half_edges_.end();
    f->id = next_id++;
    return f;
}

void HalfEdgeMesh::HalfEdgeMesh::erase_half_edge(HalfEdgeRef&& h) {
    recycled_half_edges_.emplace(h);
}

void HalfEdgeMesh::HalfEdgeMesh::erase_vertex(VertexRef&& v) {
    recycled_vertices_.emplace(v);
}

void HalfEdgeMesh::HalfEdgeMesh::erase_edge(EdgeRef&& e) {
    recycled_edges_.emplace(e);
}

void HalfEdgeMesh::HalfEdgeMesh::erase_face(FaceRef&& f) {
    recycled_faces_.emplace(f);
}
