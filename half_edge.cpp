#include <cmath>

#include "cut_mesh.h"

CutMesh::HalfEdgeRef CutMesh::emplace_half_edge() {
    auto h = half_edges_.emplace(end(half_edges_));
    //    if (recycled_half_edges_.empty()) {
    //        h = half_edges_.emplace(half_edges_.end());
    //    } else {
    //        h = recycled_half_edges_.top();
    //        recycled_half_edges_.pop();
    //    }
    h->twin = half_edges_.end();
    h->next = half_edges_.end();
    h->edge = edges_.end();
    h->vertex = vertices_.end();
    h->face = faces_.end();
    h->id = static_cast<int>(half_edges_.size()) - 1;
    return h;
}

CutMesh::VertexRef CutMesh::emplace_vertex() {
    auto v = vertices_.emplace(end(vertices_));
    //    if (recycled_vertices_.empty()) {
    //        v = vertices_.emplace(vertices_.end());
    //    } else {
    //        v = recycled_vertices_.top();
    //        recycled_vertices_.pop();
    //    }
    v->half_edge = half_edges_.end();
    v->id = static_cast<int>(vertices_.size()) - 1;
    return v;
}

CutMesh::EdgeRef CutMesh::emplace_edge() {
    auto e = edges_.emplace(end(edges_));
    //    if (recycled_edges_.empty()) {
    //        e = edges_.emplace(edges_.end());
    //    } else {
    //        e = recycled_edges_.top();
    //        recycled_edges_.pop();
    //    }
    e->half_edge = half_edges_.end();
    e->id = static_cast<int>(edges_.size()) - 1;
    return e;
}

CutMesh::FaceRef CutMesh::emplace_face() {
    auto f = faces_.emplace(end(faces_));
    //    if (recycled_faces_.empty()) {
    //        f = faces_.emplace(faces_.end());
    //    } else {
    //        f = recycled_faces_.top();
    //        recycled_faces_.pop();
    //    }
    f->half_edge = half_edges_.end();
    f->id = static_cast<int>(faces_.size());
    return f;
}

bool CutMesh::Face::enclose(const Vec2& x) const {
    Real winding_number = 0.0;
    auto h = half_edge;
    do {
        auto v0 = h->vertex;
        auto v1 = h->twin->vertex;
        Vec2 a = v0->position - x;
        Vec2 b = v1->position - x;
        winding_number += std::atan2(a.x() * b.y() - a.y() * b.x(),
                                     a.x() * b.x() + a.y() * b.y());
        h = h->next;
    } while (h != half_edge);
    // std::cerr << winding_number << '\n';
    return winding_number > std::numeric_limits<Real>::epsilon() * 1e4;
}

// void CutMesh::CutMesh::erase_half_edge(HalfEdgeRef&& h) {
//     recycled_half_edges_.emplace(h);
// }
//
// void CutMesh::CutMesh::erase_vertex(VertexRef&& v) {
//     recycled_vertices_.emplace(v);
// }
//
// void CutMesh::CutMesh::erase_edge(EdgeRef&& e) {
//     recycled_edges_.emplace(e);
// }
//
// void CutMesh::CutMesh::erase_face(FaceRef&& f) {
//     recycled_faces_.emplace(f);
// }
