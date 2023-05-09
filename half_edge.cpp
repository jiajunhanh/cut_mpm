#include "half_edge.h"
HalfEdgeMesh::HalfEdgeRef HalfEdgeMesh::emplace_half_edge() {
    HalfEdgeRef h;
    if (recycled_half_edges.empty()) {
        h = half_edges.emplace(half_edges.end());
    } else {
        h = recycled_half_edges.top();
        recycled_half_edges.pop();
    }
    h->twin = half_edges.end();
    h->next = half_edges.end();
    h->edge = edges.end();
    h->vertex = vertices.end();
    h->face = faces.end();
    h->id = next_id++;
    return h;
}

HalfEdgeMesh::VertexRef HalfEdgeMesh::emplace_vertex() {
    VertexRef v;
    if (recycled_vertices.empty()) {
        v = vertices.emplace(vertices.end());
    } else {
        v = recycled_vertices.top();
        recycled_vertices.pop();
    }
    v->half_edge = half_edges.end();
    v->id = next_id++;
    return v;
}

HalfEdgeMesh::EdgeRef HalfEdgeMesh::emplace_edge() {
    EdgeRef e;
    if (recycled_edges.empty()) {
        e = edges.emplace(edges.end());
    } else {
        e = recycled_edges.top();
        recycled_edges.pop();
    }
    e->half_edge = half_edges.end();
    e->id = next_id++;
    return e;
}

HalfEdgeMesh::FaceRef HalfEdgeMesh::emplace_face() {
    FaceRef f;
    if (recycled_faces.empty()) {
        f = faces.emplace(faces.end());
    } else {
        f = recycled_faces.top();
        recycled_faces.pop();
    }
    f->half_edge = half_edges.end();
    f->id = next_id++;
    return f;
}

void HalfEdgeMesh::HalfEdgeMesh::erase_half_edge(HalfEdgeRef h) {
    recycled_half_edges.push(h);
}

void HalfEdgeMesh::HalfEdgeMesh::erase_vertex(VertexRef v) {
    recycled_vertices.push(v);
}

void HalfEdgeMesh::HalfEdgeMesh::erase_edge(EdgeRef e) {
    recycled_edges.push(e);
}

void HalfEdgeMesh::HalfEdgeMesh::erase_face(FaceRef f) {
    recycled_faces.push(f);
}
