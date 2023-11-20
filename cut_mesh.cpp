#include "cut_mesh.h"

#include <algorithm>
#include <cmath>
#include <queue>
#include <unordered_set>

namespace {

constexpr int kGridNodesNumber = kRowSize * kRowSize;
/*
static Real random_real() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    // static std::mt19937 gen(42);
    static std::uniform_real_distribution<Real> dis(0.0, 1.0);
    return dis(gen);
}
*/

struct Vertex {
    int c[2] = {};
    Real q[2] = {};
    bool b[2] = {};

    [[nodiscard]] Real coord(int d) const {
        return static_cast<Real>(c[d]) + q[d];
    }
};

struct Edge {
   public:
    int i = 0;
    int j = 0;
    bool is_boundary = false;
    Edge(int i_, int j_, bool is_boundary_)
        : i(i_), j(j_), is_boundary(is_boundary_) {}
};

void lerp(Vertex& v, const Vertex& vi, const Vertex& vj, Real t, int d) {
    Real coord_i = vi.coord(d);
    Real coord_j = vj.coord(d);
    Real coord = coord_i + (coord_j - coord_i) * t;
    v.c[d] = static_cast<int>(std::floor(coord));
    v.q[d] = coord - std::floor(coord);
    v.b[d] = v.q[d] == 0.0f;
}

void add_grid_edges(std::vector<Vertex>& cut_vertices,
                    std::vector<Edge>& cut_edges) {
    std::vector<int> grid_cut_vertices[kGridSize * kGridSize][2];
    for (int i = kGridNodesNumber, size = static_cast<int>(cut_vertices.size());
         i < size; ++i) {
        const Vertex& v = cut_vertices[i];
        for (int d = 0; d < 2; ++d) {
            if (v.b[d]) {
                grid_cut_vertices[v.c[1] * kGridSize + v.c[0]][d ^ 1]
                    .emplace_back(i);
            }
        }
    }
    for (int r = 0; r < kGridSize; ++r) {
        for (int c = 0; c < kGridSize; ++c) {
            int node_i = r * kRowSize + c;
            for (int d = 0; d < 2; ++d) {
                auto& v = grid_cut_vertices[r * kGridSize + c][d];
                int node_j = node_i + d * kGridSize + 1;
                if (v.empty()) {
                    cut_edges.emplace_back(node_i, node_j, false);
                    continue;
                }
                std::sort(begin(v), end(v), [&](int i, int j) {
                    return cut_vertices[i].q[d] < cut_vertices[j].q[d];
                });
                cut_edges.emplace_back(node_i, v.front(), false);
                for (int i = 0, size = static_cast<int>(v.size()); i < size - 1;
                     ++i) {
                    cut_edges.emplace_back(v[i], v[i + 1], false);
                }
                cut_edges.emplace_back(v.back(), node_j, false);
            }
        }
    }
    for (int i = 0; i < kGridSize; ++i) {
        cut_edges.emplace_back(i * kRowSize + kGridSize,
                               (i + 1) * kRowSize + kGridSize, false);
        cut_edges.emplace_back(kGridSize * kRowSize + i,
                               kGridSize * kRowSize + i + 1, false);
    }
}

std::pair<std::vector<Vertex>, std::vector<Edge>>
compute_cut_vertices_and_edges(const std::vector<std::array<Real, 2>>& vertices,
                               const std::vector<std::array<int, 2>>& edges) {
    std::vector<Vertex> cut_vertices;
    std::vector<Edge> cut_edges;
    for (int r = 0; r < kRowSize; ++r) {
        for (int c = 0; c < kRowSize; ++c) {
            Vertex v{};
            v.c[0] = c;
            v.c[1] = r;
            v.b[0] = v.b[1] = true;
            cut_vertices.emplace_back(v);
        }
    }
    for (const auto& v : vertices) {
        Vertex cut_vertex;
        for (int d = 0; d < 2; ++d) {
            Real x = v[d] * kGridSize;
            cut_vertex.c[d] = static_cast<int>(std::floor(x));
            cut_vertex.q[d] = x - std::floor(x);
            cut_vertex.b[d] = (cut_vertex.q[d] == Real{0});
        }
        cut_vertices.emplace_back(cut_vertex);
    }
    for (const auto& e : edges) {
        std::vector<std::pair<Real, int>> intersection_points;
        const Vertex vi = cut_vertices[e[0] + kGridNodesNumber];
        const Vertex vj = cut_vertices[e[1] + kGridNodesNumber];
        for (int d = 0; d < 2; ++d) {
            int z_min;
            int z_max;
            if (vi.c[d] < vj.c[d] ||
                (vi.c[d] == vj.c[d] && vi.q[d] < vj.q[d])) {
                z_min = vi.b[d] ? vi.c[d] : vi.c[d] + 1;
                z_max = vj.c[d];
            } else {
                z_min = vj.b[d] ? vj.c[d] : vj.c[d] + 1;
                z_max = vi.c[d];
            }
            Real distance =
                static_cast<Real>(vj.c[d] - vi.c[d]) + (vj.q[d] - vi.q[d]);
            for (int z = z_min; z <= z_max; ++z) {
                Real t = (static_cast<Real>(z - vi.c[d]) - vi.q[d]) / distance;
                Vertex cut_v;
                lerp(cut_v, vi, vj, t, d ^ 1);
                cut_v.c[d] = z;
                cut_v.q[d] = 0;
                cut_v.b[d] = true;
                if (cut_v.b[d ^ 1]) {
                    intersection_points.emplace_back(
                        t, cut_v.c[1] * kRowSize * cut_v.c[0]);
                    continue;
                }
                cut_vertices.emplace_back(cut_v);
                intersection_points.emplace_back(
                    t, static_cast<int>(cut_vertices.size()) - 1);
            }
        }
        std::sort(
            begin(intersection_points), end(intersection_points),
            [](const auto& a, const auto& b) { return a.first < b.first; });
        if (intersection_points.empty()) {
            cut_edges.emplace_back(e[0] + kGridNodesNumber,
                                   e[1] + kGridNodesNumber, true);
            continue;
        }
        cut_edges.emplace_back(e[0] + kGridNodesNumber,
                               intersection_points.front().second, true);
        for (int i = 0, size = static_cast<int>(intersection_points.size());
             i + 1 < size; ++i) {
            if (intersection_points[i].second ==
                intersection_points[i + 1].second) {
                continue;
            }
            cut_edges.emplace_back(intersection_points[i].second,
                                   intersection_points[i + 1].second, true);
        }
        cut_edges.emplace_back(intersection_points.back().second,
                               e[1] + kGridNodesNumber, true);
    }
    add_grid_edges(cut_vertices, cut_edges);
    return {cut_vertices, cut_edges};
}

bool is_neighbor_face(CutMesh::FaceRef face0, CutMesh::FaceRef face1) {
    auto min_x = std::numeric_limits<Real>::max();
    auto max_x = std::numeric_limits<Real>::min();
    auto min_y = std::numeric_limits<Real>::max();
    auto max_y = std::numeric_limits<Real>::min();

    auto h = face0->half_edge;
    do {
        auto p = h->vertex->position;
        min_x = std::min(min_x, p.x());
        max_x = std::max(max_x, p.x());
        min_y = std::min(min_y, p.y());
        max_y = std::max(max_y, p.y());
    } while ((h = h->next) != face0->half_edge);

    min_x -= kKernelRange * kDeltaX;
    max_x += kKernelRange * kDeltaX;
    min_y -= kKernelRange * kDeltaX;
    max_y += kKernelRange * kDeltaX;

    h = face1->half_edge;
    do {
        auto p = h->vertex->position;
        if (p.x() > min_x && p.x() < max_x && p.y() > min_y && p.y() < max_y)
            return true;
    } while ((h = h->next) != face1->half_edge);
    return false;
}

}  // namespace

CutMesh construct_cut_mesh(const std::vector<std::array<Real, 2>>& vertices,
                           const std::vector<std::array<int, 2>>& edges) {
    auto [cut_vertices, cut_edges] =
        compute_cut_vertices_and_edges(vertices, edges);

    CutMesh cut_mesh;
    std::vector<CutMesh::VertexRef> mesh_vertices;
    std::vector<std::vector<CutMesh::HalfEdgeRef>> vertex_half_edges;
    cut_mesh.vertices().reserve(cut_vertices.size());
    mesh_vertices.reserve(cut_vertices.size());
    vertex_half_edges.resize(cut_vertices.size());

    for (const Vertex& cut_vertex : cut_vertices) {
        auto v = cut_mesh.emplace_vertex();
        for (int d = 0; d < 2; ++d) {
            v->position[d] =
                static_cast<Real>(cut_vertex.c[d]) + cut_vertex.q[d];
            v->on_edge[d] = cut_vertex.b[d];
        }
        v->grid_id = cut_vertex.c[1] * kGridSize + cut_vertex.c[0];
        v->position /= kGridSize;
        mesh_vertices.emplace_back(v);
    }

    cut_mesh.half_edges().reserve(2 * cut_edges.size());
    cut_mesh.edges().reserve(cut_edges.size());
    for (const Edge& cut_edge : cut_edges) {
        auto h0 = cut_mesh.emplace_half_edge();
        auto h1 = cut_mesh.emplace_half_edge();
        auto v0 = mesh_vertices[cut_edge.i];
        auto v1 = mesh_vertices[cut_edge.j];
        auto e = cut_mesh.emplace_edge();
        h0->set_tnvef(h1, h0->next, v0, e, h0->face);
        h1->set_tnvef(h0, h1->next, v1, e, h1->face);
        h0->is_boundary = h1->is_boundary = cut_edge.is_boundary;
        v0->half_edge = h0;
        v1->half_edge = h1;
        e->half_edge = h0;
        vertex_half_edges[cut_edge.i].emplace_back(h0);
        vertex_half_edges[cut_edge.j].emplace_back(h1);
    }

    for (int iv = 0, n_vertices = static_cast<int>(mesh_vertices.size());
         iv < n_vertices; ++iv) {
        auto v = mesh_vertices[iv];
        auto& half_edges = vertex_half_edges[iv];
        int n_half_edges = static_cast<int>(half_edges.size());
        std::vector<Real> angles(n_half_edges);
        for (int i = 0; i < n_half_edges; ++i) {
            auto h0 = half_edges[i];
            auto v1 = h0->twin->vertex;
            angles[i] = std::atan2(v1->position[1] - v->position[1],
                                   v1->position[0] - v->position[0]);
        }
        std::vector<int> sorted_edges(n_half_edges);
        for (int i = 0; i < n_half_edges; ++i) {
            sorted_edges[i] = i;
        }
        std::sort(
            begin(sorted_edges), end(sorted_edges),
            [&angles = angles](int i, int j) { return angles[i] < angles[j]; });
        for (int i = 0; i < n_half_edges; ++i) {
            int j = (i + 1) % n_half_edges;
            auto h0 = half_edges[sorted_edges[i]];
            auto h1 = half_edges[sorted_edges[j]];
            h1->twin->next = h0;
        }
    }

    std::unordered_set<int> visited_half_edges;
    for (auto half_edge = begin(cut_mesh.half_edges()),
              half_edges_end = end(cut_mesh.half_edges());
         half_edge != half_edges_end; ++half_edge) {
        if (visited_half_edges.count(half_edge->id)) continue;
        auto f = cut_mesh.emplace_face();
        f->half_edge = half_edge;
        auto h = half_edge;
        do {
            visited_half_edges.emplace(h->id);
        } while ((h = h->next) != half_edge);
    }

    cut_mesh.grids().resize(kGridSize * kGridSize);
    for (auto face = begin(cut_mesh.faces()), faces_end = end(cut_mesh.faces());
         face != faces_end; ++face) {
        face->calculate_center();
        Vec2 center = face->center * kGridSize;
        int r = static_cast<int>(std::floor(center.y()));
        int c = static_cast<int>(std::floor(center.x()));
        cut_mesh.grid(r, c).faces.emplace_back(face);
        auto half_edge = face->half_edge;
        auto h = half_edge;
        do {
            h->face = face;
        } while ((h = h->next) != half_edge);
    }

    for (int r = 0; r < kGridSize; ++r) {
        for (int c = 0; c < kGridSize; ++c) {
            cut_mesh.grid(r, c).vertex = mesh_vertices[r * kRowSize + c];
        }
    }

    return cut_mesh;
}

CutMesh::FaceRef CutMesh::get_enclosing_face(Vec2 x) const {
    while (true) {
        auto r = static_cast<int>(std::floor(x.y() * kGridSize));
        auto c = static_cast<int>(std::floor(x.x() * kGridSize));
        for (const auto& face : grid(r, c).faces) {
            if (face->enclose(x)) return face;
        }
        x += Vec2::Random() * 10 * std::numeric_limits<Real>::epsilon();
    }
    // assert(false && "Point must be in a face!");
    // return grid(r, c).faces[0];
}

void CutMesh::calculate_neighbor_nodes_of_faces() {
    std::vector<bool> face_visited(faces_.size());
    std::vector<bool> is_neighbor_node(vertices_.size());
    std::vector<CutMesh::FaceRef> neighbor_faces;
    std::vector<bool> neighbor_face_sides;
    std::queue<std::pair<CutMesh::FaceRef, bool>> q;
    for (auto cut_face = begin(faces_), faces_end = end(faces_);
         cut_face != faces_end; ++cut_face) {
        std::fill(begin(face_visited), end(face_visited), false);
        std::fill(begin(is_neighbor_node), end(is_neighbor_node), false);
        cut_face->neighbor_nodes.clear();
        cut_face->neighbor_node_sides.clear();
        neighbor_faces.clear();
        neighbor_face_sides.clear();
        q.emplace(cut_face, true);
        while (!q.empty()) {
            auto [face, side] = q.front();
            q.pop();
            if (face_visited[face->id]) continue;
            face_visited[face->id] = true;
            neighbor_faces.emplace_back(face);
            neighbor_face_sides.emplace_back(side);
            auto h = face->half_edge;
            do {
                auto f = h->twin->face;
                if (face_visited[f->id]) continue;
                if (!is_neighbor_face(cut_face, f)) continue;
                q.emplace(f, side ^ h->is_boundary);
            } while ((h = h->next) != face->half_edge);
        }
        for (int i = 0, n = static_cast<int>(neighbor_faces.size()); i < n;
             ++i) {
            auto f = neighbor_faces[i];
            bool s = neighbor_face_sides[i];
            auto h = f->half_edge;
            do {
                if (is_neighbor_node[h->vertex->id]) continue;
                is_neighbor_node[h->vertex->id] = true;
                cut_face->neighbor_nodes.emplace_back(h->vertex->id);
                cut_face->neighbor_node_sides.emplace_back(
                    (int{s | (i == 0)} << kValidBit));
            } while ((h = h->next) != f->half_edge);
        }
    }
}

void CutMesh::Face::calculate_center() {
    center = Vec2::Zero();
    int cnt = 0;
    auto h = half_edge;
    do {
        center += h->vertex->position;
        ++cnt;
    } while ((h = h->next) != half_edge);
    center /= cnt;
}