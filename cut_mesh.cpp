#include "cut_mesh.h"

#include <algorithm>
#include <cmath>
#include <queue>
#include <unordered_set>

namespace {

constexpr int kGridNodesNumber = kRowSize * kRowSize;

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
    v.b[d] = false;
    if (v.q[d] <= kMargin) {
        v.q[d] = 0;
        v.b[d] = true;
    }
}

void add_grid_edges(std::vector<Vertex>& cut_vertices,
                    std::vector<Edge>& cut_edges,
                    const std::unordered_set<uint64_t>& edge_hash) {
    auto add_grid_edge = [&](int i, int j) {
        if (edge_hash.count((static_cast<uint64_t>(i) << 32) |
                            static_cast<uint64_t>(j)))
            return;
        if (edge_hash.count((static_cast<uint64_t>(j) << 32) |
                            static_cast<uint64_t>(i)))
            return;
        cut_edges.emplace_back(i, j, false);
    };
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
                    add_grid_edge(node_i, node_j);
                    continue;
                }
                std::sort(begin(v), end(v), [&](int i, int j) {
                    return cut_vertices[i].q[d] < cut_vertices[j].q[d];
                });
                add_grid_edge(node_i, v.front());
                for (int i = 0, size = static_cast<int>(v.size()); i < size - 1;
                     ++i) {
                    add_grid_edge(v[i], v[i + 1]);
                }
                add_grid_edge(v.back(), node_j);
            }
        }
    }
    for (int i = 0; i < kGridSize; ++i) {
        add_grid_edge(i * kRowSize + kGridSize, (i + 1) * kRowSize + kGridSize);
        add_grid_edge(kGridSize * kRowSize + i, kGridSize * kRowSize + i + 1);
    }
}

std::pair<std::vector<Vertex>, std::vector<Edge>>
compute_cut_vertices_and_edges(
    const std::vector<std::array<Real, 2>>& vertices) {
    std::vector<Vertex> cut_vertices;
    for (int r = 0; r < kRowSize; ++r) {
        for (int c = 0; c < kRowSize; ++c) {
            Vertex v{};
            v.c[0] = c;
            v.c[1] = r;
            v.b[0] = v.b[1] = true;
            cut_vertices.emplace_back(v);
        }
    }
    std::vector<int> vertex_ids;
    for (const auto& v : vertices) {
        Vertex cut_vertex;
        for (int d = 0; d < 2; ++d) {
            Real x = v[d] * kGridSize;
            cut_vertex.c[d] = static_cast<int>(std::floor(x));
            cut_vertex.q[d] = x - std::floor(x);
            cut_vertex.b[d] = false;
            if (cut_vertex.q[d] <= kMargin) {
                cut_vertex.q[d] = 0;
                cut_vertex.b[d] = true;
            }
        }
        if (cut_vertex.b[0] && cut_vertex.b[1]) {
            vertex_ids.emplace_back(cut_vertex.c[1] * kRowSize +
                                    cut_vertex.c[0]);
            continue;
        }
        cut_vertices.emplace_back(cut_vertex);
        vertex_ids.emplace_back(static_cast<int>(cut_vertices.size()) - 1);
    }
    std::vector<Edge> cut_edges;
    std::unordered_set<uint64_t> edge_hash;
    for (int i = 0, vertex_ids_size = static_cast<int>(vertex_ids.size());
         i < vertex_ids_size; ++i) {
        auto j = (i + 1) % vertex_ids_size;
        auto id_i = vertex_ids[i];
        auto id_j = vertex_ids[j];
        std::vector<std::pair<Real, int>> intersection_points;
        const Vertex vi = cut_vertices[id_i];
        const Vertex vj = cut_vertices[id_j];
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
            if (distance == 0) continue;
            for (int z = z_min; z <= z_max; ++z) {
                Real t = (static_cast<Real>(z - vi.c[d]) - vi.q[d]) / distance;
                if (t <= 0 || t >= 1) continue;
                Vertex cut_v;
                lerp(cut_v, vi, vj, t, d ^ 1);
                cut_v.c[d] = z;
                cut_v.q[d] = 0;
                cut_v.b[d] = true;
                if (cut_v.b[d ^ 1]) {
                    auto grid_id = cut_v.c[1] * kRowSize + cut_v.c[0];
                    intersection_points.emplace_back(t, grid_id);
                    continue;
                }
                cut_vertices.emplace_back(cut_v);
                intersection_points.emplace_back(
                    t, static_cast<int>(cut_vertices.size()) - 1);
            }
        }
        intersection_points.emplace_back(Real{0}, id_i);
        intersection_points.emplace_back(Real{1}, id_j);
        std::sort(
            begin(intersection_points), end(intersection_points),
            [](const auto& a, const auto& b) { return a.first < b.first; });
        for (int k = 0, size = static_cast<int>(intersection_points.size());
             k + 1 < size; ++k) {
            auto id0 = intersection_points[k].second;
            auto id1 = intersection_points[k + 1].second;
            if (id0 == id1) {
                continue;
            }
            cut_edges.emplace_back(id0, id1, true);
            edge_hash.emplace((static_cast<uint64_t>(id0) << 32) |
                              static_cast<uint64_t>(id1));
        }
    }
    add_grid_edges(cut_vertices, cut_edges, edge_hash);
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

CutMesh::CutMesh(std::vector<std::array<Real, 2>> vertices) {
    Real orientation = 0;
    for (int i = 0, vertices_size = static_cast<int>(vertices.size());
         i < vertices_size; ++i) {
        auto j = (i + 1) % vertices_size;
        orientation += (vertices[j][0] - vertices[i][0]) *
                       (vertices[j][1] + vertices[i][1]);
    }
    if (orientation > 0) reverse(begin(vertices), end(vertices));

    auto [cut_vertices, cut_edges] = compute_cut_vertices_and_edges(vertices);

    std::vector<std::vector<CutMesh::HalfEdgeRef>> vertex_half_edges;
    vertices_.reserve(cut_vertices.size());
    vertex_half_edges.resize(cut_vertices.size());

    for (const auto& cut_vertex : cut_vertices) {
        auto v = emplace_vertex();
        for (int d = 0; d < 2; ++d) {
            v->position[d] =
                static_cast<Real>(cut_vertex.c[d]) + cut_vertex.q[d];
            // v->on_edge[d] = cut_vertex.b[d];
        }
        // v->grid_id = cut_vertex.c[1] * kGridSize + cut_vertex.c[0];
        v->position /= kGridSize;
    }

    half_edges_.reserve(2 * cut_edges.size());
    edges_.reserve(cut_edges.size());
    for (const auto& cut_edge : cut_edges) {
        auto h0 = emplace_half_edge();
        auto h1 = emplace_half_edge();
        auto v0 = begin(vertices_) + cut_edge.i;
        auto v1 = begin(vertices_) + cut_edge.j;
        auto e = emplace_edge();
        h0->set_tnvef(h1, h0->next, v0, e, h0->face);
        h1->set_tnvef(h0, h1->next, v1, e, h1->face);
        h0->is_boundary = h1->is_boundary = cut_edge.is_boundary;
        v0->half_edge = h0;
        v1->half_edge = h1;
        e->half_edge = h0;

        if (cut_edge.is_boundary) {
            v0->on_boundary = v1->on_boundary = true;
            Vec2 normal = v1->position - v0->position;
            normal = Vec2(normal.y(), -normal.x());
            normal.normalize();
            h0->normal = h1->normal = normal;
            v0->normal += normal;
            v1->normal += normal;
        }

        vertex_half_edges[cut_edge.i].emplace_back(h0);
        vertex_half_edges[cut_edge.j].emplace_back(h1);
    }

    for (auto& v : vertices_) {
        if (!v.on_boundary) continue;
        v.normal.normalize();
    }

    for (int iv = 0, n_vertices = static_cast<int>(vertices_.size());
         iv < n_vertices; ++iv) {
        auto v = begin(vertices_) + iv;
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
    for (auto half_edge = begin(half_edges_), half_edges_end = end(half_edges_);
         half_edge != half_edges_end; ++half_edge) {
        if (visited_half_edges.count(half_edge->id)) continue;
        auto f = emplace_face();
        f->half_edge = half_edge;
        auto h = half_edge;
        do {
            visited_half_edges.emplace(h->id);
        } while ((h = h->next) != half_edge);
    }

    grids_.resize(kGridSize * kGridSize);
    for (auto face = begin(faces_), faces_end = end(faces_); face != faces_end;
         ++face) {
        face->calculate_center();
        Vec2 center = face->center * kGridSize;
        int r = static_cast<int>(std::floor(center.y()));
        int c = static_cast<int>(std::floor(center.x()));
        grid(r, c).faces.emplace_back(face);
        auto half_edge = face->half_edge;
        auto h = half_edge;
        do {
            h->face = face;
        } while ((h = h->next) != half_edge);
    }

    for (int r = 0; r < kGridSize; ++r) {
        for (int c = 0; c < kGridSize; ++c) {
            grid(r, c).vertex = begin(vertices_) + (r * kRowSize + c);
        }
    }
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

void CutMesh::calculate_neighbor_nodes_and_boundaries_of_faces() {
    std::vector<bool> face_visited(faces_.size());
    std::vector<bool> is_neighbor_node(vertices_.size());
    std::vector<bool> neighbor_node_sides(vertices_.size());
    std::vector<bool> is_neighbor_boundary(half_edges_.size());
    std::vector<CutMesh::FaceRef> neighbor_faces;
    std::vector<bool> neighbor_face_sides;
    std::queue<std::pair<CutMesh::FaceRef, bool>> q;
    for (auto cut_face = begin(faces_), faces_end = end(faces_);
         cut_face != faces_end; ++cut_face) {
        std::fill(begin(face_visited), end(face_visited), false);
        std::fill(begin(is_neighbor_node), end(is_neighbor_node), false);
        std::fill(begin(neighbor_node_sides), end(neighbor_node_sides), false);
        std::fill(begin(is_neighbor_boundary), end(is_neighbor_boundary),
                  false);
        cut_face->neighbor_nodes.clear();
        cut_face->neighbor_node_sides.clear();
        cut_face->neighbor_boundaries.clear();
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
                auto id = h->vertex->id;
                if (s) neighbor_node_sides[id] = true;
                if (is_neighbor_node[id]) continue;
                is_neighbor_node[id] = true;
                cut_face->neighbor_nodes.emplace_back(id);
            } while ((h = h->next) != f->half_edge);

            h = f->half_edge;
            do {
                if (s && h->is_boundary && !is_neighbor_boundary[h->id]) {
                    is_neighbor_boundary[h->id] = true;
                    cut_face->neighbor_boundaries.emplace_back(h->id);
                }
            } while ((h = h->next) != f->half_edge);
        }
        for (const auto& i : cut_face->neighbor_nodes)
            cut_face->neighbor_node_sides.emplace_back(neighbor_node_sides[i]);
    }
}

void CutMesh::calculate_node_normals() {
    for (auto& vertex : vertices_) {
        if (vertex.on_boundary) continue;
        vertex.normal = Vec2::Zero();
        auto face = vertex.half_edge->face;
        const auto& neighbor_nodes = face->neighbor_nodes;
        const auto& neighbor_node_sides = face->neighbor_node_sides;
        Real w_sum = 0;
        for (int i = 0,
                 n_neighbor_nodes = static_cast<int>(neighbor_nodes.size());
             i < n_neighbor_nodes; ++i) {
            if (!neighbor_node_sides[i]) continue;
            auto v = begin(vertices_) + neighbor_nodes[i];
            if (!v->on_boundary) continue;
            Vec2 d = vertex.position - v->position;
            if (d.norm() > kDeltaX * 0.5) continue;
            auto w = interpolate(d);
            vertex.normal += w * v->normal;
            w_sum += w;
        }
        if (w_sum > 0) vertex.normal /= w_sum;
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