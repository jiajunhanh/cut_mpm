#include "cut_mesh.h"

#include <algorithm>
#include <cmath>
#include <queue>
#include <unordered_set>

CutMesh::CutMesh(std::vector<std::array<Real, 2>> vertices, int quality)
    : grid_size_(static_cast<int>(std::pow(2, quality - 1)) * 8),
      row_size_(grid_size_ + 1),
      n_grid_nodes_(row_size_ * row_size_),
      delta_x_(Real{1.0} / static_cast<Real>(grid_size_)),
      margin_(delta_x_ / 32)

{
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
        v->position /= static_cast<Real>(grid_size_);
        v->convex = cut_vertex.convex;
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
        if (v.on_boundary) v.normal.normalize();
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

    grids_.resize(grid_size_ * grid_size_);
    for (auto face = begin(faces_), faces_end = end(faces_); face != faces_end;
         ++face) {
        face->calculate_center();
        Vec2 center = face->center * grid_size_;
        int r = static_cast<int>(std::floor(center.y()));
        int c = static_cast<int>(std::floor(center.x()));
        grid(r, c).faces.emplace_back(face);
        auto half_edge = face->half_edge;
        auto h = half_edge;
        do {
            h->face = face;
        } while ((h = h->next) != half_edge);
    }

    for (int r = 0; r < grid_size_; ++r) {
        for (int c = 0; c < grid_size_; ++c) {
            grid(r, c).vertex = begin(vertices_) + (r * row_size_ + c);
            if (grid(r, c).faces.size() > 1) grid(r, c).near_boundary = true;
            if (r + 1 < grid_size_ && grid(r + 1, c).faces.size() > 1)
                grid(r, c).near_boundary = true;
            if (c + 1 < grid_size_ && grid(r, c + 1).faces.size() > 1)
                grid(r, c).near_boundary = true;
            if (c + 1 < grid_size_ && r + 1 < grid_size_ &&
                grid(r + 1, c + 1).faces.size() > 1)
                grid(r, c).near_boundary = true;
        }
    }
}

CutMesh::FaceRef CutMesh::get_enclosing_face(Vec2 x) const {
    while (true) {
        auto r =
            static_cast<int>(std::floor(x.y() * static_cast<Real>(grid_size_)));
        auto c =
            static_cast<int>(std::floor(x.x() * static_cast<Real>(grid_size_)));
        if (grid(r, c).faces.size() == 1) return grid(r, c).faces[0];
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
    std::vector<bool> is_neighbor_boundary(half_edges_.size());
    std::vector<CutMesh::FaceRef> neighbor_faces;
    std::queue<CutMesh::FaceRef> q;
    for (auto cut_face = begin(faces_), faces_end = end(faces_);
         cut_face != faces_end; ++cut_face) {
        std::fill(begin(face_visited), end(face_visited), false);
        std::fill(begin(is_neighbor_node), end(is_neighbor_node), false);
        std::fill(begin(is_neighbor_boundary), end(is_neighbor_boundary),
                  false);
        cut_face->neighbor_nodes.clear();
        cut_face->neighbor_boundaries.clear();
        neighbor_faces.clear();
        q.emplace(cut_face);
        while (!q.empty()) {
            auto face = q.front();
            q.pop();
            if (face_visited[face->id]) continue;
            face_visited[face->id] = true;
            neighbor_faces.emplace_back(face);
            auto h = face->half_edge;
            do {
                if (h->is_boundary) continue;
                auto f = h->twin->face;
                if (face_visited[f->id]) continue;
                if (!is_neighbor_face(cut_face, f)) continue;
                q.emplace(f);
            } while ((h = h->next) != face->half_edge);
        }
        for (auto f : neighbor_faces) {
            auto h = f->half_edge;
            do {
                if (h->vertex->convex) cut_face->near_convex = true;
                auto id = h->vertex->id;
                if (is_neighbor_node[id]) continue;
                is_neighbor_node[id] = true;
                cut_face->neighbor_nodes.emplace_back(id);
            } while ((h = h->next) != f->half_edge);

            h = f->half_edge;
            do {
                if (h->is_boundary && !is_neighbor_boundary[h->id]) {
                    is_neighbor_boundary[h->id] = true;
                    cut_face->neighbor_boundaries.emplace_back(h->id);
                }
            } while ((h = h->next) != f->half_edge);
        }
    }
}

void CutMesh::calculate_node_normals() {
    for (auto& vertex : vertices_) {
        if (vertex.on_boundary) continue;
        vertex.normal.setZero();
        auto face = vertex.half_edge->face;
        if (face->near_convex) continue;
        // Real w_sum = 0;
        for (int id : face->neighbor_boundaries) {
            auto h = begin(half_edges_) + id;
            if (!h->is_boundary) continue;
            auto v0 = h->vertex;
            auto v1 = h->twin->vertex;
            auto p0 = v0->position;
            auto p1 = v1->position;
            Real d = (vertex.position - p0).dot(h->normal);
            if (d < 0 || d > 0.5 * delta_x_) continue;
            Vec2 tangent = p1 - p0;
            Real len = tangent.norm();
            tangent.normalize();
            Real t = tangent.dot(vertex.position - p0);
            if (t < 0 || t > len) continue;
            auto w = interpolate(d);
            vertex.normal += w * h->normal;
            // w_sum += w;
        }
        // if (w_sum > 0) vertex.normal /= w_sum;
        if (!vertex.normal.isZero()) vertex.normal.normalize();
    }
}

void CutMesh::lerp(CutVertex& v, const CutVertex& vi, const CutVertex& vj,
                   Real t, int d) const {
    Real coord_i = vi.coord(d);
    Real coord_j = vj.coord(d);
    Real coord = coord_i + (coord_j - coord_i) * t;
    v.c[d] = static_cast<int>(std::floor(coord));
    v.q[d] = coord - std::floor(coord);
    v.b[d] = false;
    if (v.q[d] <= margin_) {
        v.q[d] = 0;
        v.b[d] = true;
    }
}

void CutMesh::add_grid_edges(
    std::vector<CutVertex>& cut_vertices, std::vector<CutEdge>& cut_edges,
    const std::unordered_set<uint64_t>& edge_hash) const {
    auto add_grid_edge = [&](int i, int j) {
        if (edge_hash.count((static_cast<uint64_t>(i) << 32) |
                            static_cast<uint64_t>(j)))
            return;
        if (edge_hash.count((static_cast<uint64_t>(j) << 32) |
                            static_cast<uint64_t>(i)))
            return;
        cut_edges.emplace_back(i, j, false);
    };
    std::vector<std::vector<std::vector<int>>> grid_cut_vertices(
        grid_size_ * grid_size_, std::vector<std::vector<int>>(2));
    for (int i = n_grid_nodes_, size = static_cast<int>(cut_vertices.size());
         i < size; ++i) {
        const CutVertex& v = cut_vertices[i];
        for (int d = 0; d < 2; ++d) {
            if (v.b[d]) {
                grid_cut_vertices[v.c[1] * grid_size_ + v.c[0]][d ^ 1]
                    .emplace_back(i);
            }
        }
    }
    for (int r = 0; r < grid_size_; ++r) {
        for (int c = 0; c < grid_size_; ++c) {
            int node_i = r * row_size_ + c;
            for (int d = 0; d < 2; ++d) {
                auto& v = grid_cut_vertices[r * grid_size_ + c][d];
                int node_j = node_i + d * grid_size_ + 1;
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
    for (int i = 0; i < grid_size_; ++i) {
        add_grid_edge(i * row_size_ + grid_size_,
                      (i + 1) * row_size_ + grid_size_);
        add_grid_edge(grid_size_ * row_size_ + i,
                      grid_size_ * row_size_ + i + 1);
    }
}

std::pair<std::vector<CutMesh::CutVertex>, std::vector<CutMesh::CutEdge>>
CutMesh::compute_cut_vertices_and_edges(
    const std::vector<std::array<Real, 2>>& vertices) {
    std::vector<CutVertex> cut_vertices;
    for (int r = 0; r < row_size_; ++r) {
        for (int c = 0; c < row_size_; ++c) {
            CutVertex v{};
            v.c[0] = c;
            v.c[1] = r;
            v.b[0] = v.b[1] = true;
            cut_vertices.emplace_back(v);
        }
    }
    std::vector<int> vertex_ids;
    for (int i = 0, n_vertices = static_cast<int>(vertices.size());
         i < n_vertices; ++i) {
        int h = (i + n_vertices - 1) % n_vertices;
        int j = (i + 1) % n_vertices;
        auto& v = vertices[i];
        CutVertex cut_vertex;
        for (int d = 0; d < 2; ++d) {
            Real x = v[d] * static_cast<Real>(grid_size_);
            cut_vertex.c[d] = static_cast<int>(std::floor(x));
            cut_vertex.q[d] = x - std::floor(x);
            cut_vertex.b[d] = false;
            if (cut_vertex.q[d] <= margin_) {
                cut_vertex.q[d] = 0;
                cut_vertex.b[d] = true;
            }
        }
        auto x0 = v[0] - vertices[h][0];
        auto y0 = v[1] - vertices[h][1];
        auto x1 = vertices[j][0] - v[0];
        auto y1 = vertices[j][1] - v[1];
        auto cross = x0 * y1 - x1 * y0;
        if (cross > 0) cut_vertex.convex = true;
        if (cut_vertex.b[0] && cut_vertex.b[1]) {
            auto grid_id = cut_vertex.c[1] * row_size_ + cut_vertex.c[0];
            vertex_ids.emplace_back(grid_id);
            cut_vertices[grid_id].convex = cut_vertex.convex;
            continue;
        }
        cut_vertices.emplace_back(cut_vertex);
        vertex_ids.emplace_back(static_cast<int>(cut_vertices.size()) - 1);
    }
    std::vector<CutEdge> cut_edges;
    std::unordered_set<uint64_t> edge_hash;
    for (int i = 0, vertex_ids_size = static_cast<int>(vertex_ids.size());
         i < vertex_ids_size; ++i) {
        auto j = (i + 1) % vertex_ids_size;
        auto id_i = vertex_ids[i];
        auto id_j = vertex_ids[j];
        std::vector<std::pair<Real, int>> intersection_points;
        const CutVertex vi = cut_vertices[id_i];
        const CutVertex vj = cut_vertices[id_j];
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
                CutVertex cut_v;
                lerp(cut_v, vi, vj, t, d ^ 1);
                cut_v.c[d] = z;
                cut_v.q[d] = 0;
                cut_v.b[d] = true;
                if (cut_v.b[d ^ 1]) {
                    auto grid_id = cut_v.c[1] * row_size_ + cut_v.c[0];
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

bool CutMesh::is_neighbor_face(CutMesh::FaceRef face0,
                               CutMesh::FaceRef face1) const {
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

    min_x -= kKernelRange * delta_x_;
    max_x += kKernelRange * delta_x_;
    min_y -= kKernelRange * delta_x_;
    max_y += kKernelRange * delta_x_;

    h = face1->half_edge;
    do {
        auto p = h->vertex->position;
        if (p.x() > min_x && p.x() < max_x && p.y() > min_y && p.y() < max_y)
            return true;
    } while ((h = h->next) != face1->half_edge);
    return false;
}

void CutMesh::Face::calculate_center() {
    center.setZero();
    int cnt = 0;
    auto h = half_edge;
    do {
        center += h->vertex->position;
        ++cnt;
    } while ((h = h->next) != half_edge);
    center /= static_cast<Real>(cnt);
}

Vec2 CutMesh::Vertex::project_convex_velocity(const Vec2& x,
                                              const Vec2& v) const {
    if (!convex) return project(normal, v);
    auto h = half_edge;
    Vec2 normal0 = Vec2::Zero();
    Vec2 normal1 = Vec2::Zero();
    do {
        if (!h->is_boundary) continue;
        if (!normal0.isZero()) {
            normal1 = h->normal;
            break;
        }
        normal0 = h->normal;
    } while ((h = h->twin->next) != half_edge);
    if (normal0.dot(v) > 0 || normal1.dot(v) > 0) return v;
    Vec2 perp = normal0 + normal1;
    perp = Vec2(perp.y(), -perp.x());
    Vec2 d = x - position;
    if (perp.dot(d) * perp.dot(normal0) > 0) return project(normal0, v);
    return project(normal1, v);
}
