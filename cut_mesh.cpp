#include "cut_mesh.h"
#include <algorithm>
#include <cmath>

constexpr static int n_grid_nodes = (n_grid + 1) * (n_grid + 1);

static auto lerp(Vertex &v, const Vertex &vi, const Vertex &vj, float t,
                 int d) {
    auto coord_i = vi.coord(d);
    auto coord_j = vj.coord(d);
    auto coord = coord_i + (coord_j - coord_i) * t;
    v.c[d] = static_cast<int>(std::floor(coord));
    v.q[d] = coord - std::floor(coord);
    v.b[d] = v.q[d] == 0.0f;
}

static void add_grid_edges(std::vector<Vertex> &cut_vertices,
                           std::vector<Edge> &cut_edges) {
    std::vector<size_t> grid_cut_vertices[n_grid * n_grid][2];
    for (size_t i = n_grid_nodes; i < cut_vertices.size(); ++i) {
        const auto &v = cut_vertices[i];
        for (int d = 0; d < 2; ++d) {
            if (v.b[d]) {
                grid_cut_vertices[v.c[1] * n_grid + v.c[0]][d ^ 1].emplace_back(
                    i);
            }
        }
    }
    for (int y = 0; y < n_grid; ++y) {
        for (int x = 0; x < n_grid; ++x) {
            auto node_i = y * (n_grid + 1) + x;
            for (int d = 0; d < 2; ++d) {
                auto &v = grid_cut_vertices[y * n_grid + x][d];
                int node_j = node_i;
                if (d == 0) {
                    node_j += 1;
                } else {
                    node_j += n_grid + 1;
                }
                if (v.empty()) {
                    cut_edges.emplace_back(node_i, node_j);
                    continue;
                }
                std::sort(v.begin(), v.end(), [&](size_t i, size_t j) {
                    return cut_vertices[i].q[d] < cut_vertices[j].q[d];
                });
                cut_edges.emplace_back(node_i, v.front());
                for (size_t i = 0; i < v.size() - 1; ++i) {
                    cut_edges.emplace_back(v[i], v[i + 1]);
                }
                cut_edges.emplace_back(v.back(), node_j);
            }
        }
    }
    for (int i = 0; i < n_grid; ++i) {
        cut_edges.emplace_back(n_grid + (n_grid + 1) * i,
                               n_grid + (n_grid + 1) * (i + 1));
        cut_edges.emplace_back(n_grid * (n_grid + 1) + i,
                               n_grid * (n_grid + 1) + i + 1);
    }
}

static std::pair<std::vector<Vertex>, std::vector<Edge>>
compute_cut_vertices_and_edges(
    const std::vector<std::array<float, 2>> &vertices,
    const std::vector<std::array<size_t, 2>> &edges) {
    std::vector<Vertex> cut_vertices;
    std::vector<Edge> cut_edges;
    for (int y = 0; y <= n_grid; ++y) {
        for (int x = 0; x <= n_grid; ++x) {
            Vertex v{};
            v.c[0] = x;
            v.c[1] = y;
            v.b[0] = v.b[1] = true;
            cut_vertices.emplace_back(v);
        }
    }
    for (const auto &v : vertices) {
        Vertex cut_vertex{};
        for (int d = 0; d < 2; ++d) {
            auto x = v[d] * n_grid;
            cut_vertex.c[d] = static_cast<int>(std::floor(x));
            cut_vertex.q[d] = x - std::floor(x);
            cut_vertex.b[d] = (cut_vertex.q[d] == 0.0f);
        }
        cut_vertices.emplace_back(cut_vertex);
    }
    for (const auto &e : edges) {
        std::vector<std::pair<float, size_t>> intersection_points;
        const auto vi = cut_vertices[e[0] + n_grid_nodes];
        const auto vj = cut_vertices[e[1] + n_grid_nodes];
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
            for (int z = z_min; z <= z_max; ++z) {
                auto t = (static_cast<float>(z - vi.c[d]) - vi.q[d]) /
                         (static_cast<float>(vj.c[d] - vi.c[d]) +
                          (vj.q[d] - vi.q[d]));
                Vertex cut_v{};
                lerp(cut_v, vi, vj, t, d ^ 1);
                cut_v.c[d] = z;
                cut_v.q[d] = 0;
                cut_v.b[d] = true;
                if (cut_v.b[d ^ 1]) {
                    intersection_points.emplace_back(t, cut_v.node_id());
                    continue;
                }
                cut_vertices.emplace_back(cut_v);
                intersection_points.emplace_back(t, cut_vertices.size() - 1);
            }
        }
        std::sort(intersection_points.begin(), intersection_points.end(),
                  [](const std::pair<float, size_t> &a,
                     const std::pair<float, size_t> &b) {
                      return a.first < b.first;
                  });
        if (intersection_points.empty()) {
            cut_edges.emplace_back(e[0] + n_grid_nodes, e[1] + n_grid_nodes);
            continue;
        }
        cut_edges.emplace_back(e[0] + n_grid_nodes,
                               intersection_points.front().second);
        for (size_t i = 0; i + 1 < intersection_points.size(); ++i) {
            if (intersection_points[i].second ==
                intersection_points[i + 1].second) {
                continue;
            }
            cut_edges.emplace_back(intersection_points[i].second,
                                   intersection_points[i + 1].second);
        }
        cut_edges.emplace_back(intersection_points.back().second,
                               e[1] + n_grid_nodes);
    }
    add_grid_edges(cut_vertices, cut_edges);
    return {cut_vertices, cut_edges};
}

HalfEdgeMesh
construct_cut_mesh(const std::vector<std::array<float, 2>> &vertices,
                   const std::vector<std::array<size_t, 2>> &edges) {
    auto [cut_vertices, cut_edges] =
        compute_cut_vertices_and_edges(vertices, edges);

    HalfEdgeMesh cut_mesh;
    std::vector<HalfEdgeMesh::VertexRef> vertex_map;
    std::vector<std::vector<HalfEdgeMesh::HalfEdgeRef>> vertex_half_edges;
    vertex_map.reserve(cut_vertices.size());
    vertex_half_edges.resize(cut_vertices.size());

    for (const auto &cut_vertex : cut_vertices) {
        auto v = cut_mesh.emplace_vertex();
        for (int d = 0; d < 2; ++d) {
            v->position[d] =
                static_cast<float>(cut_vertex.c[d]) + cut_vertex.q[d];
            v->position[d] /= n_grid;
        }
        vertex_map.emplace_back(v);
    }

    for (const auto &cut_edge : cut_edges) {
        auto h0 = cut_mesh.emplace_half_edge();
        auto h1 = cut_mesh.emplace_half_edge();
        auto v0 = vertex_map[cut_edge.i];
        auto v1 = vertex_map[cut_edge.j];
        auto e = cut_mesh.emplace_edge();
        h0->set_tnvef(h1, h0->next, v0, e, h0->face);
        h1->set_tnvef(h0, h1->next, v1, e, h1->face);
        v0->half_edge = h0;
        v1->half_edge = h1;
        e->half_edge = h0;
        vertex_half_edges[v0->id].emplace_back(h0);
        vertex_half_edges[v1->id].emplace_back(h1);
    }

    for (const auto &v : cut_mesh.vertices) {
        auto &half_edges = vertex_half_edges[v.id];
        std::vector<float> angles(half_edges.size());
        for (size_t i = 0; i < half_edges.size(); ++i) {
            auto h0 = half_edges[i];
            auto h1 = h0->twin;
            auto v1 = h1->vertex;
            angles[i] = std::atan2(v1->position[1] - v.position[1],
                                   v1->position[0] - v.position[0]);
        }
        std::vector<size_t> sorted_idx(half_edges.size());
        for (size_t i = 0; i < half_edges.size(); ++i) {
            sorted_idx[i] = i;
        }
        std::sort(sorted_idx.begin(), sorted_idx.end(),
                  [&angles = angles](size_t i, size_t j) {
                      return angles[i] < angles[j];
                  });
        for (size_t i = 0; i < half_edges.size(); ++i) {
            size_t j = (i + 1) % half_edges.size();
            auto h0 = half_edges[sorted_idx[i]];
            auto h1 = half_edges[sorted_idx[j]];
            h1->twin->next = h0;
        }
    }

    for (auto half_edge = cut_mesh.half_edges.begin();
         half_edge != cut_mesh.half_edges.end(); ++half_edge) {
        if (half_edge->face != cut_mesh.faces.end()) {
            continue;
        }
        auto f = cut_mesh.emplace_face();
        f->half_edge = half_edge;
        auto h = half_edge;
        do {
            h->face = f;
            h = h->next;
        } while (h != half_edge);
    }

    return cut_mesh;
}
