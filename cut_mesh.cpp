#include "cut_mesh.h"

#include <algorithm>
#include <cmath>

constexpr static int kNumberOfGridNodes = (kGridSize + 1) * (kGridSize + 1);

static void lerp(Vertex &v, const Vertex &vi, const Vertex &vj, float t,
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
    std::vector<int> grid_cut_vertices[kGridSize * kGridSize][2];
    for (int i = kNumberOfGridNodes,
             size = static_cast<int>(cut_vertices.size());
         i < size; ++i) {
        const auto &v = cut_vertices[i];
        for (int d = 0; d < 2; ++d) {
            if (v.b[d]) {
                grid_cut_vertices[v.c[1] * kGridSize + v.c[0]][d ^ 1]
                    .emplace_back(i);
            }
        }
    }
    for (int y = 0; y < kGridSize; ++y) {
        for (int x = 0; x < kGridSize; ++x) {
            auto node_i = y * (kGridSize + 1) + x;
            for (int d = 0; d < 2; ++d) {
                auto &v = grid_cut_vertices[y * kGridSize + x][d];
                int node_j = node_i + d * kGridSize + 1;
                if (v.empty()) {
                    cut_edges.emplace_back(node_i, node_j);
                    continue;
                }
                std::sort(begin(v), end(v), [&](int i, int j) {
                    return cut_vertices[i].q[d] < cut_vertices[j].q[d];
                });
                cut_edges.emplace_back(node_i, v.front());
                for (int i = 0, size = static_cast<int>(v.size()); i < size - 1;
                     ++i) {
                    cut_edges.emplace_back(v[i], v[i + 1]);
                }
                cut_edges.emplace_back(v.back(), node_j);
            }
        }
    }
    for (int i = 0; i < kGridSize; ++i) {
        cut_edges.emplace_back(i * (kGridSize + 1) + kGridSize,
                               (i + 1) * (kGridSize + 1) + kGridSize);
        cut_edges.emplace_back(kGridSize * (kGridSize + 1) + i,
                               kGridSize * (kGridSize + 1) + i + 1);
    }
}

static std::pair<std::vector<Vertex>, std::vector<Edge>>
compute_cut_vertices_and_edges(
    const std::vector<std::array<float, 2>> &vertices,
    const std::vector<std::array<int, 2>> &edges) {
    std::vector<Vertex> cut_vertices;
    std::vector<Edge> cut_edges;
    for (int y = 0; y <= kGridSize; ++y) {
        for (int x = 0; x <= kGridSize; ++x) {
            Vertex v{};
            v.c[0] = x;
            v.c[1] = y;
            v.b[0] = v.b[1] = true;
            cut_vertices.emplace_back(v);
        }
    }
    for (const auto &v : vertices) {
        Vertex cut_vertex;
        for (int d = 0; d < 2; ++d) {
            auto x = v[d] * kGridSize;
            cut_vertex.c[d] = static_cast<int>(std::floor(x));
            cut_vertex.q[d] = x - std::floor(x);
            cut_vertex.b[d] = (cut_vertex.q[d] == 0.0f);
        }
        cut_vertices.emplace_back(cut_vertex);
    }
    for (const auto &e : edges) {
        std::vector<std::pair<float, int>> intersection_points;
        const auto vi = cut_vertices[e[0] + kNumberOfGridNodes];
        const auto vj = cut_vertices[e[1] + kNumberOfGridNodes];
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
            auto distance =
                static_cast<float>(vj.c[d] - vi.c[d]) + (vj.q[d] - vi.q[d]);
            for (int z = z_min; z <= z_max; ++z) {
                auto t = (static_cast<float>(z - vi.c[d]) - vi.q[d]) / distance;
                Vertex cut_v;
                lerp(cut_v, vi, vj, t, d ^ 1);
                cut_v.c[d] = z;
                cut_v.q[d] = 0;
                cut_v.b[d] = true;
                if (cut_v.b[d ^ 1]) {
                    intersection_points.emplace_back(t, cut_v.node_id());
                    continue;
                }
                cut_vertices.emplace_back(cut_v);
                intersection_points.emplace_back(
                    t, static_cast<int>(cut_vertices.size()) - 1);
            }
        }
        std::sort(
            begin(intersection_points), end(intersection_points),
            [](const auto &a, const auto &b) { return a.first < b.first; });
        if (intersection_points.empty()) {
            cut_edges.emplace_back(e[0] + kNumberOfGridNodes,
                                   e[1] + kNumberOfGridNodes);
            continue;
        }
        cut_edges.emplace_back(e[0] + kNumberOfGridNodes,
                               intersection_points.front().second);
        for (int i = 0, size = static_cast<int>(intersection_points.size());
             i + 1 < size; ++i) {
            if (intersection_points[i].second ==
                intersection_points[i + 1].second) {
                continue;
            }
            cut_edges.emplace_back(intersection_points[i].second,
                                   intersection_points[i + 1].second);
        }
        cut_edges.emplace_back(intersection_points.back().second,
                               e[1] + kNumberOfGridNodes);
    }
    add_grid_edges(cut_vertices, cut_edges);
    return {cut_vertices, cut_edges};
}

HalfEdgeMesh construct_cut_mesh(
    const std::vector<std::array<float, 2>> &vertices,
    const std::vector<std::array<int, 2>> &edges) {
    auto [cut_vertices, cut_edges] =
        compute_cut_vertices_and_edges(vertices, edges);

    HalfEdgeMesh cut_mesh;
    std::vector<HalfEdgeMesh::VertexRef> mesh_vertices;
    std::vector<std::vector<HalfEdgeMesh::HalfEdgeRef>> vertex_half_edges;
    mesh_vertices.reserve(cut_vertices.size());
    vertex_half_edges.resize(cut_vertices.size());

    for (const auto &cut_vertex : cut_vertices) {
        auto v = cut_mesh.emplace_vertex();
        for (int d = 0; d < 2; ++d)
            v->position[d] =
                static_cast<float>(cut_vertex.c[d]) + cut_vertex.q[d];
        v->position /= kGridSize;
        mesh_vertices.emplace_back(v);
    }

    for (const auto &cut_edge : cut_edges) {
        auto h0 = cut_mesh.emplace_half_edge();
        auto h1 = cut_mesh.emplace_half_edge();
        auto v0 = mesh_vertices[cut_edge.i];
        auto v1 = mesh_vertices[cut_edge.j];
        auto e = cut_mesh.emplace_edge();
        h0->set_tnvef(h1, h0->next, v0, e, h0->face);
        h1->set_tnvef(h0, h1->next, v1, e, h1->face);
        v0->half_edge = h0;
        v1->half_edge = h1;
        e->half_edge = h0;
        vertex_half_edges[cut_edge.i].emplace_back(h0);
        vertex_half_edges[cut_edge.j].emplace_back(h1);
    }

    for (int iv = 0, n_vertices = static_cast<int>(mesh_vertices.size());
         iv < n_vertices; ++iv) {
        auto v = mesh_vertices[iv];
        auto &half_edges = vertex_half_edges[iv];
        int n_half_edges = static_cast<int>(half_edges.size());
        std::vector<float> angles(n_half_edges);
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

    for (auto half_edge = begin(cut_mesh.half_edges()),
              half_edges_end = end(cut_mesh.half_edges());
         half_edge != half_edges_end; ++half_edge) {
        if (half_edge->face != cut_mesh.faces().end()) {
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
