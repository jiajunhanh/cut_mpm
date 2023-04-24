#include "half_edge.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

constexpr int n_grid = 3;
constexpr int n_grid_nodes = (n_grid + 1) * (n_grid + 1);
constexpr float dx = 1.0f / n_grid;

class Vertex {
  public:
    std::array<int, 2> c{};
    std::array<float, 2> q{};
    std::array<bool, 2> b{};

    std::vector<size_t> edges;

    [[nodiscard]] float coord(int d) const {
        return static_cast<float>(c[d]) + q[d];
    }

    [[nodiscard]] int node_id() const { return c[1] * (n_grid + 1) + c[0]; }
};

class Edge {
  public:
    size_t i{}; // point id
    size_t j{}; // point id
};

auto lerp(Vertex &v, const Vertex &vi, const Vertex &vj, float t, int d) {
    auto coord_i = vi.coord(d);
    auto coord_j = vj.coord(d);
    auto coord = coord_i + (coord_j - coord_i) * t;
    v.c[d] = static_cast<int>(std::floor(coord));
    v.q[d] = coord - std::floor(coord);
    v.b[d] = v.q[d] == 0.0f;
}

void add_grid_edges(std::vector<Vertex> &cut_vertices,
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
    for (auto v : cut_vertices) {
        v.edges.clear();
    }
    for (size_t i = 0; i < cut_edges.size(); ++i) {
        const auto &e = cut_edges[i];
        cut_vertices[e.i].edges.push_back(i);
        cut_vertices[e.j].edges.push_back(i);
    }
}

std::pair<std::vector<Vertex>, std::vector<Edge>>
compute_cut_vertices_and_edges(const std::vector<Vertex> &vertices,
                               const std::vector<Edge> &edges) {
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
    cut_vertices.insert(cut_vertices.end(), vertices.cbegin(), vertices.cend());
    for (auto &e : edges) {
        std::vector<std::pair<float, size_t>> intersection_points;
        const auto &vi = vertices[e.i];
        const auto &vj = vertices[e.j];
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
            cut_edges.emplace_back(e.i + n_grid_nodes, e.j + n_grid_nodes);
            continue;
        }
        cut_edges.emplace_back(e.i + n_grid_nodes,
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
                               e.j + n_grid_nodes);
    }
    add_grid_edges(cut_vertices, cut_edges);
    return {cut_vertices, cut_edges};
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        return 0;
    }
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::ifstream ifs(argv[1]);
    size_t n, m;
    ifs >> n >> m;
    while (n--) {
        float x[2]{};
        Vertex v{};
        for (int d = 0; d < 2; ++d) {
            ifs >> x[d];
            assert(x[d] >= 0.0f && x[d] < 1.0f);
            x[d] *= n_grid;
            v.c[d] = static_cast<int>(std::floor(x[d]));
            v.q[d] = x[d] - std::floor(x[d]);
            v.b[d] = (v.q[d] == 0.0f);
            // std::cout << v.c[d] << ' ' << v.q[d] << ' ' << v.b[d] << '\n';
        }
        vertices.emplace_back(v);
    }
    while (m--) {
        size_t i, j;
        ifs >> i >> j;
        edges.emplace_back(i, j);
    }
    auto [cut_vertices, cut_edges] =
        compute_cut_vertices_and_edges(vertices, edges);
    std::cout << cut_vertices.size() << ' ' << cut_edges.size() << '\n';
    for (size_t i = 0; i < cut_vertices.size(); ++i) {
        const auto &v = cut_vertices[i];
        std::cout << i << ":\n";
        for (int d = 0; d < 2; ++d) {
            std::cout << v.c[d] << ' ' << v.q[d] << ' ' << v.b[d] << '\n';
        }
    }
    for (const auto &e : cut_edges) {
        std::cout << e.i << ' ' << e.j << '\n';
    }
    return 0;
}
