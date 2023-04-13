#include "half_edge.h"
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

constexpr int dimension = 2;
constexpr float dx = 1.0f;

class Vertex {
  public:
    std::array<int, dimension> c{};
    std::array<float, dimension> q{};
    std::array<bool, dimension> b{};

    [[nodiscard]] float coord(int d) const {
        return static_cast<float>(c[d]) + q[d];
    }
};

class Edge {
  public:
    size_t i{}; // point id
    size_t j{}; // point id
};

Vertex lerp(const Vertex &vi, const Vertex &vj, float t, int d_skip = -1) {
    Vertex v;
    for (int d = 0; d < dimension; ++d) {
        if (d == d_skip) {
            continue;
        }
        auto coord_i = vi.coord(d);
        auto coord_j = vj.coord(d);
        auto coord = coord_i + (coord_j - coord_i) * t;
        v.c[d] = static_cast<int>(std::floor(coord));
        v.q[d] = coord - std::floor(coord);
        v.b[d] = v.q[d] == 0.0f;
    }
    return v;
}

std::vector<Edge> compute_cut_edges(std::vector<Vertex> &vertices,
                                    const std::vector<Edge> &edges) {
    std::vector<Edge> cut_edges;
    for (auto &e : edges) {
        std::vector<std::pair<float, size_t>> intersection_points;
        const auto &vi = vertices[e.i];
        const auto &vj = vertices[e.j];
        for (int d = 0; d < dimension; ++d) {
            int z_min;
            int z_max;
            if (vi.c[d] < vj.c[d] || vi.q[d] < vj.q[d]) {
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
                intersection_points.emplace_back(t, vertices.size());
                auto cut_v = lerp(vi, vj, t, d);
                cut_v.c[d] = z;
                cut_v.q[d] = 0;
                cut_v.b[d] = true;
                vertices.emplace_back(cut_v);
            }
        }
        std::sort(intersection_points.begin(), intersection_points.end(),
                  [](const std::pair<float, size_t> &a,
                     const std::pair<float, size_t> &b) {
                      return a.first < b.first;
                  });
        if (intersection_points.empty()) {
            cut_edges.emplace_back(e.i, e.j);
            continue;
        }
        cut_edges.emplace_back(e.i, intersection_points.front().second);
        for (size_t i = 0; i + 1 < intersection_points.size(); ++i) {
            cut_edges.emplace_back(intersection_points[i].second,
                                   intersection_points[i + 1].second);
        }
        cut_edges.emplace_back(intersection_points.back().second, e.j);
    }
    return cut_edges;
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
        std::array<float, dimension> x{};
        Vertex v;
        for (int d = 0; d < dimension; ++d) {
            ifs >> x[d];
            x[d] *= dx;
            v.c[d] = static_cast<int>(std::floor(x[d]));
            v.q[d] = x[d] - std::floor(x[d]);
            v.b[d] = v.q[d] == 0.0f;
        }
        vertices.emplace_back(v);
    }
    while (m--) {
        size_t i, j;
        ifs >> i >> j;
        edges.emplace_back(i, j);
    }
    auto cut_edges = compute_cut_edges(vertices, edges);

    std::cout << vertices.size() << ' ' << cut_edges.size() << '\n';
    for (size_t i = 0; i < vertices.size(); ++i) {
        const auto &v = vertices[i];
        std::cout << i << ":\n";
        for (int d = 0; d < dimension; ++d) {
            std::cout << v.c[d] << ' ' << v.q[d] << ' ' << v.b[d] << '\n';
        }
    }
    for (const auto &e : cut_edges) {
        std::cout << e.i << ' ' << e.j << '\n';
    }
    return 0;
}
