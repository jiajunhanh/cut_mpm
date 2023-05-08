#pragma once

#include <iostream>
#include <vector>

constexpr int n_grid = 8;
constexpr int n_grid_nodes = (n_grid + 1) * (n_grid + 1);
constexpr float dx = 1.0f / n_grid;

class Vertex {
  public:
    int c[2];
    float q[2];
    bool b[2];

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

std::pair<std::vector<Vertex>, std::vector<Edge>>
compute_cut_vertices_and_edges(const std::vector<Vertex> &vertices,
                               const std::vector<Edge> &edges);