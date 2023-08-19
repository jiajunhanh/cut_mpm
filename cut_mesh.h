#pragma once

#include <iostream>
#include <vector>

#include "half_edge.h"
#include "mpm_config.h"

class Vertex {
   public:
    int c[2];
    float q[2];
    bool b[2];

    [[nodiscard]] float coord(int d) const {
        return static_cast<float>(c[d]) + q[d];
    }

    [[nodiscard]] int node_id() const { return c[1] * (n_grid + 1) + c[0]; }
};

class Edge {
   public:
    size_t i{};  // point id
    size_t j{};  // point id
};

HalfEdgeMesh construct_cut_mesh(
    const std::vector<std::array<float, 2>> &vertices,
    const std::vector<std::array<size_t, 2>> &edges);