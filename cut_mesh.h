#pragma once

#include <iostream>
#include <vector>

#include "half_edge.h"
#include "mpm_config.h"

struct Vertex {
    int c[2] = {};
    float q[2] = {};
    bool b[2] = {};

    [[nodiscard]] float coord(int d) const {
        return static_cast<float>(c[d]) + q[d];
    }

    [[nodiscard]] int node_id() const { return c[1] * (kGridSize + 1) + c[0]; }
};

struct Edge {
   public:
    int i = 0;  // point id
    int j = 0;  // point id
    Edge(int i_, int j_) : i(i_), j(j_) {}
};

HalfEdgeMesh construct_cut_mesh(
    const std::vector<std::array<float, 2>> &vertices,
    const std::vector<std::array<int, 2>> &edges);