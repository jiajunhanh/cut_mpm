#include "mpm.h"

#include <Eigen/SVD>
#include <random>

constexpr Real kParticleVolume = kDeltaX * kDeltaX * 0.25;
constexpr Real kParticleDensity = 1.0;
constexpr Real kParticleMass = kParticleVolume * kParticleDensity;
constexpr Real kYoungModulus = 1e3;
constexpr Real kPoissonRatio = 0.2;
constexpr Real kMu0 = kYoungModulus / (2.0 * (1.0 + kPoissonRatio));
constexpr Real kLambda0 = kYoungModulus * kPoissonRatio /
                          ((1.0 + kPoissonRatio) * (1.0 - 2.0 * kPoissonRatio));
const Vec2 kGravity(0.0, 1.0);

static Real random_real() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    // static std::mt19937 gen(42);
    static std::uniform_real_distribution<Real> dis(0.0, 1.0);
    return dis(gen);
}

static Real interpolate(Real x) {
    x = std::abs(x);
    /*if (x < 0.5) return 0.75 - x * x;
    if (x < 1.5) return 0.5 * (1.5 - x) * (1.5 - x);
    return 0.0;*/
    return std::max(Real{0.0}, Real{0.75} - Real{0.5} * x);
}

static Real interpolate(Vec2 x) {
    return interpolate(x.x()) * interpolate(x.y());
}

static auto svd(const Mat2& m) {
    Eigen::JacobiSVD<Mat2, Eigen::ComputeFullU | Eigen::ComputeFullV> svd(m);
    return std::tuple{svd.matrixU(), svd.singularValues(), svd.matrixV()};
}

MPM::MPM(const std::shared_ptr<CutMesh>& cut_mesh_)
    : cut_mesh_(cut_mesh_),
      face_visited_(cut_mesh_->faces().size()),
      is_neighbor_node_(cut_mesh_->vertices().size()) {
    int n_vertices = static_cast<int>(cut_mesh_->vertices().size());
    nodes_.resize(n_vertices);
    for (auto [i, v] = std::tuple(0, begin(cut_mesh_->vertices()));
         i < n_vertices; ++i, ++v) {
        nodes_[i].vertex = v;
        node_of_vertex_[v->id] = i;
    }
}

void MPM::initialize() {
    particles_.clear();
    particles_.resize(kParticleNumber);
    for (Particle& p : particles_) {
        p.x.x() = random_real() * Real{0.35} + Real{0.2};
        p.x.y() = random_real() * Real{0.35} + Real{0.2};
    }
}

void MPM::update() {
    for (GridNode& g : nodes_) {
        g.m = 0.0;
        g.v.setZero();
    }
    // P2G
    std::vector<int> neighbor_nodes;
    for (Particle& p : particles_) {
        p.F = (Mat2::Identity() + kDeltaT * p.C.block<2, 2>(0, 1)) * p.F;
        Real hardening_coefficient = 1.0;
        Real mu = kMu0 * hardening_coefficient;
        Real lambda = kLambda0 * hardening_coefficient;
        auto [U, sig, V] = svd(p.F);
        Real J = p.F.determinant();
        Mat2 PF = 2.0 * mu * (p.F - U * V.transpose()) * p.F.transpose() +
                  Mat2::Identity() * lambda * J * (J - 1.0);
        Mat3 M_inv = p.M_inv;
        Mat23 affine = Mat23::Zero();
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j)
                affine.row(i) += M_inv.row(j + 1) * PF(i, j);
        }
        affine *= -kDeltaT * kParticleVolume;
        affine += kParticleMass * p.C;
        // auto neighbor_nodes = get_neighbor_nodes(p.x);
        get_neighbor_nodes_in_place(p.x, neighbor_nodes);
        int n_neighbor_nodes = static_cast<int>(neighbor_nodes.size());
        Real weight_sum = 0.0;
        std::vector<Real> weights(n_neighbor_nodes);
        for (int i = 0; i < n_neighbor_nodes; ++i) {
            weights[i] = interpolate(
                (nodes_[neighbor_nodes[i]].vertex->position - p.x) * kGridSize);
            weight_sum += weights[i];
        }
        for (int i = 0; i < n_neighbor_nodes; ++i) {
            auto& g = nodes_[neighbor_nodes[i]];
            auto node_position = g.vertex->position;
            Vec3 distance(0.0, node_position.x() - p.x.x(),
                          node_position.y() - p.x.y());
            distance *= kGridSize;
            distance.x() = 1.0;
            Real weight = weights[i] / weight_sum;
            distance *= kDeltaX;
            g.v += weight * affine * distance;
            g.m += weight * kParticleMass;
        }
    }
    // Grid update
    for (GridNode& g : nodes_) {
        int id = g.vertex->id;
        int x = id % kRowSize;
        int y = id / kRowSize;
        if (g.m <= 0) continue;
        g.v /= g.m;
        g.v += kDeltaT * kGravity * 30;
        if (g.vertex->id >= kRowSize * kRowSize) continue;
        // if (g.vertex->id >= kRowSize * kRowSize)
        //    g.v[0] = std::max(Real{0}, g.v[0]);
        if (x < 3 && g.v[0] < 0) g.v[0] = 0;
        if (x > kGridSize - 3 && g.v[0] > 0) g.v[0] = 0;
        if (y < 3 && g.v[1] < 0) g.v[1] = 0;
        if (y > kGridSize - 3 && g.v[1] > 0) g.v[1] = 0;
    }
    // G2P
    for (Particle& p : particles_) {
        // auto neighbor_nodes = get_neighbor_nodes(p.x);
        get_neighbor_nodes_in_place(p.x, neighbor_nodes);
        int n_neighbor_nodes = static_cast<int>(neighbor_nodes.size());
        Real weight_sum = 0.0;
        std::vector<Real> weights(n_neighbor_nodes);
        for (int i = 0; i < n_neighbor_nodes; ++i) {
            weights[i] = interpolate(
                (nodes_[neighbor_nodes[i]].vertex->position - p.x) * kGridSize);
            weight_sum += weights[i];
        }
        Vec2 new_v = Vec2::Zero();
        Mat23 new_C = Mat23::Zero();
        Mat3 new_M = Mat3::Zero();
        for (int i = 0; i < n_neighbor_nodes; ++i) {
            auto& g = nodes_[neighbor_nodes[i]];
            auto node_position = g.vertex->position;
            Vec3 distance(0.0, node_position.x() - p.x.x(),
                          node_position.y() - p.x.y());
            distance *= kGridSize;
            distance.x() = 1.0;
            Real weight = weights[i] / weight_sum;
            distance *= kDeltaX;
            new_v += weight * g.v;
            new_C += weight * (g.v * distance.transpose());
            new_M += weight * distance * distance.transpose();
        }
        p.v = new_v;
        p.M_inv = new_M.inverse();
        p.C = new_C * new_M.inverse();
        p.x += kDeltaT * p.v;
    }
}

std::vector<int> MPM::get_neighbor_nodes(const Vec2& x) {
    std::vector<int> neighbor_nodes;
    get_neighbor_nodes_in_place(x, neighbor_nodes);
    return neighbor_nodes;
}

void MPM::get_neighbor_nodes_in_place(const Vec2& x,
                                      std::vector<int>& neighbor_nodes) {
    neighbor_nodes.clear();
    auto enclosing_face = cut_mesh_->get_enclosing_face(x);
    face_queue_.emplace(enclosing_face);
    while (!face_queue_.empty()) {
        auto face = face_queue_.front();
        face_queue_.pop();
        if (face_visited_[face->id]) continue;
        face_visited_[face->id] = true;
        neighbor_faces_.emplace_back(face);
        auto h = face->half_edge;
        do {
            auto f = h->twin->face;
            if (face_visited_[f->id]) continue;
            if (interpolate((f->center() - x) * kGridSize) <= 0) continue;
            face_queue_.emplace(f);
        } while ((h = h->next) != face->half_edge);
    }
    for (auto f : neighbor_faces_) {
        auto h = f->half_edge;
        do {
            if (is_neighbor_node_[h->vertex->id]) continue;
            if (interpolate(((h->vertex->position - x) * kGridSize)) <= 0)
                continue;
            is_neighbor_node_[h->vertex->id] = true;
            neighbor_nodes.emplace_back(node_of_vertex_.at(h->vertex->id));
        } while ((h = h->next) != f->half_edge);
    }
    std::fill(begin(face_visited_), end(face_visited_), false);
    std::fill(begin(is_neighbor_node_), end(is_neighbor_node_), false);
    neighbor_faces_.clear();
}
