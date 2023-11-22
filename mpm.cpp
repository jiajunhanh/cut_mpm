#include "mpm.h"

#include <Eigen/SVD>
#include <random>

namespace {

constexpr Real kParticleVolume = kDeltaX * kDeltaX * 0.25;
constexpr Real kParticleDensity = 1.0;
constexpr Real kParticleMass = kParticleVolume * kParticleDensity;
constexpr Real kYoungModulus = 1e3;
constexpr Real kPoissonRatio = 0.2;
constexpr Real kMu0 = kYoungModulus / (2.0 * (1.0 + kPoissonRatio));
constexpr Real kLambda0 = kYoungModulus * kPoissonRatio /
                          ((1.0 + kPoissonRatio) * (1.0 - 2.0 * kPoissonRatio));
const Vec2 kGravity(0.0, 1.0);

Real random_real() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    // static std::mt19937 gen(42);
    static std::uniform_real_distribution<Real> dis(0.0, 1.0);
    return dis(gen);
}

auto svd(const Mat2& m) {
    Eigen::JacobiSVD<Mat2, Eigen::ComputeFullU | Eigen::ComputeFullV> svd(m);
    return std::tuple{svd.matrixU(), svd.singularValues(), svd.matrixV()};
}

bool collision(const Vec2& px, Vec2& pv, const Vec2& o, const Vec2& d) {
    auto denominator = d.x() * pv.y() - d.y() * pv.x();
    if (denominator == 0) return false;
    Real t1 =
        (px.x() * pv.y() - px.y() * pv.x() + o.y() * pv.x() - o.x() * pv.y()) /
        denominator;
    if (t1 < -kMargin / kDeltaX || t1 > 1 + kMargin / kDeltaX) return false;
    Real t0 =
        (o.y() * d.x() - o.x() * d.y() + px.x() * d.y() - px.y() * d.x()) /
        denominator;
    if (t0 < 0 || (t0 - kDeltaT) * pv.norm() > kMargin) return false;
    pv *= std::max(Real{0}, t0 - kMargin / pv.norm()) / kDeltaT;
    return true;
}

}  // namespace

MPM::MPM(const std::shared_ptr<CutMesh>& cut_mesh_) : cut_mesh_(cut_mesh_) {
    int n_vertices = static_cast<int>(cut_mesh_->vertices().size());
    nodes_.resize(n_vertices);
    for (auto [i, v] = std::tuple(0, begin(cut_mesh_->vertices()));
         i < n_vertices; ++i, ++v)
        nodes_[i].vertex = v;
    cut_mesh_->calculate_neighbor_nodes_and_boundaries_of_faces();
    cut_mesh_->calculate_node_normals();
}

void MPM::initialize() {
    particles_.clear();
    particles_.resize(kParticleNumber);
    for (Particle& p : particles_) {
        p.x.x() = random_real() * Real{0.35} + Real{0.2};
        p.x.y() = random_real() * Real{0.35} + Real{0.02};
    }
}

void MPM::update() {
    for (GridNode& g : nodes_) {
        g.m = 0.0;
        g.v.setZero();
    }
    // P2G
    for (Particle& p : particles_) {
        p.F = (Mat2::Identity() + kDeltaT * p.C.block<2, 2>(0, 1)) * p.F;
        Real hardening_coefficient = 0.5;
        Real mu = kMu0 * hardening_coefficient;
        Real lambda = kLambda0 * hardening_coefficient;
        auto [U, sig, V] = svd(p.F);
        Real J = p.F.determinant();
        mu = 0;
        p.F = Mat2::Identity() * std::sqrt(J);
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
        auto enclosing_face = cut_mesh_->get_enclosing_face(p.x);
        const auto& neighbor_nodes = enclosing_face->neighbor_nodes;
        const auto& neighbor_node_sides = enclosing_face->neighbor_node_sides;
        int n_neighbor_nodes = static_cast<int>(neighbor_nodes.size());
        Real weight_sum = 0.0;
        std::vector<Real> weights(n_neighbor_nodes);
        for (int i = 0; i < n_neighbor_nodes; ++i) {
            if (!neighbor_node_sides[i]) continue;
            weights[i] = interpolate(
                (nodes_[neighbor_nodes[i]].vertex->position - p.x) * kGridSize);
            weight_sum += weights[i];
        }
        for (int i = 0; i < n_neighbor_nodes; ++i) {
            if (!neighbor_node_sides[i]) continue;
            if (weights[i] <= 0) continue;
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
        if (g.m <= 0) continue;
        g.v /= g.m;
        g.v += kDeltaT * kGravity * 30;
        if (g.vertex->normal != Vec2::Zero()) {
            auto dot = g.v.dot(g.vertex->normal);
            dot = std::min(dot, Real{0});
            g.v = g.v - dot * g.vertex->normal;
        }
        if (g.vertex->on_boundary) continue;
        auto id = g.vertex->id;
        int x = id % kRowSize;
        int y = id / kRowSize;
        if (x < 3 && g.v[0] < 0) g.v[0] = 0;
        if (x > kGridSize - 3 && g.v[0] > 0) g.v[0] = 0;
        if (y < 3 && g.v[1] < 0) g.v[1] = 0;
        if (y > kGridSize - 3 && g.v[1] > 0) g.v[1] = 0;
    }
    // G2P
    for (Particle& p : particles_) {
        auto enclosing_face = cut_mesh_->get_enclosing_face(p.x);
        const auto& neighbor_nodes = enclosing_face->neighbor_nodes;
        const auto& neighbor_node_sides = enclosing_face->neighbor_node_sides;
        int n_neighbor_nodes = static_cast<int>(neighbor_nodes.size());
        Real weight_sum = 0.0;
        std::vector<Real> weights(n_neighbor_nodes);
        for (int i = 0; i < n_neighbor_nodes; ++i) {
            if (!neighbor_node_sides[i]) continue;
            weights[i] = interpolate(
                (nodes_[neighbor_nodes[i]].vertex->position - p.x) * kGridSize);
            weight_sum += weights[i];
        }
        Vec2 new_v = Vec2::Zero();
        Mat23 new_C = Mat23::Zero();
        Mat3 new_M = Mat3::Zero();
        for (int i = 0; i < n_neighbor_nodes; ++i) {
            if (!neighbor_node_sides[i]) continue;
            if (weights[i] <= 0) continue;
            const auto& g = nodes_[neighbor_nodes[i]];
            Vec3 distance(0.0, g.vertex->position.x() - p.x.x(),
                          g.vertex->position.y() - p.x.y());
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

        Vec2 old_x = p.x;
        for (int i = 0; i < 8; ++i) {
            bool colliding = false;
            for (int id : enclosing_face->neighbor_boundaries) {
                auto h = begin(cut_mesh_->half_edges()) + id;
                if (p.v.dot(h->normal) > 0) continue;
                auto p0 = h->vertex->position;
                auto p1 = h->twin->vertex->position;
                if (collision(p.x, p.v, p0, p1 - p0)) {
                    colliding = true;
                    p.x += p.v * kDeltaT;
                    p.v = 2 * kMargin * h->normal / kDeltaT;
                }
            }
            if (!colliding) break;
        }
        p.x += p.v * kDeltaT;
        p.v = (p.x - old_x) / kDeltaT;
    }
}
