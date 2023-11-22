#include "mpm.h"

#include <Eigen/SVD>
#include <random>

namespace {

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

Real collision_time(const Vec2& px, const Vec2& pv, const Vec2& o,
                    const Vec2& d) {
    constexpr Real kInf = std::numeric_limits<Real>::infinity();
    auto denominator = d.x() * pv.y() - d.y() * pv.x();
    if (denominator == 0) return kInf;
    Real t1 =
        (px.x() * pv.y() - px.y() * pv.x() + o.y() * pv.x() - o.x() * pv.y()) /
        denominator;
    if (t1 < 0 || t1 > 1) return kInf;
    Real t0 =
        (o.y() * d.x() - o.x() * d.y() + px.x() * d.y() - px.y() * d.x()) /
        denominator;
    if (t0 < 0) return kInf;
    return t0;
}

}  // namespace

MPM::MPM(const std::shared_ptr<CutMesh>& cut_mesh, int quality, int material)
    : quality_(static_cast<int>(std::pow(2, quality - 1))),
      grid_size_(8 * quality_),
      row_size_(grid_size_ + 1),
      n_particles_(32 * quality_ * quality_),
      material_(material),
      delta_t_(Real{2e-3} / static_cast<Real>(quality_)),
      delta_x_(Real{1.0} / static_cast<Real>(grid_size_)),
      margin_(delta_x_ / 32),
      particle_volume_(delta_x_ * delta_x_ * Real{0.25}),
      particle_density_(1.0),
      k_particle_mass_(particle_volume_ * particle_density_),
      cut_mesh_(cut_mesh) {
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
    particles_.resize(n_particles_);
    for (Particle& p : particles_) {
        p.x.x() = random_real() * Real{0.35} + Real{0.2};
        p.x.y() = random_real() * Real{0.35} + Real{0.02};
    }
}

void MPM::update() {
#pragma omp parallel for default(none), shared(std::cout)
    for (GridNode& g : nodes_) {
        g.m = 0.0;
        g.v.setZero();
    }
    // P2G
#pragma omp parallel for default(none)
    for (Particle& p : particles_) {
        p.F = (Mat2::Identity() + delta_t_ * p.C.block<2, 2>(0, 1)) * p.F;
        Real hardening_coefficient = 0.5;
        Real mu = kMu0 * hardening_coefficient;
        Real lambda = kLambda0 * hardening_coefficient;
        auto [U, sig, V] = svd(p.F);
        Real J = p.F.determinant();
        if (material_ == 0) {
            mu = 0;
            p.F = Mat2::Identity() * std::sqrt(J);
        }
        Mat2 PF = 2.0 * mu * (p.F - U * V.transpose()) * p.F.transpose() +
                  Mat2::Identity() * lambda * J * (J - 1.0);
        Mat3 M_inv = p.M_inv;
        Mat23 affine = Mat23::Zero();
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j)
                affine.row(i) += M_inv.row(j + 1) * PF(i, j);
        }
        affine *= -delta_t_ * particle_volume_;
        affine += k_particle_mass_ * p.C;
        auto enclosing_face = cut_mesh_->get_enclosing_face(p.x);
        const auto& neighbor_nodes = enclosing_face->neighbor_nodes;
        const auto& neighbor_node_sides = enclosing_face->neighbor_node_sides;
        int n_neighbor_nodes = static_cast<int>(neighbor_nodes.size());
        Real weight_sum = 0.0;
        std::vector<Real> weights(n_neighbor_nodes);
        for (int i = 0; i < n_neighbor_nodes; ++i) {
            if (!neighbor_node_sides[i]) continue;
            weights[i] =
                interpolate((nodes_[neighbor_nodes[i]].vertex->position - p.x) *
                            grid_size_);
            weight_sum += weights[i];
        }
        for (int i = 0; i < n_neighbor_nodes; ++i) {
            if (!neighbor_node_sides[i]) continue;
            if (weights[i] <= 0) continue;
            auto& g = nodes_[neighbor_nodes[i]];
            auto node_position = g.vertex->position;
            Vec3 distance(0.0, node_position.x() - p.x.x(),
                          node_position.y() - p.x.y());
            distance *= static_cast<Real>(grid_size_);
            distance.x() = 1.0;
            Real weight = weights[i] / weight_sum;
            distance *= delta_x_;
            Vec2 v_add = weight * affine * distance;
#pragma omp atomic
            g.v[0] += v_add[0];
#pragma omp atomic
            g.v[1] += v_add[1];
#pragma omp atomic
            g.m += weight * k_particle_mass_;
        }
    }
    // Grid update
#pragma omp parallel for default(none) shared(kGravity)
    for (GridNode& g : nodes_) {
        if (g.m <= 0) continue;
        g.v /= g.m;
        g.v += delta_t_ * kGravity * 30;
        if (g.vertex->normal != Vec2::Zero()) {
            auto dot = g.v.dot(g.vertex->normal);
            dot = std::min(dot, Real{0});
            g.v = g.v - dot * g.vertex->normal;
        }
        if (g.vertex->on_boundary) continue;
        auto id = g.vertex->id;
        int x = id % row_size_;
        int y = id / row_size_;
        if (x < 3 && g.v[0] < 0) g.v[0] = 0;
        if (x > grid_size_ - 3 && g.v[0] > 0) g.v[0] = 0;
        if (y < 3 && g.v[1] < 0) g.v[1] = 0;
        if (y > grid_size_ - 3 && g.v[1] > 0) g.v[1] = 0;
    }
    // G2P
#pragma omp parallel for default(none)
    for (Particle& p : particles_) {
        auto enclosing_face = cut_mesh_->get_enclosing_face(p.x);
        const auto& neighbor_nodes = enclosing_face->neighbor_nodes;
        const auto& neighbor_node_sides = enclosing_face->neighbor_node_sides;
        int n_neighbor_nodes = static_cast<int>(neighbor_nodes.size());
        Real weight_sum = 0.0;
        std::vector<Real> weights(n_neighbor_nodes);
        for (int i = 0; i < n_neighbor_nodes; ++i) {
            if (!neighbor_node_sides[i]) continue;
            weights[i] =
                interpolate((nodes_[neighbor_nodes[i]].vertex->position - p.x) *
                            grid_size_);
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
            distance *= static_cast<Real>(grid_size_);
            distance.x() = 1.0;
            Real weight = weights[i] / weight_sum;
            distance *= delta_x_;
            new_v += weight * g.v;
            new_C += weight * (g.v * distance.transpose());
            new_M += weight * distance * distance.transpose();
        }
        p.v = new_v;
        p.M_inv = new_M.inverse();
        p.C = new_C * new_M.inverse();

        bool colliding = true;
        Vec2 old_x = p.x;
        for (int i = 0; i < 8; ++i) {
            Real min_t = std::numeric_limits<Real>::infinity();
            Vec2 normal = Vec2::Zero();
            Real distance = std::numeric_limits<Real>::infinity();
            for (int id : enclosing_face->neighbor_boundaries) {
                auto h = begin(cut_mesh_->half_edges()) + id;
                if (p.v.dot(h->normal) > 0) continue;
                auto p0 = h->vertex->position;
                auto p1 = h->twin->vertex->position;
                auto t = collision_time(p.x, p.v, p0, p1 - p0);
                if (t < min_t) {
                    min_t = t;
                    normal = h->normal;
                    distance = (p.x - p0).dot(normal);
                }
            }
            if (min_t > delta_t_) {
                colliding = false;
                break;
            }
            p.x +=
                p.v * min_t * std::max(Real{0}, distance - margin_) / distance;
            p.v = margin_ / delta_t_ * normal;
        }
        if (!colliding) {
            p.x += p.v * delta_t_;
        }
        p.v = (p.x - old_x) / delta_t_;
    }
}
