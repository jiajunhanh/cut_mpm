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
    return std::tuple{svd.matrixU(), svd.matrixV()};
}

Real collision_time(const Vec2& px, const Vec2& pv, const Vec2& o,
                    const Vec2& d) {
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
      delta_t_(Real{1.6e-3} / static_cast<Real>(quality_)),
      delta_x_(Real{1.0} / static_cast<Real>(grid_size_)),
      margin_(delta_x_ / 128),
      particle_volume_(delta_x_ * delta_x_ * Real{0.25}),
      particle_density_(1.0),
      particle_mass_(particle_volume_ * particle_density_),
      inv_delta_x_(Real{1.0} / delta_x_),
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
#pragma omp parallel for default(none)
    for (GridNode& g : nodes_) {
        g.m = 0.0;
        g.v.setZero();
    }
    // P2G
#pragma omp parallel for default(none)
    for (Particle& p : particles_) {
        p.F = (Mat2::Identity() + delta_t_ * p.C.block<2, 2>(0, 1)) * p.F;
        Real J = p.F.determinant();
        Mat2 PF;
        if (material_ != 0) {
            Real hardening_coefficient = 0.3;
            Real lambda = kLambda0 * hardening_coefficient;
            auto [U, V] = svd(p.F);
            Real mu = kMu0 * hardening_coefficient;
            PF = 2 * mu * (p.F - U * V.transpose()) * p.F.transpose() +
                 Mat2::Identity() * lambda * J * (J - 1.0);
        } else {
            Real hardening_coefficient = 1;
            Real lambda = kLambda0 * hardening_coefficient;
            PF = Mat2::Identity() * lambda * J * (J - 1.0);
            p.F = Mat2::Identity() * std::sqrt(J);
        }
        Mat3 M_inv = p.M_inv;
        Mat23 affine = Mat23::Zero();
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j)
                affine.row(i) += M_inv.row(j + 1) * PF(i, j);
        }
        affine *= -delta_t_ * particle_volume_;
        affine += particle_mass_ * p.C;
        Vec2i base = ((p.x * inv_delta_x_).array() - 0.5).cast<int>();
        if (!cut_mesh_->grid(base[1], base[0]).near_boundary) {
            Vec2 fx = p.x * inv_delta_x_ - base.cast<Real>();
            Vec2 weights[3] = {0.5 * (1.5 - fx.array()).square(),
                               0.75 - (fx.array() - 1.0).square(),
                               0.5 * (fx.array() - 0.5).square()};
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    auto offset = Vec2i(i, j);
                    auto weight = weights[i][0] * weights[j][1];
                    Vec2i grid_idx = base + offset;
                    auto& g = grid(grid_idx[0], grid_idx[1]);
                    auto node_position = g.vertex->position;
                    Vec3 distance(0.0, node_position.x() - p.x.x(),
                                  node_position.y() - p.x.y());
                    distance *= static_cast<Real>(grid_size_);
                    distance.x() = 1.0;
                    distance *= delta_x_;
                    Vec2 v_add = weight * affine * distance;
#pragma omp atomic
                    g.v[0] += v_add[0];
#pragma omp atomic
                    g.v[1] += v_add[1];
#pragma omp atomic
                    g.m += weight * particle_mass_;
                }
            }
            continue;
        }

        auto enclosing_face = cut_mesh_->get_enclosing_face(p.x);
        const auto& neighbor_nodes = enclosing_face->neighbor_nodes;
        int n_neighbor_nodes = static_cast<int>(neighbor_nodes.size());
        p.weight_sum = 0;
        p.weights.assign(n_neighbor_nodes, 0);
        for (int i = 0; i < n_neighbor_nodes; ++i) {
            auto v = nodes_[neighbor_nodes[i]].vertex;
            Vec2 d = v->position - p.x;
            auto w = interpolate(d * grid_size_);
            if (w <= 0) continue;
            if (!enclosing_face->near_convex) {
                p.weights[i] = w;
                p.weight_sum += w;
                continue;
            }
            auto min_t = kInf;
            for (auto id : enclosing_face->neighbor_boundaries) {
                auto h = begin(cut_mesh_->half_edges()) + id;
                if (h->vertex == v || h->twin->vertex == v) continue;
                auto p0 = h->vertex->position;
                auto p1 = h->twin->vertex->position;
                auto t = collision_time(p.x, d, p0, p1 - p0);
                min_t = std::min(min_t, t);
            }
            if (min_t < 1) continue;
            p.weights[i] = w;
            p.weight_sum += w;
        }
        for (int i = 0; i < n_neighbor_nodes; ++i) {
            if (p.weights[i] <= 0) continue;
            auto& g = nodes_[neighbor_nodes[i]];
            auto node_position = g.vertex->position;
            Vec3 distance(0.0, node_position.x() - p.x.x(),
                          node_position.y() - p.x.y());
            distance *= static_cast<Real>(grid_size_);
            distance.x() = 1.0;
            Real weight = p.weights[i] / p.weight_sum;
            distance *= delta_x_;
            Vec2 v_add = weight * affine * distance;
#pragma omp atomic
            g.v[0] += v_add[0];
#pragma omp atomic
            g.v[1] += v_add[1];
#pragma omp atomic
            g.m += weight * particle_mass_;
        }
    }
    // Grid update
#pragma omp parallel for default(none) shared(kGravity)
    for (GridNode& g : nodes_) {
        if (g.m <= 0) continue;
        g.v /= g.m;
        g.v += delta_t_ * kGravity * 30;
        if (!g.vertex->normal.isZero() && !g.vertex->convex)
            g.v = project(g.vertex->normal, g.v);
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
        Vec2i base = ((p.x * inv_delta_x_).array() - 0.5f).cast<int>();
        Vec2 new_v = Vec2::Zero();
        Mat23 new_C = Mat23::Zero();
        Mat3 new_M = Mat3::Zero();
        if (!cut_mesh_->grid(base[1], base[0]).near_boundary) {
            Vec2 fx = p.x * inv_delta_x_ - base.cast<Real>();
            Vec2 weights[3] = {0.5 * (1.5 - fx.array()).square(),
                               0.75 - (fx.array() - 1.0).square(),
                               0.5 * (fx.array() - 0.5).square()};
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    auto offset = Vec2i(i, j);
                    auto weight = weights[i][0] * weights[j][1];
                    Vec2i grid_idx = base + offset;
                    auto& g = grid(grid_idx[0], grid_idx[1]);
                    Vec3 distance(0.0, g.vertex->position.x() - p.x.x(),
                                  g.vertex->position.y() - p.x.y());
                    distance *= static_cast<Real>(grid_size_);
                    distance.x() = 1.0;
                    distance *= delta_x_;
                    new_v += weight * g.v;
                    new_C += weight * (g.v * distance.transpose());
                    new_M += weight * distance * distance.transpose();
                }
            }
        } else {
            const auto& neighbor_nodes = enclosing_face->neighbor_nodes;
            int n_neighbor_nodes = static_cast<int>(neighbor_nodes.size());
            for (int i = 0; i < n_neighbor_nodes; ++i) {
                if (p.weights[i] <= 0) continue;
                const auto& g = nodes_[neighbor_nodes[i]];
                Vec3 distance(0.0, g.vertex->position.x() - p.x.x(),
                              g.vertex->position.y() - p.x.y());
                distance *= static_cast<Real>(grid_size_);
                distance.x() = 1.0;
                Real weight = p.weights[i] / p.weight_sum;
                distance *= delta_x_;
                auto v = g.v;
                if (g.vertex->convex)
                    v = g.vertex->project_convex_velocity(p.x, g.v);
                new_v += weight * v;
                new_C += weight * (v * distance.transpose());
                new_M += weight * distance * distance.transpose();
            }
        }
        p.v = new_v;
        p.M_inv = new_M.inverse();
        p.C = new_C * new_M.inverse();
        if (enclosing_face->neighbor_boundaries.empty()) {
            p.x += p.v * delta_t_;
            continue;
        }

        Vec2 old_x = p.x;
        for (int i = 0; i < 4; ++i) {
            Real min_t = kInf;
            Vec2 normal = Vec2::Zero();
            Real distance = kInf;
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
            if (min_t == kInf) {
                p.x += p.v * delta_t_;
                break;
            }
            min_t *= std::max(Real{0}, (distance - margin_) / distance);
            if (min_t > delta_t_) {
                p.x += p.v * delta_t_;
                break;
            }
            p.x += p.v * min_t;
            p.v *= (delta_t_ - min_t) / delta_t_;
            p.v = project(normal, p.v);
        }
        p.v = (p.x - old_x) / delta_t_;
    }
}
