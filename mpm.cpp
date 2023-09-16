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
const MPM::Vec2 kGravity(0.0, 1.0);

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

static auto svd(const MPM::Mat2& m) {
    Eigen::JacobiSVD<MPM::Mat2, Eigen::ComputeFullU | Eigen::ComputeFullV> svd(
        m);
    return std::tuple{svd.matrixU(), svd.singularValues(), svd.matrixV()};
}

MPM::MPM(const std::shared_ptr<CutMesh>& cut_mesh_) : cut_mesh_(cut_mesh_) {
    int n_vertices = static_cast<int>(cut_mesh_->vertices().size());
    grid_nodes_.resize(n_vertices);
    for (auto [i, v] = std::tuple(0, begin(cut_mesh_->vertices()));
         i < n_vertices; ++i, ++v) {
        grid_nodes_[i].vertex = v;
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
    for (GridNode& g : grid_nodes_) {
        g.m = 0.0;
        g.v.setZero();
    }
    // P2G
    for (Particle& p : particles_) {
        Vec2i base = ((p.x * kInvDeltaX).array() - 0.5).cast<int>();
        Vec2 fx = p.x * kInvDeltaX - base.cast<Real>();
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
        Real weight_sum = 0.0;
        Real weights[3][3];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                Vec2i offset(i, j);
                Vec2 distance = (offset.cast<Real>() - fx);
                weights[i][j] =
                    interpolate(distance.x()) * interpolate(distance.y());
                weight_sum += weights[i][j];
            }
        }
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                Vec3 distance(1.0, static_cast<Real>(i) - fx.x(),
                              static_cast<Real>(j) - fx.y());
                Real weight = weights[i][j] / weight_sum;
                distance *= kDeltaX;
                int x = base.x() + i;
                int y = base.y() + j;
                GridNode& g = grid_nodes_[y * kGridRowSize + x];
                g.v += weight * affine * distance;
                g.m += weight * kParticleMass;
            }
        }
    }
    // Grid update
    for (GridNode& g : grid_nodes_) {
        int id = g.vertex->id;
        int x = id % kGridRowSize;
        int y = id / kGridRowSize;
        if (g.m <= 0.0) continue;
        g.v /= g.m;
        g.v += kDeltaT * kGravity * 30.0;
        if (x < 3 && g.v[0] < 0.0) g.v[0] = 0.0;
        if (x > kGridSize - 3 && g.v[0] > 0.0) g.v[0] = 0.0;
        if (y < 3 && g.v[1] < 0.0) g.v[1] = 0.0;
        if (y > kGridSize - 3 && g.v[1] > 0.0) g.v[1] = 0.0;
    }
    // G2P
    for (Particle& p : particles_) {
        Vec2i base = ((p.x * kInvDeltaX).array() - 0.5).cast<int>();
        Vec2 fx = p.x * kInvDeltaX - base.cast<Real>();
        Real weight_sum = 0.0;
        Real weights[3][3];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                Vec2i offset(i, j);
                Vec2 distance = (offset.cast<Real>() - fx);
                weights[i][j] =
                    interpolate(distance.x()) * interpolate(distance.y());
                weight_sum += weights[i][j];
            }
        }
        Vec2 new_v = Vec2::Zero();
        Mat23 new_C = Mat23::Zero();
        Mat3 new_M = Mat3::Zero();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                Vec3 distance(1.0, static_cast<Real>(i) - fx.x(),
                              static_cast<Real>(j) - fx.y());
                Real weight = weights[i][j] / weight_sum;
                distance *= kDeltaX;
                int x = base.x() + i;
                int y = base.y() + j;
                const GridNode& g = grid_nodes_[y * kGridRowSize + x];
                new_v += weight * g.v;
                new_C += weight * (g.v * distance.transpose());
                new_M += weight * distance * distance.transpose();
            }
        }
        p.v = new_v;
        p.M_inv = new_M.inverse();
        p.C = new_C * new_M.inverse();
        p.x += kDeltaT * p.v;
    }
}
