#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <random>

namespace aggregate {

/**
 * Portable uniform sampling from a `std::mt19937` engine.
 *
 * `std::uniform_real_distribution` is not guaranteed to produce identical
 * sequences across implementations; the engine alone is. By computing
 * the [0, 1) sample directly from the raw uint32 output, the placer's
 * behaviour stays bit-reproducible across compilers — required for the
 * diff-based ctest workflow.
 */
inline double uniform(std::mt19937 &rng)
{
    constexpr double divisor = static_cast<double>(std::mt19937::max()) + 1.0;
    return static_cast<double>(rng()) / divisor;
}

inline double uniform(std::mt19937 &rng, double lo, double hi)
{
    return lo + ( hi - lo ) * uniform(rng);
}

/// Standard normal sample via Box-Muller (only one of the two outputs used,
/// for simplicity — performance is not critical in this code path).
inline double normal(std::mt19937 &rng)
{
    double u1;
    do {
        u1 = uniform(rng);
    } while ( u1 == 0.0 );
    const double u2 = uniform(rng);
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
}

/// Uniform random direction on the unit sphere (Marsaglia: three independent
/// Gaussians, normalised). Used for fibre orientation, where the line is
/// undirected so the orientation lives on S² rather than SO(3).
inline Eigen::Vector3d uniformOnSphere(std::mt19937 &rng)
{
    Eigen::Vector3d v;
    double squaredNorm;
    do {
        v(0) = normal(rng);
        v(1) = normal(rng);
        v(2) = normal(rng);
        squaredNorm = v.squaredNorm();
    } while ( squaredNorm < 1e-30 );
    return v / std::sqrt(squaredNorm);
}

} // namespace aggregate
