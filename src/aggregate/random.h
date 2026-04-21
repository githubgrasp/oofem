#pragma once

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

} // namespace aggregate
