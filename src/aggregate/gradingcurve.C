#include "gradingcurve.h"

#include <algorithm>
#include <cmath>

#include "aggregateerror.h"
#include "random.h"

namespace aggregate {

namespace {

double ellipsoidVolume(const Eigen::Vector3d &s)
{
    return ( 4.0 / 3.0 ) * M_PI * s(0) * s(1) * s(2);
}

/// Sample one ordered triple (sx, sy, sz) with sx ≤ sy ≤ sz subject to the
/// sieve constraints for the interval [lower, upper]. Mirrors the inner
/// sampling logic of `aggregate_generation_periodic.m`.
Eigen::Vector3d sampleEllipsoid(double lower, double upper, std::mt19937 &rng)
{
    constexpr double rLower = 0.5;
    constexpr double rUpper = 1.0;
    const double r = uniform(rng, rLower, rUpper);
    const double sxLow = 0.5 * std::sqrt(2.0 / ( r * r + 1.0 )) * r * lower;
    const double sxHigh = 0.5 * upper;
    const double sx = uniform(rng, sxLow, sxHigh);
    const double sz = sx / r;

    const double sieveBound = std::sqrt(std::max(0.0, 0.5 * upper * upper - sx * sx));
    const double syHigh = std::min(sz, sieveBound);
    double syLow;
    if ( sx >= lower / std::sqrt(2.0) ) {
        syLow = sx;
    } else {
        const double minSieveBound = std::sqrt(std::max(0.0, 0.5 * lower * lower - sx * sx));
        syLow = std::max(sx, minSieveBound);
    }
    if ( syHigh < syLow ) {
        // Numerical edge case: sieve constraint impossible for this (r, sx).
        // Collapse to syLow so the sample is at least well-defined; the
        // outer loop will simply re-roll on the next iteration if the
        // resulting volume exceeds the budget.
        return Eigen::Vector3d(sx, syLow, sz);
    }
    const double sy = uniform(rng, syLow, syHigh);
    return Eigen::Vector3d(sx, sy, sz);
}

/// Sample a sphere radius uniformly in the half-sieve range [lower/2, upper/2]
/// and return it as a degenerate (sx = sy = sz = r) ellipsoid triple.
Eigen::Vector3d sampleSphere(double lower, double upper, std::mt19937 &rng)
{
    const double r = uniform(rng, 0.5 * lower, 0.5 * upper);
    return Eigen::Vector3d(r, r, r);
}

} // namespace

GradingCurve::GradingCurve(double minSieve,
                           double maxSieve,
                           double aggregateFraction,
                           double rveVol,
                           GradingShape shapeChoice)
    : dmin(minSieve),
      dmax(maxSieve),
      volumeFraction(aggregateFraction),
      rveVolume(rveVol),
      shape(shapeChoice)
{
    if ( minSieve <= 0.0 || maxSieve <= 0.0 || rveVol <= 0.0 ) {
        errorf("GradingCurve: non-positive input (dmin=%.3e dmax=%.3e V=%.3e)",
               minSieve, maxSieve, rveVol);
    }
    if ( maxSieve < 2.0 * minSieve ) {
        errorf("GradingCurve: dmax (%.3e) must be at least 2·dmin (%.3e)",
               maxSieve, minSieve);
    }
    if ( aggregateFraction < 0.0 || aggregateFraction > 1.0 ) {
        errorf("GradingCurve: aggregateFraction out of [0,1] (%.3e)", aggregateFraction);
    }
}

double GradingCurve::wriggersIntervalVolume(double lower, double upper) const
{
    const double pUpper = std::sqrt(upper / dmax);
    const double pLower = std::sqrt(lower / dmax);
    const double pmax = 1.0;
    const double pmin = std::sqrt(dmin / dmax);
    const double fraction = ( pUpper - pLower ) / ( pmax - pmin );
    return fraction * ( 1.0 - pmin ) * volumeFraction * rveVolume;
}

std::vector<Eigen::Vector3d> GradingCurve::generate(std::mt19937 &rng) const
{
    std::vector<Eigen::Vector3d> sizes;
    if ( volumeFraction == 0.0 ) {
        return sizes;
    }

    double carryover = 0.0;
    double n = 0.5 * dmax;
    while ( n >= dmin ) {
        const double m = 2.0 * n;
        const double budget = wriggersIntervalVolume(n, m) + carryover;
        double placed = 0.0;
        while ( true ) {
            Eigen::Vector3d s = ( shape == GradingShape::Sphere )
                                ? sampleSphere(n, m, rng)
                                : sampleEllipsoid(n, m, rng);
            const double v = ellipsoidVolume(s);
            if ( placed + v > budget ) {
                carryover = budget - placed;
                break;
            }
            sizes.push_back(s);
            placed += v;
        }
        n *= 0.5;
    }

    std::sort(sizes.begin(), sizes.end(),
              [](const Eigen::Vector3d &a, const Eigen::Vector3d &b) {
                  return ellipsoidVolume(a) > ellipsoidVolume(b);
              });
    return sizes;
}

} // namespace aggregate
