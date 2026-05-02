#pragma once

#include <Eigen/Dense>
#include <random>
#include <vector>

namespace aggregate {

/**
 * Generates a sieve-graded list of ellipsoid semi-axes for the placer to
 * consume. Implements Vandewalle's adaptation of Fuller's grading curve
 * with Wriggers' interval-volume formula:
 *
 *   - sieves halve from dmax/2 down to dmin
 *   - per-sieve target volume from `wriggersIntervalVolume`
 *   - within each sieve, random aggregate sizes are sampled with the
 *     constraints from `placement_aggregates_periodic.m` (sx ≤ sy ≤ sz,
 *     aspect ratio sx/sz in [0.5, 1])
 *   - residual volume below the per-sieve budget carries over to the
 *     next (smaller) sieve
 *
 * The output is sorted by volume in descending order — the placer
 * processes the large aggregates first, optionally splitting at the
 * threshold returned by `dlim()` so the smallest sit at the end.
 */
/// Aggregate shape selector. Spherical aggregates collapse `sx = sy = sz`,
/// which downstream tools (notably the generator) can handle without an
/// arbitrary-orientation ellipsoid implementation.
enum class GradingShape
{
    Ellipsoid,
    Sphere
};

class GradingCurve
{
public:
    /// Construct a grading curve targeting `aggregateFraction` of `rveVolume`,
    /// with sieve diameters bounded below by `minSieve` and above by `maxSieve`.
    /// `maxSieve` must be at least `2·minSieve` (one full sieve octave).
    /// `shape` selects between general anisotropic ellipsoids (default) and
    /// degenerate spheres (`sx = sy = sz = r`).
    GradingCurve(double minSieve,
                 double maxSieve,
                 double aggregateFraction,
                 double rveVolume,
                 GradingShape shape = GradingShape::Ellipsoid);

    /// Sample sizes (sx, sy, sz) using the supplied engine, sorted by volume desc.
    std::vector<Eigen::Vector3d> generate(std::mt19937 &rng) const;

    /// 2D variant: sample disk radii sorted by area desc. Interprets the
    /// constructor's `rveVolume` argument as RVE area; same Fuller-curve
    /// sieve logic, with disk area `π·r²` substituted for ellipsoid volume.
    std::vector<double> generateDisks(std::mt19937 &rng) const;

    /// Volume associated with the sieve interval (lower, upper) — Wriggers' formula.
    double wriggersIntervalVolume(double lower, double upper) const;

private:
    double dmin;
    double dmax;
    double volumeFraction;
    double rveVolume;
    GradingShape shape;
};

} // namespace aggregate
