#pragma once

#include <iosfwd>
#include <random>

namespace aggregate {

class Box;

/// Print placement statistics for the inclusions stored on `box`:
///   - Per-type achieved volume (3D) / area (2D) fractions vs. target.
///   - A variance-vs-scale curve at subcube grids N = 2, 4, 8, 16 (3D)
///     or N = 2, 4, 8, 16 (2D), with the per-subcube local density
///     estimated by stratified Monte-Carlo and reported as
///     `mean`, `std` (bias-corrected for MC noise), and `std/std_Poisson`.
///
/// The Poisson reference assumes inclusion centres are placed independently
/// uniformly in the cell; it gives the std a packing of the same volume
/// fraction would have under a non-interacting random process. A ratio
/// below 1 means the packing is more uniform than random (typical of
/// hard-core packings); above 1 means clustered.
///
/// `rng` is consumed for the MC sampling — pass an independent seed from
/// the placer's RNG so the stats output is reproducible without affecting
/// placement.
void printStats(const Box &box, std::mt19937 &rng, std::ostream &os);

} // namespace aggregate
