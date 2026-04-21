#pragma once

namespace aggregate {

/**
 * Diameter threshold separating aggregates placed before and after fibres.
 *
 * Replaces the symbolic Matlab solver from Vandewalle (`d_lim.m`) with a
 * companion-matrix root finder using Eigen. Solves
 *
 *   ( F·x³ + G·x² − A )² = ( A · √(x/dmax) )²
 *
 * with F = fibreFraction / (S·fibreLength), G = fibreFraction / S,
 * A = aggregateFraction, S = π·fibreDiameter² / 4. Substituting
 * u = √(x/dmax) yields two degree-6 polynomials whose smallest positive
 * real root determines x = u²·dmax.
 *
 * When fibreFraction == 0 the polynomial degenerates and the threshold is
 * simply dmax (i.e. group splitting is moot — every aggregate ends up in
 * the post-fibre group, but no fibres are placed in any case).
 */
double dlim(double fibreDiameter,
            double dmax,
            double fibreFraction,
            double fibreLength,
            double aggregateFraction);

} // namespace aggregate
