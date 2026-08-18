#pragma once

namespace aggregate {

class Ellipsoid;
class Fibre;

/**
 * Pairwise overlap predicates for the placer.
 *
 * The ellipsoid-ellipsoid test uses the Alfano-Greer (2003) characterisation:
 * compute the four eigenvalues of `A^{-1}·B`; the ellipsoids are separated
 * iff exactly two of them are real, negative, and distinct.
 *
 * Touching configurations (repeated eigenvalues, on-surface points) are
 * conservatively reported as overlapping — appropriate for a sequential
 * placer that rejects any candidate not strictly separated from existing
 * inclusions.
 */

/// True iff the two ellipsoids share any volume (touching counts as overlap).
bool ellipsoidsIntersect(const Ellipsoid &a, const Ellipsoid &b);

/// True iff the ellipsoid touches or crosses the axis-aligned plane.
/// `axis` is 0/1/2 for x/y/z; `position` is the constant coordinate.
bool ellipsoidIntersectsAxisPlane(const Ellipsoid &e, int axis, double position);

/// True iff the line segment AB of the fibre passes through the ellipsoid's
/// interior or surface. The fibre's diameter is not considered here — it is
/// the placer's responsibility to apply any safety margin via inflation.
bool fibreIntersectsEllipsoid(const Fibre &f, const Ellipsoid &e);

} // namespace aggregate
