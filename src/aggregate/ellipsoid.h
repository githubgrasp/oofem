#pragma once

#include <Eigen/Dense>
#include <iosfwd>
#include <sstream>
#include <string>

#include "inclusion.h"

namespace aggregate {

/**
 * Ellipsoid with arbitrary orientation, used as the primary inclusion type.
 *
 * Parameterised by centre, three Euler angles (Rz·Ry·Rx convention, matching
 * the original Matlab implementation), and three semi-axes. The 4×4 Alfano-
 * Greer matrix form `X·A·X^T = 0` is computed on demand for intersection
 * tests.
 */
class Ellipsoid : public Inclusion
{
public:
    /// Construct a fully-specified ellipsoid. `centre` is in world
    /// coordinates, `eulerAngles` are the (phix, phiy, phiz) triple for
    /// the Rz·Ry·Rx convention, `semiAxes` are positive (sx, sy, sz).
    Ellipsoid(int n,
              const Eigen::Vector3d &centre,
              const Eigen::Vector3d &eulerAngles,
              const Eigen::Vector3d &semiAxes);

    /// Construct an uninitialised ellipsoid; populated via `initializeFromTokens`.
    explicit Ellipsoid(int n);

    /// Parse "centre 3 cx cy cz angles 3 phix phiy phiz radii 3 sx sy sz".
    void initializeFromTokens(std::istringstream &iss);

    std::string typeName() const override { return "ellipsoid"; }
    void writeTo(std::ostream &os) const override;

    /// World-coordinate centre of the ellipsoid.
    const Eigen::Vector3d &giveCentre() const { return centre; }

    /// Euler angles (phix, phiy, phiz) for the Rz·Ry·Rx rotation convention.
    const Eigen::Vector3d &giveEulerAngles() const { return eulerAngles; }

    /// Principal semi-axis lengths (sx, sy, sz) in the ellipsoid's local frame.
    const Eigen::Vector3d &giveSemiAxes() const { return semiAxes; }

    /// Largest semi-axis — used as a cheap separation pre-check.
    double giveBoundingRadius() const { return semiAxes.maxCoeff(); }

    /// Smallest semi-axis — used to bound the centre-distance below which overlap is unavoidable.
    double giveInscribedRadius() const { return semiAxes.minCoeff(); }

    /// 3×3 rotation matrix R = Rz·Ry·Rx applied to the principal axes.
    Eigen::Matrix3d giveRotation() const;

    /// 4×4 quadric form A in the Alfano-Greer convention; `X·A·X^T = 0` on the surface.
    Eigen::Matrix4d giveQuadricForm() const;

private:
    Eigen::Vector3d centre;
    Eigen::Vector3d eulerAngles;
    Eigen::Vector3d semiAxes;
};

} // namespace aggregate
