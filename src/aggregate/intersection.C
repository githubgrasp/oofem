#include "intersection.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <complex>

#include "aggregateerror.h"
#include "ellipsoid.h"
#include "fibre.h"

namespace aggregate {

namespace {

// Tolerances for classifying eigenvalues from the Alfano-Greer test.
// The 4×4 problem is small so EigenSolver returns sub-µ accurate roots in
// well-conditioned cases; these tolerances guard the boundary cases (near-
// touching, near-degenerate axes).
constexpr double kImagTol = 1e-10;
constexpr double kReal0Tol = 1e-12;
constexpr double kDistinctTol = 1e-10;

bool isRealNegative(const std::complex<double> &z)
{
    const double mag = std::max(1.0, std::abs(z));
    return std::abs(z.imag()) <= kImagTol * mag && z.real() < -kReal0Tol * mag;
}

bool sameRealValue(const std::complex<double> &a, const std::complex<double> &b)
{
    const double mag = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a.real() - b.real()) <= kDistinctTol * mag;
}

} // namespace

bool ellipsoidsIntersect(const Ellipsoid &a, const Ellipsoid &b)
{
    const Eigen::Matrix4d A = a.giveQuadricForm();
    const Eigen::Matrix4d B = b.giveQuadricForm();
    const Eigen::Matrix4d M = A.inverse() * B;

    Eigen::EigenSolver<Eigen::Matrix4d> solver(M, /*computeEigenvectors*/ false);
    if ( solver.info() != Eigen::Success ) {
        error("ellipsoidsIntersect: EigenSolver failed to converge");
    }
    const Eigen::Vector4cd eigenvalues = solver.eigenvalues();

    int negativeRealCount = 0;
    std::complex<double> firstNegative;
    for ( int i = 0; i < 4; ++i ) {
        const auto &z = eigenvalues(i);
        if ( !isRealNegative(z) ) {
            continue;
        }
        if ( negativeRealCount == 0 ) {
            firstNegative = z;
            ++negativeRealCount;
        } else if ( !sameRealValue(z, firstNegative) ) {
            ++negativeRealCount;
        }
    }
    return negativeRealCount != 2;
}

bool ellipsoidIntersectsAxisPlane(const Ellipsoid &e, int axis, double position)
{
    if ( axis < 0 || axis > 2 ) {
        errorf("ellipsoidIntersectsAxisPlane: axis must be 0/1/2, got %d", axis);
    }
    const Eigen::Matrix3d R = e.giveRotation();
    const Eigen::Vector3d &s = e.giveSemiAxes();
    Eigen::Matrix3d S2 = Eigen::Matrix3d::Zero();
    S2(0, 0) = s(0) * s(0);
    S2(1, 1) = s(1) * s(1);
    S2(2, 2) = s(2) * s(2);
    const Eigen::Matrix3d Q = R * S2 * R.transpose();
    const double halfExtent = std::sqrt(Q(axis, axis));
    return std::abs(position - e.giveCentre()(axis)) <= halfExtent;
}

bool fibreIntersectsEllipsoid(const Fibre &f, const Ellipsoid &e)
{
    const Eigen::Matrix4d M = e.giveQuadricForm();
    Eigen::Vector4d origin;
    origin << f.giveEndpointA(), 1.0;
    Eigen::Vector4d direction;
    direction << ( f.giveEndpointB() - f.giveEndpointA() ), 0.0;

    const double a = direction.transpose() * M * direction;
    const double b = 2.0 * ( origin.transpose() * M * direction )( 0, 0 );
    const double c = origin.transpose() * M * origin;

    if ( a <= 0.0 ) {
        error("fibreIntersectsEllipsoid: degenerate quadratic (zero direction or non-positive-definite form)");
    }

    const double discr = b * b - 4.0 * a * c;
    if ( discr < 0.0 ) {
        // Line never crosses the boundary — segment is entirely inside iff endpoint A is inside.
        return c < 0.0;
    }
    const double sqrtDiscr = std::sqrt(discr);
    const double t1 = ( -b - sqrtDiscr ) / ( 2.0 * a );
    const double t2 = ( -b + sqrtDiscr ) / ( 2.0 * a );
    return t1 < 1.0 && t2 > 0.0;
}

} // namespace aggregate
