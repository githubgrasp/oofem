#include "dlim.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

#include "aggregateerror.h"

namespace aggregate {

namespace {

constexpr double kImagTol = 1e-10;
constexpr double kPositiveTol = 1e-12;

/// Smallest real-positive root of the polynomial whose coefficients are
/// stored as { c0, c1, ..., cn } with cn != 0.
double smallestPositiveRoot(const Eigen::VectorXd &coeffs)
{
    const int degree = static_cast<int>(coeffs.size()) - 1;
    if ( degree <= 0 || coeffs(degree) == 0.0 ) {
        error("dlim: degenerate polynomial passed to smallestPositiveRoot");
    }

    Eigen::MatrixXd companion = Eigen::MatrixXd::Zero(degree, degree);
    for ( int i = 1; i < degree; ++i ) {
        companion(i, i - 1) = 1.0;
    }
    for ( int i = 0; i < degree; ++i ) {
        companion(i, degree - 1) = -coeffs(i) / coeffs(degree);
    }

    Eigen::EigenSolver<Eigen::MatrixXd> solver(companion, /*computeEigenvectors*/ false);
    if ( solver.info() != Eigen::Success ) {
        error("dlim: eigenvalue solver failed on companion matrix");
    }

    double best = std::numeric_limits<double>::infinity();
    const auto eigs = solver.eigenvalues();
    for ( int i = 0; i < degree; ++i ) {
        const auto &z = eigs(i);
        const double mag = std::max(1.0, std::abs(z));
        if ( std::abs(z.imag()) > kImagTol * mag ) {
            continue;
        }
        if ( z.real() <= kPositiveTol * mag ) {
            continue;
        }
        if ( z.real() < best ) {
            best = z.real();
        }
    }
    return best;
}

} // namespace

double dlim(double fibreDiameter,
            double dmax,
            double fibreFraction,
            double fibreLength,
            double aggregateFraction)
{
    if ( fibreFraction == 0.0 ) {
        return dmax;
    }
    if ( fibreDiameter <= 0.0 || dmax <= 0.0 || fibreLength <= 0.0 || aggregateFraction <= 0.0 ) {
        errorf("dlim: non-positive input (df=%.3e dmax=%.3e L=%.3e A=%.3e)",
               fibreDiameter, dmax, fibreLength, aggregateFraction);
    }

    const double S = M_PI * fibreDiameter * fibreDiameter / 4.0;
    const double F = fibreFraction / ( S * fibreLength );
    const double G = fibreFraction / S;
    const double A = aggregateFraction;

    // Substitute u = sqrt(x/dmax). Both branches are degree-6 polynomials in u.
    // Coefficients are stored low-order first: {c0, c1, c2, c3, c4, c5, c6}.
    // P_+ : F·dmax³·u⁶ + G·dmax²·u⁴ − A·u − A
    // P_− : F·dmax³·u⁶ + G·dmax²·u⁴ + A·u − A
    Eigen::VectorXd plus(7);
    plus << -A, -A, 0.0, 0.0, G * dmax * dmax, 0.0, F * dmax * dmax * dmax;
    Eigen::VectorXd minus(7);
    minus << -A,  A, 0.0, 0.0, G * dmax * dmax, 0.0, F * dmax * dmax * dmax;

    const double uPlus = smallestPositiveRoot(plus);
    const double uMinus = smallestPositiveRoot(minus);
    const double u = std::min(uPlus, uMinus);

    if ( !std::isfinite(u) ) {
        return dmax;
    }
    return u * u * dmax;
}

} // namespace aggregate
