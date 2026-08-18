#include "stats.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <vector>

#include <Eigen/Dense>

#include "box.h"
#include "disk.h"
#include "ellipsoid.h"
#include "fibre.h"
#include "inclusion.h"

namespace aggregate {

namespace {

constexpr int kSamplesPerSubcube = 2048;
constexpr std::array<int, 4> kScales = { 2, 4, 8, 16 };

// ---- analytic volumes / areas ---------------------------------------------

double ellipsoidVolume(const Ellipsoid &e)
{
    const auto &s = e.giveSemiAxes();
    return ( 4.0 / 3.0 ) * M_PI * s(0) * s(1) * s(2);
}

double fibreVolume(const Fibre &f)
{
    const double r = 0.5 * f.giveDiameter();
    const double L = ( f.giveEndpointB() - f.giveEndpointA() ).norm();
    return M_PI * r * r * L;
}

double diskArea(const Disk &d)
{
    return M_PI * d.giveRadius() * d.giveRadius();
}

// ---- point-in-shape tests --------------------------------------------------

bool ellipsoidContains(const Ellipsoid &e, const Eigen::Vector3d &p)
{
    const Eigen::Vector3d local = e.giveRotation().transpose() * ( p - e.giveCentre() );
    const auto &s = e.giveSemiAxes();
    const double q = ( local(0) / s(0) ) * ( local(0) / s(0) )
                   + ( local(1) / s(1) ) * ( local(1) / s(1) )
                   + ( local(2) / s(2) ) * ( local(2) / s(2) );
    return q <= 1.0;
}

bool fibreContains(const Fibre &f, const Eigen::Vector3d &p)
{
    const Eigen::Vector3d &A = f.giveEndpointA();
    const Eigen::Vector3d d = f.giveEndpointB() - A;
    const Eigen::Vector3d ap = p - A;
    const double t = ap.dot(d) / d.squaredNorm();
    if ( t < 0.0 || t > 1.0 ) {
        return false;
    }
    const double r = 0.5 * f.giveDiameter();
    return ( ap - t * d ).squaredNorm() <= r * r;
}

bool diskContains(const Disk &d, const Eigen::Vector2d &p)
{
    const double r = d.giveRadius();
    return ( p - d.giveCentre() ).squaredNorm() <= r * r;
}

// ---- AABB pre-filter -------------------------------------------------------

struct Aabb3d { Eigen::Vector3d lo, hi; };
struct Aabb2d { Eigen::Vector2d lo, hi; };

bool aabbContains(const Aabb3d &b, const Eigen::Vector3d &p)
{
    return ( p.array() >= b.lo.array() ).all() && ( p.array() <= b.hi.array() ).all();
}

bool aabbContains(const Aabb2d &b, const Eigen::Vector2d &p)
{
    return ( p.array() >= b.lo.array() ).all() && ( p.array() <= b.hi.array() ).all();
}

Aabb3d ellipsoidAabb(const Ellipsoid &e)
{
    const Eigen::Matrix3d R = e.giveRotation();
    const auto &s = e.giveSemiAxes();
    Eigen::Vector3d ext;
    for ( int i = 0; i < 3; ++i ) {
        ext(i) = std::sqrt(
            ( R(i, 0) * s(0) ) * ( R(i, 0) * s(0) ) +
            ( R(i, 1) * s(1) ) * ( R(i, 1) * s(1) ) +
            ( R(i, 2) * s(2) ) * ( R(i, 2) * s(2) ) );
    }
    return { e.giveCentre() - ext, e.giveCentre() + ext };
}

Aabb3d fibreAabb(const Fibre &f)
{
    const auto &A = f.giveEndpointA();
    const auto &B = f.giveEndpointB();
    const Eigen::Vector3d r = Eigen::Vector3d::Constant(0.5 * f.giveDiameter());
    return { A.cwiseMin(B) - r, A.cwiseMax(B) + r };
}

Aabb2d diskAabb(const Disk &d)
{
    const Eigen::Vector2d r = Eigen::Vector2d::Constant(d.giveRadius());
    return { d.giveCentre() - r, d.giveCentre() + r };
}

// ---- entry types for combined real + ghost iteration -----------------------

template<class Shape, class Aabb>
struct ShapeEntry {
    const Shape *shape;
    Aabb bbox;
};

using EllipsoidEntry = ShapeEntry<Ellipsoid, Aabb3d>;
using FibreEntry     = ShapeEntry<Fibre,     Aabb3d>;
using DiskEntry      = ShapeEntry<Disk,      Aabb2d>;

template<class Entry, class Point, class Predicate>
bool anyContains(const std::vector<Entry> &entries, const Point &p, Predicate inside)
{
    for ( const auto &en : entries ) {
        if ( !aabbContains(en.bbox, p) ) {
            continue;
        }
        if ( inside(*en.shape, p) ) {
            return true;
        }
    }
    return false;
}

void collectEllipsoidsAndFibres(
    const std::vector<std::unique_ptr<Inclusion>> &list,
    std::vector<EllipsoidEntry> &ellipsoids,
    std::vector<FibreEntry> &fibres)
{
    for ( const auto &incPtr : list ) {
        if ( const auto *e = dynamic_cast<const Ellipsoid *>( incPtr.get() ) ) {
            ellipsoids.push_back({ e, ellipsoidAabb(*e) });
        } else if ( const auto *f = dynamic_cast<const Fibre *>( incPtr.get() ) ) {
            fibres.push_back({ f, fibreAabb(*f) });
        }
    }
}

void collectDisks(
    const std::vector<std::unique_ptr<Inclusion>> &list,
    std::vector<DiskEntry> &disks)
{
    for ( const auto &incPtr : list ) {
        if ( const auto *d = dynamic_cast<const Disk *>( incPtr.get() ) ) {
            disks.push_back({ d, diskAabb(*d) });
        }
    }
}

// ---- statistics ------------------------------------------------------------

struct ScaleRow {
    int N;
    long subcubes;
    double mean;
    double stdSpatial;   // MC-noise-corrected
    double stdPoisson;
};

double poissonStd(double sumVi2, double cellMeasure, double subMeasure)
{
    if ( sumVi2 <= 0.0 || cellMeasure <= 0.0 || subMeasure <= 0.0 ) {
        return 0.0;
    }
    const double fSub = subMeasure / cellMeasure;
    return std::sqrt( ( 1.0 - fSub ) * sumVi2 / ( cellMeasure * subMeasure ) );
}

// 3D Monte-Carlo: returns mean and bias-corrected std across N^3 subcubes,
// using K stratified samples per subcube. `density` is queried once per
// sample and returns 1 if the point falls inside any inclusion of the type.
template<class Density>
std::pair<double, double> mcScaleStats3d(int N, const Eigen::Vector3d &cellDims,
                                         std::mt19937 &rng, Density density)
{
    const long subcubes = static_cast<long>( N ) * N * N;
    const Eigen::Vector3d step( cellDims(0) / N, cellDims(1) / N, cellDims(2) / N );
    std::uniform_real_distribution<double> u01(0.0, 1.0);

    std::vector<double> phi;
    phi.reserve(subcubes);

    for ( int ix = 0; ix < N; ++ix ) {
        for ( int iy = 0; iy < N; ++iy ) {
            for ( int iz = 0; iz < N; ++iz ) {
                const Eigen::Vector3d lo( ix * step(0), iy * step(1), iz * step(2) );
                int hits = 0;
                for ( int k = 0; k < kSamplesPerSubcube; ++k ) {
                    const Eigen::Vector3d p( lo(0) + u01(rng) * step(0),
                                             lo(1) + u01(rng) * step(1),
                                             lo(2) + u01(rng) * step(2) );
                    if ( density(p) ) {
                        ++hits;
                    }
                }
                phi.push_back(static_cast<double>(hits) / kSamplesPerSubcube);
            }
        }
    }

    double sum = 0.0;
    for ( double x : phi ) { sum += x; }
    const double mean = sum / phi.size();

    double sumSq = 0.0;
    for ( double x : phi ) { sumSq += ( x - mean ) * ( x - mean ); }
    const double varObs = sumSq / phi.size();

    // MC noise variance for binomial sampling at K per subcube, with success
    // rate ≈ mean. Subtract to isolate the spatial component; clamp at zero.
    const double varNoise = mean * ( 1.0 - mean ) / kSamplesPerSubcube;
    const double varSpatial = std::max(0.0, varObs - varNoise);
    return { mean, std::sqrt(varSpatial) };
}

template<class Density>
std::pair<double, double> mcScaleStats2d(int N, const Eigen::Vector3d &cellDims,
                                         std::mt19937 &rng, Density density)
{
    const long subcubes = static_cast<long>( N ) * N;
    const Eigen::Vector2d step( cellDims(0) / N, cellDims(1) / N );
    std::uniform_real_distribution<double> u01(0.0, 1.0);

    std::vector<double> phi;
    phi.reserve(subcubes);

    for ( int ix = 0; ix < N; ++ix ) {
        for ( int iy = 0; iy < N; ++iy ) {
            const Eigen::Vector2d lo( ix * step(0), iy * step(1) );
            int hits = 0;
            for ( int k = 0; k < kSamplesPerSubcube; ++k ) {
                const Eigen::Vector2d p( lo(0) + u01(rng) * step(0),
                                         lo(1) + u01(rng) * step(1) );
                if ( density(p) ) {
                    ++hits;
                }
            }
            phi.push_back(static_cast<double>(hits) / kSamplesPerSubcube);
        }
    }

    double sum = 0.0;
    for ( double x : phi ) { sum += x; }
    const double mean = sum / phi.size();

    double sumSq = 0.0;
    for ( double x : phi ) { sumSq += ( x - mean ) * ( x - mean ); }
    const double varObs = sumSq / phi.size();

    const double varNoise = mean * ( 1.0 - mean ) / kSamplesPerSubcube;
    const double varSpatial = std::max(0.0, varObs - varNoise);
    return { mean, std::sqrt(varSpatial) };
}

void printRow(std::ostream &os, const ScaleRow &r)
{
    const double ratio = ( r.stdPoisson > 0.0 ) ? r.stdSpatial / r.stdPoisson : 0.0;
    os << "   " << std::setw(3) << r.N << "^3"
       << std::setw(10) << r.subcubes
       << std::setw(14) << std::fixed << std::setprecision(5) << r.stdSpatial
       << std::setw(14) << std::setprecision(3) << ratio
       << '\n';
}

void printRow2d(std::ostream &os, const ScaleRow &r)
{
    const double ratio = ( r.stdPoisson > 0.0 ) ? r.stdSpatial / r.stdPoisson : 0.0;
    os << "   " << std::setw(3) << r.N << "^2"
       << std::setw(10) << r.subcubes
       << std::setw(14) << std::fixed << std::setprecision(5) << r.stdSpatial
       << std::setw(14) << std::setprecision(3) << ratio
       << '\n';
}

void runStats3d(const Box &box, std::mt19937 &rng, std::ostream &os)
{
    std::vector<EllipsoidEntry> ellipsoids;
    std::vector<FibreEntry> fibres;
    collectEllipsoidsAndFibres(box.giveRealInclusions(), ellipsoids, fibres);

    // Σ V_i² over REAL inclusions only — used in the Poisson reference.
    double sumVeSq = 0.0;
    double sumVf = 0.0;
    double sumVfSq = 0.0;
    int nEllipsoidsReal = 0;
    int nFibresReal = 0;
    double sumVe = 0.0;

    for ( const auto &en : ellipsoids ) {
        const double v = ellipsoidVolume(*en.shape);
        sumVe   += v;
        sumVeSq += v * v;
        ++nEllipsoidsReal;
    }
    for ( const auto &en : fibres ) {
        const double v = fibreVolume(*en.shape);
        sumVf   += v;
        sumVfSq += v * v;
        ++nFibresReal;
    }

    // Now extend the entry lists with ghosts (for MC point-in-shape tests).
    collectEllipsoidsAndFibres(box.giveGhostInclusions(), ellipsoids, fibres);

    const Eigen::Vector3d cellDims = box.giveDimensions();
    const double cellVol = box.giveMeasure();
    auto subVolAt = [&]( int N ) { return cellVol / ( long(N) * N * N ); };

    const auto &grading = box.giveGradingParameters();
    const auto &fibreParams = box.giveFibreParameters();
    const double targetVe = grading.present ? grading.aggregateFraction : 0.0;
    const double targetVf = fibreParams.present ? fibreParams.fibreFraction : 0.0;
    const double actualVe = sumVe / cellVol;
    const double actualVf = sumVf / cellVol;

    os << "\n=== Placement statistics ===\n";
    os << "Achieved volume fractions:\n";
    if ( grading.present ) {
        os << "  ellipsoids: " << std::setw(4) << nEllipsoidsReal
           << " placed, target Vf = " << std::fixed << std::setprecision(4) << targetVe
           << ", actual Vf = " << actualVe;
        if ( targetVe > 0.0 ) {
            os << " (" << std::setprecision(1) << ( 100.0 * actualVe / targetVe ) << "%)";
        }
        os << '\n';
    }
    if ( fibreParams.present ) {
        os << "  fibres:     " << std::setw(4) << nFibresReal
           << " placed, target Vf = " << std::fixed << std::setprecision(4) << targetVf
           << ", actual Vf = " << actualVf;
        if ( targetVf > 0.0 ) {
            os << " (" << std::setprecision(1) << ( 100.0 * actualVf / targetVf ) << "%)";
        }
        os << '\n';
    }

    if ( ellipsoids.empty() && fibres.empty() ) {
        os << "(no inclusions placed; skipping uniformity stats)\n";
        return;
    }

    os << "\nSpatial uniformity (Monte-Carlo, K=" << kSamplesPerSubcube
       << " points/subcube; noise-corrected std):\n"
       << "   scale  subcubes           std    std/std_Poisson\n"
       << "   -------------------------------------------------\n";

    auto ellipsoidDensity = [&]( const Eigen::Vector3d &p ) {
        return anyContains(ellipsoids, p, ellipsoidContains);
    };
    auto fibreDensity = [&]( const Eigen::Vector3d &p ) {
        return anyContains(fibres, p, fibreContains);
    };

    if ( grading.present && !ellipsoids.empty() ) {
        os << "  ellipsoids (Vf = " << std::fixed << std::setprecision(5) << actualVe << ")\n";
        for ( int N : kScales ) {
            const auto [mean, stdSpatial] = mcScaleStats3d(N, cellDims, rng, ellipsoidDensity);
            ScaleRow r{ N, long(N) * N * N, mean, stdSpatial,
                        poissonStd(sumVeSq, cellVol, subVolAt(N)) };
            printRow(os, r);
        }
    }
    if ( fibreParams.present && !fibres.empty() ) {
        os << "  fibres (Vf = " << std::fixed << std::setprecision(5) << actualVf << ")\n";
        for ( int N : kScales ) {
            const auto [mean, stdSpatial] = mcScaleStats3d(N, cellDims, rng, fibreDensity);
            ScaleRow r{ N, long(N) * N * N, mean, stdSpatial,
                        poissonStd(sumVfSq, cellVol, subVolAt(N)) };
            printRow(os, r);
        }
    }
    os << "  (std/std_Poisson < 1 = more uniform than random; > 1 = clustered)\n";
}

void runStats2d(const Box &box, std::mt19937 &rng, std::ostream &os)
{
    std::vector<DiskEntry> disks;
    collectDisks(box.giveRealInclusions(), disks);

    double sumAd = 0.0;
    double sumAdSq = 0.0;
    int nDisksReal = 0;
    for ( const auto &en : disks ) {
        const double a = diskArea(*en.shape);
        sumAd   += a;
        sumAdSq += a * a;
        ++nDisksReal;
    }
    collectDisks(box.giveGhostInclusions(), disks);

    const Eigen::Vector3d cellDims = box.giveDimensions();
    const double cellArea = box.giveMeasure();
    const auto &grading = box.giveGradingParameters();
    const double targetAf = grading.present ? grading.aggregateFraction : 0.0;
    const double actualAf = sumAd / cellArea;

    os << "\n=== Placement statistics ===\n";
    os << "Achieved area fractions:\n";
    if ( grading.present ) {
        os << "  disks: " << std::setw(4) << nDisksReal
           << " placed, target Af = " << std::fixed << std::setprecision(4) << targetAf
           << ", actual Af = " << actualAf;
        if ( targetAf > 0.0 ) {
            os << " (" << std::setprecision(1) << ( 100.0 * actualAf / targetAf ) << "%)";
        }
        os << '\n';
    }

    if ( disks.empty() ) {
        os << "(no inclusions placed; skipping uniformity stats)\n";
        return;
    }

    os << "\nSpatial uniformity (Monte-Carlo, K=" << kSamplesPerSubcube
       << " points/subcell; noise-corrected std):\n"
       << "   scale  subcells           std    std/std_Poisson\n"
       << "   -------------------------------------------------\n";

    auto diskDensity = [&]( const Eigen::Vector2d &p ) {
        return anyContains(disks, p, diskContains);
    };

    os << "  disks (Af = " << std::fixed << std::setprecision(5) << actualAf << ")\n";
    for ( int N : kScales ) {
        const double subArea = cellArea / ( long(N) * N );
        const auto [mean, stdSpatial] = mcScaleStats2d(N, cellDims, rng, diskDensity);
        ScaleRow r{ N, long(N) * N, mean, stdSpatial,
                    poissonStd(sumAdSq, cellArea, subArea) };
        printRow2d(os, r);
    }
    os << "  (std/std_Poisson < 1 = more uniform than random; > 1 = clustered)\n";
}

} // namespace

void printStats(const Box &box, std::mt19937 &rng, std::ostream &os)
{
    if ( box.giveDim() == 2 ) {
        runStats2d(box, rng, os);
    } else {
        runStats3d(box, rng, os);
    }
}

} // namespace aggregate
