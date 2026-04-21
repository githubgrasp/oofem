#include "placer.h"

#include <Eigen/Geometry>
#include <cmath>
#include <memory>

#include "box.h"
#include "ellipsoid.h"
#include "fibre.h"
#include "intersection.h"
#include "random.h"

namespace aggregate {

namespace {

/// Cheap "definitely overlapping" pre-test using inscribed-radius spheres.
/// If the centres are closer than the sum of inscribed radii, the
/// ellipsoids' interiors must overlap regardless of orientation.
bool definitelyOverlap(const Ellipsoid &a, const Ellipsoid &b)
{
    const double dist = ( a.giveCentre() - b.giveCentre() ).norm();
    return dist < a.giveInscribedRadius() + b.giveInscribedRadius();
}

/// Cheap "possibly overlapping" pre-test using bounding-radius spheres.
/// If the centres are farther apart than the sum of bounding radii, the
/// ellipsoids cannot overlap.
bool possiblyOverlap(const Ellipsoid &a, const Ellipsoid &b)
{
    const double dist = ( a.giveCentre() - b.giveCentre() ).norm();
    return dist <= a.giveBoundingRadius() + b.giveBoundingRadius();
}

bool ellipsoidsOverlap(const Ellipsoid &a, const Ellipsoid &b)
{
    if ( definitelyOverlap(a, b) ) {
        return true;
    }
    if ( !possiblyOverlap(a, b) ) {
        return false;
    }
    return ellipsoidsIntersect(a, b);
}

/// Shoemake's uniform-on-SO(3) sampler (1992): three independent uniform
/// [0, 1) draws produce a unit quaternion uniformly distributed on the
/// 3-sphere, equivalent to a uniform random rotation. Returned as the
/// Z-Y-X Euler triple (phix, phiy, phiz) used by `Ellipsoid`.
Eigen::Vector3d sampleOrientation(std::mt19937 &rng)
{
    const double u1 = uniform(rng);
    const double u2 = uniform(rng);
    const double u3 = uniform(rng);
    const double sqrt1 = std::sqrt(1.0 - u1);
    const double sqrt2 = std::sqrt(u1);
    const Eigen::Quaterniond q(sqrt1 * std::cos(2.0 * M_PI * u2),
                               sqrt1 * std::sin(2.0 * M_PI * u2),
                               sqrt2 * std::sin(2.0 * M_PI * u3),
                               sqrt2 * std::cos(2.0 * M_PI * u3));
    const Eigen::Vector3d zyx = q.toRotationMatrix().eulerAngles(2, 1, 0);
    return Eigen::Vector3d(zyx(2), zyx(1), zyx(0));
}

} // namespace

std::vector<Eigen::Vector3d> ghostShifts(const Eigen::Vector3d &boxDims, int bitmask)
{
    std::vector<Eigen::Vector3d> shifts;
    const int nx = ( bitmask & 1 ) ? 2 : 1;
    const int ny = ( bitmask & 2 ) ? 2 : 1;
    const int nz = ( bitmask & 4 ) ? 2 : 1;
    for ( int kx = 0; kx < nx; ++kx ) {
        for ( int ky = 0; ky < ny; ++ky ) {
            for ( int kz = 0; kz < nz; ++kz ) {
                if ( kx == 0 && ky == 0 && kz == 0 ) {
                    continue;
                }
                shifts.emplace_back(kx * boxDims(0), ky * boxDims(1), kz * boxDims(2));
            }
        }
    }
    return shifts;
}

Placer::Placer(Box &boxRef, std::mt19937 &rngRef)
    : box(boxRef), rng(rngRef)
{}

bool Placer::detectBoundaries(Eigen::Vector3d &centre,
                              const Eigen::Vector3d &angles,
                              const Eigen::Vector3d &semiAxes,
                              int &bitmask) const
{
    bitmask = 0;
    const Eigen::Vector3d &boxDims = box.giveDimensions();
    const Eigen::Vector3i &periodic = box.givePeriodicityFlag();
    for ( int axis = 0; axis < 3; ++axis ) {
        Ellipsoid temp(0, centre, angles, semiAxes);
        if ( ellipsoidIntersectsAxisPlane(temp, axis, 0.0) ) {
            if ( periodic(axis) == 0 ) {
                return false;
            }
            bitmask |= ( 1 << axis );
        } else if ( ellipsoidIntersectsAxisPlane(temp, axis, boxDims(axis)) ) {
            if ( periodic(axis) == 0 ) {
                return false;
            }
            bitmask |= ( 1 << axis );
            centre(axis) -= boxDims(axis);
        }
    }
    return true;
}

bool Placer::overlapsAnyExisting(const Ellipsoid &candidate,
                                 const std::vector<Eigen::Vector3d> &shifts) const
{
    auto checkAgainst = [&]( const std::vector<std::unique_ptr<Inclusion>> &list ) {
        for ( const auto &incPtr : list ) {
            const auto *e = dynamic_cast<const Ellipsoid *>( incPtr.get() );
            if ( !e ) {
                continue;
            }
            if ( ellipsoidsOverlap(candidate, *e) ) {
                return true;
            }
            for ( const auto &shift : shifts ) {
                Ellipsoid copy(0,
                               candidate.giveCentre() + shift,
                               candidate.giveEulerAngles(),
                               candidate.giveSemiAxes());
                if ( ellipsoidsOverlap(copy, *e) ) {
                    return true;
                }
            }
        }
        return false;
    };

    if ( checkAgainst(box.giveRealInclusions()) ) {
        return true;
    }
    return checkAgainst(box.giveGhostInclusions());
}

bool Placer::placeOne(const Eigen::Vector3d &semiAxes)
{
    const Eigen::Vector3d &boxDims = box.giveDimensions();
    const int maxIter = box.giveMaximumIterations();

    for ( int iter = 0; iter < maxIter; ++iter ) {
        Eigen::Vector3d centre(uniform(rng, 0.0, boxDims(0)),
                               uniform(rng, 0.0, boxDims(1)),
                               uniform(rng, 0.0, boxDims(2)));
        Eigen::Vector3d angles = sampleOrientation(rng);

        int bitmask = 0;
        if ( !detectBoundaries(centre, angles, semiAxes, bitmask) ) {
            continue;
        }

        Ellipsoid candidate(0, centre, angles, semiAxes);
        const std::vector<Eigen::Vector3d> shifts = ghostShifts(boxDims, bitmask);

        if ( overlapsAnyExisting(candidate, shifts) ) {
            continue;
        }

        const int id = static_cast<int>( box.giveRealInclusions().size() ) + 1;
        box.addReal(std::make_unique<Ellipsoid>(id, centre, angles, semiAxes));
        for ( const auto &shift : shifts ) {
            box.addGhost(std::make_unique<Ellipsoid>(id, centre + shift, angles, semiAxes));
        }
        return true;
    }
    return false;
}

int Placer::placeBatch(const std::vector<Eigen::Vector3d> &sizes)
{
    int failed = 0;
    for ( const auto &s : sizes ) {
        if ( !placeOne(s) ) {
            ++failed;
        }
    }
    return failed;
}

bool Placer::detectFibreBoundaries(Eigen::Vector3d &centre,
                                   Eigen::Vector3d &endpointA,
                                   Eigen::Vector3d &endpointB,
                                   int &bitmask) const
{
    bitmask = 0;
    const Eigen::Vector3d &boxDims = box.giveDimensions();
    const Eigen::Vector3i &periodic = box.givePeriodicityFlag();
    for ( int axis = 0; axis < 3; ++axis ) {
        const double a = endpointA(axis);
        const double b = endpointB(axis);
        if ( a < 0.0 || b < 0.0 ) {
            if ( periodic(axis) == 0 ) {
                return false;
            }
            bitmask |= ( 1 << axis );
        } else if ( a > boxDims(axis) || b > boxDims(axis) ) {
            if ( periodic(axis) == 0 ) {
                return false;
            }
            bitmask |= ( 1 << axis );
            centre(axis) -= boxDims(axis);
            endpointA(axis) -= boxDims(axis);
            endpointB(axis) -= boxDims(axis);
        }
    }
    return true;
}

bool Placer::fibreOverlapsAnyEllipsoid(const Fibre &candidate,
                                       const std::vector<Eigen::Vector3d> &shifts) const
{
    auto fibreEllipsoidOverlap = [](const Fibre &f, const Ellipsoid &e) {
        const double dist = ( f.giveCentre() - e.giveCentre() ).norm();
        // Fibre centre inside the ellipsoid's inscribed sphere — must overlap.
        if ( dist < e.giveInscribedRadius() ) {
            return true;
        }
        // Fibre centre too far from ellipsoid for any reach to matter.
        if ( dist > f.giveHalfLength() + e.giveBoundingRadius() ) {
            return false;
        }
        return fibreIntersectsEllipsoid(f, e);
    };

    auto checkAgainst = [&]( const std::vector<std::unique_ptr<Inclusion>> &list ) {
        for ( const auto &incPtr : list ) {
            const auto *e = dynamic_cast<const Ellipsoid *>( incPtr.get() );
            if ( !e ) {
                continue;
            }
            if ( fibreEllipsoidOverlap(candidate, *e) ) {
                return true;
            }
            for ( const auto &shift : shifts ) {
                Fibre copy(0,
                           candidate.giveEndpointA() + shift,
                           candidate.giveEndpointB() + shift,
                           candidate.giveDiameter());
                if ( fibreEllipsoidOverlap(copy, *e) ) {
                    return true;
                }
            }
        }
        return false;
    };

    if ( checkAgainst(box.giveRealInclusions()) ) {
        return true;
    }
    return checkAgainst(box.giveGhostInclusions());
}

bool Placer::placeFibre(double length, double diameter)
{
    if ( length <= 0.0 || diameter <= 0.0 ) {
        return false;
    }
    const Eigen::Vector3d &boxDims = box.giveDimensions();
    const int maxIter = box.giveMaximumIterations();
    const double halfLength = 0.5 * length;

    for ( int iter = 0; iter < maxIter; ++iter ) {
        Eigen::Vector3d centre(uniform(rng, 0.0, boxDims(0)),
                               uniform(rng, 0.0, boxDims(1)),
                               uniform(rng, 0.0, boxDims(2)));
        const Eigen::Vector3d direction = uniformOnSphere(rng);
        Eigen::Vector3d endpointA = centre + halfLength * direction;
        Eigen::Vector3d endpointB = centre - halfLength * direction;

        int bitmask = 0;
        if ( !detectFibreBoundaries(centre, endpointA, endpointB, bitmask) ) {
            continue;
        }

        Fibre candidate(0, endpointA, endpointB, diameter);
        const std::vector<Eigen::Vector3d> shifts = ghostShifts(boxDims, bitmask);

        if ( fibreOverlapsAnyEllipsoid(candidate, shifts) ) {
            continue;
        }

        int existingFibres = 0;
        for ( const auto &inc : box.giveRealInclusions() ) {
            if ( dynamic_cast<const Fibre *>( inc.get() ) ) {
                ++existingFibres;
            }
        }
        const int id = existingFibres + 1;
        box.addReal(std::make_unique<Fibre>(id, endpointA, endpointB, diameter));
        for ( const auto &shift : shifts ) {
            box.addGhost(std::make_unique<Fibre>(id,
                                                 endpointA + shift,
                                                 endpointB + shift,
                                                 diameter));
        }
        return true;
    }
    return false;
}

int Placer::placeFibreBatch(int count, double length, double diameter)
{
    int failed = 0;
    for ( int i = 0; i < count; ++i ) {
        if ( !placeFibre(length, diameter) ) {
            ++failed;
        }
    }
    return failed;
}

} // namespace aggregate
