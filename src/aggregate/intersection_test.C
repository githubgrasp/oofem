#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "ellipsoid.h"
#include "fibre.h"
#include "intersection.h"

namespace {

int failures = 0;

#define EXPECT(cond, label)                                              \
    do {                                                                  \
        if ( !( cond ) ) {                                                \
            std::cerr << "FAIL " << __FILE__ << ':' << __LINE__           \
                      << "  " << ( label ) << '\n';                       \
            ++failures;                                                   \
        }                                                                 \
    } while ( 0 )

aggregate::Ellipsoid sphereAt(int id, double x, double y, double z, double r)
{
    return aggregate::Ellipsoid(id,
                                Eigen::Vector3d(x, y, z),
                                Eigen::Vector3d(0., 0., 0.),
                                Eigen::Vector3d(r, r, r));
}

aggregate::Ellipsoid axisAlignedEllipsoid(int id, double x, double y, double z,
                                          double sx, double sy, double sz)
{
    return aggregate::Ellipsoid(id,
                                Eigen::Vector3d(x, y, z),
                                Eigen::Vector3d(0., 0., 0.),
                                Eigen::Vector3d(sx, sy, sz));
}

void testEllipsoidEllipsoid()
{
    // Two unit spheres at distance 3 — clearly separated.
    EXPECT(!aggregate::ellipsoidsIntersect(sphereAt(1, 0, 0, 0, 1.0),
                                           sphereAt(2, 3, 0, 0, 1.0)),
           "spheres at distance 3 should not intersect");

    // Two unit spheres at distance 1 — clearly overlapping.
    EXPECT(aggregate::ellipsoidsIntersect(sphereAt(1, 0, 0, 0, 1.0),
                                          sphereAt(2, 1, 0, 0, 1.0)),
           "spheres at distance 1 should intersect");

    // Two unit spheres at distance 2 — touching, treated as intersecting.
    EXPECT(aggregate::ellipsoidsIntersect(sphereAt(1, 0, 0, 0, 1.0),
                                          sphereAt(2, 2, 0, 0, 1.0)),
           "spheres tangent at distance 2 should intersect (touching counts)");

    // Identical ellipsoids — must overlap.
    aggregate::Ellipsoid e(1,
                           Eigen::Vector3d(0.05, 0.05, 0.05),
                           Eigen::Vector3d(0.1, 0.2, 0.3),
                           Eigen::Vector3d(0.008, 0.006, 0.004));
    EXPECT(aggregate::ellipsoidsIntersect(e, e),
           "identical ellipsoids should overlap");

    // Alfano-Greer paper example (semi-axes form): primary (2,1,1) at origin,
    // secondary (3,2,4) at (7,0,0). Closest face-to-face distance along x is
    // 7 - 2 - 3 = 2 — separated.
    aggregate::Ellipsoid primary = axisAlignedEllipsoid(1, 0, 0, 0, 2., 1., 1.);
    aggregate::Ellipsoid secondary = axisAlignedEllipsoid(2, 7, 0, 0, 3., 2., 4.);
    EXPECT(!aggregate::ellipsoidsIntersect(primary, secondary),
           "Alfano-Greer scale-1 example should be separated");

    // Scale primary up to (4,2,2) — the right edge reaches x=4, secondary's
    // left edge is at x=4: tangent — counts as intersecting.
    aggregate::Ellipsoid primaryScaled = axisAlignedEllipsoid(1, 0, 0, 0, 4., 2., 2.);
    EXPECT(aggregate::ellipsoidsIntersect(primaryScaled, secondary),
           "Alfano-Greer scale-2 example should be touching (counts as overlap)");

    // Scale primary further to (5,2.5,2.5) — clearly overlapping.
    aggregate::Ellipsoid primaryBig = axisAlignedEllipsoid(1, 0, 0, 0, 5., 2.5, 2.5);
    EXPECT(aggregate::ellipsoidsIntersect(primaryBig, secondary),
           "Alfano-Greer overlap case should be detected");

    // Far-apart oriented ellipsoids — separated regardless of orientation.
    aggregate::Ellipsoid orientedA(1,
                                   Eigen::Vector3d(0., 0., 0.),
                                   Eigen::Vector3d(0.5, 0.7, 1.1),
                                   Eigen::Vector3d(1., 0.5, 0.3));
    aggregate::Ellipsoid orientedB(2,
                                   Eigen::Vector3d(10., 10., 10.),
                                   Eigen::Vector3d(0.2, -0.4, 0.9),
                                   Eigen::Vector3d(1., 0.5, 0.3));
    EXPECT(!aggregate::ellipsoidsIntersect(orientedA, orientedB),
           "far-apart oriented ellipsoids should not intersect");
}

void testEllipsoidPlane()
{
    aggregate::Ellipsoid e = axisAlignedEllipsoid(1, 0, 0, 0, 1., 2., 3.);

    EXPECT(aggregate::ellipsoidIntersectsAxisPlane(e, 0, 0.5),
           "x=0.5 should cut a (1,2,3) ellipsoid centred at origin");
    EXPECT(!aggregate::ellipsoidIntersectsAxisPlane(e, 0, 2.0),
           "x=2.0 should miss a (1,2,3) ellipsoid centred at origin");
    EXPECT(aggregate::ellipsoidIntersectsAxisPlane(e, 0, 1.0),
           "x=1.0 should be tangent to ellipsoid (counts as intersect)");

    EXPECT(aggregate::ellipsoidIntersectsAxisPlane(e, 1, 1.5),
           "y=1.5 should cut a (1,2,3) ellipsoid centred at origin");
    EXPECT(!aggregate::ellipsoidIntersectsAxisPlane(e, 1, 3.0),
           "y=3.0 should miss a (1,2,3) ellipsoid");

    EXPECT(aggregate::ellipsoidIntersectsAxisPlane(e, 2, 2.99),
           "z=2.99 should cut a (1,2,3) ellipsoid");
    EXPECT(!aggregate::ellipsoidIntersectsAxisPlane(e, 2, 3.5),
           "z=3.5 should miss a (1,2,3) ellipsoid");
}

void testFibreEllipsoid()
{
    aggregate::Ellipsoid sphere = sphereAt(1, 0, 0, 0, 1.0);

    aggregate::Fibre crossing(1,
                              Eigen::Vector3d(-2., 0., 0.),
                              Eigen::Vector3d(2., 0., 0.),
                              0.01);
    EXPECT(aggregate::fibreIntersectsEllipsoid(crossing, sphere),
           "fibre crossing unit sphere should intersect");

    aggregate::Fibre outside(1,
                             Eigen::Vector3d(-3., 0., 0.),
                             Eigen::Vector3d(-2.5, 0., 0.),
                             0.01);
    EXPECT(!aggregate::fibreIntersectsEllipsoid(outside, sphere),
           "fibre fully outside unit sphere should not intersect");

    aggregate::Fibre fullyInside(1,
                                 Eigen::Vector3d(-0.3, 0., 0.),
                                 Eigen::Vector3d(0.3, 0., 0.),
                                 0.01);
    EXPECT(aggregate::fibreIntersectsEllipsoid(fullyInside, sphere),
           "fibre fully inside unit sphere should intersect");

    aggregate::Fibre tangent(1,
                             Eigen::Vector3d(-1., 1., 0.),
                             Eigen::Vector3d(1., 1., 0.),
                             0.01);
    EXPECT(aggregate::fibreIntersectsEllipsoid(tangent, sphere),
           "fibre tangent to unit sphere should intersect (touching counts)");

    aggregate::Fibre missing(1,
                             Eigen::Vector3d(-1., 1.1, 0.),
                             Eigen::Vector3d(1., 1.1, 0.),
                             0.01);
    EXPECT(!aggregate::fibreIntersectsEllipsoid(missing, sphere),
           "fibre passing 0.1 above unit sphere should not intersect");
}

} // namespace

int main()
{
    testEllipsoidEllipsoid();
    testEllipsoidPlane();
    testFibreEllipsoid();

    if ( failures == 0 ) {
        std::cout << "intersection unit tests passed\n";
        return EXIT_SUCCESS;
    }
    std::cerr << failures << " intersection unit test(s) failed\n";
    return EXIT_FAILURE;
}
