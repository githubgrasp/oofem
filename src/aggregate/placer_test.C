#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>

#include "box.h"
#include "ellipsoid.h"
#include "fibre.h"
#include "intersection.h"
#include "placer.h"

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

/// Build a Box for testing without going through the file parser.
class TestBox : public aggregate::Box {
public:
    TestBox(double lx, double ly, double lz, int periodicX, int periodicY, int periodicZ,
            int maxIter)
    {
        std::ofstream tmp("/tmp/aggregate_placer_test_box.in");
        tmp << "#@output /tmp/aggregate_placer_test_packing.dat\n"
            << "#@box 3 " << lx << ' ' << ly << ' ' << lz << '\n'
            << "#@periodicity 3 " << periodicX << ' ' << periodicY << ' ' << periodicZ << '\n'
            << "#@maxiter " << maxIter << '\n';
        tmp.close();
        readControlRecords("/tmp/aggregate_placer_test_box.in");
    }
};

void testGhostShifts()
{
    Eigen::Vector3d boxDims(1.0, 2.0, 3.0);

    // No bits set → no ghosts.
    EXPECT(aggregate::ghostShifts(boxDims, 0).empty(),
           "bitmask 0 must produce no ghost shifts");

    // One bit → 1 ghost.
    EXPECT(aggregate::ghostShifts(boxDims, 1).size() == 1,
           "single-axis crossing must produce 1 ghost");
    EXPECT(aggregate::ghostShifts(boxDims, 2).size() == 1,
           "single-axis crossing must produce 1 ghost (y)");
    EXPECT(aggregate::ghostShifts(boxDims, 4).size() == 1,
           "single-axis crossing must produce 1 ghost (z)");

    // Two bits → 3 ghosts.
    EXPECT(aggregate::ghostShifts(boxDims, 3).size() == 3,
           "two-axis crossing must produce 3 ghosts");
    EXPECT(aggregate::ghostShifts(boxDims, 5).size() == 3,
           "two-axis crossing must produce 3 ghosts (xz)");

    // All three bits → 7 ghosts.
    EXPECT(aggregate::ghostShifts(boxDims, 7).size() == 7,
           "three-axis crossing must produce 7 ghosts");
}

void testPlaceSparse()
{
    TestBox box(1.0, 1.0, 1.0, /*periodic=*/1, 1, 1, /*maxIter=*/1000);
    std::mt19937 rng(7);
    aggregate::Placer placer(box, rng);

    // Place 5 small spheres — easy to fit.
    int placed = 0;
    for ( int i = 0; i < 5; ++i ) {
        if ( placer.placeOne(Eigen::Vector3d(0.05, 0.05, 0.05)) ) {
            ++placed;
        }
    }
    EXPECT(placed == 5, "5 small spheres should fit in unit box");
    EXPECT(box.giveRealInclusions().size() == 5, "Box should hold 5 real inclusions");
}

void testPlaceTooTight()
{
    // Tiny box, big aggregate → cannot fit.
    TestBox box(0.1, 0.1, 0.1, 0, 0, 0, /*maxIter=*/200);
    std::mt19937 rng(7);
    aggregate::Placer placer(box, rng);

    bool result = placer.placeOne(Eigen::Vector3d(0.2, 0.2, 0.2));
    EXPECT(!result, "oversize aggregate in non-periodic small box must fail to place");
    EXPECT(box.giveRealInclusions().empty(), "no inclusion should be added on failure");
}

void testDeterministic()
{
    auto run = []() {
        TestBox box(1.0, 1.0, 1.0, 1, 1, 1, 1000);
        std::mt19937 rng(1729);
        aggregate::Placer placer(box, rng);
        std::vector<Eigen::Vector3d> sizes;
        for ( int i = 0; i < 8; ++i ) {
            sizes.emplace_back(0.06, 0.08, 0.10);
        }
        placer.placeBatch(sizes);
        std::vector<Eigen::Vector3d> centres;
        for ( const auto &inc : box.giveRealInclusions() ) {
            const auto *e = dynamic_cast<const aggregate::Ellipsoid *>( inc.get() );
            if ( e ) centres.push_back(e->giveCentre());
        }
        return centres;
    };
    auto c1 = run();
    auto c2 = run();
    EXPECT(c1.size() == c2.size(), "deterministic seed must place same count");
    bool same = ( c1.size() == c2.size() );
    for ( size_t i = 0; same && i < c1.size(); ++i ) {
        if ( !c1[i].isApprox(c2[i]) ) same = false;
    }
    EXPECT(same, "deterministic seed must produce identical centres");
}

void testNoOverlapsAfterPlacement()
{
    TestBox box(1.0, 1.0, 1.0, 1, 1, 1, 1000);
    std::mt19937 rng(99);
    aggregate::Placer placer(box, rng);

    for ( int i = 0; i < 6; ++i ) {
        placer.placeOne(Eigen::Vector3d(0.08, 0.09, 0.10));
    }

    // Real inclusions among themselves: no pairwise overlap.
    const auto &reals = box.giveRealInclusions();
    bool ok = true;
    for ( size_t i = 0; ok && i < reals.size(); ++i ) {
        const auto *a = dynamic_cast<const aggregate::Ellipsoid *>( reals[i].get() );
        for ( size_t j = i + 1; ok && j < reals.size(); ++j ) {
            const auto *b = dynamic_cast<const aggregate::Ellipsoid *>( reals[j].get() );
            if ( a && b && aggregate::ellipsoidsIntersect(*a, *b) ) {
                ok = false;
            }
        }
    }
    EXPECT(ok, "placed real inclusions must be pairwise non-overlapping");
}

} // namespace

void testFibrePlaceSparse()
{
    TestBox box(1.0, 1.0, 1.0, 1, 1, 1, 1000);
    std::mt19937 rng(11);
    aggregate::Placer placer(box, rng);

    int placed = 0;
    for ( int i = 0; i < 5; ++i ) {
        if ( placer.placeFibre(0.1, 0.005) ) {
            ++placed;
        }
    }
    EXPECT(placed == 5, "5 fibres should fit in unit box with no obstacles");

    int fibreCount = 0;
    for ( const auto &inc : box.giveRealInclusions() ) {
        if ( dynamic_cast<const aggregate::Fibre *>( inc.get() ) ) {
            ++fibreCount;
        }
    }
    EXPECT(fibreCount == 5, "Box should hold 5 real fibre inclusions");
}

void testFibreAvoidsEllipsoid()
{
    TestBox box(1.0, 1.0, 1.0, 0, 0, 0, 500);
    std::mt19937 rng(17);
    aggregate::Placer placer(box, rng);

    // Plant a single fat ellipsoid in the centre. Place a number of fibres;
    // none should pierce it.
    placer.placeOne(Eigen::Vector3d(0.4, 0.4, 0.4));
    placer.placeFibreBatch(20, 0.05, 0.005);

    const aggregate::Ellipsoid *plantedEllipsoid = nullptr;
    for ( const auto &inc : box.giveRealInclusions() ) {
        if ( const auto *e = dynamic_cast<const aggregate::Ellipsoid *>( inc.get() ) ) {
            plantedEllipsoid = e;
        }
    }
    EXPECT(plantedEllipsoid != nullptr, "ellipsoid must have been planted");

    bool anyFibreOverlaps = false;
    if ( plantedEllipsoid ) {
        for ( const auto &inc : box.giveRealInclusions() ) {
            if ( const auto *f = dynamic_cast<const aggregate::Fibre *>( inc.get() ) ) {
                if ( aggregate::fibreIntersectsEllipsoid(*f, *plantedEllipsoid) ) {
                    anyFibreOverlaps = true;
                    break;
                }
            }
        }
    }
    EXPECT(!anyFibreOverlaps, "no placed fibre may intersect the planted ellipsoid");
}

void testFibreDeterministic()
{
    auto run = []() {
        TestBox box(1.0, 1.0, 1.0, 1, 1, 1, 1000);
        std::mt19937 rng(31337);
        aggregate::Placer placer(box, rng);
        placer.placeFibreBatch(6, 0.15, 0.002);
        std::vector<Eigen::Vector3d> centres;
        for ( const auto &inc : box.giveRealInclusions() ) {
            if ( const auto *f = dynamic_cast<const aggregate::Fibre *>( inc.get() ) ) {
                centres.push_back(f->giveCentre());
            }
        }
        return centres;
    };
    auto c1 = run();
    auto c2 = run();
    EXPECT(c1.size() == c2.size(), "deterministic seed must place same fibre count");
    bool same = ( c1.size() == c2.size() );
    for ( size_t i = 0; same && i < c1.size(); ++i ) {
        if ( !c1[i].isApprox(c2[i]) ) same = false;
    }
    EXPECT(same, "deterministic seed must produce identical fibre centres");
}

int main()
{
    testGhostShifts();
    testPlaceSparse();
    testPlaceTooTight();
    testDeterministic();
    testNoOverlapsAfterPlacement();
    testFibrePlaceSparse();
    testFibreAvoidsEllipsoid();
    testFibreDeterministic();

    if ( failures == 0 ) {
        std::cout << "placer unit tests passed\n";
        return EXIT_SUCCESS;
    }
    std::cerr << failures << " placer unit test(s) failed\n";
    return EXIT_FAILURE;
}
