#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

#include "dlim.h"
#include "gradingcurve.h"

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

double ellipsoidVolume(const Eigen::Vector3d &s)
{
    return ( 4.0 / 3.0 ) * M_PI * s(0) * s(1) * s(2);
}

void testDlimNoFibres()
{
    // No fibres → threshold collapses to dmax.
    const double dmax = 0.016;
    EXPECT(aggregate::dlim(0.0002, dmax, 0.0, 0.02, 0.7) == dmax,
           "fibre_fraction=0 must give dlim=dmax");
}

void testDlimWithFibres()
{
    // With fibres present, dlim should be strictly less than dmax and positive.
    const double dmax = 0.016;
    const double d = aggregate::dlim(0.0002, dmax, 0.005, 0.02, 0.4);
    EXPECT(d > 0.0 && d < dmax,
           "dlim with fibres must lie in (0, dmax)");

    // Larger fibre fraction shifts the threshold lower (smaller aggregates
    // are more likely to be placed before fibres).
    const double dHigher = aggregate::dlim(0.0002, dmax, 0.02, 0.02, 0.4);
    EXPECT(dHigher <= d,
           "larger fibre fraction should not increase dlim");
}

void testGradingEmpty()
{
    aggregate::GradingCurve gc(0.004, 0.016, 0.0, 1e-3);
    std::mt19937 rng(42);
    auto sizes = gc.generate(rng);
    EXPECT(sizes.empty(), "vol_frac=0 must give empty list");
}

void testGradingOrdered()
{
    aggregate::GradingCurve gc(0.004, 0.016, 0.5, 1e-3);
    std::mt19937 rng(42);
    auto sizes = gc.generate(rng);

    EXPECT(!sizes.empty(), "non-zero vol_frac must produce some aggregates");

    bool sortedDesc = true;
    bool sxLeqSyLeqSz = true;
    for ( size_t i = 0; i < sizes.size(); ++i ) {
        if ( !( sizes[i](0) <= sizes[i](1) && sizes[i](1) <= sizes[i](2) ) ) {
            sxLeqSyLeqSz = false;
        }
        if ( i > 0 && ellipsoidVolume(sizes[i - 1]) < ellipsoidVolume(sizes[i]) ) {
            sortedDesc = false;
        }
    }
    EXPECT(sortedDesc, "GradingCurve output must be sorted by volume desc");
    EXPECT(sxLeqSyLeqSz, "every triple must satisfy sx ≤ sy ≤ sz");
}

void testGradingDeterministic()
{
    aggregate::GradingCurve gc(0.004, 0.016, 0.5, 1e-3);
    std::mt19937 rng1(123);
    std::mt19937 rng2(123);
    auto a = gc.generate(rng1);
    auto b = gc.generate(rng2);
    EXPECT(a.size() == b.size(), "deterministic seed must give identical count");
    bool identical = ( a.size() == b.size() );
    for ( size_t i = 0; identical && i < a.size(); ++i ) {
        if ( !a[i].isApprox(b[i]) ) {
            identical = false;
        }
    }
    EXPECT(identical, "deterministic seed must give identical samples");
}

void testGradingVolumeWithinTarget()
{
    const double rveVol = 1e-3;
    const double frac = 0.4;
    const double dmin = 0.004;
    const double dmax = 0.016;
    aggregate::GradingCurve gc(dmin, dmax, frac, rveVol);
    std::mt19937 rng(42);
    auto sizes = gc.generate(rng);
    double total = 0.0;
    for ( const auto &s : sizes ) {
        total += ellipsoidVolume(s);
    }
    // Wriggers' formula explicitly truncates the Fuller curve below dmin
    // via the (1 − √(dmin/dmax)) factor — that volume isn't modelled and
    // never gets placed, so the modelled target is smaller than the
    // nominal aggregateFraction × rveVol.
    const double modelledTarget =
        frac * rveVol * ( 1.0 - std::sqrt(dmin / dmax) );
    EXPECT(total > 0.0 && total <= modelledTarget,
           "total placed volume must be positive and not exceed modelled target");
    EXPECT(total >= 0.7 * modelledTarget,
           "total placed volume should reach at least 70% of modelled target");
}

} // namespace

int main()
{
    testDlimNoFibres();
    testDlimWithFibres();
    testGradingEmpty();
    testGradingOrdered();
    testGradingDeterministic();
    testGradingVolumeWithinTarget();

    if ( failures == 0 ) {
        std::cout << "grading + dlim unit tests passed\n";
        return EXIT_SUCCESS;
    }
    std::cerr << failures << " grading/dlim unit test(s) failed\n";
    return EXIT_FAILURE;
}
