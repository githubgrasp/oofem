#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

#include <Eigen/Dense>

#include "box.h"
#include "dlim.h"
#include "gradingcurve.h"
#include "placer.h"

namespace {

double ellipsoidVolume(const Eigen::Vector3d &s)
{
    return ( 4.0 / 3.0 ) * M_PI * s(0) * s(1) * s(2);
}

/// Splits the (volume-descending) `sizes` list at the first index whose
/// volume is not greater than `eqVlim`. Front partition (group 1) is the
/// large aggregates placed before fibres; back partition (group 2) is the
/// small aggregates placed after.
std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>>
splitByVolume(const std::vector<Eigen::Vector3d> &sizes, double eqVlim)
{
    std::vector<Eigen::Vector3d> group1, group2;
    for ( const auto &s : sizes ) {
        if ( ellipsoidVolume(s) > eqVlim ) {
            group1.push_back(s);
        } else {
            group2.push_back(s);
        }
    }
    return { group1, group2 };
}

int fibreCountFromFraction(const aggregate::FibreParameters &fp,
                           double rveVolume)
{
    if ( fp.fibreFraction <= 0.0 ) {
        return 0;
    }
    const double radius = 0.5 * fp.fibreDiameter;
    const double singleFibreVolume = M_PI * radius * radius * fp.fibreLength;
    return static_cast<int>( rveVolume * fp.fibreFraction / singleFibreVolume );
}

} // namespace

int main(int argc, char *argv[])
{
    if ( argc != 2 ) {
        std::cerr << "Usage: aggregate <control-file>\n";
        return EXIT_FAILURE;
    }

    aggregate::Box box;
    box.readControlRecords(argv[1]);

    const Eigen::Vector3d &dims = box.giveDimensions();
    const double rveVolume = dims(0) * dims(1) * dims(2);
    std::mt19937 rng(box.giveRandomSeed());
    aggregate::Placer placer(box, rng);

    const auto &grading = box.giveGradingParameters();
    const auto &fibres = box.giveFibreParameters();

    std::vector<Eigen::Vector3d> group1, group2;
    if ( grading.present ) {
        aggregate::GradingCurve gc(grading.dmin,
                                   grading.dmax,
                                   grading.aggregateFraction,
                                   rveVolume);
        const auto sizes = gc.generate(rng);

        const double dlimValue =
            aggregate::dlim(fibres.fibreDiameter,
                            grading.dmax,
                            fibres.fibreFraction,
                            fibres.fibreLength,
                            grading.aggregateFraction);
        constexpr double alpha = 1.5;
        const double eqVlim =
            ( 4.0 / 3.0 ) * M_PI * std::pow(0.5 * alpha * dlimValue, 3);
        std::tie(group1, group2) = splitByVolume(sizes, eqVlim);
    }

    if ( !group1.empty() ) {
        const int failed = placer.placeBatch(group1);
        std::cout << "Group 1 ellipsoids: placed " << ( group1.size() - failed )
                  << " / " << group1.size() << "\n";
    }

    if ( fibres.present ) {
        const int count = fibreCountFromFraction(fibres, rveVolume);
        if ( count > 0 ) {
            const int failed =
                placer.placeFibreBatch(count, fibres.fibreLength, fibres.fibreDiameter);
            std::cout << "Fibres: placed " << ( count - failed )
                      << " / " << count << "\n";
        }
    }

    if ( !group2.empty() ) {
        const int failed = placer.placeBatch(group2);
        std::cout << "Group 2 ellipsoids: placed " << ( group2.size() - failed )
                  << " / " << group2.size() << "\n";
    }

    box.writePackingFile();
    std::cout << "Packing written to " << box.giveOutputFileName() << "\n";
    return EXIT_SUCCESS;
}
