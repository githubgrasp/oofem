#pragma once

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

#include "inclusion.h"

namespace aggregate {

/// Parameters captured by the optional `#@grading` directive.
struct GradingParameters
{
    bool present = false;
    double dmin = 0.0;
    double dmax = 0.0;
    double aggregateFraction = 0.0;
};

/// Parameters captured by the optional `#@fibres` directive.
struct FibreParameters
{
    bool present = false;
    double fibreFraction = 0.0;
    double fibreLength = 0.0;
    double fibreDiameter = 0.0;
};

/**
 * Periodic representative volume element (RVE) used as the placement domain.
 *
 * Holds the bounding box dimensions, per-axis periodicity flags, the random
 * number generator seed, the iteration cap for the trial-and-error placer,
 * and the lists of inclusions placed inside.
 *
 * Configured by parsing a control file of leading-`#@` directives, mirroring
 * the convention used in `src/generator/` and `src/converter/`.
 */
class Box
{
public:
    Box();

    /// Parse a control file of `#@`-directives and populate the box.
    void readControlRecords(const std::string &fileName);

    /// Write the inclusion packing to disk in the shared packing format.
    void writePackingFile() const;

    /// Path of the packing file declared by `#@output`.
    const std::string &giveOutputFileName() const { return outputFileName; }

    /// RVE side lengths in x, y, z.
    const Eigen::Vector3d &giveDimensions() const { return dimensions; }

    /// Per-axis periodicity flag (1 = periodic, 0 = real boundary).
    const Eigen::Vector3i &givePeriodicityFlag() const { return periodicityFlag; }

    /// Maximum trial-and-error attempts for placing one inclusion.
    int giveMaximumIterations() const { return maximumIterations; }

    /// Seed for the std::mt19937 RNG used by the grading curve and placer.
    unsigned int giveRandomSeed() const { return randomSeed; }

    /// Take ownership of an inclusion that should appear in the packing file.
    void addReal(std::unique_ptr<Inclusion> inclusion);

    /// Take ownership of a periodic-image inclusion used for collision tests only.
    void addGhost(std::unique_ptr<Inclusion> inclusion);

    /// Real inclusions accumulated by the placer.
    const std::vector<std::unique_ptr<Inclusion>> &giveRealInclusions() const
    { return realInclusions; }

    /// Periodic-image inclusions accumulated by the placer (not written to disk).
    const std::vector<std::unique_ptr<Inclusion>> &giveGhostInclusions() const
    { return ghostInclusions; }

    /// Aggregate grading parameters (may be absent, see `present`).
    const GradingParameters &giveGradingParameters() const { return grading; }

    /// Fibre placement parameters (may be absent, see `present`).
    const FibreParameters &giveFibreParameters() const { return fibres; }

    /// Optional VTU output path (empty if `#@vtu` was not specified).
    const std::string &giveVtuFileName() const { return vtuFileName; }

private:
    /// Apply one `#@`-directive line. Splits off the keyword and dispatches.
    void applyDirective(const std::string &line);

    std::string outputFileName;
    Eigen::Vector3d dimensions;
    Eigen::Vector3i periodicityFlag;
    unsigned int randomSeed = 1;
    int maximumIterations = 10000;

    GradingParameters grading;
    FibreParameters fibres;
    std::string vtuFileName;

    std::vector<std::unique_ptr<Inclusion>> realInclusions;
    std::vector<std::unique_ptr<Inclusion>> ghostInclusions;
};

} // namespace aggregate
