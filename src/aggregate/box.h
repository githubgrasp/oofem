#pragma once

#include <array>
#include <cstdio>
#include <string>

#include "floatarray.h"
#include "intarray.h"

namespace aggregate {

/**
 * Periodic representative volume element (RVE) used as the placement domain.
 *
 * Holds the bounding box dimensions, per-axis periodicity flags, the random
 * number generator seed, the iteration cap for the trial-and-error placer,
 * and (later) the lists of inclusions placed inside.
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

private:
    /// Apply one `#@`-directive line. Splits off the keyword and dispatches.
    void applyDirective(const std::string &line);

    std::string outputFileName;

    /// RVE side lengths in x, y, z.
    oofem::FloatArray dimensions;

    /// Per-axis periodicity flag (1 = periodic, 0 = real boundary).
    oofem::IntArray periodicityFlag;

    /// Numerical Recipes ran1 seed (negative integer triggers re-init).
    int randomSeed = -1;

    /// Maximum trial-and-error attempts for placing one inclusion.
    int maximumIterations = 10000;
};

} // namespace aggregate
