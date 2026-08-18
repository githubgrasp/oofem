#pragma once

#include <Eigen/Dense>
#include <iosfwd>
#include <sstream>
#include <string>

#include "inclusion.h"

namespace aggregate {

/**
 * Straight cylindrical fibre, parameterised by its two endpoints and a
 * single diameter. Treated as a line segment (the cylindrical body is
 * accounted for via the diameter when checking overlap with other
 * inclusions).
 */
class Fibre : public Inclusion
{
public:
    /// Construct a fully-specified fibre between two world-coordinate endpoints.
    Fibre(int n,
          const Eigen::Vector3d &endpointA,
          const Eigen::Vector3d &endpointB,
          double fibreDiameter);

    /// Construct an uninitialised fibre; populated via `initializeFromTokens`.
    explicit Fibre(int n);

    /// Parse "endpointcoords 6 ax ay az bx by bz diameter d".
    void initializeFromTokens(std::istringstream &iss);

    std::string typeName() const override { return "fibre"; }
    void writeTo(std::ostream &os) const override;

    /// First endpoint (world coordinates).
    const Eigen::Vector3d &giveEndpointA() const { return endpointA; }

    /// Second endpoint (world coordinates).
    const Eigen::Vector3d &giveEndpointB() const { return endpointB; }

    /// Midpoint of the line segment AB.
    Eigen::Vector3d giveCentre() const { return 0.5 * ( endpointA + endpointB ); }

    /// Cylindrical fibre diameter (used for volume budgets, not for geometric overlap).
    double giveDiameter() const { return diameter; }

    /// Half the segment length |AB|/2 — used as a cheap pre-test radius.
    double giveHalfLength() const { return 0.5 * ( endpointB - endpointA ).norm(); }

private:
    Eigen::Vector3d endpointA;
    Eigen::Vector3d endpointB;
    double diameter = 0.0;
};

} // namespace aggregate
