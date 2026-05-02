#pragma once

#include <Eigen/Dense>

#include "inclusion.h"

namespace aggregate {

/// 2D circular inclusion — the 2D analog of the degenerate-sphere case of
/// `Ellipsoid`. Centre is 2D, with a single radius. Used both as the
/// natural disk for textile-fibre cross-sections and as the Stage 4
/// stand-in for any 2D aggregate (we may add a general `Ellipse2D` later).
class Disk : public Inclusion
{
public:
    /// Construct a disk with given id, centre (cx, cy), and radius.
    Disk(int n, const Eigen::Vector2d &centre, double radius);

    /// Sequential id, written to the packing file.
    std::string typeName() const override { return "disk"; }

    /// Append `disk <id> centre 2 cx cy radius r` to the stream.
    void writeTo(std::ostream &os) const override;

    /// Centre of the disk.
    const Eigen::Vector2d &giveCentre() const { return centre; }

    /// Radius of the disk.
    double giveRadius() const { return radius; }

private:
    Eigen::Vector2d centre;
    double radius;
};

} // namespace aggregate
