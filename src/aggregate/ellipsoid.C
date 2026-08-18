#include "ellipsoid.h"

#include <cmath>
#include <ostream>

#include "aggregateerror.h"

namespace aggregate {

Ellipsoid::Ellipsoid(int n,
                     const Eigen::Vector3d &c,
                     const Eigen::Vector3d &angles,
                     const Eigen::Vector3d &axes)
    : Inclusion(n), centre(c), eulerAngles(angles), semiAxes(axes)
{}

Ellipsoid::Ellipsoid(int n)
    : Inclusion(n),
      centre(Eigen::Vector3d::Zero()),
      eulerAngles(Eigen::Vector3d::Zero()),
      semiAxes(Eigen::Vector3d::Zero())
{}

void Ellipsoid::initializeFromTokens(std::istringstream &iss)
{
    std::string tok;
    bool gotCentre = false, gotAngles = false, gotRadii = false;
    while ( iss >> tok ) {
        if ( tok == "centre" ) {
            int n;
            iss >> n;
            if ( n != 3 ) {
                errorf("Ellipsoid::initializeFromTokens: 'centre' expects 3 values, got %d", n);
            }
            iss >> centre(0) >> centre(1) >> centre(2);
            gotCentre = true;
        } else if ( tok == "angles" ) {
            int n;
            iss >> n;
            if ( n != 3 ) {
                errorf("Ellipsoid::initializeFromTokens: 'angles' expects 3 values, got %d", n);
            }
            iss >> eulerAngles(0) >> eulerAngles(1) >> eulerAngles(2);
            gotAngles = true;
        } else if ( tok == "radii" ) {
            int n;
            iss >> n;
            if ( n != 3 ) {
                errorf("Ellipsoid::initializeFromTokens: 'radii' expects 3 values, got %d", n);
            }
            iss >> semiAxes(0) >> semiAxes(1) >> semiAxes(2);
            gotRadii = true;
        } else {
            errorf("Ellipsoid::initializeFromTokens: unknown keyword '%s'", tok.c_str());
        }
    }
    if ( !gotCentre || !gotAngles || !gotRadii ) {
        error("Ellipsoid::initializeFromTokens: missing centre/angles/radii");
    }
    if ( ( semiAxes.array() <= 0.0 ).any() ) {
        errorf("Ellipsoid::initializeFromTokens: non-positive semi-axis (%.3e %.3e %.3e)",
               semiAxes(0), semiAxes(1), semiAxes(2));
    }
}

void Ellipsoid::writeTo(std::ostream &os) const
{
    // Degenerate ellipsoid (sx == sy == sz) is a sphere — emit the simpler
    // form so downstream tools that only know how to seed sphere surfaces
    // (e.g. the generator's InterfaceSphere) can consume it directly.
    if ( semiAxes(0) == semiAxes(1) && semiAxes(1) == semiAxes(2) ) {
        os << "sphere " << number
           << " centre 3 " << centre(0) << ' ' << centre(1) << ' ' << centre(2)
           << " radius " << semiAxes(0)
           << '\n';
        return;
    }
    os << "ellipsoid " << number
       << " centre 3 " << centre(0) << ' ' << centre(1) << ' ' << centre(2)
       << " angles 3 " << eulerAngles(0) << ' ' << eulerAngles(1) << ' ' << eulerAngles(2)
       << " radii 3 " << semiAxes(0) << ' ' << semiAxes(1) << ' ' << semiAxes(2)
       << '\n';
}

Eigen::Matrix3d Ellipsoid::giveRotation() const
{
    const double cx = std::cos(eulerAngles(0)), sx = std::sin(eulerAngles(0));
    const double cy = std::cos(eulerAngles(1)), sy = std::sin(eulerAngles(1));
    const double cz = std::cos(eulerAngles(2)), sz = std::sin(eulerAngles(2));

    Eigen::Matrix3d Rx;
    Rx << 1, 0, 0,
          0, cx, -sx,
          0, sx, cx;
    Eigen::Matrix3d Ry;
    Ry << cy, 0, sy,
          0, 1, 0,
          -sy, 0, cy;
    Eigen::Matrix3d Rz;
    Rz << cz, -sz, 0,
          sz, cz, 0,
          0, 0, 1;
    return Rz * Ry * Rx;
}

Eigen::Matrix4d Ellipsoid::giveQuadricForm() const
{
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T(3, 0) = -centre(0);
    T(3, 1) = -centre(1);
    T(3, 2) = -centre(2);

    Eigen::Matrix4d R = Eigen::Matrix4d::Identity();
    R.topLeftCorner<3, 3>() = giveRotation();

    Eigen::Matrix4d B = Eigen::Matrix4d::Zero();
    B(0, 0) = 1.0 / ( semiAxes(0) * semiAxes(0) );
    B(1, 1) = 1.0 / ( semiAxes(1) * semiAxes(1) );
    B(2, 2) = 1.0 / ( semiAxes(2) * semiAxes(2) );
    B(3, 3) = -1.0;

    return T * R * B * R.transpose() * T.transpose();
}

} // namespace aggregate
