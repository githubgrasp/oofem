#include "fibre.h"

#include <ostream>

#include "aggregateerror.h"

namespace aggregate {

Fibre::Fibre(int n,
             const Eigen::Vector3d &a,
             const Eigen::Vector3d &b,
             double d)
    : Inclusion(n), endpointA(a), endpointB(b), diameter(d)
{}

Fibre::Fibre(int n)
    : Inclusion(n),
      endpointA(Eigen::Vector3d::Zero()),
      endpointB(Eigen::Vector3d::Zero())
{}

void Fibre::initializeFromTokens(std::istringstream &iss)
{
    std::string tok;
    bool gotEndpoints = false, gotDiameter = false;
    while ( iss >> tok ) {
        if ( tok == "endpointcoords" ) {
            int n;
            iss >> n;
            if ( n != 6 ) {
                errorf("Fibre::initializeFromTokens: 'endpointcoords' expects 6 values, got %d", n);
            }
            iss >> endpointA(0) >> endpointA(1) >> endpointA(2)
                >> endpointB(0) >> endpointB(1) >> endpointB(2);
            gotEndpoints = true;
        } else if ( tok == "diameter" ) {
            iss >> diameter;
            gotDiameter = true;
        } else {
            errorf("Fibre::initializeFromTokens: unknown keyword '%s'", tok.c_str());
        }
    }
    if ( !gotEndpoints ) {
        error("Fibre::initializeFromTokens: missing 'endpointcoords'");
    }
    if ( !gotDiameter || diameter <= 0.0 ) {
        errorf("Fibre::initializeFromTokens: missing or non-positive 'diameter' (%.3e)", diameter);
    }
    if ( ( endpointA - endpointB ).squaredNorm() == 0.0 ) {
        error("Fibre::initializeFromTokens: zero-length fibre");
    }
}

void Fibre::writeTo(std::ostream &os) const
{
    os << "fibre " << number
       << " endpointcoords 6 "
       << endpointA(0) << ' ' << endpointA(1) << ' ' << endpointA(2) << ' '
       << endpointB(0) << ' ' << endpointB(1) << ' ' << endpointB(2)
       << " diameter " << diameter
       << '\n';
}

} // namespace aggregate
