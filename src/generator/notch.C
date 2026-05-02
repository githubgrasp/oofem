#include "notch.h"
#include "generatorerror.h"


Notch::Notch(int n, Grid *aGrid) : GridComponent(n, aGrid) {}

Notch::~Notch() {}


void Notch::initializeFromTokens(std::istringstream &iss)
{
    bool gotBox = false;
    std::string tok;
    while ( iss >> tok ) {
        if ( tok == "box" ) {
            int n;
            iss >> n;
            box.resize(n);
            for ( int i = 1; i <= n; ++i ) {
                iss >> box.at(i);
            }
            gotBox = true;
        } else if ( tok == "edgerefine" ) {
            iss >> edgeRefine;
        } else if ( tok == "surfacerefine" ) {
            iss >> surfaceRefine;
        } else {
            generator::errorf("Notch::initializeFromTokens: unknown keyword '%s'", tok.c_str() );
        }
    }
    if ( !gotBox ) {
        generator::error("Notch::initializeFromTokens: missing 'box' keyword");
    }
    if ( box.giveSize() != 4 && box.giveSize() != 6 ) {
        generator::errorf("Notch::initializeFromTokens: 'box' must have 4 (2D) or 6 (3D) entries, got %d",
                          box.giveSize() );
    }
}


bool Notch::containsStrictly(const oofem::FloatArray &coord, double tol) const
{
    if ( box.giveSize() == 6 ) {
        return coord.at(1) > box.at(1) + tol && coord.at(1) < box.at(4) - tol &&
               coord.at(2) > box.at(2) + tol && coord.at(2) < box.at(5) - tol &&
               coord.at(3) > box.at(3) + tol && coord.at(3) < box.at(6) - tol;
    }
    // 2D box: xmin ymin xmax ymax
    return coord.at(1) > box.at(1) + tol && coord.at(1) < box.at(3) - tol &&
           coord.at(2) > box.at(2) + tol && coord.at(2) < box.at(4) - tol;
}
