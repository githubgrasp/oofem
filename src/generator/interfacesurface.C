#include "interfacesurface.h"
#include "curve.h"
#include "vertex.h"
#include "generatorerror.h"


#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include "octreegridlocalizer.h"
#endif

InterfaceSurface::InterfaceSurface(int n, Grid *aGrid) : Inclusion(n, aGrid) //, coordinates()
{
    this->number = n;
}

InterfaceSurface::~InterfaceSurface()
// Destructor.
{}

int InterfaceSurface::generatePoints()
{
    return 1;
}

void InterfaceSurface::initializeFromTokens(std::istringstream &iss)
{
    refinement = 1.;

    bool gotCurves = false;
    std::string tok;
    while ( iss >> tok ) {
        if ( tok == "curves" ) {
            int n;
            iss >> n;
            curves.resize(n);
            for ( int i = 1; i <= n; ++i ) {
                iss >> curves.at(i);
            }
            gotCurves = true;
        } else if ( tok == "refine" ) {
            iss >> refinement;
        } else {
            generator::errorf("InterfaceSurface::initializeFromTokens: unknown keyword '%s'", tok.c_str());
        }
    }
    if ( !gotCurves ) {
        generator::error("InterfaceSurface::initializeFromTokens: missing 'curves' keyword");
    }
}


InterfaceSurface *InterfaceSurface::ofType()
{
    InterfaceSurface *surface;

    surface = new InterfaceSurface(number, grid);

    return surface;
}
