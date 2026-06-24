#include "holedisk.h"
#include "vertex.h"
#include "generatorerror.h"

#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include "octreegridlocalizer.h"
#endif


HoleDisk::HoleDisk(int n, Grid *aGrid) : Inclusion(n, aGrid) {}

HoleDisk::~HoleDisk() {}


void HoleDisk::initializeFromTokens(std::istringstream &iss)
{
    refinement = 1.;

    bool gotCentre = false, gotRadius = false;
    std::string tok;
    while ( iss >> tok ) {
        if ( tok == "centre" ) {
            int n;
            iss >> n;
            if ( n != 2 ) {
                generator::error("HoleDisk::initializeFromTokens: 'centre' must have 2 entries (cx cy)");
            }
            centre.resize(2);
            for ( int i = 1; i <= 2; ++i ) {
                iss >> centre.at(i);
            }
            gotCentre = true;
        } else if ( tok == "radius" ) {
            iss >> radius;
            gotRadius = true;
        } else if ( tok == "refine" ) {
            iss >> refinement;
        } else {
            generator::errorf("HoleDisk::initializeFromTokens: unknown keyword '%s'", tok.c_str() );
        }
    }
    if ( !gotCentre || !gotRadius ) {
        generator::error("HoleDisk::initializeFromTokens: 'centre' and 'radius' are required");
    }
}


bool HoleDisk::containsStrictly(const oofem::FloatArray &coord, double tol) const
{
    const double dx = coord.at(1) - centre.at(1);
    const double dy = coord.at(2) - centre.at(2);
    return std::sqrt(dx * dx + dy * dy) < radius - tol;
}


int HoleDisk::generatePoints()
{
    const double cx = centre.at(1), cy = centre.at(2);
    const double maxIter = grid->giveMaximumIterations();
    const double diam = grid->diameter;

    int randomIntegerOne = grid->giveRandomInteger() - 1;
    int randomIntegerTwo = grid->giveRandomInteger() - 2;

    oofem::FloatArray random(3);
    random.at(3) = 0.;

    printf("Points on hole circumference\n");

    // Rim ring on the hole/matrix boundary — these are the mechanical nodes
    // that stay ON the boundary (no converter shift needed for them).
    for ( int i = 0; i < maxIter; ++i ) {
        const double theta = 2. * M_PI * grid->ran1(& randomIntegerOne);
        random.at(1) = cx + radius * cos(theta);
        random.at(2) = cy + radius * sin(theta);

        if ( !grid->isPointInsideMeshingRegion(random) ) {
            continue;
        }
        const int flag = grid->giveGridLocalizer()->checkNodesWithinBox(random,
                             refinement * grid->giveDiameter(random) );
        if ( flag == 0 && grid->addVertex(random) ) {
            i = 0;

            // Paired sacrificial point one diameter inside the rim at the same
            // angle (mirror of the annulus inner-mirror): guarantees every rim
            // node has an inward neighbour, so the rim Delaunay comes out as a
            // clean polygon (no rim node connecting across the void).
            const double rIn = radius - diam;
            if ( rIn > 0. ) {
                oofem::FloatArray inner(3);
                inner.at(1) = cx + rIn * cos(theta);
                inner.at(2) = cy + rIn * sin(theta);
                inner.at(3) = 0.;
                grid->addVertex(inner, /*allowInsideHole*/ true);
            }
        }
    }

    printf("Sacrificial points filling hole interior\n");

    // Sacrificial interior fill (area-uniform in [0, radius)). These are placed
    // strictly inside the hole, so they bypass the matrix hole-rejection; the
    // converter drops them. They keep the rim Delaunay clean (rim nodes connect
    // inward to these instead of chording across the void).
    const double rMax = radius - diam;
    for ( int i = 0; i < maxIter && rMax > 0.; ++i ) {
        const double theta = 2. * M_PI * grid->ran1(& randomIntegerOne);
        const double u = grid->ran1(& randomIntegerTwo);
        const double r = sqrt(u) * rMax;
        random.at(1) = cx + r * cos(theta);
        random.at(2) = cy + r * sin(theta);

        const int flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, diam);
        if ( flag == 0 && grid->addVertex(random, /*allowInsideHole*/ true) ) {
            i = 0;
        }
    }

    printf("Completed hole disk\n");
    return 1;
}
