#include "interfacedisk.h"
#include "vertex.h"
#include "generatorerror.h"

#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include "octreegridlocalizer.h"
#endif


InterfaceDisk::InterfaceDisk(int n, Grid *aGrid) : Inclusion(n, aGrid) {}

InterfaceDisk::~InterfaceDisk() {}


void InterfaceDisk::initializeFromTokens(std::istringstream &iss)
{
    refinement = 1.;

    bool gotCentre = false, gotRadius = false;
    bool gotItz = false;
    std::string tok;
    while ( iss >> tok ) {
        if ( tok == "centre" ) {
            int n;
            iss >> n;
            if ( n != 2 ) {
                generator::error("InterfaceDisk::initializeFromTokens: 'centre' must have 2 entries (cx cy)");
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
        } else if ( tok == "itz" ) {
            iss >> itzThickness;
            gotItz = true;
        } else {
            generator::errorf("InterfaceDisk::initializeFromTokens: unknown keyword '%s'", tok.c_str() );
        }
    }
    if ( !gotCentre || !gotRadius ) {
        generator::error("InterfaceDisk::initializeFromTokens: 'centre' and 'radius' are required");
    }
    diameter = 2. * radius;
    if ( !gotItz ) {
        itzThickness = refinement * diameter;
    }
}


int InterfaceDisk::generatePoints()
{
    const double cx = centre.at(1), cy = centre.at(2);
    const double maxIter = grid->giveMaximumIterations();

    int randomIntegerOne = grid->giveRandomInteger() - 1;

    oofem::FloatArray random(3);
    random.at(3) = 0.;

    printf("Points on inclusion disk circumference\n");

    for ( int i = 0; i < maxIter; ++i ) {
        const double theta = 2. * M_PI * grid->ran1(& randomIntegerOne);

        random.at(1) = cx + radius * cos(theta);
        random.at(2) = cy + radius * sin(theta);

        // Skip perimeter points that fall outside the meshing region — a
        // disk centred near a non-periodic face may extend past it, and
        // those slices have no mesh to discretise. The boundary-crossing
        // anchors (placed earlier by Grid::generateInclusionBoundaryAnchors)
        // already cover the surviving part of the disk surface at the face.
        if ( !grid->isPointInsideMeshingRegion(random) ) {
            continue;
        }

        const int flag = grid->giveGridLocalizer()->checkNodesWithinBox(random,
                             refinement * grid->giveDiameter(random) );
        if ( flag == 0 && grid->addVertex(random) ) {
            i = 0;
            const int innerId = grid->giveNumberOfVertices();

            // Sister point on the outer ITZ ring. The inner-ring point and
            // sister are intentionally `itz` apart so qhull's Voronoi places
            // a cell wall exactly on the disk surface (the perpendicular
            // bisector of the close pair). Distance-checking the sister
            // against the inner ring would always reject it for small `itz`,
            // destroying that mechanism — so the inner ring is excluded from
            // the sister's neighbour check while every other vertex (notably
            // region-boundary points placed earlier in the pipeline) still
            // counts. This stops the sister from punching through the
            // boundary discretisation when an inclusion sits close to a
            // rect/prism face.
            random.at(1) = cx + ( radius + itzThickness ) * cos(theta);
            random.at(2) = cy + ( radius + itzThickness ) * sin(theta);
            GridLocalizer::nodeContainerType nbrs;
            grid->giveGridLocalizer()->giveAllNodesWithinBox(nbrs, random,
                                          refinement * grid->giveDiameter(random) );
            bool hasOtherNeighbour = false;
            for ( int id : nbrs ) {
                if ( id != innerId ) { hasOtherNeighbour = true; break; }
            }
            if ( !hasOtherNeighbour ) {
                grid->addVertex(random);
            }
        }
    }

    printf("Completed inclusion disk circumference\n");
    return 1;
}
