#include "disk.h"
#include "vertex.h"
#include "generatorerror.h"

#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include "octreegridlocalizer.h"
 #include <iostream>
#endif


Disk::Disk(int n, Grid *aGrid) : Region(n, aGrid) {}

Disk::~Disk() {}


void Disk::initializeFromTokens(std::istringstream &iss)
{
    refinement = 1.;

    bool gotCentre = false, gotRadius = false;
    std::string tok;
    while ( iss >> tok ) {
        if ( tok == "centre" ) {
            int n;
            iss >> n;
            if ( n != 2 ) {
                generator::error("Disk::initializeFromTokens: 'centre' must have 2 entries (cx cy)");
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
            generator::errorf("Disk::initializeFromTokens: unknown keyword '%s'", tok.c_str() );
        }
    }
    if ( !gotCentre || !gotRadius ) {
        generator::error("Disk::initializeFromTokens: 'centre' and 'radius' are required");
    }
    xlength = 2. * radius;
    ylength = 2. * radius;
    zlength = 0.;
}


int Disk::generatePoints()
{
    const double cx = centre.at(1), cy = centre.at(2);
    const double diam = grid->diameter;
    const double maxIter = grid->giveMaximumIterations();

    int randomIntegerOne = grid->giveRandomInteger() - 1;
    int randomIntegerTwo = grid->giveRandomInteger() - 2;

    oofem::FloatArray random(3);
    random.at(3) = 0.;

    // Centre
    random.at(1) = cx;
    random.at(2) = cy;
    grid->addVertex(random);

    printf("Points on disk circumference\n");

    // Random ring on circumference (Bolander-style: random angles, distance check).
    for ( int i = 0; i < maxIter; ++i ) {
        const double theta = 2. * M_PI * grid->ran1(& randomIntegerOne);
        random.at(1) = cx + radius * cos(theta);
        random.at(2) = cy + radius * sin(theta);

        const int flag = grid->giveGridLocalizer()->checkNodesWithinBox(random,
                             refinement * grid->giveDiameter(random) );
        if ( flag == 0 && grid->addVertex(random) ) {
            i = 0;
        }
    }

    printf("Points in disk interior\n");

    // Random interior fill (sample (theta, sqrt(u)*R) for uniform area).
    for ( int i = 0; i < maxIter; ++i ) {
        const double theta = 2. * M_PI * grid->ran1(& randomIntegerOne);
        const double r = ( radius - diam ) * sqrt(grid->ran1(& randomIntegerTwo) );

        random.at(1) = cx + r * cos(theta);
        random.at(2) = cy + r * sin(theta);

        const int flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, diam);
        if ( flag == 0 && grid->addVertex(random) ) {
            i = 0;

            // Mirror across the disk boundary so the dual cells partition
            // cleanly along the circle.
            const double rMirror = 2. * radius - r;
            random.at(1) = cx + rMirror * cos(theta);
            random.at(2) = cy + rMirror * sin(theta);
            grid->addVertex(random);
        }
    }

    return 1;
}
