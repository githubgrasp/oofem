#include "rect.h"
#include "vertex.h"
#include "generatorerror.h"

#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include <iostream>
 #include "octreegridlocalizer.h"
#endif


Rect::Rect(int n, Grid *aGrid) : Region(n, aGrid) {}

Rect::~Rect() {}


void Rect::initializeFromTokens(std::istringstream &iss)
{
    bool gotBox = false;
    refinement = 1.;
    std::string tok;
    while ( iss >> tok ) {
        if ( tok == "box" ) {
            int n;
            iss >> n;
            if ( n != 4 ) {
                generator::error("Rect::initializeFromTokens: 'box' must have 4 entries (xmin ymin xmax ymax)");
            }
            box.resize(4);
            for ( int i = 1; i <= 4; ++i ) {
                iss >> box.at(i);
            }
            gotBox = true;
        } else if ( tok == "refine" ) {
            iss >> refinement;
        } else if ( tok == "edgerefine" ) {
            iss >> edgeRefine;
        } else if ( tok == "regionrefine" ) {
            iss >> regionRefine;
        } else {
            generator::errorf("Rect::initializeFromTokens: unknown keyword '%s'", tok.c_str() );
        }
    }
    if ( !gotBox ) {
        generator::error("Rect::initializeFromTokens: missing 'box' keyword");
    }
    xlength = box.at(3) - box.at(1);
    ylength = box.at(4) - box.at(2);
    zlength = 0.;
}


void Rect::defineBoundaries(oofem::FloatArray &boundaries)
{
    boundaries.resize(6);
    boundaries.at(1) = box.at(1);  // xmin
    boundaries.at(2) = box.at(3);  // xmax
    boundaries.at(3) = box.at(2);  // ymin
    boundaries.at(4) = box.at(4);  // ymax
    boundaries.at(5) = 0.;         // zmin (degenerate)
    boundaries.at(6) = 0.;         // zmax
}


int Rect::generatePoints()
{
    printf("Generating points for rect\n");

    const double xmin = box.at(1), ymin = box.at(2);
    const double xmax = box.at(3), ymax = box.at(4);
    const double dx = xmax - xmin, dy = ymax - ymin;

    int randomIntegerOne = grid->giveRandomInteger() - 1;
    int randomIntegerTwo = grid->giveRandomInteger() - 2;

    oofem::IntArray periodicityFlag;
    grid->givePeriodicityFlag(periodicityFlag);
    // Pad to 3 entries (always treat z as non-periodic in 2D).
    if ( periodicityFlag.giveSize() < 3 ) {
        periodicityFlag.resize(3);
        periodicityFlag.at(3) = 0;
    }

    const int randomFlag = grid->giveRandomFlag();
    const double diam = grid->diameter;

    // Border margin: ranflag 0 keeps vertices off the boundary (TOL margin),
    // ranflag 1 leaves a full edge-spacing buffer so the regular edge points
    // fall on the boundary cleanly.
    double borderX = 0., borderY = 0.;
    if ( randomFlag == 0 ) {
        borderX = grid->TOL;
        borderY = grid->TOL;
    } else {
        borderX = ( periodicityFlag.at(1) == 1 ) ? grid->TOL : edgeRefine * diam;
        borderY = ( periodicityFlag.at(2) == 1 ) ? grid->TOL : edgeRefine * diam;
    }

    oofem::FloatArray random(3);
    random.at(3) = 0.;
    oofem::FloatArray newRandom(3);
    newRandom.at(3) = 0.;

    // ---- Edge placement (regular spacing along each edge) ----
    // Edges parallel to x — bottom (y = ymin) and top (y = ymax).
    if ( periodicityFlag.at(2) == 0 ) {
        int nTarget = std::max(1, (int) std::ceil(dx / ( edgeRefine * diam )));
        double newDiameter = dx / nTarget;
        if ( periodicityFlag.at(1) == 1 ) {
            nTarget--;
        }
        for ( int iy = 0; iy < 2; ++iy ) {
            random.at(2) = ( iy == 0 ) ? ymin : ymax;
            for ( int k = 0; k <= nTarget; ++k ) {
                if ( periodicityFlag.at(1) == 1 && k == 0 ) {
                    random.at(1) = xmin + 0.5 * newDiameter;
                } else {
                    random.at(1) = xmin + k * newDiameter;
                }
                int flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, 0.99 * newDiameter);
                if ( flag != 0 ) continue;
                if ( !grid->addVertex(random) ) continue;

                // Mirror across this edge (or periodic shift if x-axis periodic).
                for ( int xs = -1; xs <= 1; ++xs ) {
                    if ( xs == 0 ) continue;
                    if ( periodicityFlag.at(1) == 1 ) {
                        newRandom.at(1) = random.at(1) + xs * dx;
                    } else if ( xs == -1 ) {
                        newRandom.at(1) = 2. * xmin - random.at(1);
                    } else {
                        newRandom.at(1) = 2. * xmax - random.at(1);
                    }
                    newRandom.at(2) = random.at(2);
                    if ( grid->giveGridLocalizer()->checkNodesWithinBox(newRandom, 0.99 * newDiameter) == 0 ) {
                        grid->addVertex(newRandom);
                    }
                }
            }
        }
    }

    // Edges parallel to y — left (x = xmin) and right (x = xmax).
    if ( periodicityFlag.at(1) == 0 ) {
        int nTarget = std::max(1, (int) std::ceil(dy / ( edgeRefine * diam )));
        double newDiameter = dy / nTarget;
        if ( periodicityFlag.at(2) == 1 ) {
            nTarget--;
        }
        for ( int ix = 0; ix < 2; ++ix ) {
            random.at(1) = ( ix == 0 ) ? xmin : xmax;
            for ( int k = 0; k <= nTarget; ++k ) {
                if ( periodicityFlag.at(2) == 1 && k == 0 ) {
                    random.at(2) = ymin + 0.5 * newDiameter;
                } else {
                    random.at(2) = ymin + k * newDiameter;
                }
                int flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, 0.99 * newDiameter);
                if ( flag != 0 ) continue;
                if ( !grid->addVertex(random) ) continue;

                for ( int ys = -1; ys <= 1; ++ys ) {
                    if ( ys == 0 ) continue;
                    if ( periodicityFlag.at(2) == 1 ) {
                        newRandom.at(2) = random.at(2) + ys * dy;
                    } else if ( ys == -1 ) {
                        newRandom.at(2) = 2. * ymin - random.at(2);
                    } else {
                        newRandom.at(2) = 2. * ymax - random.at(2);
                    }
                    newRandom.at(1) = random.at(1);
                    if ( grid->giveGridLocalizer()->checkNodesWithinBox(newRandom, 0.99 * newDiameter) == 0 ) {
                        grid->addVertex(newRandom);
                    }
                }
            }
        }
    }

    printf("Finished edges. Current number of points are %d\n", grid->giveNumberOfVertices() );

    // ---- Interior random fill ----
    const double maxIter = grid->giveMaximumIterations();
    int mult = 0;
    for ( int i = 0; i < maxIter; ++i ) {
        random.at(1) = xmin + borderX + grid->ran1(& randomIntegerOne) * ( dx - 2. * borderX );
        random.at(2) = ymin + borderY + grid->ran1(& randomIntegerTwo) * ( dy - 2. * borderY );

        int flag = grid->giveGridLocalizer()->checkNodesWithinBox(random,
                       regionRefine * grid->giveDiameter(random) );

        // addVertex may silently reject a point inside a #@notch box; gating
        // `i = 0` on its return value prevents the loop from spinning.
        if ( flag == 0 && grid->addVertex(random) ) {
            i = 0;

            // Mirror / periodic-shift to all 8 neighbours in the 2D plane.
            for ( int xs = -1; xs <= 1; ++xs ) {
                for ( int ys = -1; ys <= 1; ++ys ) {
                    if ( xs == 0 && ys == 0 ) continue;
                    if ( periodicityFlag.at(1) == 1 ) {
                        newRandom.at(1) = random.at(1) + xs * dx;
                    } else if ( xs == -1 ) {
                        newRandom.at(1) = 2. * xmin - random.at(1);
                    } else if ( xs == 0 ) {
                        newRandom.at(1) = random.at(1);
                    } else {
                        newRandom.at(1) = 2. * xmax - random.at(1);
                    }
                    if ( periodicityFlag.at(2) == 1 ) {
                        newRandom.at(2) = random.at(2) + ys * dy;
                    } else if ( ys == -1 ) {
                        newRandom.at(2) = 2. * ymin - random.at(2);
                    } else if ( ys == 0 ) {
                        newRandom.at(2) = random.at(2);
                    } else {
                        newRandom.at(2) = 2. * ymax - random.at(2);
                    }
                    grid->addVertex(newRandom);
                }
            }
        }

        if ( i > mult * 10000 && i < ( mult + 1 ) * 10000 ) {
            std::cout << "Placed " << grid->giveNumberOfVertices() << " points in rect " << this->giveNumber() << "\n";
            std::cout << "Current iter " << i << " / max " << grid->giveMaximumIterations() << "\n";
            mult++;
        }
    }
    printf("Finished region. Current number of points are %d\n", grid->giveNumberOfVertices() );

    return 1;
}
