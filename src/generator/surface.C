#include "generatorlistutils.h"
#include "surface.h"
#include "curve.h"
#include "vertex.h"
#include <iostream>

#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include "octreegridlocalizer.h"
#endif

Surface::Surface(int n, Grid *aGrid) : GridComponent(n, aGrid)   //, coordinates()
{
    this->number = n;
    this->boundaryFlag = 0;
    boundaryShift.zero();
}

Surface::~Surface()
// Destructor.
{}

int
Surface::giveLocalCurve(int i)
// Returns the i-th coordinate of the receiver.
{
    if ( i > curves.giveSize() ) {
        return 0.;
    }

    return curves.at(i);
}

int Surface::giveNumberOfLocalCurves()
{
    int number = this->curves.giveSize();
    return number;
}


//generatePoints Works only for cartesian coordinate system: Check normal and then generate random points. Depending on periodicity flag, the points are mirrored or periodically shifted.

int
Surface::generatePoints()
{
    oofem::IntArray periodicityFlag(3);
    grid->givePeriodicityFlag(periodicityFlag);

    int randomFlag = grid->giveRandomFlag();

    oofem::FloatArray boundaries;
    this->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    //Get global dimensions
    oofem::FloatArray globalBoundaries;
    oofem::FloatArray globalDimension(3);
    grid->defineBoundaries(globalBoundaries);
    globalDimension.at(1) = globalBoundaries.at(2) - globalBoundaries.at(1);
    globalDimension.at(2) = globalBoundaries.at(4) - globalBoundaries.at(3);
    globalDimension.at(3) = globalBoundaries.at(6) - globalBoundaries.at(5);

    oofem::FloatArray normal(3);
    this->giveNormal(normal);


    int randomIntegerTwo = grid->giveRandomInteger() - 2;
    int randomIntegerThree = grid->giveRandomInteger() - 3;

    oofem::FloatArray random(3), newRandom(3);

    int flag = 0;

    double shift = 0.;

    double boundaryFactor = this->refinement;

    double maxIter = grid->giveMaximumIterations();


    double border = grid->diameter;

    for ( int i = 0; i < maxIter; i++ ) {
        if ( normal.at(1) == 1 && normal.at(2) == 0 && normal.at(3) == 0 ) {//y-z coordinate system
            random.at(1) = boundaries.at(1);
            random.at(2) = boundaries.at(3) + border + grid->ran1(& randomIntegerTwo) * ( boundaries.at(4) - boundaries.at(3) - 2. * border );
            random.at(3) = boundaries.at(5) + border + grid->ran1(& randomIntegerThree) * ( boundaries.at(6) - boundaries.at(5) - 2. * border );
        } else if ( normal.at(1) == 0 && normal.at(2) == 1 && normal.at(3) == 0 )      {//x-z coordinate system
            random.at(1) = boundaries.at(1) + border + grid->ran1(& randomIntegerTwo) * ( boundaries.at(2) - boundaries.at(1) - 2. * border );
            random.at(2) = boundaries.at(3);
            random.at(3) = boundaries.at(5) + border + grid->ran1(& randomIntegerThree) * ( boundaries.at(6) - boundaries.at(5) - 2. * border );
        } else if ( normal.at(1) == 0 && normal.at(2) == 0 && normal.at(3) == 1 )      {//x-y coordinate system
            random.at(1) = boundaries.at(1) + border + grid->ran1(& randomIntegerTwo) * ( boundaries.at(2) - boundaries.at(1) - 2. * border );
            random.at(2) = boundaries.at(3) + border + grid->ran1(& randomIntegerTwo) * ( boundaries.at(4) - boundaries.at(3) - 2. * border );
            random.at(3) = boundaries.at(5);
        } else   {
            printf("Will not work for this normal. Must be aligned with cartesian coordinate system.\n");
        }

        //Check if this is far enough from the others
        flag = 0;
        flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, boundaryFactor * grid->giveDiameter(random) );

        if ( flag == 0 ) {
            grid->addVertex(random);

            i = 0;

            mirrorShift(random, normal, specimenDimension, boundaries, periodicityFlag);


            //Option random=2, for which periodic points are shifted across the specimen.
            if ( randomFlag == 2 ) {
                if ( normal.at(1) == 1 && normal.at(2) == 0 && normal.at(3) == 0 ) {//x-direction
                    if ( boundaries.at(1) == globalBoundaries.at(1) ) {
                        shift = 1;
                    } else if ( boundaries.at(1) == globalBoundaries.at(2) )      {
                        shift = -1;
                    } else   {
                        printf("Surface not part of global boundary\n");
                        exit(1);
                    }

                    newRandom.at(1) = random.at(1) + shift * globalDimension.at(1);
                    newRandom.at(2) = random.at(2);
                    newRandom.at(3) = random.at(3);


                    grid->addVertex(newRandom);

                    mirrorShift(newRandom, normal, specimenDimension, boundaries, periodicityFlag);
                } else if ( normal.at(1) == 0 && normal.at(2) == 1 && normal.at(3) == 0 )     {//y-direction
                    if ( boundaries.at(3) == globalBoundaries.at(3) ) {
                        shift = 1;
                    } else if ( boundaries.at(3) == globalBoundaries.at(4) )      {
                        shift = -1;
                    } else   {
                        printf("Surface not part of global boundary\n");
                        exit(1);
                    }

                    //Shift in y-direction
                    newRandom.at(1) = random.at(1);
                    newRandom.at(2) = random.at(2) + shift * globalDimension.at(2);
                    newRandom.at(3) = random.at(3);

                    grid->addVertex(newRandom);

                    mirrorShift(newRandom, normal, specimenDimension, boundaries, periodicityFlag);
                } else if ( normal.at(1) == 0 && normal.at(2) == 0 && normal.at(3) == 1 )     {//z-direction
                    if ( boundaries.at(5) == globalBoundaries.at(5) ) {
                        shift = 1;
                    } else if ( boundaries.at(5) == globalBoundaries.at(6) )      {
                        shift = -1;
                    } else   {
                        printf("Surface not part of global boundary\n");
                        exit(1);
                    }

                    //Shift in z-direction
                    newRandom.at(1) = random.at(1);
                    newRandom.at(2) = random.at(2);
                    newRandom.at(3) = random.at(3) + shift * globalDimension.at(3);

                    grid->addVertex(newRandom);


                    mirrorShift(newRandom, normal, specimenDimension, boundaries, periodicityFlag);
                } else   {
                    printf("Error: Unknown direction.\n");
                    exit(1);
                }
            }
        }
    }

    return 1;
}


void Surface::mirrorShift(oofem::FloatArray &random, oofem::FloatArray &normal, oofem::FloatArray &specimenDimension, oofem::FloatArray &boundaries, oofem::IntArray &periodicityFlag)
{
    //Mirror (or periodic shift) with respect to two of the three axis.

    oofem::FloatArray newRandom(3);

    if ( normal.at(1) == 1 && normal.at(2) == 0 && normal.at(3) == 0 ) {//y-z coordinate system
        for ( int y = -1; y < 2; y++ ) {
            for ( int z = -1; z < 2; z++ ) {
                if ( !( z == 0 && y == 0 ) ) {
                    newRandom.at(1) = random.at(1);
                    //y-direction
                    if ( periodicityFlag.at(2) == 1 ) {
                        newRandom.at(2) = random.at(2) + y * specimenDimension.at(2);
                    } else   {
                        if ( y == -1 ) {
                            newRandom.at(2) = random.at(2) - 2. * ( random.at(2) - boundaries.at(3) );
                        } else if ( y == 0 )      {
                            newRandom.at(2) = random.at(2);
                        } else if ( y == 1 )      {
                            newRandom.at(2) = random.at(2) - 2. * ( random.at(2) - boundaries.at(4) );
                        }
                    }
                    // z-direction
                    if ( periodicityFlag.at(3) == 1 ) {
                        newRandom.at(3) = random.at(3) + z * specimenDimension.at(3);
                    } else   {
                        if ( z == -1 ) {
                            newRandom.at(3) = random.at(3) - 2. * ( random.at(3) - boundaries.at(5) );
                        } else if ( z == 0 )      {
                            newRandom.at(3) = random.at(3);
                        } else if ( z == 1 )      {
                            newRandom.at(3) = random.at(3) - 2. * ( random.at(3) - boundaries.at(6) );
                        }
                    }

                    grid->addVertex(newRandom);
                }
            }
        }
    } else if ( normal.at(1) == 0 && normal.at(2) == 1 && normal.at(3) == 0 )      {//x-z coordinate system
        for ( int x = -1; x < 2; x++ ) {
            for ( int z = -1; z < 2; z++ ) {
                if ( !( x == 0 && z == 0 ) ) {
                    //x-direction
                    if ( periodicityFlag.at(1) == 1 ) {
                        newRandom.at(1) = random.at(1) + x * specimenDimension.at(1);
                    } else   {
                        if ( x == -1 ) {
                            newRandom.at(1) = random.at(1) - 2. * ( random.at(1) - boundaries.at(1) );
                        } else if ( x == 0 )      {
                            newRandom.at(1) = random.at(1);
                        } else if ( x == 1 )      {
                            newRandom.at(1) = random.at(1) - 2. * ( random.at(1) - boundaries.at(2) );
                        }
                    }
                    //y-direction
                    newRandom.at(2) = random.at(2);
                    //z-direction
                    if ( periodicityFlag.at(3) == 1 ) {
                        newRandom.at(3) = random.at(3) + z * specimenDimension.at(3);
                    } else   {
                        if ( z == -1 ) {
                            newRandom.at(3) = random.at(3) - 2. * ( random.at(3) - boundaries.at(5) );
                        } else if ( z == 0 )      {
                            newRandom.at(3) = random.at(3);
                        } else if ( z == 1 )      {
                            newRandom.at(3) = random.at(3) - 2. * ( random.at(3) - boundaries.at(6) );
                        }
                    }

                    grid->addVertex(newRandom);
                }
            }
        }
    }
    if ( normal.at(1) == 0 && normal.at(2) == 0 && normal.at(3) == 1 ) {//x-y coordinate system
        for ( int x = -1; x < 2; x++ ) {
            for ( int y = -1; y < 2; y++ ) {
                if ( !( x == 0 && y == 0 ) ) {
                    //x-direction
                    if ( periodicityFlag.at(1) == 1 ) {
                        newRandom.at(1) = random.at(1) + x * specimenDimension.at(1);
                    } else   {
                        if ( x == -1 ) {
                            newRandom.at(1) = random.at(1) - 2. * ( random.at(1) - boundaries.at(1) );
                        } else if ( x == 0 )      {
                            newRandom.at(1) = random.at(1);
                        } else if ( x == 1 )      {
                            newRandom.at(1) = random.at(1) - 2. * ( random.at(1) - boundaries.at(2) );
                        }
                    }
                    //y-direction
                    if ( periodicityFlag.at(2) == 1 ) {
                        newRandom.at(2) = random.at(2) + y * specimenDimension.at(2);
                    } else   {
                        if ( y == -1 ) {
                            newRandom.at(2) = random.at(2) - 2. * ( random.at(2) - boundaries.at(3) );
                        } else if ( y == 0 )      {
                            newRandom.at(2) = random.at(2);
                        } else if ( y == 1 )      {
                            newRandom.at(2) = random.at(2) - 2. * ( random.at(2) - boundaries.at(4) );
                        }
                    }
                    //z-direction
                    newRandom.at(3) = random.at(3);

                    grid->addVertex(newRandom);
                }
            }
        }
    }
    return;
}

void Surface::defineBoundaries(oofem::FloatArray &boundaries)
//Determine the boundaries of the domain
{
    Curve *curve;
    oofem::IntArray localVertices;
    oofem::IntArray vertexFlag(2);

    int localCurve, localVertex;

    double x, y, z;

    boundaries.resize(6);
    boundaries.at(1) = 1.e16;
    boundaries.at(2) = -1.e16;
    boundaries.at(3) = 1.e16;
    boundaries.at(4) = -1.e16;
    boundaries.at(5) = 1.e16;
    boundaries.at(6) = -1.e16;
    for ( int i = 0; i < this->giveNumberOfLocalCurves(); i++ ) {
        localCurve = this->giveLocalCurve(i + 1);
        curve = grid->giveCurve(localCurve);

        for ( int k = 0; k < 2; k++ ) {
            localVertex = curve->giveLocalVertex(k + 1);
            x = ( grid->giveInputVertex(localVertex) )->giveCoordinate(1);
            y = ( grid->giveInputVertex(localVertex) )->giveCoordinate(2);
            z = ( grid->giveInputVertex(localVertex) )->giveCoordinate(3);

            //assign min and max values

            if ( x < boundaries.at(1) ) {
                boundaries.at(1) = x;
            }
            if ( x > boundaries.at(2) ) {
                boundaries.at(2) = x;
            }

            if ( y < boundaries.at(3) ) {
                boundaries.at(3) = y;
            }
            if ( y > boundaries.at(4) ) {
                boundaries.at(4) = y;
            }

            if ( z < boundaries.at(5) ) {
                boundaries.at(5) = z;
            }
            if ( z > boundaries.at(6) ) {
                boundaries.at(6) = z;
            }
        }
    }

    return;
}



void Surface::initializeFrom(GeneratorInputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{
    IR_GIVE_FIELD(ir, curves, _IFT_Surface_curves); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, _IFT_Surface_refine); // Macro

    normal.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, normal, _IFT_Surface_normal); //Macro

    IR_GIVE_OPTIONAL_FIELD(ir, boundaryFlag, _IFT_Surface_boundaryflag); // Macro
    this->boundaryShift.resize(3);
    if ( this->boundaryFlag == 1 ) {
        IR_GIVE_FIELD(ir, boundaryShift, _IFT_Surface_boundaryshift); // Macro
    }

    return;
}
