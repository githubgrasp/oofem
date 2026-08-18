#include "interfaceplane.h"
#include "curve.h"
#include "vertex.h"
#include "generatorerror.h"


#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include "octreegridlocalizer.h"
#endif

InterfacePlane::InterfacePlane(int n, Grid *aGrid) : Inclusion(n, aGrid) //, coordinates()
{
}

InterfacePlane::~InterfacePlane()
// Destructor.
{}


int InterfacePlane::generatePoints()
{
    //minimum points
    int minimumNumberOfPoints = 8;

    oofem::FloatArray boundaries;
    grid->defineBoundaries(boundaries);

    int randomIntegerOne = grid->giveRandomInteger() - 1;
    int randomIntegerTwo = grid->giveRandomInteger() - 2;

    double myPi = 3.14159265;


    //   int randomIntegerThree= grid->giveRandomInteger()-3;

    oofem::FloatArray random(3);
    int flag;

    double boundaryFactor = this->refinement;

    Vertex *vertex;

    double distance, newDistance, shift;
    double x, y, z;
    double maxIter = grid->giveMaximumIterations();
    oofem::FloatArray mirroredRandom(3);

    double randomAngleOne;
    double randomAngleTwo;
    double randomRadius;

    distance =  boundaryFactor * grid->diameter;

    double pointEstimate = myPi * diameter / distance;

    int numberOfIntervalsOne, numberOfIntervalsTwo, newIntervalsTwo;

    if ( pointEstimate > minimumNumberOfPoints ) {
        numberOfIntervalsTwo = ceil(pointEstimate);
    } else   {
        numberOfIntervalsTwo = minimumNumberOfPoints;
        distance = myPi * diameter / minimumNumberOfPoints;
    }

    double length = sqrt(pow(this->line.at(1) - this->line.at(4), 2.) );

    numberOfIntervalsOne = ceil(numberOfIntervalsTwo * length / ( myPi * diameter ) );


    //This implementation is hardcoding the orientation of the plane

    //Two loops
    for (int i = 0; i <= numberOfIntervalsOne; i++) {
        random.at(1) = line.at(1) + i * length / numberOfIntervalsOne;
        //     printf("random.at(1) = %e\n", random.at(1));

        //Place Midpoint
        random.at(2) = this->line.at(2);
        random.at(3) = this->line.at(3);

        grid->addVertex(random);

        for (int k = 0; k < numberOfIntervalsTwo; k++) {
            //Generate random point

            randomAngleTwo = k * 2. * myPi / numberOfIntervalsTwo;

            for (int m = 0; m < 2; m++) {
                random.at(2) = this->line.at(2) + ( 0.5 * this->diameter + ( 1. - m ) * this->itzThickness ) * cos(randomAngleTwo);
                random.at(3) = this->line.at(3) + ( 0.5 * this->diameter + ( 1. - m ) * this->itzThickness ) * sin(randomAngleTwo);

                //Check if this is far enough from the others
                flag = 0;

                grid->addVertex(random);


                if ( i == 0 || i == numberOfIntervalsOne ) {
                    printf("VertexNumbers for surface elements of cylindrical inclusion = %d\n", grid->giveNumberOfVertices() );
                }
            }//end of two layers
        }//end of iterations
    }

    printf("Completed inclusion loops\n");

    return 1;
}


void InterfacePlane::initializeFromTokens(std::istringstream &iss)
{
    refinement = 1.;

    bool gotLine = false, gotDiameter = false;
    bool gotItz = false;
    std::string tok;
    while ( iss >> tok ) {
        if ( tok == "line" ) {
            int n;
            iss >> n;
            line.resize(n);
            for ( int i = 1; i <= n; ++i ) {
                iss >> line.at(i);
            }
            gotLine = true;
        } else if ( tok == "diameter" ) {
            iss >> diameter;
            gotDiameter = true;
        } else if ( tok == "refine" ) {
            iss >> refinement;
        } else if ( tok == "itz" ) {
            iss >> itzThickness;
            gotItz = true;
        } else {
            generator::errorf("InterfacePlane::initializeFromTokens: unknown keyword '%s'", tok.c_str());
        }
    }
    if ( !gotLine || !gotDiameter ) {
        generator::error("InterfacePlane::initializeFromTokens: 'line' and 'diameter' are required");
    }
    if ( !gotItz ) {
        itzThickness = refinement * diameter;
    }
}

