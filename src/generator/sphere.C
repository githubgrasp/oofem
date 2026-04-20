#include "sphere.h"
#include "curve.h"
#include "vertex.h"
#include "generatorerror.h"

#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include "octreegridlocalizer.h"
 #include <iostream>
#endif

Sphere::Sphere(int n, Grid *aGrid) : Region(n, aGrid)
{
    this->number = n;
}

Sphere::~Sphere()
{}

int Sphere::generatePoints()
{
    //minimum points
    int randomIntegerOne = grid->giveRandomInteger() - 1;
    int randomIntegerTwo = grid->giveRandomInteger() - 2;
    int randomIntegerThree = grid->giveRandomInteger() - 3;

    double myPi = 3.14159265;

    oofem::FloatArray random(3);

    int flag;

    Vertex *vertex;

    double maxIter = grid->giveMaximumIterations();

    double randomTheta;
    double randomPhi;
    double randomRadius;

    //Place vertices

    //A) Centre
    random.at(1) = this->centre.at(1);
    random.at(2) = this->centre.at(2);
    random.at(3) = this->centre.at(3);

    grid->addVertex(random);

    printf("Points on sphere surface\n");

    for (int i = 0; i < maxIter; i++) {
        randomTheta = 2 * myPi * grid->ran1(& randomIntegerOne);
        randomPhi = acos(2 * grid->ran1(& randomIntegerTwo) - 1);

        random.at(1) = this->centre.at(1) + this->radius * cos(randomTheta) * sin(randomPhi);
        random.at(2) = this->centre.at(2) + this->radius * sin(randomTheta) * sin(randomPhi);
        random.at(3) = this->centre.at(3) + this->radius * cos(randomPhi);

        //Check if this is far enough from the others.

        flag = 0;
        flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, refinement * grid->diameter);



        if ( flag == 0 ) {
            i = 0;

            grid->addVertex(random);
        }
    }

    printf("Completed sphere surface\n");
    printf("number of vertices=%d\n", grid->giveNumberOfVertices() );

    printf("Points in sphere\n");

    for (int i = 0; i < maxIter; i++) {
        randomTheta = 2 * myPi * grid->ran1(& randomIntegerOne);
        randomPhi = acos(2 * grid->ran1(& randomIntegerTwo) - 1);
        randomRadius = ( this->radius - grid->diameter ) * pow(grid->ran1(& randomIntegerThree), 1. / 3.);

        random.at(1) = this->centre.at(1) + randomRadius * cos(randomTheta) * sin(randomPhi);
        random.at(2) = this->centre.at(2) + randomRadius * sin(randomTheta) * sin(randomPhi);
        random.at(3) = this->centre.at(3) + randomRadius * cos(randomPhi);

        //Check if this is far enough from the others.

        flag = 0;
        flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, grid->diameter);

        if ( flag == 0 ) {
            i = 0;

            grid->addVertex(random);

            //Mirror vertices with respect to the  sphere surface
            random.at(1) = this->centre.at(1) + ( 2. * this->radius - randomRadius ) * cos(randomTheta) * sin(randomPhi);
            random.at(2) = this->centre.at(2) + ( 2. * this->radius - randomRadius ) * sin(randomTheta) * sin(randomPhi);
            random.at(3) = this->centre.at(3) + ( 2. * this->radius - randomRadius ) * cos(randomPhi);

            grid->addVertex(random);
        }
    }


    return 1;
}



int Sphere::generatePeriodicPoints()
{
    return 1;
}



void Sphere::initializeFromTokens(std::istringstream &iss)
{
    refinement = 1.;

    bool gotCentre = false, gotRadius = false;
    std::string tok;
    while ( iss >> tok ) {
        if ( tok == "centre" ) {
            int n;
            iss >> n;
            centre.resize(n);
            for ( int i = 1; i <= n; ++i ) {
                iss >> centre.at(i);
            }
            gotCentre = true;
        } else if ( tok == "radius" ) {
            iss >> radius;
            gotRadius = true;
        } else if ( tok == "refine" ) {
            iss >> refinement;
        } else {
            generator::errorf("Sphere::initializeFromTokens: unknown keyword '%s'", tok.c_str());
        }
    }
    if ( !gotCentre || !gotRadius ) {
        generator::error("Sphere::initializeFromTokens: 'centre' and 'radius' are required");
    }
}
