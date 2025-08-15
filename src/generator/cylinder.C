#include "cylinder.h"
#include "curve.h"
#include "vertex.h"
#include "record.h"


#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include "octreegridlocalizer.h"
#endif

Cylinder::Cylinder(int n, Grid *aGrid) : Region(n, aGrid) //, coordinates()
{
    this->number = n;
}

Cylinder::~Cylinder()
// Destructor.
{}

int Cylinder::generatePoints()
{
    int randomIntegerOne = grid->giveRandomInteger() - 1;
    int randomIntegerTwo = grid->giveRandomInteger() - 2;
    int randomIntegerThree = grid->giveRandomInteger() - 3;

    double myPi = 3.14159265;


    oofem::FloatArray random(3);
    int flag;

    Vertex *vertex;

    double distance;
    double maxIter = grid->giveMaximumIterations();

    double randomTheta;
    double randomRadius;
    double randomXCoord;
    double randomYCoord;
    double randomZCoord;


    printf("Generate centres of end circles\n");
    printf("First point\n");
    random.at(1) = this->line.at(1);
    random.at(2) = this->line.at(2);
    random.at(3) = this->line.at(3);

    flag = 0;
    flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, grid->giveDiameter(random) );

    if ( flag == 0 ) {
        grid->addVertex(random);
    }
    printf("Second point\n");
    random.at(1) = this->line.at(4);
    random.at(2) = this->line.at(5);
    random.at(3) = this->line.at(6);


    flag = 0;
    flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, grid->giveDiameter(random) );

    if ( flag == 0 ) {
        grid->addVertex(random);
    }

    printf("Placed vertices at the end cylinders in the centre of circles. Vertex number = %d\n", grid->giveNumberOfVertices() );

    printf("Start to generate points on circles\n");

    printf("Start with circle one\n");

    double spacing = refinement * grid->diameter;
    int intervals = ceil(2 * myPi * this->radius / spacing);
    random.at(1) = this->line.at(1);
    for (int i = 0; i < intervals; i++) {
        randomTheta = 2. * myPi / intervals * i;
        random.at(2) = this->line.at(2) + this->radius * cos(randomTheta);
        random.at(3) = this->line.at(3) + this->radius * sin(randomTheta);

        flag = 0;
        flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, 0.9 * spacing);

        if ( flag == 0 ) {
            grid->addVertex(random);
        }
    }

    printf("Placed vertices for circle one. Vertex number = %d\n", grid->giveNumberOfVertices() );


    printf("Start with circle two\n");
    spacing = refinement * grid->diameter;
    intervals = ceil(2 * myPi * this->radius / spacing);
    random.at(1) = this->line.at(4);
    for (int i = 0; i < intervals; i++) {
        randomTheta = 2. * myPi / intervals * i;
        random.at(2) = this->line.at(5) + this->radius * cos(randomTheta);
        random.at(3) = this->line.at(6) + this->radius * sin(randomTheta);

        flag = 0;
        flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, 0.9 * spacing);

        if ( flag == 0 ) {
            grid->addVertex(random);
        }
    }

    printf("Placed vertices on circle two. Vertex number = %d\n", grid->giveNumberOfVertices() );

    printf("Start with disk one\n");


    spacing = 0. * refinement * grid->diameter;
    random.at(1) = this->line.at(1);
    for (int i = 0; i < maxIter; i++) {
        randomTheta = 2 * myPi * grid->ran1(& randomIntegerOne);
        randomRadius = ( this->radius - spacing ) * grid->ran1(& randomIntegerTwo);
        random.at(2) = this->line.at(2) + randomRadius * cos(randomTheta);
        random.at(3) = this->line.at(3) + randomRadius * sin(randomTheta);

        //Check if this is far enough from the others.

        flag = 0;
        flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, refinement * grid->giveDiameter(random) );

        if ( flag == 0 ) {
            i = 0;

            grid->addVertex(random);

            //Mirror point with respect to circle
            random.at(2) = this->line.at(2) + ( 2. * this->radius - randomRadius ) * cos(randomTheta);
            random.at(3) = this->line.at(3) + ( 2. * this->radius - randomRadius ) * sin(randomTheta);
            grid->addVertex(random);
        }
    }

    printf("Placed vertices on disk one. Vertex number = %d\n", grid->giveNumberOfVertices() );

    printf("Start with disk two\n");


    spacing = 0.0 * refinement * grid->diameter;
    random.at(1) = this->line.at(4);
    for (int i = 0; i < maxIter; i++) {
        randomTheta = 2 * myPi * grid->ran1(& randomIntegerOne);
        randomRadius = ( this->radius - spacing ) * grid->ran1(& randomIntegerTwo);
        random.at(2) = this->line.at(5) + randomRadius * cos(randomTheta);
        random.at(3) = this->line.at(6) + randomRadius * sin(randomTheta);

        //Check if this is far enough from the others.

        flag = 0;
        flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, refinement * grid->giveDiameter(random) );

        if ( flag == 0 ) {
            i = 0;

            grid->addVertex(random);

            //Mirror point with respect to circle
            random.at(2) = this->line.at(5) + ( 2. * this->radius - randomRadius ) * cos(randomTheta);
            random.at(3) = this->line.at(6) + ( 2. * this->radius - randomRadius ) * sin(randomTheta);
            grid->addVertex(random);
        }
    }
    printf("Placed vertices on disk two. Vertex number = %d\n", grid->giveNumberOfVertices() );

    printf("now the surface of the cylinder\n");
    spacing = 0. * refinement * grid->diameter;
    for (int i = 0; i < maxIter; i++) {
        randomXCoord = this->line.at(1) + spacing + ( line.at(4) - line.at(1) - 2. * spacing ) *
                       grid->ran1(& randomIntegerOne);
        random.at(1) = randomXCoord;
        randomTheta = 2 * myPi * grid->ran1(& randomIntegerOne);
        random.at(2) = this->line.at(2) + this->radius * cos(randomTheta);
        random.at(3) = this->line.at(3) + this->radius * sin(randomTheta);

        //Check if this is far enough from the others.

        flag = 0;
        flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, refinement * grid->giveDiameter(random) );

        if ( flag == 0 ) {
            i = 0;

            grid->addVertex(random);

            //Mirror point with respect to first end circles
            random.at(1) = 2. * line.at(1) - randomXCoord;

            grid->addVertex(random);

            //Mirror point with respect to second end circle
            random.at(1) = 2. * line.at(4) - randomXCoord;

            grid->addVertex(random);
        }
    }
    printf("Surface is completed. Vertex number = %d\n", grid->giveNumberOfVertices() );

    printf("now the inside of the cylinder\n");

    spacing = refinement * grid->diameter;

    for (int i = 0; i < maxIter; i++) {
        randomXCoord = this->line.at(1) + spacing +
                       ( line.at(4) - line.at(1) - 2. * spacing ) * grid->ran1(& randomIntegerOne);

        randomRadius = ( this->radius - spacing ) * grid->ran1(& randomIntegerTwo);
        randomTheta = 2 * myPi * grid->ran1(& randomIntegerThree);

        randomYCoord = this->line.at(2) + randomRadius * cos(randomTheta);
        randomZCoord = this->line.at(3) + randomRadius * sin(randomTheta);

        random.at(1) = randomXCoord;
        random.at(2) = randomYCoord;
        random.at(3) = randomZCoord;

        //Check if this is far enough from the others.

        flag = 0;
        flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, grid->giveDiameter(random) );

        if ( flag == 0 ) {
            i = 0;

            grid->addVertex(random);

            //Mirror original point with respect to cylinder surface
            random.at(1) = randomXCoord;
            random.at(2) = this->line.at(2) + ( 2. * this->radius - randomRadius ) * cos(randomTheta);
            random.at(3) = this->line.at(3) + ( 2. * this->radius - randomRadius ) * sin(randomTheta);

            grid->addVertex(random);

            //Mirror original point with respect to first end circles
            random.at(1) = 2. * line.at(1) - randomXCoord;
            random.at(2) = randomYCoord;
            random.at(3) = randomZCoord;

            grid->addVertex(random);

            //Mirror mirrored point (first end circle) with respect to cylinder surface
            random.at(1) = 2. * line.at(1) - randomXCoord;
            random.at(2) = this->line.at(2) + ( 2. * this->radius - randomRadius ) * cos(randomTheta);
            random.at(3) = this->line.at(3) + ( 2. * this->radius - randomRadius ) * sin(randomTheta);
            grid->addVertex(random);

            //Mirror original point with respect to second end circle
            random.at(1) = 2. * line.at(4) - randomXCoord;
            random.at(2) = randomYCoord;
            random.at(3) = randomZCoord;

            grid->addVertex(random);

            //Mirror mirrored point (second end circle) with respect to cylinder surface
            random.at(1) = 2. * line.at(4) - randomXCoord;
            random.at(2) = this->line.at(2) + ( 2. * this->radius - randomRadius ) * cos(randomTheta);
            random.at(3) = this->line.at(3) + ( 2. * this->radius - randomRadius ) * sin(randomTheta);
            grid->addVertex(random);
        }
    }

    printf("Inside is completed. Vertex number = %d\n", grid->giveNumberOfVertices() );

    printf("Completed cylinder point generation.\n");


    return 1;
}

void
Cylinder::initializeFrom(GeneratorInputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{
    // IntArray *dofIDArry;
    IR_GIVE_FIELD(ir, line, _IFT_Cylinder_line);
    IR_GIVE_FIELD(ir, radius, _IFT_Cylinder_radius); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, _IFT_Cylinder_refine); // Macro
    return;
}



Cylinder *Cylinder::ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    Cylinder *cylinder;

    cylinder = new Cylinder(number, grid);

    return cylinder;
}
