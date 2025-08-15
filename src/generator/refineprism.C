#include "refineprism.h"
#include "curve.h"
#include "vertex.h"
#include "record.h"


#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include "octreegridlocalizer.h"
#endif


RefinePrism::RefinePrism(int n, Grid *aGrid) : Refinement(n, aGrid) //, coordinates()
{
    this->number = n;
}

RefinePrism::~RefinePrism()
{}

void
RefinePrism::initializeFrom(GeneratorInputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{
    IR_GIVE_FIELD(ir, this->box, _IFT_RefinePrism_box);

    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, _IFT_RefinePrism_refine); // Macro

    return;
}

double
RefinePrism::giveDiameter(oofem::FloatArray &coord) {
    //check of coord is in prism
    double diam = 0.;

    if ( coord.at(1) >= box.at(1) && coord.at(2) >= box.at(2) && coord.at(3) >= box.at(3) && coord.at(1) <= box.at(4) && coord.at(2) <= box.at(5) && coord.at(3) <= box.at(6) ) {//inside the box
        diam = this->refinement * grid->diameter;
    } else   {
        diam = grid->diameter;
    }

    return diam;
}

RefinePrism *
RefinePrism::ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    RefinePrism *prism;

    prism = new RefinePrism(number, grid);

    return prism;
}
