#include "interfacecylinder.h"
#include "curve.h"
#include "vertex.h"


#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#endif

InterfaceCylinder :: InterfaceCylinder(int n,Grid* aGrid) : Inclusion(n,aGrid) //, coordinates()
{
  this->number = n;
}


void
InterfaceCylinder :: initializeFrom(ConverterInputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{
    IR_GIVE_FIELD(ir, line, _IFT_InterfaceCylinder_line);
    IR_GIVE_FIELD(ir, radius, _IFT_InterfaceCylinder_radius); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, _IFT_InterfaceCylinder_refine); // Macro
    itzThickness = refinement*grid->giveDiameter();
    IR_GIVE_OPTIONAL_FIELD(ir, itzThickness, _IFT_InterfaceCylinder_itz); // Macro
    return;
}


InterfaceCylinder *InterfaceCylinder :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    InterfaceCylinder *cylinder;

    cylinder = new InterfaceCylinder(number,grid);

    return cylinder;
}
