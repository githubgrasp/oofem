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

InterfaceCylinder :: ~InterfaceCylinder()
// Destructor.
{
}


IRResultType
InterfaceCylinder :: initializeFrom(InputRecord *ir)
// Gets from the source line from the data file all the data of the receiver.
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, line, IFT_InterfaceCylinder_line, "line");
    IR_GIVE_FIELD(ir, radius, IFT_InterfaceCylinder_radius, "radius"); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, IFT_InterfaceCylinder_refine, "refine"); // Macro
    itzThickness = refinement*grid->giveDiameter();
    IR_GIVE_OPTIONAL_FIELD(ir, itzThickness, IFT_InterfaceCylinder_itz, "itz"); // Macro
    return IRRT_OK;

}


InterfaceCylinder *InterfaceCylinder :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    InterfaceCylinder *cylinder;

    cylinder = new InterfaceCylinder(number,grid);

    return cylinder;
}
