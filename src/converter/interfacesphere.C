#include "interfacesphere.h"
#include "curve.h"
#include "vertex.h"


#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#endif

InterfaceSphere :: InterfaceSphere(int n,Grid* aGrid) : Inclusion(n,aGrid) //, coordinates()
{
  this->number = n;  
}

InterfaceSphere :: ~InterfaceSphere()
// Destructor.
{
}


void
InterfaceSphere :: initializeFrom(ConverterInputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{
    IR_GIVE_FIELD(ir, centre, _IFT_InterfaceSphere_centre);    
    IR_GIVE_FIELD(ir, radius, _IFT_InterfaceSphere_diameter); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, _IFT_InterfaceSphere_refine); // Macro
    itzThickness = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, itzThickness, _IFT_InterfaceSphere_itz); // Macro 
   return;
}


InterfaceSphere *InterfaceSphere :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    InterfaceSphere *sphere;
    
    sphere = new InterfaceSphere(number,grid);

    return sphere;
}
