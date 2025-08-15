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


IRResultType
InterfaceSphere :: initializeFrom(InputRecord *ir)
// Gets from the source line from the data file all the data of the receiver.
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // int j, size;
    // FloatArray vertices;
    // IntArray *dofIDArry;
    IR_GIVE_FIELD(ir, centre, IFT_InterfaceSphere_centre, "centre");    
    IR_GIVE_FIELD(ir, radius, IFT_InterfaceSphere_diameter, "radius"); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, IFT_InterfaceSphere_refine, "refine"); // Macro
    itzThickness = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, itzThickness, IFT_InterfaceSphere_itz, "itz"); // Macro 
   return IRRT_OK;
   
   //Print command for double
   //   printf("diameter = %e\n", diameter);
   //Print command for integer
   //   printf("number = %d\n", number);
   //Print command for array
   //   centre.printYourself();
}


InterfaceSphere *InterfaceSphere :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    InterfaceSphere *sphere;
    
    sphere = new InterfaceSphere(number,grid);

    return sphere;
}
