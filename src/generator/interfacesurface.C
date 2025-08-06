#include "interfacesurface.h"
#include "curve.h"
#include "vertex.h"
#include "record.h"


#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#include "octreegridlocalizer.h"
#endif

InterfaceSurface :: InterfaceSurface(int n,Grid* aGrid) : Inclusion(n,aGrid) //, coordinates()
{
  this->number = n;  
}

InterfaceSurface :: ~InterfaceSurface()
// Destructor.
{
}

int InterfaceSurface :: generatePoints()
{
  return 1;
}

IRResultType
InterfaceSurface :: initializeFrom(InputRecord *ir)
// Gets from the source line from the data file all the data of the receiver.
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    
    IR_GIVE_FIELD(ir, curves, IFT_InterfaceSurface_curves, "curves"); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, IFT_InterfaceSurface_refine, "refine"); // Macro
    return IRRT_OK;

}

InterfaceSurface *InterfaceSurface :: ofType()
{
    InterfaceSurface *surface;
    
    surface = new InterfaceSurface(number,grid);

    return surface;
}
