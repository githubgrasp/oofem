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

void
InterfaceSurface :: initializeFrom(GeneratorInputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{
     
    IR_GIVE_FIELD(ir, curves, _IFT_InterfaceSurface_curves); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, _IFT_InterfaceSurface_refine); // Macro
    return;

}

InterfaceSurface *InterfaceSurface :: ofType()
{
    InterfaceSurface *surface;
    
    surface = new InterfaceSurface(number,grid);

    return surface;
}
