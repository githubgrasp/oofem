#include "curve.h"
#include "vertex.h"
#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#endif

Curve :: Curve(int n,Grid* aGrid) : GridComponent(n,aGrid) //, coordinates()
{
  this->number =n;
}

Curve :: ~Curve()
// Destructor.
{
}


int 
Curve :: giveLocalVertex(int i)
// Returns the i-th coordinate of the receiver.
{
  if ( i > vertices.giveSize() ) {
    return 0.;
  }
  
  return vertices.at(i);
}

void
Curve :: initializeFrom(ConverterInputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{
    IR_GIVE_FIELD(ir, vertices, _IFT_Curve_vertices); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, _IFT_Curve_refine); // Macro
    return;
}
