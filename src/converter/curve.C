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


IRResultType
Curve :: initializeFrom(InputRecord *ir)
// Gets from the source line from the data file all the data of the receiver.
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, vertices, IFT_Curve_vertices, "vertices"); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, IFT_Curve_refine, "refine"); // Macro
    return IRRT_OK;
}



Curve *Curve :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, LineSide,..).
{
    Curve *curve;
    
    curve = new Curve(number,grid);

    return curve;
}
