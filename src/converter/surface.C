#include "surface.h"
#include "curve.h"
#include "vertex.h"

#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#endif

Surface :: Surface(int n,Grid* aGrid) : GridComponent(n,aGrid) //, coordinates()
{
  this->number = n;  
}


Surface :: ~Surface()
// Destructor.
{
}


int 
Surface :: giveLocalCurve(int i)
// Returns the i-th coordinate of the receiver.
{
  if ( i > curves.giveSize() ) {
    return 0.;
  }
  
  return curves.at(i);
}

int Surface :: giveNumberOfLocalCurves()
{
  int number = this->curves.giveSize(); 
  return number;
}


IRResultType
Surface :: initializeFrom(InputRecord *ir)
// Gets from the source line from the data file all the data of the receiver.
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    //    int j, size;
    //    FloatArray vertices;
    // IntArray *dofIDArry;
    
    IR_GIVE_FIELD(ir, curves, IFT_Surface_curves, "curves"); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, IFT_Surface_refine, "refine"); // Macro
    return IRRT_OK;
}



Surface *Surface :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, LineSide,..).
{
    Surface *surface;
    
    surface = new Surface(number,grid);

    return surface;
}
