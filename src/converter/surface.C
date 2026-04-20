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





Surface *Surface :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, LineSide,..).
{
    Surface *surface;
    
    surface = new Surface(number,grid);

    return surface;
}
