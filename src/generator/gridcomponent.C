#include "gridcomponent.h"


#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#endif

GridComponent :: GridComponent(int n, Grid* aGrid) //, coordinates()
{
  grid = aGrid;
  number = n;
}

