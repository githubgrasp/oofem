#include "region.h"
#include "curve.h"
#include "vertex.h"
#include "surface.h"

#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
#endif

Region::Region(int n, Grid *aGrid) : GridComponent(n, aGrid)
{
    this->number = n;
}
