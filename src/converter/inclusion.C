#include "inclusion.h"
#include "curve.h"
#include "vertex.h"

#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
#endif

Inclusion::Inclusion(int n, Grid *aGrid) : GridComponent(n, aGrid)
{
    this->number = n;
}
