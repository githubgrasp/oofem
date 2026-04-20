#include "curve.h"
#include "vertex.h"
#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
#endif

Curve::Curve(int n, Grid *aGrid) : GridComponent(n, aGrid) //, coordinates()
{
    this->number = n;
}


int
Curve::giveLocalVertex(int i)
// Returns the i-th coordinate of the receiver.
{
    if ( i > vertices.giveSize() ) {
        return 0.;
    }

    return vertices.at(i);
}

