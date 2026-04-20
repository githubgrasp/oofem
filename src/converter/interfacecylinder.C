#include "interfacecylinder.h"
#include "curve.h"
#include "vertex.h"


#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
#endif

InterfaceCylinder::InterfaceCylinder(int n, Grid *aGrid) : Inclusion(n, aGrid)
{
    this->number = n;
}




InterfaceCylinder *InterfaceCylinder::ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    InterfaceCylinder *cylinder;

    cylinder = new InterfaceCylinder(number, grid);

    return cylinder;
}
