#include "interfacesphere.h"
#include "curve.h"
#include "vertex.h"


#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
#endif

InterfaceSphere::InterfaceSphere(int n, Grid *aGrid) : Inclusion(n, aGrid) //, coordinates()
{
    this->number = n;
}




InterfaceSphere *InterfaceSphere::ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    InterfaceSphere *sphere;

    sphere = new InterfaceSphere(number, grid);

    return sphere;
}
