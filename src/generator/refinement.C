#include "refinement.h"
#include "curve.h"
#include "vertex.h"

#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include "octreegridlocalizer.h"
#endif

Refinement::Refinement(int n, Grid *aGrid) : GridComponent(n, aGrid) //, coordinates()
{
    this->number = n;
}

Refinement::~Refinement()
// Destructor.
{}
