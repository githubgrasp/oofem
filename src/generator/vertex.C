#include "vertex.h"
#include "gridcomponent.h"

#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
#endif

Vertex::Vertex(int n, Grid *aGrid) : GridComponent(n, aGrid)  //, coordinates()
{
    this->number = n;
}




Vertex::~Vertex()
// Destructor.
{}


double Vertex::giveCoordinate(int i)
// Returns the i-th coordinate of the receiver.
{
    if ( i > coordinates.giveSize() ) {
        return 0.;
    }

    return coordinates.at(i);
}


void
Vertex::initializeFrom(GeneratorInputRecord &ir)
{
    int j, size;
    oofem::FloatArray triplets;

    IR_GIVE_FIELD(ir, coordinates, _IFT_Vertex_coords); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, _IFT_Vertex_refine); // Macro
    radius = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, radius, _IFT_Vertex_radius); // Macro


    return;
}
