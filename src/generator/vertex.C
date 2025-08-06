#include "vertex.h"
#include "gridcomponent.h"
#include "alist.h"

#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#endif

Vertex :: Vertex(int n, Grid* aGrid) : GridComponent(n,aGrid) //, coordinates()
{
  this->number = n;  
  }




Vertex :: ~Vertex()
// Destructor.
{
}


double Vertex :: giveCoordinate(int i)
// Returns the i-th coordinate of the receiver.
{
    if ( i > coordinates.giveSize() ) {
        return 0.;
    }

    return coordinates.at(i);
}


IRResultType
Vertex :: initializeFrom(InputRecord *ir)
// Gets from the source line from the data file all the data of the receiver.
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    int j, size;
    FloatArray triplets;
    // IntArray *dofIDArry;

#  ifdef VERBOSE
    // VERBOSE_PRINT1("Instanciating node ",number)
#  endif

    IR_GIVE_FIELD(ir, coordinates, IFT_Vertex_coords, "coords"); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, IFT_Vertex_refine, "refine"); // Macro
    radius = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, radius, IFT_Vertex_radius, "radius"); // Macro

    //randomswitch = 0;
    // IR_GIVE_OPTIONAL_FIELD(ir, randomswitch, IFT_Vertex_randomswitch, "randomswitch"); // Macro

    return IRRT_OK;
}


Vertex *Vertex :: ofType(){
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).

  
  Vertex *vertex;
  
  vertex = new Vertex(number,grid);
  
  return vertex;
  
  

}
