#include "tetra.h"
#include "vertex.h"
#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include <iostream>
#endif

Tetra :: Tetra(int n, Grid *aGrid) : GridComponent(n, aGrid) {
    this->number = n;
}

Tetra :: ~Tetra()
// Destructor.
{}

int
Tetra :: giveLocalVertex(int i)
// Returns the i-th coordinate of the receiver.
{
    if ( i > vertices.giveSize() ) {
        exit(1);
    }

    return vertices.at(i);
}


Tetra *Tetra :: ofType()
// Returns a new Tetrahedra, which has the same number than the receiver.
{
    Tetra *tetra;

    tetra = new Tetra(number, grid);

    return tetra;
}


void Tetra :: giveCoordinates(oofem::FloatArray& coords)
//gives the coordinates of the nodes of the tetra in a vector
{
  coords.resize(12);
  oofem::FloatArray nodeCoords;
  for(int i = 1; i<=4; i++){    
    this->grid->giveDelaunayVertex(vertices.at(i))->giveCoordinates(nodeCoords);
    coords.at(1+(i-1)*3) = nodeCoords.at(1);
    coords.at(2+(i-1)*3) = nodeCoords.at(2);
    coords.at(3+(i-1)*3) = nodeCoords.at(3);
  }  
}
