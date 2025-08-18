#include "vertex.h"
#include "gridcomponent.h"

#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#endif

Vertex :: Vertex(int n, Grid *aGrid) : GridComponent(n, aGrid)//, coordinates()
{
  this->printFlag = 0;
  this->number = n;
    this->outsideFlag = 0;
    this->helpOutsideFlag.zero();    
    this->location = 0;
    this->boundaryFlag = -1;
    this->cellVertices.zero();
    this->cellElements.zero();
}


Vertex :: ~Vertex()
// Destructor.
{}

double Vertex :: giveCoordinate(int i)
// Returns the i-th coordinate of the receiver.
{
    if ( i > coordinates.giveSize() ) {
        return 0.;
    }

    return coordinates.at(i);
}

void
Vertex :: initializeFrom(ConverterInputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{
    IR_GIVE_FIELD(ir, coordinates, _IFT_Vertex_coords); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, _IFT_Vertex_refine); // Macro
    radius = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, radius, _IFT_Vertex_radius); // Macro

    return;
}

Vertex *Vertex :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, LineSide,..).
{
    Vertex *vertex;

    vertex = new Vertex(number, grid);

    return vertex;
}

void Vertex :: setLocalLine(int elem)
{
    int size = localLines.giveSize();
    localLines.resize(size + 1);
    localLines.at(size + 1) = elem;
    return;
}

void Vertex :: setLocalTetra(int elem)
{
    int size = localTetras.giveSize();
    localTetras.resize(size + 1);
    localTetras.at(size + 1) = elem;
    return;
}



void Vertex :: setLocalLink(int elem)
{
    int size = localLinks.giveSize();
    localLinks.resize(size + 1);
    localLinks.at(size + 1) = elem;
    return;
}


void
Vertex :: updateCellVertices(oofem::IntArray &nodes)
//Update cross-section vertices
{
    int flag;
    int size;
    for ( int m = 0; m < nodes.giveSize(); m++ ) {
        flag = 0;
        //Check if this node exist already
        for ( int i = 0; i < this->cellVertices.giveSize(); i++ ) {
            if ( this->cellVertices.at(i + 1) == nodes.at(m + 1) ) {
                flag = 1;
                break;
            }
        }
        if ( flag == 0 ) {
            size = cellVertices.giveSize();
            cellVertices.resize(size + 1);
            cellVertices.at(size + 1) = nodes.at(m + 1);
        }
    }

    return;
}


void
Vertex :: updateCellElements(oofem::IntArray &elements)
//Update cross-section vertices
{
    int flag;
    int size;
    for ( int m = 0; m < elements.giveSize(); m++ ) {
        flag = 0;
        //Check if this node exist already
        for ( int i = 0; i < this->cellElements.giveSize(); i++ ) {
            if ( this->cellElements.at(i + 1) == elements.at(m + 1) ) {
                flag = 1;
                break;
            }
        }	
        if ( flag == 0 ) {
            size = cellElements.giveSize();
            cellElements.resize(size + 1);
            cellElements.at(size + 1) = elements.at(m + 1);
        }
    }

    return;
}
