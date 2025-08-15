#include "line.h"
#include "vertex.h"
#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include <iostream>
#endif

Line :: Line(int n, Grid *aGrid) : GridComponent(n, aGrid) //, coordinates()
{
    this->number = n;
    this->periodicElement = 0;
}

Line :: ~Line()
// Destructor.
{}


int
Line :: giveLocalVertex(int i)
// Returns the i-th coordinate of the receiver.
{
    if ( i > vertices.giveSize() ) {
        exit(1);
    }

    return vertices.at(i);
}

void
Line :: updateCrossSectionVertices(oofem::IntArray &nodes)
//Update cross-section vertices
{
    int flag;
    int size;
    for ( int m = 0; m < nodes.giveSize(); m++ ) {
        flag = 0;
        //Check if this node exist already
        for ( int i = 0; i < this->crossSectionVertices.giveSize(); i++ ) {
            if ( this->crossSectionVertices.at(i + 1) == nodes.at(m + 1) ) {
                flag = 1;
                break;
            }
        }
        if ( flag == 0 ) {
            size = crossSectionVertices.giveSize();
            crossSectionVertices.resize(size + 1);
            crossSectionVertices.at(size + 1) = nodes.at(m + 1);
        }
    }

    return;
}


void
Line :: updateCrossSectionElements(oofem::IntArray &elements)
//Update cross-section vertices
{
    int flag;
    int size;
    for ( int m = 0; m < elements.giveSize(); m++ ) {
        flag = 0;
        //Check if this node exist already
        for ( int i = 0; i < this->crossSectionElements.giveSize(); i++ ) {
            if ( this->crossSectionElements.at(i + 1) == elements.at(m + 1) ) {
                flag = 1;
                break;
            }
        }
        if ( flag == 0 ) {
            size = crossSectionElements.giveSize();
            crossSectionElements.resize(size + 1);
            crossSectionElements.at(size + 1) = elements.at(m + 1);
        }
    }

    return;
}


void
Line :: updateCrossSectionElement(int element)
//Update cross-section vertices
{
    int flag = 0;
    int size;
    //Check if this node exist already
    for ( int i = 0; i < this->crossSectionElements.giveSize(); i++ ) {
        if ( this->crossSectionElements.at(i + 1) == element ) {
            flag = 1;
            break;
        }
    }
    if ( flag == 0 ) {
        size = crossSectionElements.giveSize();
        crossSectionElements.resize(size + 1);
        crossSectionElements.at(size + 1) = element;
    }
    return;
}







Line *Line :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, LineSide,..).
{
    Line *line;

    line = new Line(number, grid);

    return line;
}


// update material
void Line:: updateMaterial(int typeOfMaterial)
{
    m_typeOfMaterial=typeOfMaterial;
    return;
}

int Line:: delaunayAreaCheck()
{
    //coordinates of the two nodes
    Vertex *vertexA, *vertexB;
    oofem::FloatArray coordsA(3), coordsB(3);

    vertexA  = grid->giveDelaunayVertex(this->vertices.at(1));
    vertexB  = grid->giveDelaunayVertex(this->vertices.at(2));

    for ( int i = 0; i < 3; i++ ) {
        coordsA.at(i + 1) =  vertexA->giveCoordinate(i + 1);
        coordsB.at(i + 1) =  vertexB->giveCoordinate(i + 1);
    }    
    
    oofem::FloatArray polygonCoords(3*this->crossSectionVertices.giveSize());
    //Read all the coordinates of the cross-section nodes and store them in an array
    for ( int i = 0; i < this->crossSectionVertices.giveSize(); i++){
      for (int k = 0; k < 3; k++){
	polygonCoords.at(3*i+k+1) = (grid->giveVoronoiVertex(crossSectionVertices.at(i+1)))->giveCoordinate(k+1);
      }
    }
    
    //Calculate normal vector
    oofem::FloatArray normal(3);
    for ( int i = 0; i < 3; i++ ) {
       normal.at(i + 1) = coordsB.at(i + 1) - coordsA.at(i + 1);
    }

    double length  = sqrt(pow(normal.at(1), 2.) + pow(normal.at(2), 2.) + pow(normal.at(3), 2.) );

    // Compute midpoint
    oofem::FloatArray midPoint(3);

    for ( int i = 0; i < 3; i++ ) {
      midPoint.at(i + 1) = 0.5 * ( coordsB.at(i + 1) + coordsA.at(i + 1) );
    }

    for ( int i = 0; i < 3; i++ ) {
      normal.at(i + 1) /= length;
    }

    //Construct two perpendicular axis so that n is normal to the plane which they create
    //Check, if one of the components of the normal-direction is zero
    oofem::FloatArray s(3), t(3);
    if (normal.at(1) == 0 ) {
        s.at(1) = 0.;
        s.at(2) = normal.at(3);
        s.at(3) = -normal.at(2);
    } else if ( normal.at(2) == 0 ) {
        s.at(1) = normal.at(3);
        s.at(2) = 0.;
        s.at(3) = -normal.at(1);
    } else {
        s.at(1) = normal.at(2);
        s.at(2) = -normal.at(1);
        s.at(3) = 0.;
    }

    s.normalize();

    t.beVectorProductOf(normal, s);
    t.normalize();

    //Set up rotation matrix
    FloatMatrix lcs(3, 3);

    for ( int i = 1; i <= 3; i++ ) {
        lcs.at(1, i) = normal.at(i);
        lcs.at(2, i) = s.at(i);
        lcs.at(3, i) = t.at(i);
    }


    //Calculate the local coordinates of the polygon vertices
    oofem::FloatArray help(3), test(3);
    oofem::FloatArray lpc(3 * crossSectionVertices.giveSize());
    for ( int k = 0; k < crossSectionVertices.giveSize(); k++ ) {
        for ( int n = 0; n < 3; n++ ) {
            help(n) = polygonCoords(3 * k + n);
        }

        test.beProductOf(lcs, help);
        for ( int n = 0; n < 3; n++ ) {
            lpc(3 * k + n) = test(n);
        }
    }

    
    double area = 0.;

    for ( int k = 0; k < crossSectionVertices.giveSize(); k++ ) {
      if ( k < crossSectionVertices.giveSize() - 1 ) {
            area += lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 2) - lpc(3 * ( k + 1 ) + 1) * lpc(3 * k + 2);
        } else {   //Back to zero for n+1
            area += lpc(3 * k + 1) * lpc(2) - lpc(1) * lpc(3 * k + 2);
        }
    }

    area *= 0.5;
    
    int areaFlag = 1;

    if(sqrt(fabs(area))/grid->giveDiameter() < grid->giveTol()){
      areaFlag = 0;
    }
    else{
      areaFlag = 1;
    }
    return areaFlag;
}



