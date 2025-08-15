#include "cylinder.h"
#include "curve.h"
#include "vertex.h"
#include "line.h"

#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#endif

Cylinder :: Cylinder(int n,Grid* aGrid) : Region(n,aGrid) //, coordinates()
{
  this->number = n;  
}

Cylinder :: ~Cylinder()
// Destructor.
{
}


IRResultType
Cylinder :: initializeFrom(InputRecord *ir)
// Gets from the source line from the data file all the data of the receiver.
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // int j, size;
    // oofem::FloatArray vertices;
    // IntArray *dofIDArry;
    IR_GIVE_FIELD(ir, line, IFT_Cylinder_line, "line");    
    IR_GIVE_FIELD(ir, radius, IFT_Cylinder_radius, "radius"); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, IFT_Cylinder_refine, "refine"); // Macro
   return IRRT_OK;

}


void
Cylinder :: defineBoundaries(oofem::FloatArray &boundaries)
//Determine the boundaries of the grid
{

    boundaries.resize(6);
    boundaries.zero();
    return;
}


void Cylinder :: findOutsiders(oofem::FloatArray &boundaries)
{
    //This class finds nodes and elements which are either outside or on the boundary of the specimen.
    //It finds also periodic nodes and elements which are needed for the output of periodic cells.
    //This is only working for a rectangular region.
  

    //Values of the outside flag:
    //For nodes: 0 = inside,  1 = outside, 2 = on boundary
   //For elements: 0 = inside, 1 = completetly outside (one node could be on boundary), 2 = one node inside and one outside, 3 = on boundary 

    IntArray nodes;
    oofem::FloatArray coords, coordsOne, coordsTwo;
    int outsideFlag;
    IntArray locationArray(2);

    double newTol = 1.e3*this->grid->giveTol();

    double help;
    double help2;
    double distance=0.;
    
    //DelaunayVertices
    for ( int i = 0; i < this->grid->giveNumberOfDelaunayVertices(); i++ ) {
      outsideFlag = 0;
      this->grid->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
      //Need to check three cases, which are two boundary circles at the end and the boundary cylinder.
      //Start with the circles
      distance = sqrt(pow(coords.at(2)-this->line.at(2),2.) + pow(coords.at(3)-this->line.at(3),2.));
      if(distance - newTol >  this->radius || coords.at(1)+newTol < line.at(1) || coords.at(1)-newTol >line.at(4)){
	outsideFlag = 1;
      }
      else if( fabs(distance - this->radius)< newTol || fabs(coords.at(1)-line.at(1))<newTol || fabs(coords.at(1)-line.at(4))<newTol ){ //On the boundary
	outsideFlag = 2;
      } 
      this->grid->giveDelaunayVertex(i + 1)->setOutsideFlag(outsideFlag);       
    }


    int updateFlag = 0;
    //VoronoiVertices
    for ( int i = 0; i < this->grid->giveNumberOfVoronoiVertices(); i++ ) {
        outsideFlag = 0;
        this->grid->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
	
	//Need to check three cases, which are two boundary circles at the end and the boundary cylinder.
	//Start with the circles
	distance = sqrt(pow(coords.at(2)-this->line.at(2),2.) + pow(coords.at(3)-this->line.at(3),2.));
	if(distance - newTol >  this->radius || coords.at(1) + newTol < line.at(1) || coords.at(1) - newTol >line.at(4)){
	  outsideFlag = 1;
	}
	else if( fabs(distance - this->radius)< newTol || fabs(coords.at(1)-line.at(1))<newTol || fabs(coords.at(1)-line.at(4))<newTol ){ //On the boundary
	  outsideFlag = 2;
	} 
 	
        if ( outsideFlag == 1) {
            //Check if nodes are too far outside and shift them to the boundaries of mirrored sphere
            //This is needed to get the octree localizer working
	  updateFlag = 0;
	  if ( ( coords.at(1) + newTol ) < (2.*line.at(1) - line.at(4)) ) {
	    coords.at(1) = 2.*line.at(1)-line.at(4);
	    updateFlag = 1;
	  }
	  if ( ( coords.at(1) + newTol ) > (2.*line.at(4) - line.at(1)) ) {
	    coords.at(1) = 2.*line.at(4)-line.at(1);
	    updateFlag = 1;
	  }
	  if(distance >2.*this->radius){
	    coords.at(2) = 2.*this->radius/distance*coords.at(2);
	    coords.at(3) = 2.*this->radius/distance*coords.at(3);
	    updateFlag = 1;
	  }
	  if( updateFlag == 1 ) {
	    this->grid->giveVoronoiVertex(i + 1)->setCoordinates(coords);
	  }
        }
        this->grid->giveVoronoiVertex(i + 1)->setOutsideFlag(outsideFlag);
    }


    
    //Delaunay elements
    int boundaryNumber = -1;
    for ( int i = 0; i < this->grid->giveNumberOfDelaunayLines(); i++ ) {
        outsideFlag = 0;
        this->grid->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
        if ( (this->grid->giveDelaunayVertex( nodes.at(1) )->giveOutsideFlag() == 1 &&
              this->grid->giveDelaunayVertex( nodes.at(2) )->giveOutsideFlag() == 0 ) ||
             (this->grid->giveDelaunayVertex( nodes.at(1) )->giveOutsideFlag() == 0 &&
              this->grid->giveDelaunayVertex( nodes.at(2) )->giveOutsideFlag() == 1 ) ) {
	  outsideFlag = 2;
        } else if ((this->grid->giveDelaunayVertex( nodes.at(1) )->giveOutsideFlag() == 1 &&
		    this->grid->giveDelaunayVertex( nodes.at(2) )->giveOutsideFlag() == 1) ||
		   (this->grid->giveDelaunayVertex( nodes.at(1) )->giveOutsideFlag() == 2 &&
		    this->grid->giveDelaunayVertex( nodes.at(2) )->giveOutsideFlag() == 1) ||
		   (this->grid->giveDelaunayVertex( nodes.at(1) )->giveOutsideFlag() == 1 &&
		    this->grid->giveDelaunayVertex( nodes.at(2) )->giveOutsideFlag() == 2)) {
	  outsideFlag = 1;
        } else if (this->grid->giveDelaunayVertex( nodes.at(1) )->giveOutsideFlag() == 2 &&
                   this->grid->giveDelaunayVertex( nodes.at(2) )->giveOutsideFlag() == 2 ) {
	  outsideFlag = 3;
        }
        this->grid->giveDelaunayLine(i + 1)->setOutsideFlag(outsideFlag);
    }

    //Voronoi elements.
    for ( int i = 0; i < this->grid->giveNumberOfVoronoiLines(); i++ ) {
        outsideFlag = 0;
        this->grid->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
        if ( nodes.at(1) == 0 || nodes.at(2) == 0 ) {
            outsideFlag++;
        } else {
            if ( ( this->grid->giveVoronoiVertex( nodes.at(1) )->giveOutsideFlag() == 1 &&
                   this->grid->giveVoronoiVertex( nodes.at(2) )->giveOutsideFlag() == 0 ) ||
                 ( this->grid->giveVoronoiVertex( nodes.at(1) )->giveOutsideFlag() == 0 &&
                   this->grid->giveVoronoiVertex( nodes.at(2) )->giveOutsideFlag() == 1 ) ) {
                outsideFlag = 2;
            } else if ( this->grid->giveVoronoiVertex( nodes.at(1) )->giveOutsideFlag() == 1 ||
                        this->grid->giveVoronoiVertex( nodes.at(2) )->giveOutsideFlag() == 1 ) {
                outsideFlag = 1;
            } else if ( this->grid->giveVoronoiVertex( nodes.at(1) )->giveOutsideFlag() == 2 &&
                        this->grid->giveVoronoiVertex( nodes.at(2) )->giveOutsideFlag() == 2 ) {
                outsideFlag = 3;
            }
        }
        this->grid->giveVoronoiLine(i + 1)->setOutsideFlag(outsideFlag);
    }
    
    
    return;
}



int
Cylinder :: modifyVoronoiCrossSection(int elementNumber)
{
    IntArray crossSectionElements;
    IntArray crossSectionVertices;

    IntArray nodes(2);

    oofem::FloatArray coords(3);

    double newTol = 1.e3*this->grid->giveTol();
    
    this->grid->giveDelaunayLine(elementNumber)->giveCrossSectionElements(crossSectionElements);
    this->grid->giveDelaunayLine(elementNumber)->giveCrossSectionVertices(crossSectionVertices);

    int elementSize = crossSectionElements.giveSize();
    int vertexSize = crossSectionVertices.giveSize();
    int elementCounter = 0;
    double distance = 0.;
    for ( int i = 0; i < elementSize; i++ ) {
        this->grid->giveVoronoiLine( crossSectionElements.at(i + 1) )->giveLocalVertices(nodes);
        if ( this->grid->giveVoronoiLine( crossSectionElements.at(i + 1) )->giveOutsideFlag() == 1 || this->grid->giveVoronoiLine( crossSectionElements.at(i + 1) )->giveOutsideFlag() == 3 ) {//completely outside or on the surface
            crossSectionElements.at(i + 1) = 0;
            elementCounter++;
        } else if ( this->grid->giveVoronoiLine( crossSectionElements.at(i + 1) )->giveOutsideFlag() == 2 )         {//crossing boundary
            //Keep element but shift nodes to surface of cylinder;
            for ( int k = 0; k < 2; k++ ) {
                if ( this->grid->giveVoronoiVertex( nodes.at(k + 1) )->giveOutsideFlag() == 1 ) {
                    this->grid->giveVoronoiVertex( nodes.at(k + 1) )->giveCoordinates(coords);

		    if( coords.at(1) + newTol < line.at(1) ){
		      coords.at(1) = line.at(1);
		    }
		    if( coords.at(1) - newTol > line.at(4) ){
		      coords.at(1) = line.at(4);
		    }

		    distance = sqrt(pow(coords.at(2)-this->line.at(2),2.) + pow(coords.at(3)-this->line.at(3),2.));
		    if(distance-newTol>radius){
		      coords.at(2) = this->radius/distance*coords.at(2);
		      coords.at(3) = this->radius/distance*coords.at(3);
		    }
                    this->grid->giveVoronoiVertex( nodes.at(k + 1) )->setCoordinates(coords);
		    this->grid->giveVoronoiVertex( nodes.at(k + 1) )->setOutsideFlag(2);
                }
            }
            this->grid->giveVoronoiLine( crossSectionElements.at(i + 1) )->setOutsideFlag(0);
        }
    }

    //Resize cross-section elements
    int newSize = elementSize - elementCounter;
    IntArray modifiedCrossSectionElements(newSize);
    int help = 0;
    for ( int i = 0; i < elementSize; i++ ) {
        if ( crossSectionElements.at(i + 1) == 0 ) {
            help++;
        } else   {
            modifiedCrossSectionElements.at(i - help + 1) = crossSectionElements.at(i + 1);
        }
    }
    this->grid->giveDelaunayLine(elementNumber)->setCrossSectionElements(modifiedCrossSectionElements);

    //Determine how many vertices should be deleted;
    int nodeCounter = 0;
    for ( int i = 0; i < vertexSize; i++ ) {
        if ( this->grid->giveVoronoiVertex( crossSectionVertices.at(i + 1) )->giveOutsideFlag() == 1 ) {
            crossSectionVertices.at(i + 1) = 0;
            nodeCounter++;
        }
    }

    //Write modified vertex vector
    IntArray modifiedCrossSectionVertices(vertexSize - nodeCounter);
    help = 0;
    for ( int i = 0; i < vertexSize; i++ ) {
        if ( crossSectionVertices.at(i + 1) == 0 ) {
            help++;
        } else   {
            modifiedCrossSectionVertices.at(i - help + 1) = crossSectionVertices.at(i + 1);
        }
    }
    this->grid->giveDelaunayLine(elementNumber)->setCrossSectionVertices(modifiedCrossSectionVertices);

    int flag = 1;
    if(modifiedCrossSectionVertices.giveSize() < 3){
      flag = 0;
      this->grid->giveDelaunayLine(elementNumber)->setOutsideFlag(1);
    }

    
    return flag;
}


int
Cylinder :: areaCheck(int elementNumber)
{

    IntArray crossSectionElements;
    IntArray crossSectionVertices;

    IntArray nodes(2);

    oofem::FloatArray coords(3);

    double newTol = 1.e3*this->grid->giveTol();
    
    this->grid->giveDelaunayLine(elementNumber)->giveCrossSectionElements(crossSectionElements);
    this->grid->giveDelaunayLine(elementNumber)->giveCrossSectionVertices(crossSectionVertices);

    int elementSize = crossSectionElements.giveSize();
    int vertexSize = crossSectionVertices.giveSize();

    return 1;
}



Cylinder *Cylinder :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    Cylinder *cylinder;
    
    cylinder = new Cylinder(number,grid);

    return cylinder;
}
