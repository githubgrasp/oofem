#include "prism.h"
#include "curve.h"
#include "vertex.h"
#include "surface.h"
#include "intarray.h"
#include "line.h"
#include "fibre.h"
#include "datareader.h"
#include "octreegridlocalizer.h"
#include "tetra.h"

#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
#endif

Prism :: Prism(int n, Grid *aGrid) : Region(n, aGrid) //, coordinates()
{
    this->number = n;
}

Prism :: ~Prism()
// Destructor.
{}

int
Prism :: giveLocalSurface(int i)
// Returns the i-th coordinate of the receiver.
{
    if ( i > surfaces.giveSize() ) {
        return 0.;
    }

    return surfaces.at(i);
}

void
Prism :: initializeFrom(ConverterInputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{
    IR_GIVE_FIELD(ir, this->box, _IFT_Prism_box);    

    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, _IFT_Prism_refine); // Macro
    return;
}

Prism *Prism :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    Prism *region;

    region = new Prism(number, grid);

    return region;
}



int
Prism :: giveSwitches(oofem::IntArray &switches, const oofem::FloatArray &coords) {
  //Find location and switches for coords
  //The problem is in here. I cannot use the ordinary switch function to fix the problem with vertices on edges. Location is wrong allocated.
  oofem::IntArray periodicityFlag(3);
  grid->givePeriodicityFlag(periodicityFlag);
  
  oofem::FloatArray boundaries(3);
  this->defineBoundaries(boundaries);
  oofem::FloatArray specimenDimension(3);
  specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
  specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
  specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

  oofem::FloatArray coordsTest(3);
  int counter = 1;
  int location = 0.;
  switches.resize(3);
  switches.zero();
  for ( int x = -1; x <  2; x++ ) {
    for ( int y = -1; y <  2; y++ ) {
      for ( int z = -1; z <  2; z++ ) {
        if ( !( z == 0 && y == 0 && x == 0 ) ) {
        //Idea is to shift the coordinates into the opposite direction
        coordsTest.at(1) = coords.at(1) - x * specimenDimension.at(1);
        coordsTest.at(2) = coords.at(2) - y * specimenDimension.at(2);
        coordsTest.at(3) = coords.at(3) - z * specimenDimension.at(3);
        location++;
        //Make a check that coords are in the cell including boundaries
        if ( coordsTest.at(1) > boundaries.at(1) -this->grid->giveTol() && coordsTest.at(1) < boundaries.at(2) +this->grid->giveTol() &&
	       coordsTest.at(2) > boundaries.at(3) -this->grid->giveTol() && coordsTest.at(2) < boundaries.at(4) +this->grid->giveTol() &&
	       coordsTest.at(3) > boundaries.at(5) -this->grid->giveTol() && coordsTest.at(3) < boundaries.at(6) +this->grid->giveTol() ) {
	  if(!((periodicityFlag.at(1) == 0 && x != 0) || (periodicityFlag.at(2) == 0 && y != 0) || (periodicityFlag.at(3) == 0 && z != 0)) ){
	    switches.at(1) = x;
	    switches.at(2) = y;
	    switches.at(3) = z;
	    return location;
	  }
        }
        }
      }
    }
  }
  return location;
  //  printf("pb in giveSwitches \n");
}

int Prism::giveCornerSwitches(oofem::IntArray& switches, const oofem::FloatArray& coords)
{
    //Find location and switches for corners

    oofem::IntArray periodicityFlag(3);
    grid->givePeriodicityFlag(periodicityFlag);

    oofem::FloatArray boundaries(6);
    this->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    switches.resize(3);
    switches.zero();

    if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0 ) {
        //Periodicity in x
        if ( fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() ) {
	  switches.at(1)=1;
	  switches.at(2)=0;
	  switches.at(3)=0;
	  
            return 22;
        } else {
            return 0;
        }
    } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
        //Periodicity in y
        if ( fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() ) {
            switches.at(1)=0;
            switches.at(2)=1;
            switches.at(3)=0;

            return 16;
        } else {
            return 0;
        }
    } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
        //Periodicity in z
        if ( fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() ) {
            switches.at(1)=0;
            switches.at(2)=0;
            switches.at(3)=1;

            return 14;
        } else {
            return 0;
        }
    } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
        //Periodicity in x y
        if ( fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() &&
             ( fabs(coords.at(3) - boundaries.at(5)) < this->grid->giveTol() ||
             fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() ) &&
             fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() )  { //corner at [xmax, ymax, zmin or zmax]
            switches.at(1) = 1;
            switches.at(2) = 1;
            switches.at(3) = 0;

            return 25;
        } else if ( fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() &&
                    ( fabs(coords.at(3) - boundaries.at(5)) < this->grid->giveTol() ||
                    fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() ) &&
                    fabs(coords.at(1) - boundaries.at(1)) < this->grid->giveTol() )  { //corner at [xmin, ymax, zmin or zmax]
            switches.at(1) = 0;
            switches.at(2) = 1;
            switches.at(3) = 0;

            return 16;
        } else if ( fabs(coords.at(2) - boundaries.at(3)) < this->grid->giveTol() &&
                    ( fabs(coords.at(3) - boundaries.at(5)) < this->grid->giveTol() ||
                    fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() ) &&
                    fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() )  { //corner at [xmax, ymin, zmin or zmax]
            switches.at(1) = 1;
            switches.at(2) = 0;
            switches.at(3) = 0;

            return 22;
        } else {
            return 0;
        }
    } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
        //Periodicity in x z
        if ( fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() &&
             ( fabs(coords.at(2) - boundaries.at(3)) < this->grid->giveTol() ||
             fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() ) &&
             fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() )  { //corner at [xmax, ymin or ymax, zmax]
            switches.at(1) = 1;
            switches.at(2) = 0;
            switches.at(3) = 1;

            return 23;
        } else if ( fabs(coords.at(3) - boundaries.at(5)) < this->grid->giveTol() &&
                    ( fabs(coords.at(2) - boundaries.at(3)) < this->grid->giveTol() ||
                    fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() ) &&
                    fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() )  { //corner at [xmax, ymin or ymax, zmin]
            switches.at(1) = 1;
            switches.at(2) = 0;
            switches.at(3) = 0;

            return 22;
        } else if ( fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() &&
                    ( fabs(coords.at(2) - boundaries.at(3)) < this->grid->giveTol() ||
                    fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() ) &&
                    fabs(coords.at(1) - boundaries.at(1)) < this->grid->giveTol() )  { //corner at [xmin, ymin or ymax, zmax]
            switches.at(1) = 0;
            switches.at(2) = 0;
            switches.at(3) = 1;

            return 14;
        } else {
            return 0;
        }
    } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
        //Periodicity in y z
        if ( fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() &&
             ( fabs(coords.at(1) - boundaries.at(1)) < this->grid->giveTol() ||
             fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() ) &&
             fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() )  { //corner at [xmin or xmax, ymax, zmax]
            switches.at(1) = 0;
            switches.at(2) = 1;
            switches.at(3) = 1;

            return 17;
        } else if ( fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() &&
                    ( fabs(coords.at(1) - boundaries.at(1)) < this->grid->giveTol() ||
                    fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() ) &&
                    fabs(coords.at(2) - boundaries.at(3)) < this->grid->giveTol() )  { //corner at [xmin or xmax, ymin, zmax]
            switches.at(1) = 0;
            switches.at(2) = 0;
            switches.at(3) = 1;

            return 14;
        } else if ( fabs(coords.at(3) - boundaries.at(5)) < this->grid->giveTol() &&
                    ( fabs(coords.at(1) - boundaries.at(1)) < this->grid->giveTol() ||
                    fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() ) &&
                    fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() )  { //corner at [xmin or xmax, ymax, zmin]
            switches.at(1) = 0;
            switches.at(2) = 1;
            switches.at(3) = 0;

            return 16;
        } else {
            return 0;
        }
    } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
        //Periodicity in x y z
        //For a cube, there are 7 corners that must be mapped to the one periodic image.
        //add 3 remaining corners, this function
        if ( fabs(coords.at(2) - boundaries.at(3)) < this->grid->giveTol() &&
             fabs(coords.at(3) - boundaries.at(5)) < this->grid->giveTol() &&
             fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() )  { //corner at [xmax, ymin, zmin]
            switches.at(1) = 1;
            switches.at(2) = 0;
            switches.at(3) = 0;

            return 22;
        } else if ( fabs(coords.at(2) - boundaries.at(3)) < this->grid->giveTol() &&
                    fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() &&
                    fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() ) { //corner at [xmax, ymin, zmax]
            switches.at(1) = 1;
            switches.at(2) = 0;
            switches.at(3) = 1;

            return 23;
        } else if ( fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() &&
                    fabs(coords.at(3) - boundaries.at(5)) < this->grid->giveTol() &&
                    fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() ) { //corner at [xmax, ymax, zmin]
            switches.at(1) = 1;
            switches.at(2) = 1;
            switches.at(3) = 0;

            return 25;
        } else if ( fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() &&
                    fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() &&
                    fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() ) { //corner at [xmax, ymax, zmax]
            switches.at(1) = 1;
            switches.at(2) = 1;
            switches.at(3) = 1;

            return 26;
        } else if ( fabs(coords.at(1) - boundaries.at(1)) < this->grid->giveTol() &&
                    fabs(coords.at(3) - boundaries.at(5)) < this->grid->giveTol() &&
                    fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() ) { //corner at [xmin, ymax, zmin]
            switches.at(1)=0;
            switches.at(2)=1;
            switches.at(3)=0;

            return 16;
        } else if ( fabs(coords.at(1) - boundaries.at(1)) < this->grid->giveTol() &&
                    fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() &&
                    fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() ) { //corner at [xmin, ymax, zmax]
            switches.at(1)=0;
            switches.at(2)=1;
            switches.at(3)=1;

            return 17;
        } else if ( fabs(coords.at(2) - boundaries.at(3)) < this->grid->giveTol() &&
                    fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() &&
                    fabs(coords.at(1) - boundaries.at(1)) < this->grid->giveTol() ) { //corner at [xmin, ymin, zmax]
            switches.at(1)=0;
            switches.at(2)=0;
            switches.at(3)=1;

            return 14;
        } else {
            return 0;
        }
    } else {
      converter::error("No periodicity detected.\n");
            return 0;
    }
}

int Prism::giveEdgeSwitches(oofem::IntArray& switches, const oofem::FloatArray& coords)
{
    //Find location and switches for coords

    oofem::IntArray periodicityFlag(3);
    grid->givePeriodicityFlag(periodicityFlag);

    oofem::FloatArray boundaries(6);
    this->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    switches.resize(3);
    switches.zero();

    if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0 ) {
        //Periodicity in x
        if ( fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() ) {
            switches.at(1)=1;
            switches.at(2)=0;
            switches.at(3)=0;

            return 22;
        } else {
            return 0;
        }
    } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
        //Periodicity in y
        if ( fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() ) {
            switches.at(1)=0;
            switches.at(2)=1;
            switches.at(3)=0;

            return 16;
        } else {
            return 0;
        }
    } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
        //Periodicity in z
        if ( fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() ) {
            switches.at(1)=0;
            switches.at(2)=0;
            switches.at(3)=1;

            return 14;
        } else {
            return 0;
        }
    } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
        //Periodicity in x y
        if ( fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol()  &&
             fabs(coords.at(2) - boundaries.at(4)) > this->grid->giveTol() ) { //nodes on edges [xmax,zmin], [xmax,zmax] and [xmax,ymin]
            switches.at(1) = 1;
            switches.at(2) = 0;
            switches.at(3) = 0;

            return 22;
        } else if ( fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() &&
                    fabs(coords.at(1) - boundaries.at(2)) > this->grid->giveTol() ) { //nodes on edges [ymax,zmin], [ymax,zmax] and [xmin,ymax]
            switches.at(1) = 0;
            switches.at(2) = 1;
            switches.at(3) = 0;

            return 16;
        } else if ( fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() &&
                    fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() ) { //nodes on edge [xmax,ymax]
            switches.at(1) = 1;
            switches.at(2) = 1;
            switches.at(3) = 0;

            return 25;
        } else {
            return 0;
        }
    } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
        //Periodicity in x z
        if ( fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol()  &&
             fabs(coords.at(3) - boundaries.at(6)) > this->grid->giveTol() ) { //nodes on edges [xmax,ymin], [xmax,ymax] and [xmax,zmin]
            switches.at(1) = 1;
            switches.at(2) = 0;
            switches.at(3) = 0;

            return 22;
        } else if ( fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() &&
                    fabs(coords.at(1) - boundaries.at(2)) > this->grid->giveTol() ) { //nodes on edges [ymin,zmax], [ymax,zmax] and [xmin,zmax]
            switches.at(1) = 0;
            switches.at(2) = 0;
            switches.at(3) = 1;

            return 14;
        } else if ( fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() &&
                    fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() ) { //nodes on edge [xmax,zmax]
            switches.at(1) = 1;
            switches.at(2) = 0;
            switches.at(3) = 1;

            return 23;
        } else {
            return 0;
        }
    } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
        //Periodicity in y z
        if ( fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol()  &&
             fabs(coords.at(3) - boundaries.at(6)) > this->grid->giveTol() ) { //nodes on edges [xmin,ymax], [xmax,ymax] and [ymax,zmin]
            switches.at(1) = 0;
            switches.at(2) = 1;
            switches.at(3) = 0;

            return 16;
        } else if ( fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() &&
                    fabs(coords.at(2) - boundaries.at(4)) > this->grid->giveTol() ) { //nodes on edges [xmin,zmax], [xmax,zmax] and [ymin,zmax]
            switches.at(1) = 0;
            switches.at(2) = 0;
            switches.at(3) = 1;

            return 14;
        } else if ( fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() &&
                    fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() ) { //nodes on edge [ymax,zmax]
            switches.at(1) = 0;
            switches.at(2) = 1;
            switches.at(3) = 1;

            return 17;
        } else {
            return 0;
        }
    } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
        //Periodicity in x y z
        if ( (fabs(coords.at(2) - boundaries.at(3)) < this->grid->giveTol()  ||
             fabs(coords.at(3) - boundaries.at(5)) < this->grid->giveTol()) &&
             fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() ) { //nodes on edges [xmax,ymin] and [xmax,zmin]
            switches.at(1) = 1;
            switches.at(2) = 0;
            switches.at(3) = 0;

            return 22;
        } else if ( fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() &&
                    fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() ) { //nodes on edge [xmax,zmax]
            switches.at(1) = 1;
            switches.at(2) = 0;
            switches.at(3) = 1;

            return 23;
        } else if ( fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() &&
                    fabs(coords.at(1) - boundaries.at(2)) < this->grid->giveTol() )  { //nodes on edge [xmax,ymax]
            switches.at(1) = 1;
            switches.at(2) = 1;
            switches.at(3) = 0;

            return 25;
        } else if ( (fabs(coords.at(1) - boundaries.at(1)) < this->grid->giveTol() ||
                    fabs(coords.at(3) - boundaries.at(5)) < this->grid->giveTol()) &&
                    fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() ) { //nodes on edges [xmin,ymax] and [ymax,zmin]
            switches.at(1)=0;
            switches.at(2)=1;
            switches.at(3)=0;

            return 16;
        } else if ( fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() &&
                    fabs(coords.at(2) - boundaries.at(4)) < this->grid->giveTol() ) { //nodes on edge [ymax,zmax]
            switches.at(1)=0;
            switches.at(2)=1;
            switches.at(3)=1;

            return 17;
        } else if ( (fabs(coords.at(2) - boundaries.at(3)) < this->grid->giveTol() ||
                    fabs(coords.at(1) - boundaries.at(1)) < this->grid->giveTol()) &&
                    fabs(coords.at(3) - boundaries.at(6)) < this->grid->giveTol() ) { //nodes on edges [ymin,zmax] and [xmin,zmax]
            switches.at(1)=0;
            switches.at(2)=0;
            switches.at(3)=1;

            return 14;
        } else {
            return 0;
        }
    } else {
      converter::error("No periodicity detected\n");
        return 0;
    }
}


int
Prism :: modifyVoronoiCrossSection(int elementNumber)
{
    oofem::IntArray crossSectionElements;
    oofem::IntArray crossSectionVertices;

    oofem::IntArray nodes(2);

    oofem::FloatArray coords(3);

    oofem::FloatArray boundaries(3);
    this->defineBoundaries(boundaries);

    this->grid->giveDelaunayLine(elementNumber)->giveCrossSectionElements(crossSectionElements);
    this->grid->giveDelaunayLine(elementNumber)->giveCrossSectionVertices(crossSectionVertices);

    int elementSize = crossSectionElements.giveSize();
    int vertexSize = crossSectionVertices.giveSize();
    int elementCounter = 0;
    for ( int i = 0; i < elementSize; i++ ) {
        this->grid->giveVoronoiLine( crossSectionElements.at(i + 1) )->giveLocalVertices(nodes);
        if ( this->grid->giveVoronoiLine( crossSectionElements.at(i + 1) )->giveOutsideFlag() == 1 || this->grid->giveVoronoiLine( crossSectionElements.at(i + 1) )->giveOutsideFlag() == 3 ) {//completely outside or on the surface
	  //	  if(this->grid->giveVoronoiLine( crossSectionElements.at(i + 1) )->giveOutsideFlag() == 1){
	    this->grid->giveVoronoiLine( crossSectionElements.at(i + 1) )->setOutsideFlag(-1);//Set it negative so that is not dealt with in periodicity check.
	    //	  }
	  crossSectionElements.at(i + 1) = 0;
	  elementCounter++;
        } else if ( this->grid->giveVoronoiLine( crossSectionElements.at(i + 1) )->giveOutsideFlag() == 2 )         {//crossing boundary
            //Keep element but shift nodes;

            for ( int k = 0; k < 2; k++ ) {
                if ( this->grid->giveVoronoiVertex( nodes.at(k + 1) )->giveOutsideFlag() == 1 ) {
                    this->grid->giveVoronoiVertex( nodes.at(k + 1) )->giveCoordinates(coords);
                    for ( int n = 0; n < 3; n++ ) {
                        if ( ( coords.at(n + 1) +this->grid->giveTol() ) < boundaries.at(2 * n + 1) ) {
                            coords.at(n + 1) = boundaries.at(2 * n + 1);
                            this->grid->giveVoronoiVertex( nodes.at(k + 1) )->setOutsideFlag(2);
                            break;
                        } else if ( ( coords.at(n + 1) -this->grid->giveTol() ) > boundaries.at(2 * n + 2) )   {
                            coords.at(n + 1) = boundaries.at(2 * n + 2);
                            this->grid->giveVoronoiVertex( nodes.at(k + 1) )->setOutsideFlag(2);
                            break;
                        }
                    }
                    this->grid->giveVoronoiVertex( nodes.at(k + 1) )->setCoordinates(coords);
                }
            }
            this->grid->giveVoronoiLine( crossSectionElements.at(i + 1) )->setOutsideFlag(0);
        }
    }


    //Resize cross-section elements
    int newSize = elementSize - elementCounter;
    oofem::IntArray modifiedCrossSectionElements(newSize);
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
    oofem::IntArray modifiedCrossSectionVertices(vertexSize - nodeCounter);
    help = 0;
    for ( int i = 0; i < vertexSize; i++ ) {
        if ( crossSectionVertices.at(i + 1) == 0 ) {
            help++;
        } else   {
            modifiedCrossSectionVertices.at(i - help + 1) = crossSectionVertices.at(i + 1);
        }
    }
    this->grid->giveDelaunayLine(elementNumber)->setCrossSectionVertices(modifiedCrossSectionVertices);
    return 1;
}



void Prism :: findOutsiders(oofem::FloatArray &boundaries)
{
  //This class finds nodes and elements which are either outside or on the boundary of the specimen.
  //It finds also periodic nodes and elements which are needed for the output of periodic cells.
  //This is only working for a rectangular region.
  
  //Values of the outside flag:
  //For nodes: 0 = inside,  1 = outside, 2 = on boundary
  //For line elements: 0 = inside, 1 = completetly outside (one node could be on boundary), 2 = one node inside and one outside, 3 = on boundary
  //For tetras: 0 = inside, 1 = completely outside, 2 = some nodes intside others outside, 3 = all nodes on boundary.
  
  //Introduce a help node outside flag array, which uses the same notation, but for the three directions separately
  oofem::IntArray helpOutsideFlag(6);
  
  oofem::IntArray nodes;
  oofem::FloatArray coords(3), coordsOne(3), coordsTwo(3);
  int outsideFlag;
  oofem::IntArray locationArray(2);

  double newTol = 1.e-3 *this->grid->giveTol();

  //The implementation should depend on the type of region. Ideally findOutsiders should be implemented in the region.
  //Types of regions could be prisms (this is currently the region), cylinders, spheres (what we are currently working on), etc.

  //Determine specimen dimensions
  oofem::FloatArray specimenDimension(3);
  specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
  specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
  specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

  oofem::IntArray periodicityFlag(3);
  grid->givePeriodicityFlag(periodicityFlag);
    
  double help;
  double help2;
  //DelaunayVertices
  for ( int i = 0; i < this->grid->giveNumberOfDelaunayVertices(); i++ ) {
    helpOutsideFlag.zero();
    outsideFlag = 0;
    this->grid->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
    for ( int k = 0; k < 3; k++ ) {
      if ( ( coords.at(k + 1) + newTol ) < boundaries.at(2 * k + 1) ) {
	helpOutsideFlag.at(2*k+1) = 1;
      } else if ( fabs( coords.at(k + 1) - boundaries.at(2 * k + 1) ) < newTol ) {                 //On the boundary
	helpOutsideFlag.at(2.*k+1) = 2;
      } else if ( ( coords.at(k + 1) - newTol ) > boundaries.at(2 * k + 2) ) {
	helpOutsideFlag.at(2*k+2) = 1;
      } else if ( fabs( coords.at(k + 1) - boundaries.at(2 * k + 2) ) < newTol ) {
	helpOutsideFlag.at(2*k+2) = 2;
      }
    }

    //Assign outside flag based on helpOutsideFlag values
    for(int k = 0; k<6;k++){
      if(helpOutsideFlag.at(k+1) == 1){
	outsideFlag = 1;
	break;
      }
      else if(helpOutsideFlag.at(k+1) == 2){
	outsideFlag = 2;
      }
    }
        
    this->grid->giveDelaunayVertex(i + 1)->setOutsideFlag(outsideFlag);
    this->grid->giveDelaunayVertex(i + 1)->setHelpOutsideFlag(helpOutsideFlag);
  }

  int updateFlag = 0;
  //Voronoi vertices
  for ( int i = 0; i < this->grid->giveNumberOfVoronoiVertices(); i++ ) {
    outsideFlag = 0;
    helpOutsideFlag.zero();
    this->grid->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
    for ( int k = 0; k < 3; k++ ) {
      if ( ( coords.at(k + 1) + newTol ) < boundaries.at(2 * k + 1) ) {
	helpOutsideFlag.at(2*k+1) = 1;
      } else if ( fabs( coords.at(k + 1) - boundaries.at(2 * k + 1) ) < newTol ) {                 //On the boundary
	helpOutsideFlag.at(2*k+1) = 2;
      } else if ( ( coords.at(k + 1) - newTol ) > boundaries.at(2 * k + 2) ) {
	helpOutsideFlag.at(2*k+2) = 1;
      } else if ( fabs( coords.at(k + 1) - boundaries.at(2 * k + 2) ) < newTol ) {
	helpOutsideFlag.at(2*k+2) = 2;
      }
    }
    //assign now outside flag based on helpOutsideFlag
    //Assign outside flag based on helpOutsideFlag values
    for(int k = 0; k<6;k++){
      if(helpOutsideFlag.at(k+1) == 1){
	outsideFlag = 1;
	break;
      }
      else if(helpOutsideFlag.at(k+1) == 2){
	outsideFlag = 2;
      }
    }      
    


    if ( (outsideFlag == 1 || outsideFlag == 2) && (periodicityFlag.at(1) == 1 || periodicityFlag.at(2) == 1 || periodicityFlag.at(3) == 1) ) {
      //Could cause problems with older lattice approaches, which would then need to be resolved.
      //Allow case with nodes on boundary. 
      //Check if nodes are too far outside and shift them to the boundaries of periodic neighbours
      //This is needed to get the octree localizer working
      updateFlag = 0;
      for ( int k = 0; k < 3; k++ ) {
	if ( ( coords.at(k + 1) + newTol ) < boundaries.at(2 * k + 1) - specimenDimension.at(k + 1) ) {
	  coords.at(k + 1) = boundaries.at(2 * k + 1) - specimenDimension.at(k + 1);
	  updateFlag = 1;
	} else if ( ( coords.at(k + 1) - newTol ) > boundaries.at(2 * k + 2) + specimenDimension.at(k + 1) ) {
	  coords.at(k + 1) = boundaries.at(2 * k + 2) + specimenDimension.at(k + 1);
	}
	updateFlag = 1;
      }
      if ( updateFlag == 1 ) {
	this->grid->giveVoronoiVertex(i + 1)->setCoordinates(coords);
      }
    }

    this->grid->giveVoronoiVertex(i + 1)->setOutsideFlag(outsideFlag);
    this->grid->giveVoronoiVertex(i + 1)->setHelpOutsideFlag(helpOutsideFlag);
  }
    
  //Reinforcements nodes
  for ( int i = 0; i < this->grid->giveNumberOfReinforcementNode(); i++ ) {
    outsideFlag = 0;
    this->grid->giveReinforcementNode(i + 1)->giveCoordinates(coords);
    for ( int k = 0; k < 3; k++ ) {
      if ( ( coords.at(k + 1) + newTol ) < boundaries.at(2 * k + 1) ) {
	outsideFlag = 1;
	break;
      } else if ( fabs( coords.at(k + 1) - boundaries.at(2 * k + 1) ) < newTol ) {                 //On the boundary
	outsideFlag = 2;
      } else if ( ( coords.at(k + 1) - newTol ) > boundaries.at(2 * k + 2) ) {
	outsideFlag = 1;
	break;
      } else if ( fabs( coords.at(k + 1) - boundaries.at(2 * k + 2) ) < newTol ) {
	outsideFlag = 2;
      }
    }
    this->grid->giveReinforcementNode(i + 1)->setOutsideFlag(outsideFlag);
  }

  int nodeFlag1=0,nodeFlag2=0,nodeFlag3=0,nodeFlag4=0;
  if(grid->giveRegularFlag() == 2){
    //Tetrahedral elements
    for ( int i = 0; i < this->grid->giveNumberOfDelaunayTetras(); i++ ) {
      outsideFlag = 0;
      this->grid->giveDelaunayTetra(i + 1)->giveLocalVertices(nodes);
      nodeFlag1 = this->grid->giveDelaunayVertex( nodes.at(1) )->giveOutsideFlag();
      nodeFlag2 = this->grid->giveDelaunayVertex( nodes.at(2) )->giveOutsideFlag();
      nodeFlag3 = this->grid->giveDelaunayVertex( nodes.at(3) )->giveOutsideFlag();
      nodeFlag4 = this->grid->giveDelaunayVertex( nodes.at(4) )->giveOutsideFlag();

      if( (nodeFlag1 == 0 || nodeFlag1 == 2)  &&  (nodeFlag2 == 0 || nodeFlag2 == 2) && (nodeFlag3 == 0 || nodeFlag3 == 2) && (nodeFlag4 == 0 || nodeFlag4 == 2) ){
	outsideFlag = 0;
      }
      else if( nodeFlag1 == 0 || nodeFlag2 == 0 || nodeFlag3 == 0 || nodeFlag4 == 0){
	outsideFlag = 2;
      }
      else{
	outsideFlag = 1;
      }
      this->grid->giveDelaunayTetra(i + 1)->setOutsideFlag(outsideFlag);
    }
  }

  int boundaryNumber = -1;
  oofem::IntArray helpNodeFlag1(6);
  oofem::IntArray helpNodeFlag2(6);
  if(grid->giveRegularFlag() !=2){
    //Delaunay elements     
    for ( int i = 0; i < this->grid->giveNumberOfDelaunayLines(); i++ ) {
      outsideFlag = 0;
      this->grid->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
      nodeFlag1 = this->grid->giveDelaunayVertex( nodes.at(1) )->giveOutsideFlag();
      nodeFlag2 = this->grid->giveDelaunayVertex( nodes.at(2) )->giveOutsideFlag();
      this->grid->giveDelaunayVertex( nodes.at(1) )->giveHelpOutsideFlag(helpNodeFlag1);
      this->grid->giveDelaunayVertex( nodes.at(2) )->giveHelpOutsideFlag(helpNodeFlag2);
      if ( ( nodeFlag1 == 1 && (nodeFlag2 == 0 || nodeFlag2 == 2) ) ||
	   ((nodeFlag1  == 0 || nodeFlag1  == 2) && nodeFlag2 == 1 )) {

	//Check that the node that is outside is on the same boundary as the inside node
	outsideFlag = 2;
	if(nodeFlag1 == 2 || nodeFlag2 == 2){
	  for(int k=0;k <3;k++){
	    if((helpNodeFlag1.at(2*k+1) == 1 || helpNodeFlag1.at(2*k+2) == 1 || helpNodeFlag2.at(2*k+1) == 1 || helpNodeFlag2.at(2*k+2) == 1 ) && periodicityFlag.at(k+1) != 1){
	      outsideFlag = 1;
	    }
	  }
	}

	if(outsideFlag == 2){      
	  //One node is outside and one inside.
	  /**Store in the inside node info which boundary is crossed. */
	
	  for ( int k = 0; k < 2; k++ ) {
	    if (this->grid->giveDelaunayVertex( nodes.at(k + 1) )->giveOutsideFlag() == 1 ) {
	    this->grid->giveDelaunayVertex( nodes.at(k + 1) )->giveCoordinates(coords);
	    for ( int m = 0; m < 3; m++ ) {
	      if ( ( coords.at(m + 1) + newTol ) < boundaries.at(2 * m + 1) ) {
		boundaryNumber = 2 * m + 1;
		break;
	      } else if ( ( coords.at(m + 1) - newTol ) > boundaries.at(2 * m + 2) ) {
		boundaryNumber = 2 * m + 2;
		break;
	      }
	    }
	  }
	}
	
	for ( int k = 0; k < 2; k++ ) {
	  if (this->grid->giveDelaunayVertex( nodes.at(k + 1) )->giveOutsideFlag() == 0 && boundaryNumber != -1 ) {
	    this->grid->giveDelaunayVertex( nodes.at(k + 1) )->setBoundaryFlag(boundaryNumber);
	  }
	}
	}
      } else if (nodeFlag1 == 1 || nodeFlag2 == 1 ) {
	outsideFlag = 1;
      } else if (nodeFlag1 == 2 && nodeFlag2 == 2 ) {
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
    
    //Beam elements
    boundaryNumber = -1;
    for ( int i = 0; i < this->grid->giveNumberOfLatticeBeams(); i++ ) {
      outsideFlag = 0;
      this->grid->giveLatticeBeam(i + 1)->giveLocalVertices(nodes);
      if ( (this->grid->giveReinforcementNode( nodes.at(1) )->giveOutsideFlag() == 1 &&
	    this->grid->giveReinforcementNode( nodes.at(2) )->giveOutsideFlag() == 0 ) ||
	   (this->grid->giveReinforcementNode( nodes.at(1) )->giveOutsideFlag() == 0 &&
	    this->grid->giveReinforcementNode( nodes.at(2) )->giveOutsideFlag() == 1 ) ||
	   (this->grid->giveReinforcementNode( nodes.at(1) )->giveOutsideFlag() == 0 &&
            this->grid->giveReinforcementNode( nodes.at(2) )->giveOutsideFlag() == 2 ) ) {
	outsideFlag = 2;
          
	//One node is outside and one inside.
	/**Store in the inside node info which boundary is crossed. */
                
	for ( int k = 0; k < 2; k++ ) {
	  if (this->grid->giveReinforcementNode( nodes.at(k + 1) )->giveOutsideFlag() == 1 ) {
	    this->grid->giveReinforcementNode( nodes.at(k + 1) )->giveCoordinates(coords);
	    for ( int m = 0; m < 3; m++ ) {
	      if ( ( coords.at(m + 1) + newTol ) < boundaries.at(2 * m + 1) ) {
		boundaryNumber = 2 * m + 1;
		break;
	      } else if ( ( coords.at(m + 1) - newTol ) > boundaries.at(2 * m + 2) ) {
		boundaryNumber = 2 * m + 2;
		break;
	      }
	    }
	  }
	}
                
	//One node is inside and one is on the surface.

	for ( int k = 0; k < 2; k++ ) {
	  if (this->grid->giveReinforcementNode( nodes.at(k + 1) )->giveOutsideFlag() == 2 ) {
	    this->grid->giveReinforcementNode( nodes.at(k + 1) )->giveCoordinates(coords);
	    for ( int m = 0; m < 3; m++ ) {
	      if ( fabs( coords.at(m + 1) - boundaries.at(2 * m + 1) ) < newTol  ) {
		boundaryNumber = 2 * m + 1;
		break;
	      } else if ( fabs( coords.at(m + 1) - boundaries.at(2 * m + 2) ) < newTol  ) {
		boundaryNumber = 2 * m + 2;
		break;
	      }
	    }
	  }
	}

	for ( int k = 0; k < 2; k++ ) {
	  if (this->grid->giveReinforcementNode( nodes.at(k + 1) )->giveOutsideFlag() == 0 && boundaryNumber != -1 ) {
	    this->grid->giveReinforcementNode( nodes.at(k + 1) )->setBoundaryFlag(boundaryNumber);
	  }
	}
      } else if (this->grid->giveReinforcementNode( nodes.at(1) )->giveOutsideFlag() == 1 ||
		 this->grid->giveReinforcementNode( nodes.at(2) )->giveOutsideFlag() == 1 ) {
	outsideFlag = 1;
      } else if (this->grid->giveReinforcementNode( nodes.at(1) )->giveOutsideFlag() == 2 &&
		 this->grid->giveReinforcementNode( nodes.at(2) )->giveOutsideFlag() == 2 ) {
	outsideFlag = 3;
      }
      this->grid->giveLatticeBeam(i + 1)->setOutsideFlag(outsideFlag);
    }
    
    //link elements
    boundaryNumber = -1;
    for ( int i = 0; i < this->grid->giveNumberOfLatticeLinks(); i++ ) {
      outsideFlag = 0;
      this->grid->giveLatticeLink(i + 1)->giveLocalVertices(nodes);
      if ( (this->grid->giveReinforcementNode( nodes.at(1) )->giveOutsideFlag() == 1 &&
	    this->grid->giveDelaunayVertex( nodes.at(2) )->giveOutsideFlag() == 0 ) ||
	   (this->grid->giveReinforcementNode( nodes.at(1) )->giveOutsideFlag() == 0 &&
            this->grid->giveDelaunayVertex( nodes.at(2) )->giveOutsideFlag() == 1 ) ) {
	outsideFlag = 2;
                
	//One node is outside and one inside.
	/**Store in the inside node info which boundary is crossed. */
                                
	if (this->grid->giveReinforcementNode( nodes.at(1) )->giveOutsideFlag() == 1 ) {
	  this->grid->giveReinforcementNode( nodes.at(1) )->giveCoordinates(coords);
	  for ( int m = 0; m < 3; m++ ) {
	    if ( ( coords.at(m + 1) + newTol ) < boundaries.at(2 * m + 1) ) {
	      boundaryNumber = 2 * m + 1;
	      break;
	    } else if ( ( coords.at(m + 1) - newTol ) > boundaries.at(2 * m + 2) ) {
	      boundaryNumber = 2 * m + 2;
	      break;
	    }
	  }
	}
                
	if (this->grid->giveDelaunayVertex( nodes.at(2) )->giveOutsideFlag() == 1 ) {
	  this->grid->giveDelaunayVertex( nodes.at(2) )->giveCoordinates(coords);
	  for ( int m = 0; m < 3; m++ ) {
	    if ( ( coords.at(m + 1) + newTol ) < boundaries.at(2 * m + 1) ) {
	      boundaryNumber = 2 * m + 1;
	      break;
	    } else if ( ( coords.at(m + 1) - newTol ) > boundaries.at(2 * m + 2) ) {
	      boundaryNumber = 2 * m + 2;
	      break;
	    }
	  }
	}
                
	if (this->grid->giveReinforcementNode( nodes.at( 1) )->giveOutsideFlag() == 0 && boundaryNumber != -1 ) {
	  this->grid->giveReinforcementNode( nodes.at(1) )->setBoundaryFlag(boundaryNumber);
	}
                
	if (this->grid->giveDelaunayVertex( nodes.at(2) )->giveOutsideFlag() == 0 && boundaryNumber != -1 ) {
	  this->grid->giveDelaunayVertex( nodes.at(2) )->setBoundaryFlag(boundaryNumber);
	}
                
      } else if (grid->giveReinforcementNode( nodes.at(1) )->giveOutsideFlag() == 1 ||
		 this->grid->giveDelaunayVertex( nodes.at(2) )->giveOutsideFlag() == 1 ) {
	outsideFlag = 1;
      } else if (this->grid->giveReinforcementNode( nodes.at(1) )->giveOutsideFlag() == 2 &&
		 this->grid->giveDelaunayVertex( nodes.at(2) )->giveOutsideFlag() == 2 ) {
	outsideFlag = 3;
      }
      this->grid->giveLatticeLink(i + 1)->setOutsideFlag(outsideFlag);
    }
  }
    
  //Modify cross-sections of elements. This should fix the problem in the subsequent part. However, for elements on the boundary whichc cross one other boundary, this does not fix it since those elements are given a outside flag = 2
  for ( int i = 0; i < this->grid->giveNumberOfDelaunayLines(); i++ ) {
    if ( this->grid->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3) {
      this->modifyVoronoiCrossSection(i + 1);
    }
    // else if(this->grid->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2){
    //   //Check if this is one of the elemnets that this on the boundary    
    // }
  }
    
  //Find periodic nodes and elements.
  nodeContainerType nodeSet;
  nodeContainerType :: iterator pos;
  oofem::IntArray switches(3), crossSectionElements, newNodes(2), lineCandidates;
  oofem::FloatArray newCoords(3);
  int location, periodicNode, flag;
  oofem::FloatArray coordsInside(3), mirrorCoords(3);
  if ( periodicityFlag.at(1) == 1 || periodicityFlag.at(2) == 1 || periodicityFlag.at(3) == 1) {
    //Delaunay Lines
    for ( int i = 0; i < this->grid->giveNumberOfDelaunayLines(); i++ ) {
      if (this->grid->giveDelaunayLine(i + 1)->giveOutsideFlag() != 1 ) {
	if (this->grid->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
	  this->grid->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
	  location = 0;
	  for ( int m = 0; m < 2; m++ ) {
	    if (this->grid->giveDelaunayVertex( nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
	      this->grid->giveDelaunayVertex( nodes.at(m + 1) )->giveCoordinates(coords);
	      location = this->giveSwitches(switches, coords);
	      //new coords
	      for ( int n = 0; n < 3; n++ ) {
		newCoords.at(n + 1) = coords.at(n + 1) - switches.at(n + 1) * specimenDimension.at(n + 1);
	      }
	      nodeSet.clear();
	      grid->delaunayLocalizer->giveAllNodesWithinBox(nodeSet, newCoords,this->grid->giveTol(), 0);
	      if ( nodeSet.size() == 1 ) { //Found nodes
		periodicNode = * nodeSet.begin();
	      }else {
		newCoords.printYourself();
		converter::errorf( "Could not find periodic node for delaunay elements. Node set is %d", nodeSet.size() );
	      }
	      this->grid->giveDelaunayVertex( nodes.at(m + 1) )->setPeriodicNode(periodicNode);
	      this->grid->giveDelaunayVertex( nodes.at(m + 1) )->setLocation(location);
	      //used for the Delaunay transport model
	      if ( m == 0 ) {
		this->grid->giveDelaunayVertex( nodes.at(2) )->giveCoordinates(coordsInside);
	      } else {
		this->grid->giveDelaunayVertex( nodes.at(1) )->giveCoordinates(coordsInside);
	      }
	      this->grid->giveDelaunayVertex(periodicNode)->giveLocalLines(lineCandidates);
	      oofem::IntArray mirrorNodes(2);
	      for ( int k = 0; k < lineCandidates.giveSize(); k++ ) {
		this->grid->giveDelaunayLine( lineCandidates.at(k + 1) )->giveLocalVertices(mirrorNodes);
		for ( int z = 0; z < 2; z++ ) {
		  this->grid->giveDelaunayVertex( mirrorNodes.at(z + 1) )->giveCoordinates(mirrorCoords);
		  location = this->giveSwitches(switches, mirrorCoords);
		  for ( int n = 0; n < 3; n++ ) {
		    mirrorCoords.at(n + 1) = mirrorCoords.at(n + 1) - switches.at(n + 1) * specimenDimension.at(n + 1);
		  }
		  if ( fabs( mirrorCoords.at(1) - coordsInside.at(1) ) <this->grid->giveTol() &&
		       fabs( mirrorCoords.at(2) - coordsInside.at(2) ) <this->grid->giveTol() &&
		       fabs( mirrorCoords.at(3) - coordsInside.at(3) ) <this->grid->giveTol() ) {
		    //Found mirror element      
		    this->grid->giveDelaunayLine(i + 1)->setPeriodicElement( lineCandidates.at(k + 1) );
		    break;
		  }
		}
	      }//until here
	    }
	  }
	}
	    
		
	//Check if one is outside and mirror. Peter: I do not understand what the point about this is.
	if(periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(2) == 1){ 
	  this->grid->giveDelaunayLine(i + 1)->giveCrossSectionElements(crossSectionElements);
	  for ( int k = 0; k < crossSectionElements.giveSize(); k++ ) {
	    if ( this->grid->giveVoronoiLine( crossSectionElements.at(k + 1) )->giveOutsideFlag() == 1 ) {
	      this->grid->giveVoronoiLine( crossSectionElements.at(k + 1) )->giveLocalVertices(nodes);
		      
	      //Loop to check if the cross-section line nodes are on different locations
	      //If they do, one of the switches array from the nodes is used to shift the line
	      for ( int m = 0; m < 2; m++ ) {
		this->grid->giveVoronoiVertex( nodes.at(m + 1) )->giveCoordinates(coords);
		locationArray.at(m + 1) = this->giveSwitches(switches, coords);
	      }

	      for ( int m = 0; m < 2; m++ ) {
		this->grid->giveVoronoiVertex( nodes.at(m + 1) )->giveCoordinates(coords);

		if ( locationArray.at(1) == locationArray.at(2) || m == 0 ) {
		  location = this->giveSwitches(switches, coords);
		}

		//new coords
		for ( int n = 0; n < 3; n++ ) {
		  newCoords.at(n + 1) = coords.at(n + 1) - switches.at(n + 1) * specimenDimension.at(n + 1);
		}

		nodeSet.clear();
		grid->voronoiLocalizer->giveAllNodesWithinBox(nodeSet, newCoords,this->grid->giveTol(), 1);
		if ( nodeSet.size() == 1 ) { //Found nodes
		  newNodes.at(m + 1) = * nodeSet.begin();
		} else {
		  newCoords.printYourself();
		  newNodes.at(1) = 0;
		  converter::errorf( "Could not find periodic node of cross-section elements of delaunay elements. Node set is %d", nodeSet.size() );
		  break;
		}
	      }
	      if(newNodes.at(1) != 0){
		this->grid->giveVoronoiVertex( newNodes.at(1) )->giveLocalLines(lineCandidates);

		for ( int m = 0; m < lineCandidates.giveSize(); m++ ) {
		  this->grid->giveVoronoiLine( lineCandidates.at(m + 1) )->giveLocalVertices(nodes);
		  if ( ( nodes.at(1) == newNodes.at(1) && nodes.at(2) == newNodes.at(2) ) ||
		       ( nodes.at(1) == newNodes.at(2) && nodes.at(2) == newNodes.at(1) ) ) {
		    this->grid->giveVoronoiLine( crossSectionElements.at(k + 1) )->setPeriodicElement( lineCandidates.at(m + 1) );
		    break;
		  }
		}
	      }
	    }
	  } //End of cross-section check
	}
      }
    }

    printf("Finished finding periodic info for Delaunay elements\n");
	
    //Beam elements
    periodicNode = 0;
    for ( int i = 0; i < this->grid->giveNumberOfLatticeBeams(); i++ ) {
      if (this->grid->giveLatticeBeam(i + 1)->giveOutsideFlag() != 1 ) {
	if (this->grid->giveLatticeBeam(i + 1)->giveOutsideFlag() == 2 ) {
	  this->grid->giveLatticeBeam(i + 1)->giveLocalVertices(nodes);
	  location = 0;
	  for ( int m = 0; m < 2; m++ ) {
	    if (this->grid->giveReinforcementNode( nodes.at(m + 1) )->giveOutsideFlag() == 1 ||
		this->grid->giveReinforcementNode( nodes.at(m + 1) )->giveOutsideFlag() == 2  ) {
	      this->grid->giveReinforcementNode( nodes.at(m + 1) )->giveCoordinates(coords);
	      location = this->giveSwitches(switches, coords);
	      //new coords
	      for ( int n = 0; n < 3; n++ ) {
		newCoords.at(n + 1) = coords.at(n + 1) - switches.at(n + 1) * specimenDimension.at(n + 1);
	      }
            
	      periodicNode = (Fibre::reinforcementLocalizer(newCoords,grid, 2*10*grid->giveTol()))->giveNumber();
	      this->grid->giveReinforcementNode( nodes.at(m + 1) )->setPeriodicNode(periodicNode);
	      this->grid->giveReinforcementNode( nodes.at(m + 1) )->setLocation(location);
	      //used for the Delaunay transport model
	      if ( m == 0 ) {
		this->grid->giveReinforcementNode( nodes.at(2) )->giveCoordinates(coordsInside);
	      } else {
		this->grid->giveReinforcementNode( nodes.at(1) )->giveCoordinates(coordsInside);
	      }
                            
	      this->grid->giveReinforcementNode(periodicNode)->giveLocalLines(lineCandidates);
	      oofem::IntArray mirrorNodes(2);
	      for ( int k = 0; k < lineCandidates.giveSize(); k++ ) {
		this->grid->giveLatticeBeam( lineCandidates.at(k + 1) )->giveLocalVertices(mirrorNodes);
		for ( int z = 0; z < 2; z++ ) {
		  this->grid->giveReinforcementNode( mirrorNodes.at(z + 1) )->giveCoordinates(mirrorCoords);
		  location = this->giveSwitches(switches, mirrorCoords);
		  for ( int n = 0; n < 3; n++ ) {
		    mirrorCoords.at(n + 1) = mirrorCoords.at(n + 1) - switches.at(n + 1) * specimenDimension.at(n + 1);
		  }
		  if ( fabs( mirrorCoords.at(1) - coordsInside.at(1) ) <this->grid->giveTol() &&
		       fabs( mirrorCoords.at(2) - coordsInside.at(2) ) <this->grid->giveTol() &&
		       fabs( mirrorCoords.at(3) - coordsInside.at(3) ) <this->grid->giveTol() ) {
		    //Found mirror element
		    this->grid->giveLatticeBeam(i + 1)->setPeriodicElement( lineCandidates.at(k + 1) );
		    break;
		  }
		}
	      }//until here
	    }
	  }
	}
                
      }
    }

    printf("Finished finding periodic info for Beam elements\n");
        
                
    //Link elements
    for ( int i = 0; i < this->grid->giveNumberOfLatticeLinks(); i++ ) {
      if (this->grid->giveLatticeLink(i + 1)->giveOutsideFlag() != 1 ) {
	if (this->grid->giveLatticeLink(i + 1)->giveOutsideFlag() == 2 ) {
	  this->grid->giveLatticeLink(i + 1)->giveLocalVertices(nodes);
	  location = 0;
                    
	  int m=0;
	  if (this->grid->giveReinforcementNode( nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
	    this->grid->giveReinforcementNode( nodes.at(m + 1) )->giveCoordinates(coords);
	    location = this->giveSwitches(switches, coords);
	    //new coords
	    for ( int n = 0; n < 3; n++ ) {
	      newCoords.at(n + 1) = coords.at(n + 1) - switches.at(n + 1) * specimenDimension.at(n + 1);
	    }

                        
	    periodicNode = (Fibre::reinforcementLocalizer(newCoords,grid,2*10*grid->giveTol()))->giveNumber();
	    this->grid->giveReinforcementNode( nodes.at(m + 1) )->setPeriodicNode(periodicNode);
	    this->grid->giveReinforcementNode( nodes.at(m + 1) )->setLocation(location);
	    //used for the Delaunay transport model
                       
	    this->grid->giveDelaunayVertex( nodes.at(2) )->giveCoordinates(coordsInside);
                        
                        
	    this->grid->giveReinforcementNode(periodicNode)->giveLocalLinks(lineCandidates);
	    oofem::IntArray mirrorNodes(2);
	    for ( int k = 0; k < lineCandidates.giveSize(); k++ ) {
	      this->grid->giveLatticeLink( lineCandidates.at(k + 1) )->giveLocalVertices(mirrorNodes);
	      for ( int z = 0; z < 2; z++ ) {
		if (z==0) {
		  this->grid->giveReinforcementNode( mirrorNodes.at(z + 1) )->giveCoordinates(mirrorCoords);}
		else {grid->giveDelaunayVertex( mirrorNodes.at(z + 1) )->giveCoordinates(mirrorCoords);}
		location = this->giveSwitches(switches, mirrorCoords);
		for ( int n = 0; n < 3; n++ ) {
		  mirrorCoords.at(n + 1) = mirrorCoords.at(n + 1) - switches.at(n + 1) * specimenDimension.at(n + 1);
		}
		if ( fabs( mirrorCoords.at(1) - coordsInside.at(1) ) <this->grid->giveTol() &&
		     fabs( mirrorCoords.at(2) - coordsInside.at(2) ) <this->grid->giveTol() &&
		     fabs( mirrorCoords.at(3) - coordsInside.at(3) ) <this->grid->giveTol() ) {
		  //Found mirror element
		  this->grid->giveLatticeLink(i + 1)->setPeriodicElement( lineCandidates.at(k + 1) );
		  break;
		}
	      }
	    }//until here
	  }
                    
	  m=1;
	  if (this->grid->giveDelaunayVertex( nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
	    this->grid->giveDelaunayVertex( nodes.at(m + 1) )->giveCoordinates(coords);
	    location = this->giveSwitches(switches, coords);
	    //new coords
	    for ( int n = 0; n < 3; n++ ) {
	      newCoords.at(n + 1) = coords.at(n + 1) - switches.at(n + 1) * specimenDimension.at(n + 1);
	    }
	    nodeSet.clear();
	    grid->delaunayLocalizer->giveAllNodesWithinBox(nodeSet, newCoords,this->grid->giveTol(), 0);
	    if ( nodeSet.size() == 1 ) { //Found nodes
	      periodicNode = * nodeSet.begin();
	    } else {
	      newCoords.printYourself();
	      converter::errorf( "Could not find periodic node for link elments. Node set is %d", nodeSet.size() );
	    }
	    this->grid->giveDelaunayVertex( nodes.at(m + 1) )->setPeriodicNode(periodicNode);
	    this->grid->giveDelaunayVertex( nodes.at(m + 1) )->setLocation(location);
	    //used for the Delaunay transport model
                        
	    this->grid->giveReinforcementNode( nodes.at(1) )->giveCoordinates(coordsInside);
                        
                        
	    this->grid->giveDelaunayVertex(periodicNode)->giveLocalLinks(lineCandidates);
	    oofem::IntArray mirrorNodes(2);
	    for ( int k = 0; k < lineCandidates.giveSize(); k++ ) {
	      this->grid->giveLatticeLink( lineCandidates.at(k + 1) )->giveLocalVertices(mirrorNodes);
	      for ( int z = 0; z < 2; z++ ) {
		if (z==0) {
		  this->grid->giveReinforcementNode( mirrorNodes.at(z + 1) )->giveCoordinates(mirrorCoords);}
		else {grid->giveDelaunayVertex( mirrorNodes.at(z + 1) )->giveCoordinates(mirrorCoords);}
		location = this->giveSwitches(switches, mirrorCoords);
		for ( int n = 0; n < 3; n++ ) {
		  mirrorCoords.at(n + 1) = mirrorCoords.at(n + 1) - switches.at(n + 1) * specimenDimension.at(n + 1);
		}
		if ( fabs( mirrorCoords.at(1) - coordsInside.at(1) ) <this->grid->giveTol() &&
		     fabs( mirrorCoords.at(2) - coordsInside.at(2) ) <this->grid->giveTol() &&
		     fabs( mirrorCoords.at(3) - coordsInside.at(3) ) < this->grid->giveTol() ) {
		  //Found mirror element
		  this->grid->giveLatticeLink(i + 1)->setPeriodicElement( lineCandidates.at(k + 1) );
		  break;
		}
	      }
	    }//until here
	  }
                    
	}
                
      }
    }
        
    printf("Finished finding periodic info for link elements\n");

	
    //Find periodic info for tetrahedral elements. Assume that nodes are on the surface.
    //Extend this later to take into account the case of crossing surfaces.
    oofem::FloatArray coords(3);
	
    if(grid->giveRegularFlag() == 2){//Only implemented for nodes on surfaces
      for ( int i = 0; i < this->grid->giveNumberOfDelaunayTetras(); i++ ) {

	if(this->grid->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 0){
	
	  this->grid->giveDelaunayTetra(i + 1)->giveLocalVertices(nodes);
	    
	  for(int k = 0; k < nodes.giveSize(); k++){
	    location=0;
	    if(grid->giveDelaunayVertex(nodes.at(k+1))->giveOutsideFlag() == 2){//node on surface
	      grid->giveDelaunayVertex(nodes.at(k+1))->giveCoordinates(coords);
	      //compute number of coordinates at boundaries
	      int nBoundCoord=0;
	      for ( int l=1; l <=3; l++ ) {
		if ( fabs(coords.at(l) - boundaries.at(2*l-1)) < grid->giveTol() ) {
		  nBoundCoord++;
		} else if ( fabs(coords.at(l) - boundaries.at(2*l)) < grid->giveTol() ) {
		  nBoundCoord++;
		}
	      }

	      if ( nBoundCoord == 1 ) {
		//node located at face
		location = this->giveSwitches(switches, coords);
	      } else if ( nBoundCoord == 2 ) {
		//node located on edge
		location = this->giveEdgeSwitches(switches, coords);
	      } else if ( nBoundCoord == 3 ) {
		//corner
		location = this->giveCornerSwitches(switches, coords);
	      }

	      //new coords
	      for ( int n = 0; n < 3; n++ ) {
		newCoords.at(n + 1) = coords.at(n + 1) - switches.at(n + 1) * specimenDimension.at(n + 1);
	      }
	      nodeSet.clear();
	      grid->delaunayLocalizer->giveAllNodesWithinBox(nodeSet, newCoords,this->grid->giveTol(), 0);
	      if ( nodeSet.size() == 1 ) { //Found nodes
		periodicNode = * nodeSet.begin();
	      } else {
		converter::errorf( "Could not find periodic node. Node set is %d", nodeSet.size() );
	      }
	      this->grid->giveDelaunayVertex( nodes.at(k + 1) )->setPeriodicNode(periodicNode);
	      this->grid->giveDelaunayVertex( nodes.at(k + 1) )->setLocation(location);
	    }
	  }
	}
      }
    }//end tetra
	
	
	
    //Voronoi Lines
    if(grid->giveRegularFlag() !=2 && (periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1) ){//Does not work if nodes are on surface	
      for ( int i = 0; i < this->grid->giveNumberOfVoronoiLines(); i++ ) {
	if ( this->grid->giveVoronoiLine(i + 1)->giveOutsideFlag() != 1 ) {
	  //Do this for all elements which are not completely outside
	  if ( this->grid->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {
	    this->grid->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
	    location = 0;
	    for ( int m = 0; m < 2; m++ ) {
	      if ( this->grid->giveVoronoiVertex( nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
		this->grid->giveVoronoiVertex( nodes.at(m + 1) )->giveCoordinates(coords);
		location = this->giveSwitches(switches, coords);
		//new coords
		for ( int n = 0; n < 3; n++ ) {
		  newCoords.at(n + 1) = coords.at(n + 1) - switches.at(n + 1) * specimenDimension.at(n + 1);
		}
		nodeSet.clear();
		grid->voronoiLocalizer->giveAllNodesWithinBox(nodeSet, newCoords,this->grid->giveTol(), 1);
		if ( nodeSet.size() == 1 ) { //Found nodes
		  periodicNode = * nodeSet.begin();
		} else {
		  converter::errorf( "Could not find periodic node for voronoi lines. Node set is %d", nodeSet.size() );
		}
		this->grid->giveVoronoiVertex( nodes.at(m + 1) )->setPeriodicNode(periodicNode);
		this->grid->giveVoronoiVertex( nodes.at(m + 1) )->setLocation(location);

		//Find the periodic voronoi line crossing the boundaries
		//This information is used for periodic random fields

		if ( m == 0 ) {
		  this->grid->giveVoronoiVertex( nodes.at(2) )->giveCoordinates(coordsInside);
		} else {
		  this->grid->giveVoronoiVertex( nodes.at(1) )->giveCoordinates(coordsInside);
		}

		this->grid->giveVoronoiVertex(periodicNode)->giveLocalLines(lineCandidates);
		oofem::IntArray mirrorNodes(2);
		for ( int k = 0; k < lineCandidates.giveSize(); k++ ) {
		  this->grid->giveVoronoiLine( lineCandidates.at(k + 1) )->giveLocalVertices(mirrorNodes);
		  for ( int z = 0; z < 2; z++ ) {
		    this->grid->giveVoronoiVertex( mirrorNodes.at(z + 1) )->giveCoordinates(mirrorCoords);
		    location = this->giveSwitches(switches, mirrorCoords);
		    for ( int n = 0; n < 3; n++ ) {
		      mirrorCoords.at(n + 1) = mirrorCoords.at(n + 1) - switches.at(n + 1) * specimenDimension.at(n + 1);
		    }
		    if ( fabs( mirrorCoords.at(1) - coordsInside.at(1) ) <this->grid->giveTol() &&
			 fabs( mirrorCoords.at(2) - coordsInside.at(2) ) <this->grid->giveTol() &&
			 fabs( mirrorCoords.at(3) - coordsInside.at(3) ) <this->grid->giveTol() ) {
		      //Found mirror element
		      this->grid->giveVoronoiLine(i + 1)->setPeriodicElement( lineCandidates.at(k + 1) );
		      break;
		    }
		  }
		}
	      }
	    }
	  }

	  //Check if crosssection elements are outside and mirror
	  this->grid->giveVoronoiLine(i + 1)->giveCrossSectionElements(crossSectionElements);
	  for ( int k = 0; k < crossSectionElements.giveSize(); k++ ) {
	    if ( this->grid->giveDelaunayLine( crossSectionElements.at(k + 1) )->giveOutsideFlag() == 1 ) {
	      this->grid->giveDelaunayLine( crossSectionElements.at(k + 1) )->giveLocalVertices(nodes);

	      //Loop to check if the cross-section line nodes are on different locations
	      //If they do, one of the switches array from the nodes is used to shift the line
	      for ( int m = 0; m < 2; m++ ) {
		this->grid->giveDelaunayVertex( nodes.at(m + 1) )->giveCoordinates(coords);
		locationArray.at(m + 1) = this->giveSwitches(switches, coords);
	      }

	      for ( int m = 0; m < 2; m++ ) {
		this->grid->giveDelaunayVertex( nodes.at(m + 1) )->giveCoordinates(coords);


		if ( locationArray.at(1) == locationArray.at(2) || m == 0 ) {
		  location = this->giveSwitches(switches, coords);
		}

		//new coords
		for ( int n = 0; n < 3; n++ ) {
		  newCoords.at(n + 1) = coords.at(n + 1) - switches.at(n + 1) * specimenDimension.at(n + 1);
		}
		nodeSet.clear();
		grid->delaunayLocalizer->giveAllNodesWithinBox(nodeSet, newCoords,this->grid->giveTol(), 0);
		if ( nodeSet.size() == 1 ) { //Found nodes
		  newNodes.at(m + 1) = * nodeSet.begin();
		} else {
		  converter::errorf( "Could not find periodic node. Node set is %d", nodeSet.size() );
		}
	      }
	      this->grid->giveDelaunayVertex( newNodes.at(1) )->giveLocalLines(lineCandidates);
	      for ( int m = 0; m < lineCandidates.giveSize(); m++ ) {
		this->grid->giveDelaunayLine( lineCandidates.at(m + 1) )->giveLocalVertices(nodes);
		if ( ( nodes.at(1) == newNodes.at(1) && nodes.at(2) == newNodes.at(2) ) ||
		     ( nodes.at(1) == newNodes.at(2) && nodes.at(2) == newNodes.at(1) ) ) {
		  //Found mirror element
		  this->grid->giveDelaunayLine( crossSectionElements.at(k + 1) )->setPeriodicElement( lineCandidates.at(m + 1) );
		  break;
		}
	      }
	    }
	  } //End of cross-section check
	}
      }
	
    }//end of check of regular flag
    printf("Finished finding periodic info for voronoi elements\n");
  }  

	
  return;
  }

void Prism :: defineBoundaries(oofem::FloatArray &boundaries)
//Determine the boundaries of the grid
{
  boundaries.resize(6);
  
  boundaries.at(1) = this->box.at(1);
  boundaries.at(2) = this->box.at(4);
  boundaries.at(3) = this->box.at(2);
  boundaries.at(4) = this->box.at(5);
  boundaries.at(5) = this->box.at(3);
  boundaries.at(6) = this->box.at(6);
  
  return;
}

