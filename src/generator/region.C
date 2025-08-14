#include "generatorlistutils.h"
#include "generatorerror.h"
#include "region.h"
#include "curve.h"
#include "vertex.h"
#include "surface.h"
#include "intarray.h"


#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include "octreegridlocalizer.h"
 #include <iostream>
#endif

Region :: Region(int n, Grid *aGrid) : GridComponent(n, aGrid)
{
    this->number = n;
}

Region :: ~Region()
// Destructor.
{

}


int
Region :: giveLocalSurface(int i)
// Returns the i-th coordinate of the receiver.
{
    if ( i > surfaces.giveSize() ) {
        return 0.;
    }

    return surfaces.at(i);
}


int Region :: generatePoints()
{


    oofem::IntArray curves;

    oofem::FloatArray boundaries;
    this->defineBoundaries(boundaries);

    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);
    
    int randomIntegerOne = grid->giveRandomInteger() - 1;
    int randomIntegerTwo = grid->giveRandomInteger() - 2;
    int randomIntegerThree = grid->giveRandomInteger() - 3;

    oofem::FloatArray random(3);
    int flag;



    oofem::IntArray periodicityFlag(3);
    grid->givePeriodicityFlag(periodicityFlag);
    
    double maxIter = grid->giveMaximumIterations();
    oofem::FloatArray newRandom(3);

    int randomFlag = grid->giveRandomFlag();

    double border=0.;
    if(randomFlag == 0){
      border = grid->TOL;
    }
    else{
      border = grid->diameter;
    }
    
    int vertexNumber = generator::size1(grid->vertexList);

    int tempSize = 1e9;
    generator::ensure_size1(grid->vertexList,tempSize);
    
    //Generate random vertices
    for ( int i = 0; i < maxIter; i++ ) {

      random.at(1) = boundaries.at(1) + border + grid->ran1(& randomIntegerOne) * ( boundaries.at(2) - boundaries.at(1) -  2.*border );
      random.at(2) = boundaries.at(3) + border + grid->ran1(& randomIntegerTwo) * ( boundaries.at(4) - boundaries.at(3) - 2.*border );
      random.at(3) = boundaries.at(5) + border + grid->ran1(& randomIntegerThree) * ( boundaries.at(6) - boundaries.at(5) - 2*border );

      //Check if this is far enough from the others
        flag = 0;
        flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, this->refinement * grid->giveDiameter(random) );

        if ( flag == 0 ) {

	  /* auto *v = new Vertex(vertexNumber+1, grid); */
	  /* v->setCoordinates(random); */
	  /* grid->setVertex(vertexNumber+1, v); */
	  /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1, random); */
	  /* vertexNumber++; */

	  grid->addVertex(random);
	  
            /* vertex = ( Vertex * ) ( Vertex(vertexNumber + 1, grid).ofType() ); */
            /* vertex->setCoordinates(random); */
            /* grid->setVertex(vertexNumber + 1, vertex); */
            /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber + 1, random); */
            /* vertexNumber++; */

	    i = 0;

            //Do now the mirroring and shifting.
	    
            //Do the periodic shift
            for ( int x = -1; x < 2; x++ ) {
	      for ( int y = -1; y < 2; y++ ) {
		for ( int z = -1; z < 2; z++ ) {
		  if ( !( x == 0 && y == 0 && z == 0 ) ) {
		    if(periodicityFlag.at(1) == 1){//periodic shift
		      newRandom.at(1) = random.at(1) + x * specimenDimension.at(1);
		    }
		    else{
		      if( x == -1 ){//mirorroring
			newRandom.at(1) = random.at(1) - 2.*(random.at(1) - boundaries.at(1));
		      }
		      else if ( x == 0 ){
			newRandom.at(1) = random.at(1);
		      }
		      else if( x == 1 ){
			newRandom.at(1) = random.at(1) - 2.*(random.at(1) - boundaries.at(2));
		      }
		    }

		    if(periodicityFlag.at(2) == 1){//periodic shift
		      newRandom.at(2) = random.at(2) + y * specimenDimension.at(2);
		    }

		    else{
		      if( y == -1 ){//mirorroring
			newRandom.at(2) = random.at(2) - 2.*(random.at(2) - boundaries.at(3));
		      }
		      else if ( y == 0 ){
			newRandom.at(2) = random.at(2);
		      }
		      else if( y == 1 ){
			newRandom.at(2) = random.at(2) - 2.*(random.at(2) - boundaries.at(4));
		      }
		    }

		    if(periodicityFlag.at(3) == 1){//periodic shift
		      newRandom.at(3) = random.at(3) + z * specimenDimension.at(3);
		    }

		    else{
		      if( z == -1 ){//mirorroring
			newRandom.at(3) = random.at(3) - 2.*(random.at(3) - boundaries.at(5));
		      }
		      else if ( z == 0 ){
			newRandom.at(3) = random.at(3);
		      }
		      else if( z == 1 ){
			newRandom.at(3) = random.at(3) - 2.*(random.at(3) - boundaries.at(6));
		      }
		    }


		    /* auto *v = new Vertex(vertexNumber+1, grid); */
		    /* v->setCoordinates(newRandom); */
		    /* grid->setVertex(vertexNumber+1, v); */
		    /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1, newRandom); */
		    /* vertexNumber++; */

		    grid->addVertex(newRandom);
		    
		    /* vertex = ( Vertex * ) ( Vertex(vertexNumber + 1, grid).ofType() ); */
		    /* vertex->setCoordinates(newRandom); */
		    /* grid->setVertex(vertexNumber + 1, vertex); */
		    /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber + 1, newRandom); */
		    /* vertexNumber++; */


		    
		  }		  
		}
	      }
	    }
	}
    }
    generator::ensure_size1(grid->vertexList,vertexNumber);

    return 1;
}


void Region :: defineBoundaries(oofem::FloatArray &boundaries)
//Determine the boundaries of the domain
{
    Surface *surface;
    Curve *curve;
    oofem::IntArray localVertices;
    oofem::IntArray vertexFlag(2);

    int localSurface, localCurve, localVertex;
  
    oofem::IntArray surfaces;
    this->giveLocalSurfaces(surfaces);

    double x, y, z;

    boundaries.resize(6);
    boundaries.at(1) = 1.e16;
    boundaries.at(2) = -1.e16;
    boundaries.at(3) = 1.e16;
    boundaries.at(4) = -1.e16;
    boundaries.at(5) = 1.e16;
    boundaries.at(6) = -1.e16;
    for ( int m = 0; m < surfaces.giveSize(); m++ ) {
        localSurface = this->giveLocalSurface(m + 1);
        surface = grid->giveSurface(localSurface);
        for ( int i = 0; i < surface->giveNumberOfLocalCurves(); i++ ) {
            localCurve = surface->giveLocalCurve(i + 1);
            curve = grid->giveCurve(localCurve);

            for ( int k = 0; k < 2; k++ ) {
                localVertex = curve->giveLocalVertex(k + 1);
                x = ( grid->giveInputVertex(localVertex) )->giveCoordinate(1);
                y = ( grid->giveInputVertex(localVertex) )->giveCoordinate(2);
                z = ( grid->giveInputVertex(localVertex) )->giveCoordinate(3);

                //assign min and max values

                if ( x < boundaries.at(1) ) {
                    boundaries.at(1) = x;
                }
                if ( x > boundaries.at(2) ) {
                    boundaries.at(2) = x;
                }

                if ( y < boundaries.at(3) ) {
                    boundaries.at(3) = y;
                }
                if ( y > boundaries.at(4) ) {
                    boundaries.at(4) = y;
                }

                if ( z < boundaries.at(5) ) {
                    boundaries.at(5) = z;
                }
                if ( z > boundaries.at(6) ) {
                    boundaries.at(6) = z;
                }
            }
        }
    }

    return;
}


void
Region :: initializeFrom(GeneratorInputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{

    IR_GIVE_FIELD(ir, surfaces, _IFT_Region_surfaces); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, _IFT_Region_refine); // Macro

    return;
}



/* Region *Region :: ofType() */
/* // Returns a new DofManager, which has the same number than the receiver, */
/* // but belongs to aClass (Node, ElementSide,..). */
/* { */
/*     Region *region; */

/*     region = new Region(number, grid); */

/*     return region; */
/* } */


int Region :: generateRegularPoints1()
{
    Curve *curve;

    oofem::IntArray curves;
    double x, y, z;

    oofem::FloatArray boundaries;
    grid->defineBoundaries(boundaries);
    oofem::FloatArray random(3);

    int n1edges = grid->xyzEdges.at(1);
    int n2edges = grid->xyzEdges.at(2);
    int n3edges = grid->xyzEdges.at(3);
    double n1length = fabs( boundaries.at(1) - boundaries.at(2) );
    double n2length = fabs( boundaries.at(3) - boundaries.at(4) );
    double n3length = fabs( boundaries.at(5) - boundaries.at(6) );

    int vertexNumber = generator::size1(grid->vertexList);
    int tempSize = 1e9;
    generator::ensure_size1(grid->vertexList,tempSize);

    for ( int i = 0; i < 2 * n3edges - 1; i++ ) {
        random.at(3) = boundaries.at(5) + ( i + 1 ) * n3length / ( 2 * ( double ) n3edges );

        //Generate precise points
        if ( i % 2 ) { // B types odd
            for ( int j = 0; j < n2edges - 1; j++ ) {
                for ( int l = 0; l < n1edges - 1; l++ ) {
                    random.at(1) = boundaries.at(1) + ( 1 + l ) * n1length / ( ( double ) n1edges );
                    random.at(2) = boundaries.at(3) + ( 1 + j ) * n2length / ( ( double ) n2edges );

                    if ( random.at(1) + grid->TOL < boundaries.at(1) ||
                         random.at(1) - grid->TOL > boundaries.at(2) ||
                         random.at(2) + grid->TOL < boundaries.at(3) ||
                         random.at(2) - grid->TOL > boundaries.at(4) ||
                         random.at(3) + grid->TOL < boundaries.at(5) ||
                         random.at(3) - grid->TOL > boundaries.at(6) ) {
		      generator::error("The nodes are located in the wrong position");
                    }

                    /* vertex = ( Vertex * ) ( Vertex(vertexNumber + 1, grid).ofType() ); */
                    /* vertex->setCoordinates(random); */
                    /* grid->setVertex(vertexNumber + 1, vertex); */
                    /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber + 1, random); */
                    /* vertexNumber++; */

		    /* auto *v = new Vertex(vertexNumber+1, grid); */
		    /* v->setCoordinates(random); */
		    /* grid->setVertex(vertexNumber+1, v); */
		    /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1, random); */
		    /* vertexNumber++; */

		    grid->addVertex(random);
		    
                }
            }
        }  else {
            for ( int k = 0; k < n2edges; k++ ) {
                for ( int n = 0; n < n1edges; n++ ) {
                    random.at(1) = boundaries.at(1) + ( 0.5 + n ) * n1length / ( ( double ) n1edges );
                    random.at(2) = boundaries.at(3) + ( 0.5 + k ) * n2length / ( ( double ) n2edges );

                    if ( random.at(1) + grid->TOL < boundaries.at(1) ||
                         random.at(1) - grid->TOL > boundaries.at(2) ||
                         random.at(2) + grid->TOL < boundaries.at(3) ||
                         random.at(2) - grid->TOL > boundaries.at(4) ||
                         random.at(3) + grid->TOL < boundaries.at(5) ||
                         random.at(3) - grid->TOL > boundaries.at(6) ) {
		      generator::error("The nodes are located in the wrong position");
                    }

		    /* auto *v = new Vertex(vertexNumber+1, grid); */
		    /* v->setCoordinates(random); */
		    /* grid->setVertex(vertexNumber+1, v); */
		    /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1, random); */
		    /* vertexNumber++; */
		    
		    grid->addVertex(random);
		    
		    /* vertex = ( Vertex * ) ( Vertex(vertexNumber + 1, grid).ofType() ); */
                    /* vertex->setCoordinates(random); */
                    /* grid->setVertex(vertexNumber + 1, vertex); */
                    /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber + 1, random); */
                    /* vertexNumber++; */
                }
            }
        }
    }

    generator::ensure_size1(grid->vertexList,vertexNumber);

    return 1;
}



int Region :: generateRegularPoints2()
{
    Surface *surface;
    Curve *curve;
    Vertex *vertex;

    //  int firstVertex = giveLocalVertex(1);
    int localSurface;
    int localCurve;
    int localVertex;
    oofem::IntArray curves;
    double x, y, z;

    oofem::FloatArray boundaries;
    grid->defineBoundaries(boundaries);

    oofem::FloatArray random(3);

    int n1edges = grid->xyzEdges.at(1);
    int n2edges = grid->xyzEdges.at(2);
    int n3edges = grid->xyzEdges.at(3);
    double n1length = fabs( boundaries.at(1) - boundaries.at(2) );
    double n2length = fabs( boundaries.at(3) - boundaries.at(4) );
    double n3length = fabs( boundaries.at(5) - boundaries.at(6) );


    int vertexNumber = generator::size1(grid->vertexList);
    int tempSize = 1e9;
    generator::ensure_size1(grid->vertexList,tempSize);
    for ( int i = 0; i < n3edges * 2 - 1; i++ ) {
        random.at(3) = boundaries.at(5) + ( i + 1 ) * n3length / ( ( double ) n3edges * 2 );

        //Generate precise points
        if ( i % 2 ) { // B types odd
            for ( int j = 0; j < n2edges * 2 - 1; j++ ) {
                if ( j % 2 ) {
                    for ( int l = 0; l < n1edges - 1; l++ ) {
                        random.at(1) = boundaries.at(1) + ( 1 + l ) * n1length / ( ( double ) n1edges );
                        random.at(2) = boundaries.at(3) + ( 1 + j ) * n2length / ( ( double ) n2edges * 2 );


		    /* auto *v = new Vertex(vertexNumber+1, grid); */
		    /* v->setCoordinates(random); */
		    /* grid->setVertex(vertexNumber+1, v); */
		    /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1, random); */
		    /* vertexNumber++; */

		    		    grid->addVertex(random);
			
                        /* vertex = ( Vertex * ) ( Vertex(vertexNumber + 1, grid).ofType() ); */
                        /* vertex->setCoordinates(random); */
                        /* grid->setVertex(vertexNumber + 1, vertex); */
                        /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber + 1, random); */
                        /* vertexNumber++; */

                        if ( random.at(1) + grid->TOL <  boundaries.at(1) ||
                             random.at(1) - grid->TOL > boundaries.at(2) ||
                             random.at(2) + grid->TOL <  boundaries.at(3) ||
                             random.at(2) - grid->TOL > boundaries.at(4) ||
                             random.at(3) + grid->TOL <  boundaries.at(5) ||
                             random.at(3) - grid->TOL > boundaries.at(6) ) {

			  generator::error("The nodes are located in the wrong position");			  
                        }
                    }
                } else   {
                    random.at(2) = boundaries.at(3) + ( 1 + j ) * n2length / ( n2edges * 2 );
                    int k = 0;
                    while ( k < n1edges ) {
                        random.at(1) = boundaries.at(1) + ( k + 0.5 ) * n1length / ( ( double ) n1edges );
                        k++;

			
                        /* vertex = ( Vertex * ) ( Vertex(vertexNumber + 1, grid).ofType() ); */
                        /* vertex->setCoordinates(random); */
                        /* grid->setVertex(vertexNumber + 1, vertex); */
                        /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber + 1, random); */
                        /* vertexNumber++; */

		    /* auto *v = new Vertex(vertexNumber+1, grid); */
		    /* v->setCoordinates(random); */
		    /* grid->setVertex(vertexNumber+1, v); */
		    /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1, random); */
		    /* vertexNumber++; */

		    grid->addVertex(random);
			
                        if ( random.at(1) + grid->TOL < boundaries.at(1) ||
                             random.at(1) - grid->TOL > boundaries.at(2) ||
                             random.at(2) + grid->TOL < boundaries.at(3) ||
                             random.at(2) - grid->TOL > boundaries.at(4) ||
                             random.at(3) + grid->TOL < boundaries.at(5) ||
                             random.at(3) - grid->TOL > boundaries.at(6) ) {
			  generator::error("The nodes are located in the wrong position");
                        }
                    }
                }
            }
        } else   { //i is a multiple of three
            for ( int n = 0; n < n2edges * 2 - 1; n++ ) {
                if ( n % 2 ) {
                    for ( int m = 0; m < n1edges; m++ ) {
                        random.at(1) = boundaries.at(1) + ( 0.5 + m ) * n1length / ( ( double ) n1edges );
                        random.at(2) = boundaries.at(3) + ( n + 1 ) * n2length / ( ( double ) n2edges * 2 );

                        /* vertex = ( Vertex * ) ( Vertex(vertexNumber + 1, grid).ofType() ); */
                        /* vertex->setCoordinates(random); */
                        /* grid->setVertex(vertexNumber + 1, vertex); */
                        /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber + 1, random); */
                        /* vertexNumber++; */

		    /* auto *v = new Vertex(vertexNumber+1, grid); */
		    /* v->setCoordinates(random); */
		    /* grid->setVertex(vertexNumber+1, v); */
		    /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1, random); */
		    /* vertexNumber++; */
		    grid->addVertex(random);
			

                        if ( random.at(1) + grid->TOL < boundaries.at(1) ||
                             random.at(1) - grid->TOL > boundaries.at(2) ||
                             random.at(2) + grid->TOL < boundaries.at(3) ||
                             random.at(2) - grid->TOL > boundaries.at(4) ||
                             random.at(3) + grid->TOL < boundaries.at(5) ||
                             random.at(3) - grid->TOL > boundaries.at(6) ) {
			  generator::error("The nodes are located in the wrong position");
                        }
                    }
                } else   {
                    for ( int s = 0; s < n1edges - 1; s++ ) {
                        random.at(1) =  boundaries.at(1) + ( s + 1 ) * n1length / ( ( double ) n1edges );
                        random.at(2) =  boundaries.at(3) + ( n + 1 ) * n2length / ( ( double ) n2edges * 2 );


		    /* auto *v = new Vertex(vertexNumber+1, grid); */
		    /* v->setCoordinates(random); */
		    /* grid->setVertex(vertexNumber+1, v); */
		    /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1, random); */
		    /* vertexNumber++; */

		    		    grid->addVertex(random);
			
                        /* vertex = ( Vertex * ) ( Vertex(vertexNumber + 1, grid).ofType() ); */
                        /* vertex->setCoordinates(random); */
                        /* grid->setVertex(vertexNumber + 1, vertex); */
                        /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber + 1, random); */
                        /* vertexNumber++; */

                        if ( random.at(1) + grid->TOL < boundaries.at(1) ||
                             random.at(1) - grid->TOL > boundaries.at(2) ||
                             random.at(2) + grid->TOL < boundaries.at(3) ||
                             random.at(2) - grid->TOL > boundaries.at(4) ||
                             random.at(3) + grid->TOL < boundaries.at(5) ||
                             random.at(3) - grid->TOL > boundaries.at(6) ) {

			  generator::error("The nodes are located in the wrong position\n");
                        }
                    }
                }
            }
        }
    }

    generator::ensure_size1(grid->vertexList,vertexNumber);

    return 1;
}




int Region :: generatePeriodicPoints()
{
    Vertex *vertex;

    oofem::FloatArray boundaries;
    grid->defineBoundaries(boundaries);

    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int randomIntegerOne = grid->giveRandomInteger() - 1;
    int randomIntegerTwo = grid->giveRandomInteger() - 2;
    int randomIntegerThree = grid->giveRandomInteger() - 3;

    oofem::FloatArray random(3);
    oofem::FloatArray coordsAtPeriodicity(3);
    int flag;


    oofem::IntArray periodicityFlag;
    grid->givePeriodicityFlag(periodicityFlag);
    
    double maxIter = grid->giveMaximumIterations();
    oofem::FloatArray mirroredRandom(3);

    int vertexNumber = generator::size1(grid->vertexList);
    int tempSize = 1e9;
    generator::ensure_size1(grid->vertexList,tempSize);

    int mult = 1;
    int tempIter = 0;
    oofem::FloatArray randomPeriodic(3);
    for ( int i = 0; i < maxIter; i++ ) {
        //Generate random point

        random.at(1) = boundaries.at(1) + grid->TOL + grid->ran1(& randomIntegerOne) * ( boundaries.at(2) - boundaries.at(1) - 2 * grid->TOL );
        random.at(2) = boundaries.at(3) + grid->TOL + grid->ran1(& randomIntegerTwo) * ( boundaries.at(4) - boundaries.at(3) - 2. * grid->TOL );
        random.at(3) = boundaries.at(5) + grid->TOL + grid->ran1(& randomIntegerThree) * ( boundaries.at(6) - boundaries.at(5) - 2. * grid->TOL );

        //Check if this is far enough from the others
        flag = 0;
        flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, this->refinement * grid->giveDiameter(random) );

        if ( flag == 0 ) {

		    /* auto *v = new Vertex(vertexNumber+1, grid); */
		    /* v->setCoordinates(random); */
		    /* grid->setVertex(vertexNumber+1, v); */
		    /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1, random); */
		    /* vertexNumber++; */

		    grid->addVertex(random);
		    
	  
            /* vertex = ( Vertex * ) ( Vertex(vertexNumber + 1, grid).ofType() ); */
            /* vertex->setCoordinates(random); */
            /* grid->setVertex(vertexNumber + 1, vertex); */
            /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber + 1, random); */
            /* vertexNumber++; */

            //Do the periodic shift
            for ( int x = -1; x < 2; x++ ) {
                for ( int y = -1; y < 2; y++ ) {
                    for ( int z = -1; z < 2; z++ ) {
                        if ( !( z == 0 && y == 0 && x == 0 ) ) {
                            randomPeriodic.at(1) = random.at(1) + x * specimenDimension.at(1);
                            randomPeriodic.at(2) = random.at(2) + y * specimenDimension.at(2);
                            randomPeriodic.at(3) = random.at(3) + z * specimenDimension.at(3);

		    /* auto *v = new Vertex(vertexNumber+1, grid); */
		    /* v->setCoordinates(randomPeriodic); */
		    /* grid->setVertex(vertexNumber+1, v); */
		    /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1, randomPeriodic); */
		    /* vertexNumber++; */

		    grid->addVertex(randomPeriodic);
		    
                            /* vertex = ( Vertex * ) ( Vertex(vertexNumber + 1, grid).ofType() ); */
                            /* vertex->setCoordinates(randomPeriodic); */
                            /* grid->setVertex(vertexNumber + 1, vertex); */
                            /* grid->giveGridLocalizer()->insertSequentialNode(vertexNumber + 1, randomPeriodic); */
                            /* vertexNumber++; */

			    //Print every 1m node message. Debug. 
			    if(tempIter > mult*10) {
			      std :: cout << "Placed " << vertexNumber << " points in region "<< this->giveNumber() <<"\n";
			      std :: cout << "Current greatest iterator is " << tempIter << ". Maximum allowed iterator is "<< grid->giveMaximumIterations() <<".\n";
			      mult++;
			    }
                        }
                    }
                }
            }

	    if(i>tempIter){
	      tempIter = i;
	    }
            i = 0;
        }
    }

    generator::ensure_size1(grid->vertexList,vertexNumber);


    return 1;
}
