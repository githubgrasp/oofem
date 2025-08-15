#include "generatorlistutils.h"
#include "curve.h"
#include "vertex.h"
#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include "octreegridlocalizer.h"
 #include <iostream>
#endif

Curve :: Curve(int n, Grid *aGrid) : GridComponent(n, aGrid) //, coordinates()
{
    this->number = n;
}

Curve :: ~Curve()
// Destructor.
{}


int
Curve :: giveLocalVertex(int i)
// Returns the i-th coordinate of the receiver.
{
    if ( i > vertices.giveSize() ) {
        return 0.;
    }

    return vertices.at(i);
}


void
Curve :: initializeFrom(GeneratorInputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{

    IR_GIVE_FIELD(ir, vertices, _IFT_Curve_vertices); // Macro
    refinement = 1.;
    IR_GIVE_FIELD(ir, refinement, _IFT_Curve_refine); // Macro
    normal.zero();
    IR_GIVE_FIELD(ir, normal, _IFT_Curve_refine); // Macro
    
    return;
}



/* Curve *Curve :: ofType() */
/* // Returns a new DofManager, which has the same number than the receiver, */
/* // but belongs to aClass (Node, ElementSide,..). */
/* { */
/*     Curve *curve; */

/*     curve = new Curve(number, grid); */

/*     return curve; */
/* } */


int Curve :: generatePoints()
{

  double xOne = ( grid->giveInputVertex( giveLocalVertex(1) ) )->giveCoordinate(1);
  double yOne = ( grid->giveInputVertex( giveLocalVertex(1) ) )->giveCoordinate(2);
  double zOne = ( grid->giveInputVertex( giveLocalVertex(1) ) )->giveCoordinate(3);

  double xTwo = ( grid->giveInputVertex( giveLocalVertex(2) ) )->giveCoordinate(1);
  double yTwo = ( grid->giveInputVertex( giveLocalVertex(2) ) )->giveCoordinate(2);
  double zTwo = ( grid->giveInputVertex( giveLocalVertex(2) ) )->giveCoordinate(3);

  int randomIntegerOne = grid->giveRandomInteger() - 1;

  oofem::FloatArray random(3);
  int flag;

  double boundaryFactor = this->refinement;
    
  Vertex *vertex;

  double shiftX = 0.;
  double shiftY = 0.;
  double shiftZ = 0.;
    
  double distance, shift;
  double x, y, z;
  double maxIter = grid->giveMaximumIterations();

  oofem::IntArray periodicityFlag(3);
  grid->givePeriodicityFlag(periodicityFlag);

  int randomFlag = grid->giveRandomFlag();
    
  oofem::FloatArray boundaries;
  grid->defineBoundaries(boundaries);
  oofem::FloatArray specimenDimension(3);
  specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
  specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
  specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

  oofem::FloatArray newRandom(3);

  for ( int i = 0; i < maxIter; i++ ) {
    //Generate random point
    random.at(1) = xOne + grid->ran1(& randomIntegerOne) * ( xTwo - xOne );
    random.at(2) = yOne + grid->ran1(& randomIntegerOne) * ( yTwo - yOne );
    random.at(3) = zOne + grid->ran1(& randomIntegerOne) * ( zTwo - zOne );

    //Check if this is far enough from the others
    flag = 0;
    flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, boundaryFactor * grid->giveDiameter(random) );

    if ( flag == 0 ) {


      		    grid->addVertex(random);

      i = 0;


      mirrorShift(random, normal,specimenDimension,boundaries,periodicityFlag);

      //Option random=2, for which periodic points are shifted across the specimen.
      if(randomFlag ==2){	      
	if( normal.at(1) == 1 && normal.at(2) == 0 && normal.at(3) == 0){//x-direction

	  if(yOne == boundaries.at(3)){//Check that this is the right line
	    shiftY = 1.;
	  }
	  else{
	    shiftY = -1.;
	  }
	  
	  if(zOne == boundaries.at(5)){//Check that this is the right line
	    shiftZ = 1.;
	  }
	  else{
	    shiftZ = -1.;
	  }
		    
	  //Shift in y-direction
	  newRandom.at(1) = random.at(1);
	  newRandom.at(2) = random.at(2) + shiftY*specimenDimension.at(2);
	  newRandom.at(3) = random.at(3);


		    grid->addVertex(newRandom);
	  
	  i = 0;
	  
	  mirrorShift(newRandom, normal,specimenDimension,boundaries,periodicityFlag);

	  //Shift in z-direction
	  newRandom.at(1) = random.at(1);
	  newRandom.at(2) = random.at(2);
	  newRandom.at(3) = random.at(3) + shiftZ*specimenDimension.at(3);

		    grid->addVertex(newRandom);

	  i = 0;
	  
	  mirrorShift(newRandom, normal,specimenDimension,boundaries,periodicityFlag);

	  //Shift in z and y-direction
	  newRandom.at(1) = random.at(1);
	  newRandom.at(2) = random.at(2) + shiftY*specimenDimension.at(2);
	  newRandom.at(3) = random.at(3) + shiftZ*specimenDimension.at(3);

	  grid->addVertex(newRandom);
	  

	  i = 0;
	  
	  mirrorShift(newRandom, normal,specimenDimension,boundaries,periodicityFlag);
	}
	else if( normal.at(1) == 0 && normal.at(2) == 1 && normal.at(3) == 0){//y-direction
	  
	  if(xOne == boundaries.at(1)){//Check that this is the right line
	    shiftX = 1.;
	  }
	  else{
	    shiftX = -1;
	  }
	  
	  if(zOne == boundaries.at(5)){//Check that this is the right line
	    shiftZ = 1.;
	  }
	  else{
	    shiftZ = -1;
	  }
				
	  //Shift in x-direction
	  newRandom.at(1) = random.at(1) + shiftX*specimenDimension.at(1);
	  newRandom.at(2) = random.at(2);
	  newRandom.at(3) = random.at(3);

	  grid->addVertex(newRandom);
	  
	  i = 0;
	  
	  mirrorShift(newRandom, normal,specimenDimension,boundaries,periodicityFlag);

		
	  //Shift in z-direction
	  newRandom.at(1) = random.at(1);
	  newRandom.at(2) = random.at(2);
	  newRandom.at(3) = random.at(3) + shiftZ*specimenDimension.at(3);

	  grid->addVertex(newRandom);

	  i = 0;
		  
	  mirrorShift(newRandom, normal,specimenDimension,boundaries,periodicityFlag);


	  //Shift in x and z-direction
	  newRandom.at(1) = random.at(1) + shiftX*specimenDimension.at(1);
	  newRandom.at(2) = random.at(2);
	  newRandom.at(3) = random.at(3) + shiftZ*specimenDimension.at(3);

	  grid->addVertex(newRandom);
	  

	  i = 0;
		  
	  mirrorShift(newRandom, normal,specimenDimension,boundaries,periodicityFlag);		  
	}
	else if( normal.at(1) == 0 && normal.at(2) == 0 && normal.at(3) == 1){//z-direction
	  
	  if(xOne == boundaries.at(1)){//Check that this is the right line
	    shiftX = 1.;
	  }
	  else{
	    shiftX = -1.;
	  }

	  if(yOne == boundaries.at(3)){//Check that this is the right line
	    shiftY = 1.;
	  }
	  else{
	    shiftY = -1.;
	  }
	       	  
	  //Shift in x-direction
	  newRandom.at(1) = random.at(1) + shiftX*specimenDimension.at(1);
	  newRandom.at(2) = random.at(2);
	  newRandom.at(3) = random.at(3);



	  grid->addVertex(newRandom);

	  i = 0;
	  	  
	  mirrorShift(newRandom, normal,specimenDimension,boundaries,periodicityFlag);
	
	  //Shift in y-direction
	  newRandom.at(1) = random.at(1);
	  newRandom.at(2) = random.at(2) + shiftY*specimenDimension.at(2);
	  newRandom.at(3) = random.at(3);

	  grid->addVertex(newRandom);
	  
	  i = 0;
	  
	  mirrorShift(newRandom, normal,specimenDimension,boundaries,periodicityFlag);
	  
	  
	  //Shift in x and y-direction
	  newRandom.at(1) = random.at(1) + shiftX*specimenDimension.at(1);
	  newRandom.at(2) = random.at(2) + shiftY*specimenDimension.at(2);
	  newRandom.at(3) = random.at(3);


	  grid->addVertex(newRandom);

	  i = 0;
	  
	  mirrorShift(newRandom, normal,specimenDimension,boundaries,periodicityFlag);
	}	    
	else{
	  printf("Error: Unknown direction.\n");
	  exit(1);
	}
      }
    }
  }	    
    
    
  return 1;
}

	
void
Curve :: mirrorShift(oofem::FloatArray &random, oofem::FloatArray &normal, oofem::FloatArray &specimenDimension, oofem::FloatArray &boundaries, oofem::IntArray& periodicityFlag) {
  
  Vertex *vertex;
  oofem::FloatArray newRandom(3);
  int randomFlag = grid->giveRandomFlag();  
  
  //Do mirroring and shift
  if(normal.at(1) == 1 && normal.at(2) == 0 && normal.at(3) == 0){//x-direction
    for ( int x = -1; x< 2; x++)
      if(x != 0) {
	if(periodicityFlag.at(1) == 1){//periodic shift along line
	  newRandom.at(1) = random.at(1) + x * specimenDimension.at(1);
	}
	else if(periodicityFlag.at(1) == 0){//mirror along line
	  if (x==-1){
	    newRandom.at(1) = random.at(1) - 2. * ( random.at(1) - boundaries.at(1) );
	  }
	  else if (x == 1){
	    newRandom.at(1) = random.at(1) - 2. * ( random.at(1) - boundaries.at(2) );
	  }
	}		
	newRandom.at(2) = random.at(2);
	newRandom.at(3) = random.at(3);


	  grid->addVertex(newRandom);
      }
  }
  
  
  //Do mirroring and shift
  if(normal.at(1) == 0 && normal.at(2) == 1 && normal.at(3) == 0){//y-direction
    for ( int y = -1; y< 2; y++)
      if(y != 0) {
	if(periodicityFlag.at(2) == 1){//shift along line
	  newRandom.at(2) = random.at(2) + y * specimenDimension.at(2);
	}
	else if(periodicityFlag.at(2) == 0){//mirror along line
	  if (y == -1){
	    newRandom.at(2) = random.at(2) - 2. * ( random.at(2) - boundaries.at(3) );
	  }
	  else if (y == 1){
	    newRandom.at(2) = random.at(2) - 2. * ( random.at(2) - boundaries.at(4) );
	  }
	}		
	newRandom.at(1) = random.at(1);
	newRandom.at(3) = random.at(3);


	  grid->addVertex(newRandom);
      }
  }
  
  //Do mirroring and shift
  if(normal.at(1) == 0 && normal.at(2) == 0 && normal.at(3) == 1){//z-direction
    for ( int z = -1; z< 2; z++)
      if(z != 0) {
	if(periodicityFlag.at(3) == 1){//shift along line
	  newRandom.at(3) = random.at(3) + z * specimenDimension.at(3);
	}
	else if(periodicityFlag.at(3) == 0){//mirror along line
	  if (z == -1){
	    newRandom.at(3) = random.at(3) - 2. * ( random.at(3) - boundaries.at(5) );
	  }
	  else if (z == 1){
	    newRandom.at(3) = random.at(3) - 2. * ( random.at(3) - boundaries.at(6) );
	  }
	}		
	newRandom.at(1) = random.at(1);
	newRandom.at(2) = random.at(2);

	  grid->addVertex(newRandom);
	  
	
      }
  }
  return;
}


