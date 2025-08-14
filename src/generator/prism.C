#include "prism.h"
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

Prism :: Prism(int n, Grid *aGrid) : Region(n, aGrid)
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


int Prism :: generatePoints()
{
    oofem::IntArray curves;

    printf("Generating points for prism\n");
    
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

    double borderX=0.;
    double borderY=0.;
    double borderZ=0.;
    if(randomFlag == 0){
      borderX = grid->TOL;
      borderY = grid->TOL;
      borderZ = grid->TOL;
    }
    else{
      //x-direction
      if(periodicityFlag.at(1) == 1){ 
	borderX = grid->TOL;
      }
      else{
	borderX = grid->diameter;
      }
      //y-direction
      if(periodicityFlag.at(2) == 1){ 
	borderY = grid->TOL;
      }
      else{
	borderY = grid->diameter;
      }
      //z-direction
      if(periodicityFlag.at(3) == 1){ 
	borderZ = grid->TOL;
      }
      else{
	borderZ = grid->diameter;
      }	    
    }
    
    //Generation of vertices needs to be split into three parts. First the edges, then the surfaces and finally the region.

    printf("Start with regular points for edges.\n");
    //Generate vertices on edges regularly
    //For x-direction
    int nTarget = 0;
    double newDiameter;
    if(periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0){
      nTarget = ceil(specimenDimension.at(1)/grid->diameter);
      newDiameter = specimenDimension.at(1)/nTarget;
      if(periodicityFlag.at(1) == 1) nTarget--;
      //First of four lines
      for(int z = 0;z<2;z++){
	for(int y = 0;y<2;y++){
	  for(int k = 0;k<=nTarget;k++){
	    if(periodicityFlag.at(1) == 1 && k == 0){
	      random.at(1) = boundaries.at(1) + (0.5 + k)*newDiameter;   
	    }
	    else{
	      random.at(1) = boundaries.at(1) + k*newDiameter;
	    }
	    random.at(2) = boundaries.at(3) + y*specimenDimension.at(2);
	    random.at(3) = boundaries.at(5) + z*specimenDimension.at(3);
	    
	    //Check if this is far enough from the others
	    flag = 0;
	    flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, 0.99*newDiameter );
	    
	    if ( flag == 0 ) {


	      grid->addVertex(random);
	      
	      for ( int x = -1; x< 2; x++){
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
		  flag = 0;
		  flag = grid->giveGridLocalizer()->checkNodesWithinBox( newRandom, 0.99*newDiameter ); //so that corner points are not duplicated
		  if(flag == 0){
		    grid->addVertex(newRandom);
		    
		  }
		}
	      }//end of shifting mirroring
	    }
	  }
	}
      }
    }

    //For y-direction
    if(periodicityFlag.at(1) == 0 && periodicityFlag.at(3) == 0){
      nTarget = ceil(specimenDimension.at(2)/grid->diameter);
      newDiameter = specimenDimension.at(2)/nTarget;
      //First of four lines
      for(int z = 0;z<2;z++){
	for(int x = 0;x<2;x++){
	  for(int k = 0;k<=nTarget;k++){
	    random.at(1) = boundaries.at(1) + x*specimenDimension.at(1);
	    random.at(2) = boundaries.at(3) + k*newDiameter;
	    random.at(3) = boundaries.at(5) + z*specimenDimension.at(3);
	    
	    //Check if this is far enough from the others
	    flag = 0;
	    flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, 0.99*newDiameter );
	    
	    if ( flag == 0 ) {
	      grid->addVertex(random);

	      for ( int y = -1; y< 2; y++){
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

		  flag = 0;
		  flag = grid->giveGridLocalizer()->checkNodesWithinBox( newRandom, 0.99*newDiameter ); //so that corner points are not duplicated
		  if(flag == 0){		  
		    grid->addVertex(newRandom);


		  }
		}
	      }//end of shifting/mirroring	      
	    }
	  }
	}
      }
    }


    //For z-direction
    if(periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 0){
      nTarget = ceil(specimenDimension.at(3)/grid->diameter);
      newDiameter = specimenDimension.at(3)/nTarget;

      //First of four lines
      for(int y = 0;y<2;y++){
	for(int x = 0;x<2;x++){
	  for(int k = 0;k<=nTarget;k++){
	    random.at(1) = boundaries.at(1) + x*specimenDimension.at(1);
	    random.at(2) = boundaries.at(3) + y*specimenDimension.at(2);
	    random.at(3) = boundaries.at(5) + k*newDiameter;
	    
	    //Check if this is far enough from the others
	    flag = 0;
	    flag = grid->giveGridLocalizer()->checkNodesWithinBox( random,  0.99*newDiameter );
	    
	    if ( flag == 0 ) {
	      grid->addVertex(random);

	      for ( int z = -1; z< 2; z++){
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

		  flag = 0;
		  flag = grid->giveGridLocalizer()->checkNodesWithinBox( newRandom, 0.99*newDiameter ); //so that corner points are not duplicated
		  if(flag == 0){		  		  
		    grid->addVertex(newRandom);

		  }
		}
	      }//end of shifting/mirroring 	      
	    }
	  }
	}
      }
    }

    printf("Finished edges. Current number of points are %d\n", grid->giveNumberOfVertices());
    
    printf("Start with surfaces\n");
    oofem::FloatArray normal(3);
    //Generate random vertices on surface
    if(periodicityFlag.at(1) == 0){//surface with normals 1,0,0 are generated

      //Define normal
      normal.at(1) = 1;
      normal.at(2) = 0;
      normal.at(3) = 0;

      for ( int i = 0; i < maxIter; i++ ) {
	random.at(1) = boundaries.at(1);
	random.at(2) = boundaries.at(3) + borderY + grid->ran1(& randomIntegerTwo) * ( specimenDimension.at(2) -2.*borderY);
	random.at(3) = boundaries.at(5) + borderZ + grid->ran1(& randomIntegerThree) * ( specimenDimension.at(3) -2.*borderZ);  

	flag = 0;
	flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, grid->giveDiameter(random) );
	
	if ( flag == 0 ) {
	  grid->addVertex(random);
	  
	  i = 0;	  
	  
	  mirrorShiftSurface(random, normal,specimenDimension,boundaries,periodicityFlag);

	  //Shift node over to other side
	  newRandom.at(1) = random.at(1) + specimenDimension.at(1);
	  newRandom.at(2) = random.at(2);
	  newRandom.at(3) = random.at(3);

	  grid->addVertex(newRandom);

	  mirrorShiftSurface(newRandom, normal,specimenDimension,boundaries,periodicityFlag);
	}
      }
    }

    printf("Finished surfaces with normal 1,0,0. Current number of points are %d\n", grid->giveNumberOfVertices());
    
    if(periodicityFlag.at(2) == 0){//surface with normals 0,1,0 are generated
      //Define normal
      normal.at(1) = 0;
      normal.at(2) = 1;
      normal.at(3) = 0;

      for ( int i = 0; i < maxIter; i++ ) {
	random.at(1) = boundaries.at(1) + borderX + grid->ran1(& randomIntegerTwo) * ( specimenDimension.at(1) -2.*borderX);
	random.at(2) = boundaries.at(3);
	random.at(3) = boundaries.at(5) + borderZ + grid->ran1(& randomIntegerThree) * ( specimenDimension.at(3) -2.*borderZ);	
	flag = 0;
	flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, grid->giveDiameter(random) );
	
	if ( flag == 0 ) {
	  grid->addVertex(random);

	  i = 0;	  
	  
	  mirrorShiftSurface(random, normal,specimenDimension,boundaries,periodicityFlag);

	  //Shift node over to other side
	  newRandom.at(1) = random.at(1);
	  newRandom.at(2) = random.at(2) + specimenDimension.at(2);
	  newRandom.at(3) = random.at(3);

	  grid->addVertex(newRandom);
	  
	  mirrorShiftSurface(newRandom, normal,specimenDimension,boundaries,periodicityFlag);
	}
      }
    }

    printf("Finished surfaces with normal 0,1,0. Current number of points are %d\n", grid->giveNumberOfVertices());

    if(periodicityFlag.at(3) == 0){//surface with normals 0,0,1 are generated
      //Define normal
      normal.at(1) = 0;
      normal.at(2) = 0;
      normal.at(3) = 1;

      for ( int i = 0; i < maxIter; i++ ) {
	random.at(1) = boundaries.at(1) + borderX + grid->ran1(& randomIntegerTwo) * ( specimenDimension.at(1) -2.*borderX);
	random.at(2) = boundaries.at(3) + borderY + grid->ran1(& randomIntegerThree) * ( specimenDimension.at(2) -2.*borderY);
	random.at(3) = boundaries.at(5);

	flag = 0;
	flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, grid->giveDiameter(random) );
	
	if ( flag == 0 ) {
	  grid->addVertex(random);

	  i = 0;	  
	  
	  mirrorShiftSurface(random, normal,specimenDimension,boundaries,periodicityFlag);

	  //Shift node over to other side
	  newRandom.at(1) = random.at(1);
	  newRandom.at(2) = random.at(2);
	  newRandom.at(3) = random.at(3) + specimenDimension.at(3);
	  
	  grid->addVertex(newRandom);
	  mirrorShiftSurface(newRandom, normal,specimenDimension,boundaries,periodicityFlag);
	}
      }
    }
    printf("Finished surfaces with normal 0,0,1. Current number of points are %d\n", grid->giveNumberOfVertices());
    
     printf("Start with region.\n");
    
    //Generate random vertices within the prism
    int mult = 0;
    for ( int i = 0; i < maxIter; i++ ) {

      random.at(1) = boundaries.at(1) + borderX + grid->ran1(& randomIntegerOne) * ( specimenDimension.at(1) -  2.*borderX );
      random.at(2) = boundaries.at(3) + borderY + grid->ran1(& randomIntegerTwo) * ( specimenDimension.at(2) - 2.*borderY );
      random.at(3) = boundaries.at(5) + borderZ + grid->ran1(& randomIntegerThree) * ( specimenDimension.at(3) - 2*borderZ );

      //Check if this is far enough from the others
        flag = 0;
        flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, grid->giveDiameter(random) );

        if ( flag == 0 ) {
	  grid->addVertex(random);


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

		    grid->addVertex(newRandom);
		    
		  }		  
		}
	      }
	    }
	}

	if(i > mult*10000 && i < (mult+1)*10000) {
	  std :: cout << "Placed " << grid->giveNumberOfVertices() << " points in region "<< this->giveNumber() <<"\n";
	  std :: cout << "Current greatest iterator is " << i << ". Maximum allowed iterator is "<< grid->giveMaximumIterations() <<".\n";
	  mult++;
	}
	
	
    }

    printf("Finished region. Current number of points are %d\n", grid->giveNumberOfVertices());
    

    return 1;
}

    
void Prism :: defineBoundaries(oofem::FloatArray &boundaries)
//Determine the boundaries of the domain
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



void Prism :: mirrorShiftSurface(oofem::FloatArray& random, oofem::FloatArray& normal,oofem::FloatArray& specimenDimension,oofem::FloatArray& boundaries, oofem::IntArray& periodicityFlag)
{

  //Mirror (or periodic shift) with respect to two of the three axis.

  oofem::FloatArray newRandom(3);
  
  if(normal.at(1) == 1 && normal.at(2) == 0 && normal.at(3) == 0){//y-z coordinate system
    for ( int y = -1; y < 2; y++ ) {
      for ( int z = -1; z < 2; z++ ) {
	if ( !( z == 0 && y == 0 ) ) {
	  newRandom.at(1) = random.at(1);
	  //y-direction
	  if(periodicityFlag.at(2) == 1){
	    newRandom.at(2) = random.at(2) + y * specimenDimension.at(2);
	  }
	  else{
	    if( y==-1 ){
	      newRandom.at(2) = random.at(2) - 2. * ( random.at(2) - boundaries.at(3) );
	    }
	    else if( y==0 ){
	      newRandom.at(2) = random.at(2);
	    }
	    else if( y==1 ){
	      newRandom.at(2) = random.at(2) - 2. * ( random.at(2) - boundaries.at(4) );
	    }
	  }
	  // z-direction
	  if(periodicityFlag.at(3) == 1){
	    newRandom.at(3) = random.at(3) + z * specimenDimension.at(3);
	  }
	  else{
	    if( z== -1 ){
	      newRandom.at(3) = random.at(3) - 2. * ( random.at(3) - boundaries.at(5) );
	    }
	    else if( z==0 ){
	      newRandom.at(3) = random.at(3);
	    }
	    else if( z==1 ){
	      newRandom.at(3) = random.at(3) - 2. * ( random.at(3) - boundaries.at(6) );
	    }
	  }

	  grid->addVertex(newRandom);
	}
      }
    }	    
  }
  else if(normal.at(1) == 0 && normal.at(2) == 1 && normal.at(3) == 0){//x-z coordinate system
    for ( int x = -1; x < 2; x++ ) {
      for ( int z = -1; z < 2; z++ ) {
	if ( !( x == 0 && z == 0 ) ) {
	  //x-direction
	  if(periodicityFlag.at(1) == 1){
	    newRandom.at(1) = random.at(1) + x * specimenDimension.at(1);
	  }
	  else{
	    if( x==-1 ){
	      newRandom.at(1) = random.at(1) - 2. * ( random.at(1) - boundaries.at(1) );
	    }
	    else if( x==0 ){
	      newRandom.at(1) = random.at(1);
	    }
	    else if( x==1 ){
	      newRandom.at(1) = random.at(1) - 2. * ( random.at(1) - boundaries.at(2) );
	    }
	  }
	  //y-direction
	  newRandom.at(2) = random.at(2);
	  //z-direction
	  if(periodicityFlag.at(3) == 1){
	    newRandom.at(3) = random.at(3) + z * specimenDimension.at(3);
	  }
	  else{
	    if( z== -1 ){
	      newRandom.at(3) = random.at(3) - 2. * ( random.at(3) - boundaries.at(5) );
	    }
	    else if( z==0 ){
	      newRandom.at(3) = random.at(3);
	    }
	    else if( z==1 ){
	      newRandom.at(3) = random.at(3) - 2. * ( random.at(3) - boundaries.at(6) );
	    }
	  }
	  grid->addVertex(newRandom);

	}
      }
    }
  }
  if(normal.at(1) == 0 && normal.at(2) == 0 && normal.at(3) == 1){//x-y coordinate system
    for ( int x = -1; x < 2; x++ ) {
      for ( int y = -1; y < 2; y++ ) {
	if ( !( x == 0 && y == 0 ) ) {
	  //x-direction
	  if(periodicityFlag.at(1) == 1){
	    newRandom.at(1) = random.at(1) + x * specimenDimension.at(1);
	  }
	  else{
	    if( x==-1 ){
	      newRandom.at(1) = random.at(1) - 2. * ( random.at(1) - boundaries.at(1) );
	    }
	    else if( x==0 ){
	      newRandom.at(1) = random.at(1);
	    }
	    else if( x==1 ){
	      newRandom.at(1) = random.at(1) - 2. * ( random.at(1) - boundaries.at(2) );
	    }
	  }
	  //y-direction
	  if(periodicityFlag.at(2) == 1){
	    newRandom.at(2) = random.at(2) + y * specimenDimension.at(2);
	  }
	  else{
	    if( y==-1 ){
	      newRandom.at(2) = random.at(2) - 2. * ( random.at(2) - boundaries.at(3) );
	    }
	    else if( y==0 ){
	      newRandom.at(2) = random.at(2);
	    }
	    else if( y==1 ){
	      newRandom.at(2) = random.at(2) - 2. * ( random.at(2) - boundaries.at(4) );
	    }
	  }
	  //z-direction
	  newRandom.at(3) = random.at(3);
	  grid->addVertex(newRandom);

	}		
      }
    }
  }
  return;
}


void
Prism :: initializeFrom(GeneratorInputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{
    //two points to define prism
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


int Prism :: generateRegularPoints1()
{
    oofem::IntArray curves;

    oofem::FloatArray boundaries;
    this->defineBoundaries(boundaries);
    oofem::FloatArray random(3);

    int n1edges = grid->xyzEdges.at(1);
    int n2edges = grid->xyzEdges.at(2);
    int n3edges = grid->xyzEdges.at(3);
    double n1length = fabs( boundaries.at(1) - boundaries.at(2) );
    double n2length = fabs( boundaries.at(3) - boundaries.at(4) );
    double n3length = fabs( boundaries.at(5) - boundaries.at(6) );

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
                        std :: cout << "The nodes are located in the wrong position" << '\n';
                        exit(0);
                    }

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
                        std :: cout << "The nodes are located in the wrong position" << '\n';
                        exit(0);
                    }

		    		    grid->addVertex(random);
                }
            }
        }
    }


    return 1;
}



int Prism :: generateRegularPoints2()
{

    oofem::IntArray curves;


    oofem::FloatArray boundaries;
    this->defineBoundaries(boundaries);

    oofem::FloatArray random(3);

    int n1edges = grid->xyzEdges.at(1);
    int n2edges = grid->xyzEdges.at(2);
    int n3edges = grid->xyzEdges.at(3);
    double n1length = fabs( boundaries.at(1) - boundaries.at(2) );
    double n2length = fabs( boundaries.at(3) - boundaries.at(4) );
    double n3length = fabs( boundaries.at(5) - boundaries.at(6) );


    for ( int i = 0; i < n3edges * 2 - 1; i++ ) {
        random.at(3) = boundaries.at(5) + ( i + 1 ) * n3length / ( ( double ) n3edges * 2 );

        //Generate precise points
        if ( i % 2 ) { // B types odd
            for ( int j = 0; j < n2edges * 2 - 1; j++ ) {
                if ( j % 2 ) {
                    for ( int l = 0; l < n1edges - 1; l++ ) {
                        random.at(1) = boundaries.at(1) + ( 1 + l ) * n1length / ( ( double ) n1edges );
                        random.at(2) = boundaries.at(3) + ( 1 + j ) * n2length / ( ( double ) n2edges * 2 );

		    grid->addVertex(random);

                        if ( random.at(1) + grid->TOL <  boundaries.at(1) ||
                             random.at(1) - grid->TOL > boundaries.at(2) ||
                             random.at(2) + grid->TOL <  boundaries.at(3) ||
                             random.at(2) - grid->TOL > boundaries.at(4) ||
                             random.at(3) + grid->TOL <  boundaries.at(5) ||
                             random.at(3) - grid->TOL > boundaries.at(6) ) {
                            std :: cout << "The nodes are located in the wrong position" << '\n';
                            exit(0);
                        }
                    }
                } else   {
                    random.at(2) = boundaries.at(3) + ( 1 + j ) * n2length / ( n2edges * 2 );
                    int k = 0;
                    while ( k < n1edges ) {
                        random.at(1) = boundaries.at(1) + ( k + 0.5 ) * n1length / ( ( double ) n1edges );
                        k++;
		    grid->addVertex(random);


                        if ( random.at(1) + grid->TOL < boundaries.at(1) ||
                             random.at(1) - grid->TOL > boundaries.at(2) ||
                             random.at(2) + grid->TOL < boundaries.at(3) ||
                             random.at(2) - grid->TOL > boundaries.at(4) ||
                             random.at(3) + grid->TOL < boundaries.at(5) ||
                             random.at(3) - grid->TOL > boundaries.at(6) ) {
                            std :: cout << "The nodes are located in the wrong position" << '\n';
                            exit(0);
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

		    grid->addVertex(random);


                        if ( random.at(1) + grid->TOL < boundaries.at(1) ||
                             random.at(1) - grid->TOL > boundaries.at(2) ||
                             random.at(2) + grid->TOL < boundaries.at(3) ||
                             random.at(2) - grid->TOL > boundaries.at(4) ||
                             random.at(3) + grid->TOL < boundaries.at(5) ||
                             random.at(3) - grid->TOL > boundaries.at(6) ) {
                            std :: cout << "The nodes are located in the wrong position" << '\n';
                            exit(0);
                        }
                    }
                } else   {
                    for ( int s = 0; s < n1edges - 1; s++ ) {
                        random.at(1) =  boundaries.at(1) + ( s + 1 ) * n1length / ( ( double ) n1edges );
                        random.at(2) =  boundaries.at(3) + ( n + 1 ) * n2length / ( ( double ) n2edges * 2 );

		    grid->addVertex(random);

                        if ( random.at(1) + grid->TOL < boundaries.at(1) ||
                             random.at(1) - grid->TOL > boundaries.at(2) ||
                             random.at(2) + grid->TOL < boundaries.at(3) ||
                             random.at(2) - grid->TOL > boundaries.at(4) ||
                             random.at(3) + grid->TOL < boundaries.at(5) ||
                             random.at(3) - grid->TOL > boundaries.at(6) ) {
                            std :: cout << "The nodes are located in the wrong position \n";
                            exit(0);
                        }
                    }
                }
            }
        }
    }


    return 1;
}


int Prism :: generatePeriodicPoints()
{
  //    Vertex *vertex;

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
    oofem::FloatArray coordsAtPeriodicity(3);
    int flag;

    oofem::IntArray periodicityFlag;
    grid->givePeriodicityFlag(periodicityFlag);
    
    double maxIter = grid->giveMaximumIterations();
    oofem::FloatArray mirroredRandom(3);

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
		    grid->addVertex(random);

            //Do the periodic shift
            for ( int x = -1; x < 2; x++ ) {
                for ( int y = -1; y < 2; y++ ) {
                    for ( int z = -1; z < 2; z++ ) {
                        if ( !( z == 0 && y == 0 && x == 0 ) ) {
                            randomPeriodic.at(1) = random.at(1) + x * specimenDimension.at(1);
                            randomPeriodic.at(2) = random.at(2) + y * specimenDimension.at(2);
                            randomPeriodic.at(3) = random.at(3) + z * specimenDimension.at(3);

			    grid->addVertex(random);

			    //Print every 1m node message. Debug. 
			    if(tempIter > mult*10) {
			      std :: cout << "Placed " << grid->giveNumberOfVertices() << " points in region "<< this->giveNumber() <<"\n";
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


    return 1;
}
