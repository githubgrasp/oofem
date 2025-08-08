#include "interfacesphere.h"
#include "curve.h"
#include "vertex.h"
#include "record.h"


#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#include "octreegridlocalizer.h"
#endif

InterfaceSphere :: InterfaceSphere(int n,Grid* aGrid) : Inclusion(n,aGrid) //, coordinates()
{
  this->number = n;  
}

InterfaceSphere :: ~InterfaceSphere()
{

}


int InterfaceSphere :: generatePoints()
{
  /*
    November 7 2018.
    Generate randomly placed points on sphere surface.
    Use method that Adrien introduced for fibres.
    If finely meshes, no need really to calculate minimum distance using sphere circle. 
    However, could be done like this if necessary
   */

  
  //minimum points
  int randomIntegerOne= grid->giveRandomInteger()-1;
  int randomIntegerTwo= grid->giveRandomInteger()-2;
  int randomIntegerThree = grid->giveRandomInteger()-3;
   
  double myPi = 3.14159265;
  
  oofem::FloatArray random(3);
  
  int flag;

  Vertex *vertex;
  
   double maxIter= grid->giveMaximumIterations();
   int vertexNumber = grid->giveNumberOfVertices();
   int tempSize = 1.e9;
   generator::ensure_size1(grid->vertexList,tempSize);

   double randomTheta;
   double randomPhi;
   double randomRadius;

   //Place vertices
   
   vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
   vertex->setCoordinates(random);
   grid->setVertex(vertexNumber+1, vertex);
   grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
   
   vertexNumber++;

   printf("Points on inclusion sphere surface\n");
   //This is not the right way of doing it. Since you will have too many points close to the apex.
       
   for(int i= 0;i<maxIter;i++){   
     
     randomTheta = 2*myPi*grid->ran1(&randomIntegerOne);
     randomPhi = acos(2*grid->ran1(&randomIntegerTwo)-1);
     
     random.at(1) = this->centre.at(1) + this->radius*cos(randomTheta)*sin(randomPhi);
     random.at(2) = this->centre.at(2) + this->radius*sin(randomTheta)*sin(randomPhi);
     random.at(3) = this->centre.at(3) + this->radius*cos(randomPhi);
     
     //Check if this is far enough from the others.
     
     flag = 0;
     flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, refinement*grid->giveDiameter(random) );
     
          
     if(flag == 0){
       i = 0;
       
       vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
       vertex->setCoordinates(random);
       grid->setVertex(vertexNumber+1, vertex);
       grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
       vertexNumber++;

       //Generate siblings node separated by ITZ
       random.at(1) = this->centre.at(1) + (this->radius+this->itzThickness)*cos(randomTheta)*sin(randomPhi);
       random.at(2) = this->centre.at(2) + (this->radius+this->itzThickness)*sin(randomTheta)*sin(randomPhi);
       random.at(3) = this->centre.at(3) + (this->radius+this->itzThickness)*cos(randomPhi);
       
       vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
       vertex->setCoordinates(random);
       grid->setVertex(vertexNumber+1, vertex);
       grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
       vertexNumber++;
       
     }
   }
   
   printf("Completed inclusion sphere surface\n");

   generator::ensure_size1(grid->vertexList,vertexNumber);
    
  return 1;

}



// int InterfaceSphere :: generatePeriodicPoints()
// {
//   //minimum points
//   int minimumNumberOfPoints = 10; 

//   int numberOfSurfacePoints = 0;
  
//   FloatArray boundaries;
//    grid->defineBoundaries(boundaries);
 
//    int randomIntegerOne= grid->giveRandomInteger()-1;
//    int randomIntegerTwo= grid->giveRandomInteger()-2;
   
//    double myPi = 3.14159265;


// //   int randomIntegerThree= grid->giveRandomInteger()-3;
  
//    FloatArray random(3);
//    FloatArray randomOne(3);

//    int flag;

//    double boundaryFactor = this->refinement;
  
//    Vertex *vertex;
   
//    double distance,shift;
//    double x,y,z;
//    double maxIter= grid->giveMaximumIterations();
//    FloatArray mirroredRandom(3);
//    int vertexNumber = grid->giveNumberOfVertices();
//    int tempSize = 10000000;
//    grid->vertexList->growTo(tempSize);

//    double randomAngleOne;
//    double randomAngleTwo;

//    distance =  boundaryFactor*grid->giveDiameter();   

//    double pointEstimate = myPi*diameter/distance;
//    double deviation;
   
//    int numberOfIntervalsOne,numberOfIntervalsTwo ;
   
//    if(pointEstimate > minimumNumberOfPoints){
//      numberOfIntervalsOne = ceil(pointEstimate);
//    }
//    else{
//      numberOfIntervalsOne = minimumNumberOfPoints;
//    }
   
//    printf("numberOfIntervalsOne = %d\n", numberOfIntervalsOne);

//    int maxSurfacePoints = pow(numberOfIntervalsOne,2.);
   
//    FloatArray surfaceCoords(3*maxSurfacePoints);

//    //Place vertices
//    //Centre
//    random.at(1) = this->centre.at(1);
//    random.at(2) = this->centre.at(2);
//    random.at(3) = this->centre.at(3);
   
//    vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
//    vertex->setCoordinates(random);
//    grid->setVertex(vertexNumber+1, vertex);
//    grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
   
//    vertexNumber++;
   
//   for(int i= 0;i<maxIter;i++){   
//     if(numberOfSurfacePoints == maxSurfacePoints){
//       printf("maximum number of surfacePoints reached\n");
//       break;
//     }

//      randomAngleOne = grid->ran1(&randomIntegerOne)*2.*myPi;     
//      numberOfIntervalsTwo = ceil(sin(randomAngleOne)*numberOfIntervalsOne);
//      randomAngleTwo = grid->ran1(&randomIntegerTwo)*2.*myPi;
     
//        for(int m = 0;m<2;m++){
	 
// 	 randomOne.at(1) = this->centre.at(1) + (0.5*this->diameter + this->itzThickness)*cos(randomAngleTwo)*sin(randomAngleOne);
// 	 randomOne.at(2) = this->centre.at(2) + (0.5*this->diameter + this->itzThickness)*sin(randomAngleTwo)*sin(randomAngleOne);
// 	 randomOne.at(3) = this->centre.at(3) + (0.5*this->diameter + this->itzThickness)*cos(randomAngleOne);
	 
// 	 random.at(1) = this->centre.at(1) + 0.5*this->diameter*cos(randomAngleTwo+deviation)*sin(randomAngleOne);
// 	 random.at(2) = this->centre.at(2) + 0.5*this->diameter*sin(randomAngleTwo+deviation)*sin(randomAngleOne);
// 	 random.at(3) = this->centre.at(3) + 0.5*this->diameter*cos(randomAngleOne);
	 
// 	 //Check if this is far enough from the others.
//  	 flag = 0;

// 	 for (int i=0;i<numberOfSurfacePoints;i++){
// 	   distance = sqrt(pow(surfaceCoords.at(3*i+1)-random.at(1),2.) 
// 			   + pow(surfaceCoords.at(3*i+2)-random.at(2),2.) 
// 			   + pow(surfaceCoords.at(3*i+3)-random.at(3),2.));
	   
// 	   if(distance < boundaryFactor*grid->giveDiameter()){
// 	     flag = 1;
// 	     break;
// 	   }
// 	 }

// 	 if(flag == 0){
// 	   //	   i = 0;

// 	   vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
// 	   vertex->setCoordinates(random);
// 	   grid->setVertex(vertexNumber+1, vertex);
// 	   grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
// 	   vertexNumber++; 

// 	   //Create points
// 	   surfaceCoords.at(3*numberOfSurfacePoints+1) = random.at(1);
// 	   surfaceCoords.at(3*numberOfSurfacePoints+2) = random.at(2);
// 	   surfaceCoords.at(3*numberOfSurfacePoints+3) = random.at(3);
// 	   numberOfSurfacePoints++;
       
// 	   vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
// 	   vertex->setCoordinates(randomOne);
// 	   grid->setVertex(vertexNumber+1, vertex);
// 	   grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,randomOne);
// 	   vertexNumber++; 

// 	 }


//    }

	        
//   }
  
//   printf("Completed sphere loops\n");
  
//   grid->vertexList->growTo(vertexNumber);
  
  
//   return 1;
// }



void
InterfaceSphere :: initializeFrom(GeneratorInputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{ 
    IR_GIVE_FIELD(ir, centre, _IFT_InterfaceSphere_centre);    
    IR_GIVE_FIELD(ir, radius, _IFT_InterfaceSphere_radius); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, _IFT_InterfaceSphere_refine); // Macro
    itzThickness = refinement*diameter;
    IR_GIVE_OPTIONAL_FIELD(ir, itzThickness, _IFT_InterfaceSphere_itz); // Macro
    return;
}



InterfaceSphere *InterfaceSphere :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    InterfaceSphere *sphere;
    
    sphere = new InterfaceSphere(number,grid);

    return sphere;
}
