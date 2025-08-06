#include "sphere.h"
#include "curve.h"
#include "vertex.h"
#include "record.h"


#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#include "octreegridlocalizer.h"
 #include <iostream>
#endif

Sphere :: Sphere(int n,Grid* aGrid) : Region(n,aGrid) //, coordinates()
{
  this->number = n;  
}

Sphere :: ~Sphere()
{
}


int Sphere :: generatePoints()
{

  //minimum points
  int randomIntegerOne= grid->giveRandomInteger()-1;
  int randomIntegerTwo= grid->giveRandomInteger()-2;
  int randomIntegerThree = grid->giveRandomInteger()-3;
   
  double myPi = 3.14159265;
  
  FloatArray random(3);
  
  int flag;

  Vertex *vertex;
  
   double maxIter= grid->giveMaximumIterations();
   int vertexNumber = grid->giveNumberOfVertices();
   int tempSize = 1.e9;
   grid->vertexList->growTo(tempSize);

   double randomTheta;
   double randomPhi;
   double randomRadius;

   //Place vertices

   //A) Centre
   random.at(1) = this->centre.at(1);
   random.at(2) = this->centre.at(2);
   random.at(3) = this->centre.at(3);
   
   vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
   vertex->setCoordinates(random);
   grid->setVertex(vertexNumber+1, vertex);
   grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
   
   vertexNumber++;

   printf("Points on sphere surface\n");
       
   for(int i= 0;i<maxIter;i++){   
     
     randomTheta = 2*myPi*grid->ran1(&randomIntegerOne);
     randomPhi = acos(2*grid->ran1(&randomIntegerTwo)-1);
     
     random.at(1) = this->centre.at(1) + this->radius*cos(randomTheta)*sin(randomPhi);
     random.at(2) = this->centre.at(2) + this->radius*sin(randomTheta)*sin(randomPhi);
     random.at(3) = this->centre.at(3) + this->radius*cos(randomPhi);
     
     //Check if this is far enough from the others.
     
     flag = 0;
     flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, refinement*grid->diameter );
     
     
     
     if(flag == 0){
       i = 0;
       
       vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
       vertex->setCoordinates(random);
       grid->setVertex(vertexNumber+1, vertex);
       grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
       vertexNumber++; 
     }
   }
   
   printf("Completed sphere surface\n");
   printf("number of vertices=%d\n",vertexNumber);
   
   printf("Points in sphere\n");
   
     for(int i= 0;i<maxIter;i++){   
   
       randomTheta = 2*myPi*grid->ran1(&randomIntegerOne);
       randomPhi = acos(2*grid->ran1(&randomIntegerTwo)-1);
       randomRadius = (this->radius-grid->diameter)*pow(grid->ran1(&randomIntegerThree),1./3.);
     	 
       random.at(1) = this->centre.at(1) + randomRadius*cos(randomTheta)*sin(randomPhi);
       random.at(2) = this->centre.at(2) + randomRadius*sin(randomTheta)*sin(randomPhi);
       random.at(3) = this->centre.at(3) + randomRadius*cos(randomPhi);
	 
       //Check if this is far enough from the others.

       flag = 0;
       flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, grid->diameter );



	 if(flag == 0){
	   i = 0;

	   vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
	   vertex->setCoordinates(random);
	   grid->setVertex(vertexNumber+1, vertex);
	   grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
	   vertexNumber++;

	   //Mirror vertices with respect to the  sphere surface
	   random.at(1) = this->centre.at(1) + (2.*this->radius-randomRadius)*cos(randomTheta)*sin(randomPhi);
	   random.at(2) = this->centre.at(2) + (2.*this->radius-randomRadius)*sin(randomTheta)*sin(randomPhi);
	   random.at(3) = this->centre.at(3) + (2.*this->radius-randomRadius)*cos(randomPhi);
	   vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
	   vertex->setCoordinates(random);
	   grid->setVertex(vertexNumber+1, vertex);
	   grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
	   vertexNumber++;
	   
	 }
     }
  
  grid->vertexList->growTo(vertexNumber);
  
  return 1;
}



int Sphere :: generatePeriodicPoints()
{
  
  return 1;
}



IRResultType
Sphere :: initializeFrom(InputRecord *ir)
// Gets from the source line from the data file all the data of the receiver.
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, centre, IFT_Sphere_centre, "centre");
    IR_GIVE_FIELD(ir, radius, IFT_Sphere_radius, "radius");
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, IFT_Sphere_refine, "refine");
    
    return IRRT_OK;

}



Sphere *Sphere :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    Sphere *sphere;
    
    sphere = new Sphere(number,grid);

    return sphere;
}
