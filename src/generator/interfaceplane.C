#include "interfaceplane.h"
#include "curve.h"
#include "vertex.h"
#include "record.h"


#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#include "octreegridlocalizer.h"
#endif

InterfacePlane :: InterfacePlane(int n,Grid* aGrid) : Inclusion(n,aGrid) //, coordinates()
{
  this->number = n;  
}

InterfacePlane :: ~InterfacePlane()
// Destructor.
{
}


int InterfacePlane :: generatePoints()
{
  //minimum points
  int minimumNumberOfPoints = 8; 
  
  FloatArray boundaries;
   grid->defineBoundaries(boundaries);
 
   int randomIntegerOne= grid->giveRandomInteger()-1;
   int randomIntegerTwo= grid->giveRandomInteger()-2;
   
   double myPi = 3.14159265;


//   int randomIntegerThree= grid->giveRandomInteger()-3;
  
   FloatArray random(3);
   int flag;

   double boundaryFactor = this->refinement;
  
   Vertex *vertex;
   
   double distance, newDistance,shift;
   double x,y,z;
   double maxIter = grid->giveMaximumIterations();
   FloatArray mirroredRandom(3);
   int vertexNumber = grid->giveNumberOfVertices();
   int tempSize = 10000000;
   grid->vertexList->growTo(tempSize);

   double randomAngleOne;
   double randomAngleTwo;
   double randomRadius;

   distance =  boundaryFactor*grid->diameter;   

   double pointEstimate = myPi*diameter/distance;
   
   int numberOfIntervalsOne,numberOfIntervalsTwo, newIntervalsTwo ;

   if(pointEstimate > minimumNumberOfPoints){
     numberOfIntervalsTwo = ceil(pointEstimate);
   }
   else{
     numberOfIntervalsTwo = minimumNumberOfPoints;
     distance = myPi*diameter/minimumNumberOfPoints;
   }
   
   double length = sqrt(pow(this->line.at(1)-this->line.at(4),2.));

   numberOfIntervalsOne = ceil(numberOfIntervalsTwo*length/(myPi*diameter));

   
   //This implementation is hardcoding the orientation of the plane
  
   //Two loops
   for(int i= 0;i<=numberOfIntervalsOne;i++){
     random.at(1) = line.at(1) + i*length/numberOfIntervalsOne;     
     //     printf("random.at(1) = %e\n", random.at(1));

     //Place Midpoint
     random.at(2) = this->line.at(2);
     random.at(3) = this->line.at(3);
     
     vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
     vertex->setCoordinates(random);
     grid->setVertex(vertexNumber+1, vertex);
     grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
     
     //	 i = 0;
     vertexNumber++;
     
          
     for(int k = 0;k<numberOfIntervalsTwo;k++){
  
       //Generate random point
     
       randomAngleTwo = k*2.*myPi/numberOfIntervalsTwo;
              
       for(int m = 0;m<2;m++){
	 random.at(2) = this->line.at(2) + (0.5*this->diameter + (1.-m)*this->itzThickness)*cos(randomAngleTwo);
	 random.at(3) = this->line.at(3) + (0.5*this->diameter + (1.-m)*this->itzThickness)*sin(randomAngleTwo);
	 
	 //Check if this is far enough from the others
 	 flag = 0;
	        
	 vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
	 vertex->setCoordinates(random);
	 grid->setVertex(vertexNumber+1, vertex);
	 grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
	 
	 
	 //	 i = 0;
	 vertexNumber++;
	 if(i==0||i==numberOfIntervalsOne){
	   printf("VertexNumbers for surface elements of cylindrical inclusion = %d\n", vertexNumber);
	 }
 
       }//end of two layers
     }//end of iterations
     
   }
   
   printf("Completed inclusion loops\n");

   grid->vertexList->growTo(vertexNumber);
   
   return 1;
}


IRResultType
InterfacePlane :: initializeFrom(InputRecord *ir)
// Gets from the source line from the data file all the data of the receiver.
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // int j, size;
    // FloatArray vertices;
    // IntArray *dofIDArry;
    IR_GIVE_FIELD(ir, line, IFT_InterfacePlane_line, "line");    
    IR_GIVE_FIELD(ir, diameter, IFT_InterfacePlane_diameter, "diameter"); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, IFT_InterfacePlane_refine, "refine"); // Macro
    itzThickness = refinement*diameter;
    IR_GIVE_OPTIONAL_FIELD(ir, itzThickness, IFT_InterfacePlane_itz, "itz"); // Macro
    return IRRT_OK;

}



InterfacePlane *InterfacePlane :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    InterfacePlane *plane;
    
    plane = new InterfacePlane(number,grid);

    return plane;
}
