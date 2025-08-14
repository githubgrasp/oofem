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
   
  double myPi = 3.14159265;
  
  oofem::FloatArray random(3);
  
  int flag;

  
   double maxIter= grid->giveMaximumIterations();


   double randomTheta;
   double randomPhi;

     grid->addVertex(random);
		    

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
       
     grid->addVertex(random);
       

       //Generate siblings node separated by ITZ
       random.at(1) = this->centre.at(1) + (this->radius+this->itzThickness)*cos(randomTheta)*sin(randomPhi);
       random.at(2) = this->centre.at(2) + (this->radius+this->itzThickness)*sin(randomTheta)*sin(randomPhi);
       random.at(3) = this->centre.at(3) + (this->radius+this->itzThickness)*cos(randomPhi);


       grid->addVertex(random);
		    
              
     }
   }
   
   printf("Completed inclusion sphere surface\n");

  return 1;

}



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


