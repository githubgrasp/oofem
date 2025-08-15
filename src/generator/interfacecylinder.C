#include "interfacecylinder.h"
#include "curve.h"
#include "vertex.h"
#include "record.h"


#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#include "octreegridlocalizer.h"
#endif

InterfaceCylinder :: InterfaceCylinder(int n,Grid* aGrid) : Inclusion(n,aGrid) //, coordinates()
{
  this->number = n;  
}

InterfaceCylinder :: ~InterfaceCylinder()
// Destructor.
{
}

int InterfaceCylinder :: generatePoints()
{

  printf("In generatePoints in InterfaceCylinder\n");
  
  
  double myPi = 3.14159265;
  
  
  oofem::FloatArray random(3);
  int flag;
  

  double randomTheta;
  double randomXCoord;

   printf("In point generation of interface cylinder\n");
   printf("Vertex number at the start of interface cylinder = %d\n", grid->giveNumberOfVertices());
   printf("Generate midpoints of end circles\n");
   printf("First point\n");
   random.at(1) = this->line.at(1);
   random.at(2) = this->line.at(2);
   random.at(3) = this->line.at(3);

   flag = 0;
   flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, refinement*grid->giveDiameter(random) );

   if(flag == 0){

     grid->addVertex(random);
          
   }

   printf("Second point\n");
   random.at(1) = this->line.at(4);
   random.at(2) = this->line.at(5);
   random.at(3) = this->line.at(6);
   
   flag = 0;
   flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, refinement*grid->giveDiameter(random) );

   if(flag == 0){


     grid->addVertex(random);
     
   }

  
   printf("Generate regular points on lines\n");
   printf("Start with circe one\n");
   double spacing = refinement*grid->diameter;
   int intervals = ceil(2*myPi*this->radius/spacing);
   printf("Spacing for circe one = %e and intervals = %d\n", spacing, intervals);
   random.at(1) = this->line.at(1);
   for(int i= 0;i<intervals;i++){     
     randomTheta = 2.*myPi/intervals*i; 
     random.at(2) = this->line.at(2) + this->radius*cos(randomTheta);
     random.at(3) = this->line.at(3) + this->radius*sin(randomTheta);

     flag = 0;
     flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, 0.9*spacing );
     
     if(flag == 0){
       
       grid->addVertex(random);

       //Generate point on other side
       random.at(2) = this->line.at(2) + (this->radius+this->itzThickness)*cos(randomTheta);
       random.at(3) = this->line.at(3) + (this->radius+this->itzThickness)*sin(randomTheta);

       grid->addVertex(random);
       
     }
     else{
       printf("SHOULD NOT BE HERE\n");
     }

   }
   printf("circle one completed. Vertex number = %d\n", grid->giveNumberOfVertices());

   
   printf("Start with circle two\n");
   random.at(1) = this->line.at(4);
   spacing = refinement*grid->diameter;
   intervals = ceil(2*myPi*this->radius/spacing);
   printf("Spacing for circe one = %e and intervals = %d\n", spacing, intervals);
   for(int i= 0;i<intervals;i++){      
     randomTheta = 2*myPi/intervals*i;
     random.at(2) = this->line.at(5) + this->radius*cos(randomTheta);
     random.at(3) = this->line.at(6) + this->radius*sin(randomTheta);

     //Check if this is far enough from the others.
     
     flag = 0;
     flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, 0.9*spacing );
     
     if(flag == 0){

       grid->addVertex(random);

       //Generate point on other side
       random.at(2) = this->line.at(2) + (this->radius+itzThickness)*cos(randomTheta);
       random.at(3) = this->line.at(3) + (this->radius+itzThickness)*sin(randomTheta);

       grid->addVertex(random);

     }
     else{
       printf("SHOULD NOT BE HERE\n");
     }
   }

   printf("circle two completed. Vertex number = %d\n", grid->giveNumberOfVertices());

   
   printf("now regular points on the surface of the cylinder\n");
   spacing = refinement*grid->diameter;
   double intervalsX = ceil((line.at(4)-line.at(1))/(0.9*spacing));

   for(int i= 1;i<intervalsX;i++) {
       for ( int k = 0; k < intervals; k++ ) {
           randomXCoord   = this->line.at( 1 ) + ( line.at( 4 ) - line.at( 1 ) ) / intervalsX * i;
           random.at( 1 ) = randomXCoord;

           randomTheta = 2 * myPi / intervals * k;
           random.at( 2 ) = this->line.at( 2 ) + this->radius * cos( randomTheta );
           random.at( 3 ) = this->line.at( 3 ) + this->radius * sin( randomTheta );

           // Check if this is far enough from the others.

           flag = 0;
           flag = grid->giveGridLocalizer()->checkNodesWithinBox( random, 0.7 * spacing );

           if ( flag == 0 ) {

               grid->addVertex( random );

               random.at( 2 ) = this->line.at( 2 ) + ( this->radius + itzThickness ) * cos( randomTheta );
               random.at( 3 ) = this->line.at( 3 ) + ( this->radius + itzThickness ) * sin( randomTheta );

               grid->addVertex( random );

               if ( k == 0 ) {
                   random.at( 2 ) = this->line.at( 2 );
                   random.at( 3 ) = this->line.at( 3 );

                   grid->addVertex( random );

               }
               // Mirror point with respect to first end circles
               random.at( 1 ) = 2. * line.at( 1 ) - randomXCoord;
               random.at( 2 ) = this->line.at( 2 ) + this->radius * cos( randomTheta );
               random.at( 3 ) = this->line.at( 3 ) + this->radius * sin( randomTheta );

               grid->addVertex( random );

               random.at( 2 ) = this->line.at( 2 ) + ( this->radius + this->itzThickness ) * cos( randomTheta );
               random.at( 3 ) = this->line.at( 3 ) + ( this->radius + this->itzThickness ) * sin( randomTheta );

               grid->addVertex( random );

               random.at( 2 ) = this->line.at( 2 );
               random.at( 3 ) = this->line.at( 3 );

               grid->addVertex( random );

               // Mirror point with respect to second end circle
               random.at( 1 ) = 2. * line.at( 4 ) - randomXCoord;
               random.at( 2 ) = this->line.at( 2 ) + this->radius * cos( randomTheta );
               random.at( 3 ) = this->line.at( 3 ) + this->radius * sin( randomTheta );

               grid->addVertex( random );

               random.at( 2 ) = this->line.at( 2 ) + ( this->radius + this->itzThickness ) * cos( randomTheta );
               random.at( 3 ) = this->line.at( 3 ) + ( this->radius + this->itzThickness ) * sin( randomTheta );

               grid->addVertex( random );


               random.at( 2 ) = this->line.at( 2 );
               random.at( 3 ) = this->line.at( 3 );

               grid->addVertex( random );

           } else {
               printf( "SHOULD NOT BE HERE\n" );
           }
       }
   }

   printf("surface of interface cylinder completed. Vertex number = %d\n", grid->giveNumberOfVertices());

   printf("Completed inclusion loops\n");

   return 1;
}


void
InterfaceCylinder :: initializeFrom(GeneratorInputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{
    IR_GIVE_FIELD(ir, line, _IFT_InterfaceCylinder_line);    
    IR_GIVE_FIELD(ir, radius, _IFT_InterfaceCylinder_radius); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, _IFT_InterfaceCylinder_refine); // Macro
    itzThickness = refinement*diameter;
    IR_GIVE_OPTIONAL_FIELD(ir, itzThickness, _IFT_InterfaceCylinder_itz); // Macro

    return;
}



InterfaceCylinder *InterfaceCylinder :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    InterfaceCylinder *cylinder;
    
    cylinder = new InterfaceCylinder(number,grid);

    return cylinder;
}
