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

Region :: Region(int n,Grid* aGrid) : GridComponent(n,aGrid) //, coordinates()
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
  Surface *surface;
  Curve *curve;
  Vertex *vertex;
  
  //  int firstVertex = giveLocalVertex(1);
  int localSurface;
  int localCurve;
  int localVertex;
  IntArray curves;
  double x,y,z;

  FloatArray boundaries;
  grid->defineBoundaries(boundaries);
 
  //int randomIntegerOne= grid->giveRandomInteger()-1;
  //int randomIntegerTwo= grid->giveRandomInteger()-2;
  //int randomIntegerThree= grid->giveRandomInteger()-3;
  
  FloatArray random(3);
  // int flag;

  //double boundaryFactor = this->refinement;
   int n1edges = this->xedges;   
   int n2edges = this->yedges;   
   int n3edges = this->zedges;  
   double n1length = this ->xlength; 
   double n2length = this ->ylength; 
   double n3length = this ->zlength; 
   //std::cout << n1length << ' ' << n3length << ' ' << n3length <<'\n';
   //std::cout << n1edges << ' ' << n3edges << ' ' << n3edges <<'\n';
  double distance;
  //  double x,y,z;
  // double maxIter = grid->giveMaximumIterations();
  FloatArray mirroredRandom(3);

  int vertexNumber = grid->vertexList->giveSize();
  
  //std::cout<<"kab00m "<< n1length <<" " << n2length <<" " << n3length <<'\n';
  //Try to overcome mirroring by incrasing region size
  /* boundaries.at(1) -= 5* grid->giveDiameter();
  boundaries.at(3) -= 5* grid->giveDiameter();
  boundaries.at(5) -= 5* grid->giveDiameter();
 
  boundaries.at(2) += 5* grid->giveDiameter();
  boundaries.at(4) += 5* grid->giveDiameter();
  boundaries.at(6) += 5* grid->giveDiameter();*/

 //   FloatArray centre(3);
 //   centre.at(1) = (boundaries.at(2)-boundaries.at(1))/2.;
 //    centre.at(2) = (boundaries.at(4)-boundaries.at(3))/2.;
 //   centre.at(3) = (boundaries.at(6)-boundaries.at(5))/2.;

  int tempSize= 5000000;
  grid->vertexList->growTo(tempSize);
  for(int i= 1;i < n3edges*2 ;i++){

    random.at(3) =  i*n3length/(2*n3edges);
    
    //Generate precise points
    if (i % 2){ // i odd

      //new
 for (int j =0 ; j < n2edges ; j++){
	if(j % 2){
	  random.at(2) = (0.66666666 + j)*n2length/n2edges;
	  int k =1;
	  while( k < n1edges-1 ){

	    random.at(1) = (1.5+k)*n1length/n1edges;
	  k++;
	 
	  
	  vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
	  vertex->setCoordinates(random);
	  grid->setVertex(vertexNumber+1, vertex);
	//  
	  
	  grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
      
	  //i = 0;*/
	  vertexNumber++;
	  if ( random.at(1) <= 0. || random.at(1) >= n1length || random.at (2) <= 0. || random.at(2) >= n2length ||  random.at (3) <= 0. || random.at(3) >= n3length){
	    std::cout<<"FAIL FAIL 1"<<'\n';
	  }
	  
	  //std::cout << "x " << random.at(1) << " y "<< random.at(2) << " z "<< random.at(3)<< '\n'; 
	  }
	  
	}
	  
       	else{
	 
	  for(int l = 0 ; l < n1edges ; l++){	 
	 	   
	  random.at(1) = (0.5+l)*n1length/n1edges;
	  
	  random.at(2) = (0.66666666 + j)*n2length/n2edges;
	  //std::cout<< " show it ! " <<" y "<<random.at(2) << '\n';              	  
	  vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
	  vertex->setCoordinates(random);
	  grid->setVertex(vertexNumber+1, vertex);
	  //  
	
	  grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
	  
	  //i = 0;*/
	  vertexNumber++;
	  //std::cout<< " show it ! " << "x "<< random.at(1) <<" y "<<random.at(2) << " z " << random.at(3) <<'\n';          
	  if ( random.at(1) <= 0. || random.at(1) >= n1length  || random.at (2) <= 0. || random.at(2) >= n2length ||  random.at (3) <= 0. || random.at(3) >= n3length){
	    std::cout<<"FAIL FAIL 2"<<'\n';
	  }

    	  
	  }
	  
	}


		
      //std::cout << "x " << random.at(1) << " y "<< random.at(2) << " z "<< random.at(3)<< '\n'; 
      }
    
    }
	
   



    //end new

    else if(i % 4 == 0){ //i is a multiple of four

      for( int n = 1 ; n < n2edges ; n++){
	
	if (n % 2) { 
	  
	  for(int m = 0 ; m < n1edges ; m++){
	    //   std::cout<<"kab00m "<<m<<'\n'; 	          	
	    //std::cout<<"kab00m "<<'\n';
	    random.at(1) =  (0.5 + m)*n1length/n1edges ;
	    // std::cout<< " show it ! " << n << '\n';    
	    random.at(2) =  n*n2length/n2edges;
	   
	    vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
	    vertex->setCoordinates(random);
	    grid->setVertex(vertexNumber+1, vertex);
	    //  
	
	    grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
	
	    //i = 0;*/
	    vertexNumber++;
	    if ( random.at(1) <= 0. || random.at(1) >= 1. || random.at (2) <= 0. || random.at(2) >= 1.03923 ||  random.at (3) <= 0. || random.at(3) >= 0.979795){
	      std::cout<<"FAIL FAIL 3"<<'\n';
	    }
	    
	  }
	  
     	}
	else{
	  for(int s = 1 ; s < n1edges ; s++){
	    random.at(1) =   s*n1length/n1edges ;
	    //std::cout<<"da en " << n << '\n';
	    random.at(2) =  n*n2length/n2edges;
	    //std::cout<< " show it ! " << "x "<< random.at(1) <<" y "<<random.at(2) << '\n';     
	    vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
	    vertex->setCoordinates(random);
	    grid->setVertex(vertexNumber+1, vertex);
	  //  
	  
	  grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
	  
	  //i = 0;*/
	  vertexNumber++;
	  if ( random.at(1) <= 0. || random.at(1) >= 1. || random.at (2) <= 0. || random.at(2) >= 1.03923 ||  random.at (3) <= 0. || random.at(3) >= 0.979795){
	    std::cout<<"FAIL FAIL 4"<<'\n';
	  }
	  
	  }
	  
	}
	
 
	//maybe heere }


	
	
      //std::cout << "x " << random.at(1) << " y "<< random.at(2) << " z "<< random.at(3)<< '\n'; 
      }
    
    }
    else{

 for (int j =0 ; j < n2edges ; j++){
	if(j % 2){
	  random.at(2) = (0.333333333 + j)*n2length/n2edges;
	  int k =1;
	  while( k < n1edges ){

	  random.at(1) = k*n1length/n1edges;
	  k++;
	 
	  
	  vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
	  vertex->setCoordinates(random);
	  grid->setVertex(vertexNumber+1, vertex);
	//  
	  
	  grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
      
	  //i = 0;*/
	  vertexNumber++;
	  if ( random.at(1) <= 0. || random.at(1) >= n1length || random.at (2) <= 0. || random.at(2) >= n2length ||  random.at (3) <= 0. || random.at(3) >= n3length){
	    std::cout<<"FAIL FAIL 1"<<'\n';
	  }
	  
	  //std::cout << "x " << random.at(1) << " y "<< random.at(2) << " z "<< random.at(3)<< '\n'; 
	  }
	  
	}
	  
       	else{
	 
	  for(int l = 0 ; l < n1edges ; l++){	 
	 	   
	  random.at(1) = (0.5+l)*n1length/n1edges;
	  
	  random.at(2) = (0.333333333 + j)*n2length/n2edges;
	  //std::cout<< " show it ! " <<" y "<<random.at(2) << '\n';              	  
	  vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
	  vertex->setCoordinates(random);
	  grid->setVertex(vertexNumber+1, vertex);
	  //  
	
	  grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
	  
	  //i = 0;*/
	  vertexNumber++;
	  //std::cout<< " show it ! " << "x "<< random.at(1) <<" y "<<random.at(2) << " z " << random.at(3) <<'\n';          
	  if ( random.at(1) <= 0. || random.at(1) >= n1length  || random.at (2) <= 0. || random.at(2) >= n2length ||  random.at (3) <= 0. || random.at(3) >= n3length){
	    std::cout<<"FAIL FAIL 2"<<'\n';
	  }

    	  
	  }
	  
	}

 }

    }

   
    
      // std::cout<<"kab00m "<< random.at(1) <<' '<< random.at(2) << ' ' << random.at(3) <<' ' << vertexNumber <<'\n'; 
    //Check if this is far enough from the others
    //flag = 0;

    //     //    Debug: Introduce mesh refinement temporarily
    //     distance = sqrt(pow(centre.at(1)-random.at(1),2.)+pow(centre.at(2)-random.at(2),2.)+pow(centre.at(3)-random.at(3),2.));
    
    //     if(distance <0.1*centre.at(1)){
    //       flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, 0.5*boundaryFactor*grid->giveDiameter());
    //     }
    //     else{
    //       flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, boundaryFactor*grid->giveDiameter());
    //     }
    
    /*flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, boundaryFactor*grid->giveDiameter());

    //    flag = grid->giveGridLocalizer()->checkNodesWithinBox(random, boundaryFactor*grid->giveDiameter());
    if(flag == 0){*/

    /* vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
      vertex->setCoordinates(random);
      grid->setVertex(vertexNumber+1, vertex);
      //  
      
      grid->giveGridLocalizer()->insertSequentialNode(vertexNumber+1,random);
      
      //i = 0;
      vertexNumber++;*/
    
      // Old mirroring stuff 
//      //Mirrors in this plane (anti-clockwise)
//       // Zone 1     
//       vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
//       mirroredRandom = random;
//       mirroredRandom.at(1) = random.at(1) -2.*(random.at(1) - boundaries.at(2));
//       vertex->setCoordinates(mirroredRandom);
//       grid->setVertex(vertexNumber+1, vertex);      
//       vertexNumber++;
//       //Zone 2
//       vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
//       mirroredRandom = random;
//       mirroredRandom.at(1) = random.at(1) - 2.*(random.at(1) - boundaries.at(2));
//       mirroredRandom.at(2) = random.at(2) - 2.*(random.at(2) - boundaries.at(4));
//       vertex->setCoordinates(mirroredRandom);
//       grid->setVertex(vertexNumber+1, vertex);      
//       vertexNumber++;

//       //Zone 3
//       vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
//       mirroredRandom = random;
//       mirroredRandom.at(1) = random.at(1);
//       mirroredRandom.at(2) = random.at(2) - 2.*(random.at(2) - boundaries.at(4));
//       vertex->setCoordinates(mirroredRandom);
//       grid->setVertex(vertexNumber+1, vertex);      
//       vertexNumber++;
      
//       //Zone 4
//       vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
//       mirroredRandom = random;
//       mirroredRandom.at(1) = random.at(1) - 2.*(random.at(1) - boundaries.at(1));
//       mirroredRandom.at(2) = random.at(2) - 2.*(random.at(2) - boundaries.at(4));
//       vertex->setCoordinates(mirroredRandom);
//       grid->setVertex(vertexNumber+1, vertex);      
//       vertexNumber++;

//       //Zone 5
//       vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
//       mirroredRandom = random;
//       mirroredRandom.at(1) = random.at(1) - 2.*(random.at(1) - boundaries.at(1));
//       mirroredRandom.at(2) = random.at(2);
//       vertex->setCoordinates(mirroredRandom);
//       grid->setVertex(vertexNumber+1, vertex);      
//       vertexNumber++;

//       //Zone 6
//       vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
//       mirroredRandom = random;
//       mirroredRandom.at(1) = random.at(1) - 2.*(random.at(1) - boundaries.at(1));
//       mirroredRandom.at(2) = random.at(2) - 2.*(random.at(2) - boundaries.at(3));
//       vertex->setCoordinates(mirroredRandom);
//       grid->setVertex(vertexNumber+1, vertex);      
//       vertexNumber++;

//       //Zone 7
//       vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
//       mirroredRandom = random;
//       mirroredRandom.at(1) = random.at(1);
//       mirroredRandom.at(2) = random.at(2) - 2.*(random.at(2) - boundaries.at(3));
//       vertex->setCoordinates(mirroredRandom);
//       grid->setVertex(vertexNumber+1, vertex);      
//       vertexNumber++;

//       //Zone 8
//       vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
//       mirroredRandom = random;
//       mirroredRandom.at(1) = random.at(1) - 2.*(random.at(1) - boundaries.at(2));
//       mirroredRandom.at(2) = random.at(2) - 2.*(random.at(2) - boundaries.at(3));
//       vertex->setCoordinates(mirroredRandom);
//       grid->setVertex(vertexNumber+1, vertex);      
//       vertexNumber++;

//       //Repeat this with two different mother vertices


//       for(int k = 0;k<2;k++){
// 	//Mother vertex
// 	vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
// 	mirroredRandom = random;
// 	mirroredRandom.at(3) = random.at(3)-2.*(random.at(3) - boundaries.at(5+k));
// 	vertex->setCoordinates(mirroredRandom);
// 	grid->setVertex(vertexNumber+1, vertex);
// 	i = 0;
// 	vertexNumber++;

// 	//Mirrors in this plane (anti-clockwise)
// 	// Zone 1     
// 	vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
// 	mirroredRandom = random;
// 	mirroredRandom.at(1) = random.at(1) -2.*(random.at(1) - boundaries.at(2));
// 	mirroredRandom.at(3) = random.at(3)-2.*(random.at(3) - boundaries.at(5+k));
// 	vertex->setCoordinates(mirroredRandom);
// 	grid->setVertex(vertexNumber+1, vertex);      
// 	vertexNumber++;

// 	//Zone 2
// 	vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
// 	mirroredRandom = random;
// 	mirroredRandom.at(1) = random.at(1) - 2.*(random.at(1) - boundaries.at(2));
// 	mirroredRandom.at(2) = random.at(2) - 2.*(random.at(2) - boundaries.at(4));
// 	mirroredRandom.at(3) = random.at(3)-2.*(random.at(3) - boundaries.at(5+k));
// 	vertex->setCoordinates(mirroredRandom);
// 	grid->setVertex(vertexNumber+1, vertex);      
// 	vertexNumber++;

// 	//Zone 3
// 	vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
// 	mirroredRandom = random;
// 	mirroredRandom.at(1) = random.at(1);
// 	mirroredRandom.at(2) = random.at(2) - 2.*(random.at(2) - boundaries.at(4));
// 	mirroredRandom.at(3) = random.at(3)-2.*(random.at(3) - boundaries.at(5+k));
// 	vertex->setCoordinates(mirroredRandom);
// 	grid->setVertex(vertexNumber+1, vertex);      
// 	vertexNumber++;

// 	//Zone 4
// 	vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
// 	mirroredRandom = random;
// 	mirroredRandom.at(1) = random.at(1) - 2.*(random.at(1) - boundaries.at(1));
// 	mirroredRandom.at(2) = random.at(2) - 2.*(random.at(2) - boundaries.at(4));
// 	mirroredRandom.at(3) = random.at(3)-2.*(random.at(3) - boundaries.at(5+k));
// 	vertex->setCoordinates(mirroredRandom);
// 	grid->setVertex(vertexNumber+1, vertex);      
// 	vertexNumber++;

// 	//Zone 5
// 	vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
// 	mirroredRandom = random;
// 	mirroredRandom.at(1) = random.at(1) - 2.*(random.at(1) - boundaries.at(1));
// 	mirroredRandom.at(2) = random.at(2);
// 	mirroredRandom.at(3) = random.at(3)-2.*(random.at(3) - boundaries.at(5+k));
// 	vertex->setCoordinates(mirroredRandom);
// 	grid->setVertex(vertexNumber+1, vertex);      
// 	vertexNumber++;

// 	//Zone 6
// 	vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
// 	mirroredRandom = random;
// 	mirroredRandom.at(1) = random.at(1) - 2.*(random.at(1) - boundaries.at(1));
// 	mirroredRandom.at(2) = random.at(2) - 2.*(random.at(2) - boundaries.at(3));
// 	mirroredRandom.at(3) = random.at(3)-2.*(random.at(3) - boundaries.at(5+k));
// 	vertex->setCoordinates(mirroredRandom);
// 	grid->setVertex(vertexNumber+1, vertex);      
// 	vertexNumber++;

// 	//Zone 7
// 	vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
// 	mirroredRandom = random;
// 	mirroredRandom.at(1) = random.at(1);
// 	mirroredRandom.at(2) = random.at(2) - 2.*(random.at(2) - boundaries.at(3));
// 	mirroredRandom.at(3) = random.at(3) - 2.*(random.at(3) - boundaries.at(5+k));
// 	vertex->setCoordinates(mirroredRandom);
// 	grid->setVertex(vertexNumber+1, vertex);      
// 	vertexNumber++;

// 	//Zone 8
// 	vertex = ( Vertex * ) ( Vertex(vertexNumber+1,grid).ofType()); 
// 	mirroredRandom = random;
// 	mirroredRandom.at(1) = random.at(1) - 2.*(random.at(1) - boundaries.at(2));
// 	mirroredRandom.at(2) = random.at(2) - 2.*(random.at(2) - boundaries.at(3));
// 	mirroredRandom.at(3) = random.at(3) - 2.*(random.at(3) - boundaries.at(5+k));
// 	vertex->setCoordinates(mirroredRandom);
// 	grid->setVertex(vertexNumber+1, vertex);      
// 	vertexNumber++;

//       }
      //}



  }

   grid->vertexList->growTo(vertexNumber);
 

  //  printf("numberOfVertices for region %d = %d\n",this->number, vertexNumber);
  return 1;
}



IRResultType
Region :: initializeFrom(InputRecord *ir)
// Gets from the source line from the data file all the data of the receiver.
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    //    int j, size;
    //    FloatArray vertices;
    // IntArray *dofIDArry;
    
    IR_GIVE_FIELD(ir, surfaces, IFT_Region_surfaces, "surfaces"); // Macro
    refinement = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, IFT_Region_refine, "refine"); // Macro
   
    IR_GIVE_FIELD(ir, xedges, IFT_xedges, "xedges"); // Macro
    IR_GIVE_FIELD(ir, yedges, IFT_yedges, "yedges"); // Macro
    IR_GIVE_FIELD(ir, zedges, IFT_zedges, "zedges"); // Macro
    IR_GIVE_FIELD(ir, xlength, IFT_xlength, "xlength"); // Macro
    IR_GIVE_FIELD(ir, ylength, IFT_ylength, "ylength"); // Macro
    IR_GIVE_FIELD(ir, zlength, IFT_zlength, "zlength"); // Macro
    return IRRT_OK; 
}



Region *Region :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    Region *region;
    
    region = new Region(number,grid);

    return region;
}
