#include "grid.h"
#include "vertex.h"
#include "curve.h"
#include "surface.h"
#include "prism.h"
#include "region.h"
#include "inclusion.h"
#include "refinement.h"
#include "octreegridlocalizer.h"
//#include "oofem_limits.h"
#include "datareader.h"
#include "domain.h"
#include "sphere.h"
#include "interfacesphere.h"
#include "interfacecylinder.h"
#include "refineprism.h"
#include "cylinder.h"
#include <iostream>


#ifndef __MAKEDEPEND
 #include <string.h>
 #include <stdarg.h>
 #ifdef HAVE_STRINGS_H
 #endif
 #include <strings.h>
 #include <math.h>
 #include <stdio.h>
 #include <cstdlib>
 #include <ctype.h>
#endif

using namespace oofem;

//Definitions for random number generator
#define IA 16807
#define IM 2147483647
#define AM ( 1.0 / IM )
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV ( 1 + ( IM - 1 ) / NTAB )
#define EPS 1.2e-7
#define RNMX ( 1.0 - EPS )



Grid :: Grid(int i)
// Constructor. Creates a new domain.
{
    //  maxIter = 0;

    /* vertexList           = new AList< Vertex >(0); */
    /* inputVertexList     = new AList< Vertex >(0); */
    /* controlVertexList = new AList< Vertex >(0); */
    /* curveList        = new AList< Curve >(0); */
    /* surfaceList          = new AList< Surface >(0); */
    /* regionList                = new AList< Region >(0); */
    /* inclusionList                = new AList< Inclusion >(0); */
    /* refinementList                = new AList< Refinement >(0); */


    
    gridLocalizer      = NULL;

    TOL = 0.;
}

Grid :: ~Grid()
// Destructor.
{
    /* delete vertexList; */
    /* delete inputVertexList; */
    /* delete controlVertexList; */
    /* delete curveList; */
    /* delete surfaceList; */
    /* delete regionList; */
    /* delete inclusionList; */
    /* delete refinementList; */

    for (auto v : vertexList) delete v;
    for (auto v : inputVertexList) delete v;
    for (auto v : controlVertexList) delete v;
    for (auto c : curveList) delete c;
    for (auto s : surfaceList) delete s;
    for (auto r : regionList) delete r;
    for (auto i : inclusionList) delete i;
    for (auto r : refinementList) delete r;

    delete gridLocalizer;
    
    delete gridLocalizer;
    
}



GridLocalizer *Grid :: giveGridLocalizer()
//
// return connectivity Table - if no defined - creates new one
//
{
    if ( gridLocalizer == NULL ) {
        gridLocalizer = new OctreeGridLocalizer(1, this);
    }

    return gridLocalizer;
}



Vertex *Grid::giveVertex(int n)
{
    if (n >= 1 && n <= static_cast<int>(vertexList.size()) && vertexList[n - 1] != nullptr) {
        return vertexList[n - 1];  // 1-based to 0-based
    } else {
        printf("giveVertex: undefined vertex (%d)\n", n);
        exit(1);
    }

    return nullptr;
}


/* Vertex *Grid :: giveVertex(int n) */
/* // Returns the n-th vertex.  */
/* { */
/*     if ( vertexList->includes(n) ) { */
/*         return vertexList->at(n); */
/*     } else { */
/*         printf("giveVertex: undefined vertex (%d)", n); */
/*         exit(1); */
/*     } */

/*     return NULL; */
/* } */


Vertex *Grid :: giveInputVertex(int n)
// Returns the n-th vertex. Creates this node if it does not exist yet.
{

    if (n >= 1 && n <= static_cast<int>(inputVertexList.size()) && inputVertexList[n - 1] != nullptr) {
        return inputVertexList[n - 1];  // 1-based to 0-based
    } else {
        printf("giveInputVertex: undefined inputVertex (%d)\n", n);
        exit(1);
    }
    return nullptr;
}



Vertex *Grid :: giveControlVertex(int n)
// Returns the n-th vertex. Creates this node if it does not exist yet.
{

    if (n >= 1 && n <= static_cast<int>(controlVertexList.size()) && controlVertexList[n - 1] != nullptr) {
        return controlVertexList[n - 1];  // 1-based to 0-based
    } else {
        printf("giveControlVertex: undefined controlVertex (%d)\n", n);
        exit(1);
    }
    return nullptr;

  /* if ( controlVertexList->includes(n) ) { */
    /*     return controlVertexList->at(n); */
    /* } else { */
    /*     printf("giveControlVertex: undefined inputVertex (%d)", n); */
    /*     exit(1); */
    /* } */

    /* return NULL; */
}



Curve *Grid :: giveCurve(int n)
// Returns the n-th element. Generates error if it is not defined yet.
{

    if (n >= 1 && n <= static_cast<int>(curveList.size()) && curveList[n - 1] != nullptr) {
        return curveList[n - 1];  // 1-based to 0-based
    } else {
        printf("giveCurve: undefined curve (%d)\n", n);
        exit(1);
    }
    return nullptr;

  
    /* if ( curveList->includes(n) ) { */
    /*     return curveList->at(n); */
    /* } else { */
    /*     printf("giveCurve: undefined curve (%d)", n); */
    /*     exit(1); */
    /* } */

    /* return NULL; */
}

Surface *Grid :: giveSurface(int n)
// Returns the n-th element. Generates error if it is not defined yet.
{

    if (n >= 1 && n <= static_cast<int>(surfaceList.size()) && surfaceList[n - 1] != nullptr) {
        return surfaceList[n - 1];  // 1-based to 0-based
    } else {
        printf("giveSurface: undefined surface (%d)\n", n);
        exit(1);
    }
    return nullptr;


  /* if ( surfaceList->includes(n) ) { */
  /*       return surfaceList->at(n); */
  /*   } else { */
  /*       printf("giveSurface: undefined surface (%d)", n); */
  /*       exit(1); */
  /*   } */

  /*   return NULL; */
}

double Grid :: giveDiameter(FloatArray &coords) {
  double diam = this->diameter;
  for (int n=1;n<=this->giveNumberOfRefinements();n++){
    diam = this->giveRefinement(n)->giveDiameter(coords);
  }

  return diam;
}


Region *Grid :: giveRegion(int n)
// Returns the n-th element. Generates error if it is not defined yet.
{
    if (n >= 1 && n <= static_cast<int>(regionList.size()) && regionList[n - 1] != nullptr) {
        return regionList[n - 1];  // 1-based to 0-based
    } else {
        printf("giveRegion: undefined region (%d)\n", n);
        exit(1);
    }
    return nullptr;
}

Inclusion *Grid :: giveInclusion(int n)
// Returns the n-th element. Generates error if it is not defined yet.
{
  if (n >= 1 && n <= static_cast<int>(inclusionList.size()) && inclusionList[n - 1] != nullptr) {
        return inclusionList[n - 1];  // 1-based to 0-based
    } else {
        printf("giveInclusion: undefined inclusion (%d)\n", n);
        exit(1);
    }
    return nullptr;
}

Refinement *Grid :: giveRefinement(int n)
// Returns the n-th element. Generates error if it is not defined yet.
{
  if (n >= 1 && n <= static_cast<int>(refinementList.size()) && refinementList[n - 1] != nullptr) {
        return refinementList[n - 1];  // 1-based to 0-based
    } else {
        printf("giveRefinement: undefined refinement (%d)\n", n);
        exit(1);
    }
    return nullptr;
}

/* void Grid :: resizeVertices(int _newSize) { vertexList->growTo(_newSize); } */
/* void Grid :: resizeInputVertices(int _newSize) { inputVertexList->growTo(_newSize); } */
/* void Grid :: resizeControlVertices(int _newSize) { controlVertexList->growTo(_newSize); } */
/* void Grid :: resizeCurves(int _newSize) { curveList->growTo(_newSize); } */
/* void Grid :: resizeSurfaces(int _newSize) { surfaceList->growTo(_newSize); } */
/* void Grid :: resizeRegions(int _newSize) { regionList->growTo(_newSize); } */
/* void Grid :: resizeInclusions(int _newSize) { inclusionList->growTo(_newSize); } */
/* void Grid :: resizeRefinements(int _newSize) { refinementList->growTo(_newSize); } */

void Grid :: resizeVertices(int _newSize) { vertexList.resize(_newSize, nullptr); }
void Grid :: resizeInputVertices(int _newSize) { inputVertexList.resize(_newSize, nullptr); }
void Grid :: resizeControlVertices(int _newSize) { controlVertexList.resize(_newSize, nullptr); }
void Grid :: resizeCurves(int _newSize) { curveList.resize(_newSize, nullptr); }
void Grid :: resizeSurfaces(int _newSize) { surfaceList.resize(_newSize, nullptr); }
void Grid :: resizeRegions(int _newSize) { regionList.resize(_newSize, nullptr); }
void Grid :: resizeInclusions(int _newSize) { inclusionList.resize(_newSize, nullptr); }
void Grid :: resizeRefinements(int _newSize) { refinementList.resize(_newSize, nullptr); }


/* void Grid :: setVertex(int i, Vertex *obj) { vertexList->put(i, obj); } */
/* void Grid :: setInputVertex(int i, Vertex *obj) { inputVertexList->put(i, obj); } */
/* void Grid :: setControlVertex(int i, Vertex *obj) { controlVertexList->put(i, obj); } */
/* void Grid :: setCurve(int i, Curve *obj) { curveList->put(i, obj); } */
/* void Grid :: setSurface(int i, Surface *obj) { surfaceList->put(i, obj); } */
/* void Grid :: setRegion(int i, Region *obj) { regionList->put(i, obj); } */
/* void Grid :: setInclusion(int i, Inclusion *obj) { inclusionList->put(i, obj); } */
/* void Grid :: setRefinement(int i, Refinement *obj) { refinementList->put(i, obj); } */


void Grid :: setVertex(int i, Vertex *obj) { vertexList.at(i - 1) = obj; }
void Grid :: setInputVertex(int i, Vertex *obj) { inputVertexList.at(i - 1) = obj; }
void Grid :: setControlVertex(int i, Vertex *obj) { controlVertexList.at(i - 1) = obj; }
void Grid :: setCurve(int i, Curve *obj) { curveList.at(i - 1) = obj; }
void Grid :: setSurface(int i, Surface *obj) { surfaceList.at(i - 1) = obj; }
void Grid :: setRegion(int i, Region *obj) { regionList.at(i - 1) = obj; }
void Grid :: setInclusion(int i, Inclusion *obj) { inclusionList.at(i - 1) = obj; }
void Grid :: setRefinement(int i, Refinement *obj) { refinementList.at(i - 1) = obj; }



int Grid :: generatePoints()
{
  //@todo: Is the regular generation the same for normal and periodic?
  if ( this->regularFlag == 1 ) { //regular
    this->generateRegularPoints();
  } else { //random (irregular or periodic)
    this->generateRandomPoints();
  }
  return 1;
}

int Grid :: generateRandomPoints()
{
  //Generation of mixture of periodic and random points.
  //Consider only four cases so far
  //1) no direction periodic, 2) x-direction periodic, 3) x- and y- direction periodic, 4) all direction periodic 
  IntArray periodicityFlag;
  this->givePeriodicityFlag(periodicityFlag);
  
  FloatArray lineNormal;
  FloatArray surfaceNormal;
  // measure time consumed by point generation
  clock_t start;
  clock_t sc;
  clock_t ec;
  long nsec;
  
  start = :: clock();

  if(periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0 && randomFlag != 2){//Region without periodicity

    if(randomFlag == 1){ //Input vertices for Grassl's mesh approach
      this->generateInputPoints();
    }


    //generate inclusions
    sc = :: clock();
    //Inclusions
    //@todo: Inclusions need to have a random and regular point generation
    for ( int i = 0; i < this->giveNumberOfInclusions(); i++ ) {
      ( this->giveInclusion(i + 1) )->generatePoints();
      printf( "numberOfVertices = %d\n", this->giveNumberOfVertices() );
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    oofem::OOFEM_LOG_INFO("Points for inclusions generated in %lds \n", nsec);
    
    
    //Control vertices
    this->generateControlPoints(); //This should now include also periodic shifts
    
    if(randomFlag == 1){ //Input vertices for Grassl's mesh approach
      //Curves
      sc = :: clock();
      for (int i = 0; i < this->giveNumberOfCurves(); i++ ) {
	(this->giveCurve(i+1))->generatePoints();
	printf("numberOfVertices = %d\n", this->giveNumberOfVertices());
      }
      ec = :: clock();
      nsec = ( ec - sc ) / CLOCKS_PER_SEC;
      OOFEM_LOG_INFO("Points on curves generated in %lds \n", nsec);
      
      
      //Surfaces
      sc = :: clock();
      for (int i = 0; i < this->giveNumberOfSurfaces(); i++ ) {
	(this->giveSurface(i+1))->generatePoints();
	printf("numberOfVertices = %d\n", this->giveNumberOfVertices());
      }
      ec = :: clock();
      nsec = ( ec - sc ) / CLOCKS_PER_SEC;
      OOFEM_LOG_INFO("Points on surfaces generated in %lds \n", nsec);
    }
    
    //Regions
    sc = :: clock();
    for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
      ( this->giveRegion(i + 1) )->generatePoints();
      printf( "numberOfVertices = %d\n", this->giveNumberOfVertices() );
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    OOFEM_LOG_INFO("Points in regions generated in %lds \n", nsec);
    
    nsec = ( ec - start ) / CLOCKS_PER_SEC;
    OOFEM_LOG_INFO("All points generated in %lds \n", nsec);
    
  }
  else if( this->randomFlag == 2 ){//Generate periodic surfaces for Adam
    
    /*Each direction should be generated only once and then shifted by the specimen dimension*/
    IntArray lineCounter(3);
    FloatArray lineNormal(3);
    lineCounter.zero();
    
    this->generateInputPoints();
    
    this->generateControlPoints(); //This should now include also periodic shifts

    sc = :: clock();
    for (int i = 0; i < this->giveNumberOfCurves(); i++ ) {
      this->giveCurve(i+1)->giveNormal(lineNormal);
      if(lineNormal.at(1) == 1 && lineNormal.at(2) == 0 && lineNormal.at(3) == 0 && lineCounter.at(1) == 0){  
	(this->giveCurve(i+1))->generatePoints();
	printf("numberOfVertices = %d\n", this->giveNumberOfVertices());
	lineCounter.at(1) = 1;
      }
      else if(lineNormal.at(1) == 0 && lineNormal.at(2) == 1 && lineNormal.at(3) == 0 && lineCounter.at(2) == 0){
	(this->giveCurve(i+1))->generatePoints();
	printf("numberOfVertices = %d\n", this->giveNumberOfVertices());
	lineCounter.at(2) = 1;
      }
      else if(lineNormal.at(1) == 0 && lineNormal.at(2) == 0 && lineNormal.at(3) == 1 && lineCounter.at(3) == 0){
	(this->giveCurve(i+1))->generatePoints();
	printf("numberOfVertices = %d\n", this->giveNumberOfVertices());
	lineCounter.at(3) = 1;	
      }
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    OOFEM_LOG_INFO("Points on curves generated in %lds \n", nsec);
    
    //Surfaces
    IntArray surfaceCounter(3);
    FloatArray surfaceNormal(3);
    surfaceCounter.zero();

    sc = :: clock();
    for (int i = 0; i < this->giveNumberOfSurfaces(); i++ ) {
      this->giveSurface(i+1)->giveNormal(surfaceNormal);      
      if(surfaceNormal.at(1) == 1 && surfaceNormal.at(2) == 0 && surfaceNormal.at(3) == 0 && surfaceCounter.at(1) == 0){
	(this->giveSurface(i+1))->generatePoints();
	printf("numberOfVertices = %d\n", this->giveNumberOfVertices());
	surfaceCounter.at(1) = 1;
      }
      else if(surfaceNormal.at(1) == 0 && surfaceNormal.at(2) == 1 && surfaceNormal.at(3) == 0 && surfaceCounter.at(2) == 0){
	(this->giveSurface(i+1))->generatePoints();
	printf("numberOfVertices = %d\n", this->giveNumberOfVertices());
	surfaceCounter.at(2) = 1;       
      }
      else if(surfaceNormal.at(1) == 0 && surfaceNormal.at(2) == 0 && surfaceNormal.at(3) == 1 && surfaceCounter.at(3) == 0){
	(this->giveSurface(i+1))->generatePoints();
	printf("numberOfVertices = %d\n", this->giveNumberOfVertices());
	surfaceCounter.at(3) = 1;
      }      
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    OOFEM_LOG_INFO("Points on surfaces generated in %lds \n", nsec);

    //Regions
    sc = :: clock();
    for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
      ( this->giveRegion(i + 1) )->generatePoints();
      printf( "numberOfVertices = %d\n", this->giveNumberOfVertices() );
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    OOFEM_LOG_INFO("Points in regions generated in %lds \n", nsec);
    
    nsec = ( ec - start ) / CLOCKS_PER_SEC;
    OOFEM_LOG_INFO("All points generated in %lds \n", nsec);
        
  }  
  else if( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0){
    //This is the case of a region with periodicity in x-direction, as is the case for a periodic beam
    sc = :: clock();

    //Curves
    for (int i = 0; i < this->giveNumberOfCurves(); i++ ) {
      this->giveCurve(i+1)->giveNormal(lineNormal);
      if( lineNormal.at(1) == 1 && lineNormal.at(2) == 0 && lineNormal.at(3) == 0){// line in x-direction
	  (this->giveCurve(i+1))->generatePoints();
	  printf("numberOfVertices = %d\n", this->giveNumberOfVertices());
      }           
      ec = :: clock();
      nsec = ( ec - sc ) / CLOCKS_PER_SEC;
      OOFEM_LOG_INFO("Points on curves generated in %lds \n", nsec);
    }

    //Surfaces
    sc = :: clock();
    for (int i = 0; i < this->giveNumberOfSurfaces(); i++ ) {
       this->giveSurface(i+1)->giveNormal(surfaceNormal);
      if( !(surfaceNormal.at(1) == 1 && surfaceNormal.at(2) == 0 && surfaceNormal.at(3) == 0)){// all surfaces without normal in x-direction
	(this->giveSurface(i+1))->generatePoints();
	printf("numberOfVertices = %d\n", this->giveNumberOfVertices());
      }
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    OOFEM_LOG_INFO("Points on surfaces generated in %lds \n", nsec);
    
    //Regions
    sc = :: clock();
    for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
      ( this->giveRegion(i + 1) )->generatePoints();
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    OOFEM_LOG_INFO("Points in regions generated in %lds \n", nsec);
  }
  else if(periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0){
    //Region with periodicity in x and y directions. This could be the case of a slab.
    
    this->generateControlPoints(); //This should now include also periodic shifts
    
    //Surfaces
    sc = :: clock();
    for (int i = 0; i < this->giveNumberOfSurfaces(); i++ ) {
      this->giveSurface(i+1)->giveNormal(surfaceNormal);
      if( surfaceNormal.at(1) == 0 && surfaceNormal.at(2) == 0 && surfaceNormal.at(3) == 1){// all surfaces without normal in x-direction
	(this->giveSurface(i+1))->generatePoints();
	printf("numberOfVertices = %d\n", this->giveNumberOfVertices());
      }
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    OOFEM_LOG_INFO("Points on surfaces generated in %lds \n", nsec);
        
    //Regions
    sc = :: clock();
    for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
      ( this->giveRegion(i + 1) )->generatePoints();
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    OOFEM_LOG_INFO("Points in regions generated in %lds \n", nsec);
    
  }
  else if(periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1){
    //Region with periodicity in all directions.

    sc = :: clock();
    this->generateControlPoints();
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    OOFEM_LOG_INFO("Control points generated in %lds \n", nsec);
        
    sc = :: clock();
    //Inclusions
    //@todo: Inclusions need to have a random and regular point generation
    for ( int i = 0; i < this->giveNumberOfInclusions(); i++ ) {
        ( this->giveInclusion(i + 1) )->generatePeriodicPoints();
        printf( "numberOfVertices = %d\n", this->giveNumberOfVertices() );
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    OOFEM_LOG_INFO("Points for inclusions generated in %lds \n", nsec);

    //Regions
    sc = :: clock();
    for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
        ( this->giveRegion(i + 1) )->generatePoints();
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    OOFEM_LOG_INFO("Points in regions generated in %lds \n", nsec);    
  }
  else{//should not be here
    printf("Something wrong with mixed point generation and periodic flag in grid->generateRandomPoints\n");
  }

  return 1;
  
}

void Grid :: generateControlPoints()
{
    Vertex *vertex;

    int vertexNumber = this->vertexList.size();
    
    FloatArray boundaries;
    defineBoundaries(boundaries);

    FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);
    
    //Generate control vertices
    FloatArray coords(3), mirroredCoords(3);
    
    int controlVertexNumber = this->controlVertexList.size();
    for ( int i = 0; i < controlVertexNumber; i++ ) {
      this->giveControlVertex(i+1)->giveCoordinates(coords);
      vertex = ( Vertex * ) ( Vertex(vertexNumber + 1, this).ofType() );
      vertex->setCoordinates(coords);
      this->setVertex(vertexNumber + 1, vertex);
      this->giveGridLocalizer()->insertSequentialNode(vertexNumber + 1, coords);
      vertexNumber++;


      if(this->randomFlag == 0){ //Here it is assumed that for Grassl's meshing approach control nodes are on the surface, and for Bolander's approach they are not.
	//Do now the mirroring.
	//It works only for rectangular specimen in a cartesian coordinate system.
	//x-direction


      for ( int k = 1; k <= 2; k++ ) {
	  mirroredCoords.at(1) = coords.at(1) - 2. * ( coords.at(1) - boundaries.at(k) );
	  mirroredCoords.at(2) = coords.at(2);
	  mirroredCoords.at(3) = coords.at(3);
	  vertex = ( Vertex * ) ( Vertex(vertexNumber + 1, this).ofType() );
	  vertex->setCoordinates(mirroredCoords);
	  this->setVertex(vertexNumber + 1, vertex);
	  this->giveGridLocalizer()->insertSequentialNode(vertexNumber + 1, mirroredCoords);
	  vertexNumber++;
	}
	
	//y-direction
	for ( int k = 3; k <= 4; k++ ) {
	  mirroredCoords.at(1) = coords.at(1);
	  mirroredCoords.at(2) = coords.at(2) - 2. * ( coords.at(2) - boundaries.at(k) );
	  mirroredCoords.at(3) = coords.at(3);
	  vertex = ( Vertex * ) ( Vertex(vertexNumber + 1, this).ofType() );
	  vertex->setCoordinates(mirroredCoords);
	  this->setVertex(vertexNumber + 1, vertex);
	  this->giveGridLocalizer()->insertSequentialNode(vertexNumber + 1, mirroredCoords);
	  vertexNumber++;
	}
	
	//z-direction
	for ( int k = 5; k <= 6; k++ ) {
	  mirroredCoords.at(1) = coords.at(1);
	  mirroredCoords.at(2) = coords.at(2);
	  mirroredCoords.at(3) = coords.at(3) - 2. * ( coords.at(3) - boundaries.at(k) );
	  vertex = ( Vertex * ) ( Vertex(vertexNumber + 1, this).ofType() );
	  vertex->setCoordinates(mirroredCoords);
	  this->setVertex(vertexNumber + 1, vertex);
	  this->giveGridLocalizer()->insertSequentialNode(vertexNumber + 1, mirroredCoords);
	  vertexNumber++;
	}
      }
    }
      
  return;
}




void Grid :: generateInputPoints()
{
    Vertex *vertex;

    int vertexNumber = this->vertexList.size();
    
    //Generate control vertices
    FloatArray coords(3);
    
    int inputVertexNumber = this->inputVertexList.size();
    for ( int i = 0; i < inputVertexNumber; i++ ) {
      this->giveInputVertex(i+1)->giveCoordinates(coords);
      vertex = ( Vertex * ) ( Vertex(vertexNumber + 1, this).ofType() );
      vertex->setCoordinates(coords);
      this->setVertex(vertexNumber + 1, vertex);
      this->giveGridLocalizer()->insertSequentialNode(vertexNumber + 1, coords);
      vertexNumber++;
    }
    
  return;
}



int Grid :: generateRegularPoints()
{
    //Regular approach for both normal and periodic
    clock_t start;
    clock_t sc;
    clock_t ec;
    long nsec;

    start = :: clock();
    //The old method is reviced, the only generation needed is the region one.
    //The user has to set in the input a larger region, than the one in the converter, in order to avoid boundary discrepancies.
    //Regions
    for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
        if ( regType == 1 ) { //BCC
            ( this->giveRegion(i + 1) )->generateRegularPoints1();
        } else if ( regType == 2 )     { //FCC
            ( this->giveRegion(i + 1) )->generateRegularPoints2();
        } else   {
            printf("Error: The regular grid type was not determined\n");
            exit(0);
        }
    }
    sc = :: clock();
    return 1;
}


int Grid :: instanciateYourself(DataReader *dr)
// Creates all objects mentioned in the data file.
{
  //    const char *__keyword, *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //    IRResultType result;                            // Required by IR_GIVE_FIELD macro

    int i, num;
    std::string name;

    Vertex *vertex;
    Curve *curve;
    Surface *surface;
    Prism *prism;
    Sphere *sphere;
    Cylinder *cylinder;
    
    InterfaceSphere *intersphere;
    InterfaceCylinder *intercylinder;

    RefinePrism *refineprism;
    
    //    InputRecord *ir;

    /* read type of Grid to be solved
     * This information is currently not used, since we do not distinguish between different domain types.
     * However, later we should have different domain types which are inherited from domain.
     */
    //    ir = dr->giveInputRecord(DataReader :: IR_domainRec, 1);
    auto &ir = dr->giveInputRecord(oofem::DataReader::IR_domainRec, 1);
    IR_GIVE_FIELD(ir, name, IFT_Domain_type);
    ir.finish();

    // read
    ir = dr->giveInputRecord(DataReader :: IR_controlRec, 1);
    IR_GIVE_FIELD(ir, diameter, IFT_diam, "diam"); // Macro
    TOL = 1.e-6 * diameter;

    //    IR_GIVE_FIELD(ir, density, IFT_dens, "dens"); // Macro

    maxIter = 1.e6;
    IR_GIVE_OPTIONAL_FIELD(ir, maxIter, IFT_maxiter, "maxiter"); // Macro

    //This will need to be removed later again, once the aggregate generator works
    aggregateFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, aggregateFlag, IFT_aggflag, "aggflag"); // Macro

    if ( aggregateFlag == 1 ) {
        IR_GIVE_OPTIONAL_FIELD(ir, targetAggregates, IFT_aggflag, "target"); // Macro
    }

    // We need to write the generator so that it can be used for both periodic and normal meshes.
    // Consequently, we should introduce two flags. One regular flag and one periodicity flag.
    // The standard case should be an irregular generator without periodicity, which would correspond to
    // regularFlag = periodicityFlag = 0

    //when regularFlag=1 regular generator is used.
    this->regularFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->regularFlag, IFT_regflag, "regflag"); // Macro
    if ( regularFlag < 0 || regularFlag > 1 ) {
        printf("Error: Unknown regular flag. Should be 0 or 1.\n");
        exit(0);
    }

    //Only needed for the regular generator, but then it should be mandatory
    if ( this->regularFlag == 1 ) {
        xyzEdges.resize(0);
        IR_GIVE_FIELD(ir, xyzEdges, IFT_xyzedges, "xyzedges"); // Macro

        regType = 0;
        IR_GIVE_FIELD(ir, regType, IFT_regtype, "regtype"); // Macro
    }

    
    //when periodicityFlag=1 periodic cell generator is used.
    periodicityFlag.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, periodicityFlag, IFT_perflag, "perflag"); // Macro
    if ( periodicityFlag.size() != 3 ) {
        printf("Error: Unknown periodicity flag. Should have three components.\n");
        exit(0);
    }

    //two types of random generator
    //0 - Bolander's way of generating meshes. No nodes on boundary. This approach creates problems for coupled analyses, but works well for structural simulations.
    //1 - Grassl's way of generating meshes. Nodes on boundaries. Better for coupled analyses.
    //2 - Grassl's way of generating meshes, but periodic nodes on boundaries.
    randomFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, randomFlag, IFT_ranflag, "ranflag"); // Macro
    if ( randomFlag < 0 || randomFlag > 2 ) {
        printf("Error: Unknown random flag. Should be 0, 1 or 2.\n");
        exit(0);
    }
    

    //This is useful to be able to regenerate the same mesh.
    //Simply set it to a negative value and this is used to generate the random numbers.
    //If it is not set, then time is used to generate the interger, which is used to generate the random numbers.
    randomInteger = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, randomInteger, IFT_ranint, "ranint"); // Macro
    if ( randomInteger >= 0 ) {
        randomInteger = -time(NULL);
        printf("random integer is determined as %d\n", randomInteger);
    }

    ir->finish();

    // read domain description
    int nvertex, ncontrolvertex=0, nelem, ncurve, nsurface, nregion, ninclusion,nrefinement;
    ir = dr->giveInputRecord(DataReader :: IR_domainCompRec, 1);
    IR_GIVE_FIELD(ir, nvertex, IFT_nvertex, "nvertex"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, ncontrolvertex, IFT_ncontrolvertex, "ncontrolvertex"); // Macro
    IR_GIVE_FIELD(ir, ncurve, IFT_ncurve, "ncurve"); // Macro
    IR_GIVE_FIELD(ir, nsurface, IFT_nsurface, "nsurface"); // Macro
    IR_GIVE_FIELD(ir, nregion, IFT_nregion, "nregion"); // Macro


    IR_GIVE_FIELD(ir, ninclusion, IFT_ninclusion, "ninclusion"); // Macro
    IR_GIVE_FIELD(ir, nrefinement, IFT_nrefinement, "nrefinement"); // Macro
    ir->finish();

    //read vertices
    //    inputVertexList->growTo(nvertex);
    inputVertexList.resize(nvertex, nullptr);
    
    for ( i = 0; i < nvertex; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_vertexRec, i + 1);
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num, MAX_NAME_LENGTH);


        ( vertex = ( Vertex * ) ( Vertex(num, this).ofType() ) )->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > nvertex ) ) {
            printf("instanciateYourself: Invalid vertex number (num=%d)", num);
            exit(0);
        }



        if ( !inputVertexList->includes(num) ) {
            setInputVertex(num, vertex);
        } else {
            printf("instanciateYourself: Vertex entry already exist (num=%d)", num);
            exit(0);
        }

        ir->finish();
    }


    //read control vertices
    controlVertexList->growTo(ncontrolvertex);
    controlVertexList.resize(ncontrolvertex, nullptr);

    for ( i = 0; i < ncontrolvertex; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_vertexRec, i + 1);
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num, MAX_NAME_LENGTH);


        ( vertex = ( Vertex * ) ( Vertex(num, this).ofType() ) )->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > ncontrolvertex ) ) {
            printf("instanciateYourself: Invalid vertex number (num=%d)", num);
            exit(0);
        }



        if ( !controlVertexList->includes(num) ) {
            setControlVertex(num, vertex);
        } else {
            printf("instanciateYourself: Vertex entry already exist (num=%d)", num);
            exit(0);
        }

        ir->finish();
    }
    
    // read curves
    //    curveList->growTo(ncurve);
    curveList.resize(ncurve, nullptr);    
    for ( i = 0; i < ncurve; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_curveRec, i + 1);
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num, MAX_NAME_LENGTH);

        ( curve = ( Curve * ) ( Curve(num, this).ofType() ) )->initializeFrom(ir);

        if ( ( num < 1 ) || ( num > ncurve ) ) {
            printf("instanciateYourself: Invalid curve number (num=%d)", num);
            exit(0);
        }

        if ( !curveList->includes(num) ) {
            setCurve(num, curve);
        } else {
            printf("instanciateYourself: Curve entry already exist (num=%d)", num);
            exit(0);
        }
        ir->finish();
    }


    // read surfaces
    //    surfaceList->growTo(nsurface);
    surfaceList.resize(nsurface, nullptr);
    for ( i = 0; i < nsurface; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_surfaceRec, i + 1);
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num, MAX_NAME_LENGTH);

        ( surface = ( Surface * ) ( Surface(num, this).ofType() ) )->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > nsurface ) ) {
            printf("instanciateYourself: Invalid curve number (num=%d)", num);
            exit(0);
        }

        if ( !surfaceList->includes(num) ) {
            setSurface(num, surface);
        } else {
            printf("instanciateYourself: Surface entry already exist (num=%d)", num);
            exit(0);
        }
        ir->finish();
    }

    // read regions
    //    regionList->growTo(nregion);
    regionList.resize(nregion, nullptr);
    for ( i = 0; i < nregion; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_regionRec, i + 1);
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num, MAX_NAME_LENGTH);

        if ( !strcmp(name, "prism") ) { //Prism
            ( prism = ( Prism * ) ( Prism(num, this).ofType() ) )->initializeFrom(ir);

            if ( ( num < 1 ) || ( num > nregion ) ) {
                printf("instanciateYourself: Invalid region number (num=%d)", num);
                exit(0);
            }

            if ( !regionList->includes(num) ) {
                setRegion(num, prism);
            } else {
                printf("instanciateYourself: Region entry already exist (num=%d)", num);
                exit(0);
            }
        }
	else if ( !strcmp(name, "cylinder") ) { //Cylinder
            ( cylinder = ( Cylinder * ) ( Cylinder(num, this).ofType() ) )->initializeFrom(ir);

            if ( ( num < 1 ) || ( num > nregion ) ) {
                printf("instanciateYourself: Invalid region number (num=%d)", num);
                exit(0);
            }

            if ( !regionList->includes(num) ) {
                setRegion(num, cylinder);
            } else {
                printf("instanciateYourself: Region entry already exist (num=%d)", num);
                exit(0);
            }
        }
	else if ( !strcmp(name, "sphere") ) { //Sphere
            ( sphere = ( Sphere * ) ( Sphere(num, this).ofType() ) )->initializeFrom(ir);

            if ( ( num < 1 ) || ( num > nregion ) ) {
                printf("instanciateYourself: Invalid region number (num=%d)", num);
                exit(0);
            }

            if ( !regionList->includes(num) ) {
                setRegion(num, sphere);
            } else {
                printf("instanciateYourself: Regio entry already exist (num=%d)", num);
                exit(0);
            }
        }
	else {
            printf("Error: Unknown region type\n");
            exit(0);
	}
	
        ir->finish();
    }


    // read inclusions
    //    inclusionList->growTo(ninclusion);
    inclusionList.resize(ninclusion, nullptr);
    for ( i = 0; i < ninclusion; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_inclusionRec, i + 1);
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num, MAX_NAME_LENGTH);

        if ( !strcmp(name, "intersphere") ) { //Sphere
            ( intersphere = ( InterfaceSphere * ) ( InterfaceSphere(num, this).ofType() ) )->initializeFrom(ir);

            if ( ( num < 1 ) || ( num > ninclusion ) ) {
                printf("instanciateYourself: Invalid inclusion number (num=%d)", num);
                exit(0);
            }

            if ( !inclusionList->includes(num) ) {
                setInclusion(num, intersphere);
            } else {
                printf("instanciateYourself: Inclusion entry already exist (num=%d)", num);
                exit(0);
            }
        }
	else if ( !strcmp(name, "interfacecylinder") ) {       //Interface Cylinder
            ( intercylinder = ( InterfaceCylinder * ) ( InterfaceCylinder(num, this).ofType() ) )->initializeFrom(ir);

            if ( ( num < 1 ) || ( num > ninclusion ) ) {
                printf("instanciateYourself: Invalid inclusion number (num=%d)", num);
                exit(0);
            }

            if ( !inclusionList->includes(num) ) {
                setInclusion(num, intercylinder);
            } else {
                printf("instanciateYourself: Inclusion entry already exist (num=%d)", num);
                exit(0);
            }
        } else {
            printf("Error: Unknown inclusion type\n");
            exit(0);
        }


        ir->finish();
    }
    
    // read refinements
    //    refinementList->growTo(nrefinement);
    refinementList.resize(nrefinement, nullptr);
    for ( i = 0; i < nrefinement; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_refinementRec, i + 1);
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num, MAX_NAME_LENGTH);

        if ( !strcmp(name, "refineprism") ) { //Prism
            ( refineprism = ( RefinePrism * ) ( RefinePrism(num, this).ofType() ) )->initializeFrom(ir);

            if ( ( num < 1 ) || ( num > nrefinement ) ) {
                printf("instanciateYourself: Invalid refinement number (num=%d)", num);
                exit(0);
            }

            if ( !refinementList->includes(num) ) {
                setRefinement(num, refineprism);
            } else {
                printf("instanciateYourself: Refinement entry already exist (num=%d)", num);
                exit(0);
            }
        } else {
            printf("Error: Unknown refinement type\n");
            exit(0);
        }


        ir->finish();
    }

  
    this->giveGridLocalizer()->init(true);

    return 1;
}




void Grid :: defineBoundaries(FloatArray &boundaries)
//Determine the boundaries of the domain
{  
    Region *region;
    Surface *surface;
  
    int localSurface, localCurve, localVertex;
    if(this->giveNumberOfRegions() >0){
      if(this->giveNumberOfRegions() != 1){
	printf("Error. Cannot defined boundaries for multiple regions yet!\n");
	exit(1);
      }
      this->giveRegion(1)->defineBoundaries(boundaries);
    }
    else if(this->giveNumberOfSurfaces() >0){
      if(this->giveNumberOfSurfaces() !=1){
	printf("Error. Cannot defined boundaries for multiple surfaces yet!\n");
	exit(1);
      }
      this->giveSurface(1)->defineBoundaries(boundaries);
      return;
    }
    return;
}



void Grid :: giveOutput(FILE *outputStream)
{
    int dimension = 3;
    fprintf(outputStream, "%d\n", dimension);

    fprintf( outputStream, "%d\n", this->giveNumberOfVertices() );

    for ( int i = 0; i < this->giveNumberOfVertices(); i++ ) {
        fprintf( outputStream, "%.16e %.16e %.16e\n", ( this->giveVertex(i + 1) )->giveCoordinate(1), ( this->giveVertex(i + 1) )->giveCoordinate(2), ( this->giveVertex(i + 1) )->giveCoordinate(3) );
    }


    //temporary measure to be able to generate simple aggregate placement.
    //Needs to be moved to aggregate generator later
    if ( aggregateFlag == 1 ) {

      

      
      FloatArray boundaries;
        defineBoundaries(boundaries);
        FILE *aggregateStream;
        if ( ( aggregateStream = fopen("aggregate.dat", "w") ) == NULL ) {
            printf("Can't open output file aggregate.dat");
            exit(1);
        }

        int outsideFlag;
        FloatArray coords;
        int counter = 0;
        for ( int i = 0; i < this->giveNumberOfVertices(); i++ ) {
            if ( ( ( i ) / 27. ) >= targetAggregates ) {
                break;
            }
            outsideFlag = 0;
            this->giveVertex(i + 1)->giveCoordinates(coords);
            for ( int k = 0; k < 3; k++ ) {
                if ( ( coords.at(k + 1) + TOL + diameter / 2. ) < boundaries.at(2 * k + 1) ) {
                    outsideFlag = 1;
                    break;
                } else if ( ( coords.at(k + 1) - TOL - diameter / 2. ) > boundaries.at(2 * k + 2) ) {
                    outsideFlag = 1;
                    break;
                }
            }
            if ( outsideFlag != 1 ) {
                counter++;
                fprintf(aggregateStream, "sphere %d centre 3 %.8e %.8e %.8e diameter %.4e\n", counter, coords.at(1), coords.at(2), coords.at(3), this->diameter);
            }
        }
    }


    return;
}


double Grid :: ran1(int *idum)
{
    // random number generator from "Numerical recipes book"
    int j;
    long k;
    static long iy = 0;
    static long iv [ NTAB ];
    float temp;

    if ( * idum <= 0 || !iy ) {
        if ( -( * idum ) < 1 ) {
            * idum = 1;
        } else {
            * idum = -( * idum );
        }
        for ( j = NTAB + 7; j >= 0; j-- ) {
            k = ( * idum ) / IQ;
            * idum = IA * ( * idum - k * IQ ) - IR * k;
            if ( * idum < 0 ) {
                * idum += IM;
            }
            if ( j < NTAB ) {
                iv [ j ] = * idum;
            }
        }
        iy = iv [ 0 ];
    }
    k = ( * idum ) / IQ;
    * idum = IA * ( * idum - k * IQ ) - IR * k;
    if ( * idum < 0 ) {
        * idum += IM;
    }
    j = iy / NDIV;
    iy = iv [ j ];
    iv [ j ] = * idum;
    if ( ( temp = AM * iy ) > RNMX ) {
        return RNMX;
    } else {
        return temp;
    }
}



//#endif
