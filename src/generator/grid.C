#include "generatorlistutils.h"
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
#include "generatordatareader.h"
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
  if (n >= 1 && n <= static_cast<int>(generator::size1(vertexList)) && vertexList[n - 1] != nullptr) {
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

  if (n >= 1 && n <= static_cast<int>(generator::size1(inputVertexList)) && inputVertexList[n - 1] != nullptr) {
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

  if (n >= 1 && n <= static_cast<int>(generator::size1(controlVertexList)) && controlVertexList[n - 1] != nullptr) {
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

  if (n >= 1 && n <= static_cast<int>(generator::size1(curveList)) && curveList[n - 1] != nullptr) {
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

  if (n >= 1 && n <= static_cast<int>(generator::size1(surfaceList)) && surfaceList[n - 1] != nullptr) {
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

double Grid :: giveDiameter(oofem::FloatArray &coords) {
  double diam = this->diameter;
  for (int n=1;n<=this->giveNumberOfRefinements();n++){
    diam = this->giveRefinement(n)->giveDiameter(coords);
  }

  return diam;
}


Region *Grid :: giveRegion(int n)
// Returns the n-th element. Generates error if it is not defined yet.
{
  if (n >= 1 && n <= static_cast<int>(generator::size1(regionList)) && regionList[n - 1] != nullptr) {
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
  if (n >= 1 && n <= static_cast<int>(generator::size1(inclusionList)) && inclusionList[n - 1] != nullptr) {
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
  if (n >= 1 && n <= static_cast<int>(generator::size1(refinementList)) && refinementList[n - 1] != nullptr) {
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


void Grid :: setVertex(int i, Vertex *obj) {
  generator::put1_replace(vertexList, i, obj);
}

void Grid :: setInputVertex(int i, Vertex *obj) {
  generator::put1_replace(inputVertexList, i, obj);
}

void Grid :: setControlVertex(int i, Vertex *obj) {
  generator::put1_replace(controlVertexList, i, obj);
}

void Grid::setCurve(int i, Curve* obj){
  generator::put1_replace(curveList, i, obj);
}

void Grid :: setSurface(int i, Surface *obj) {
  generator::put1_replace(surfaceList, i, obj);
}

void Grid :: setRegion(int i, Region *obj) {
  generator::put1_replace(regionList, i, obj);
}

void Grid :: setInclusion(int i, Inclusion *obj) {
  generator::put1_replace(inclusionList, i, obj);
}

void Grid :: setRefinement(int i, Refinement *obj) {
  generator::put1_replace(refinementList, i, obj);
}

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
  oofem::IntArray periodicityFlag;
  this->givePeriodicityFlag(periodicityFlag);
  
  oofem::FloatArray lineNormal;
  oofem::FloatArray surfaceNormal;
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
    printf("Points for inclusions generated in %lds \n", nsec);
    
    
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
      printf("Points on curves generated in %lds \n", nsec);
      
      
      //Surfaces
      sc = :: clock();
      for (int i = 0; i < this->giveNumberOfSurfaces(); i++ ) {
	(this->giveSurface(i+1))->generatePoints();
	printf("numberOfVertices = %d\n", this->giveNumberOfVertices());
      }
      ec = :: clock();
      nsec = ( ec - sc ) / CLOCKS_PER_SEC;
      printf("Points on surfaces generated in %lds \n", nsec);
    }
    
    //Regions
    sc = :: clock();
    for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
      ( this->giveRegion(i + 1) )->generatePoints();
      printf( "numberOfVertices = %d\n", this->giveNumberOfVertices() );
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    printf("Points in regions generated in %lds \n", nsec);
    
    nsec = ( ec - start ) / CLOCKS_PER_SEC;
    printf("All points generated in %lds \n", nsec);
    
  }
  else if( this->randomFlag == 2 ){//Generate periodic surfaces for Adam
    
    /*Each direction should be generated only once and then shifted by the specimen dimension*/
    oofem::IntArray lineCounter(3);
    oofem::FloatArray lineNormal(3);
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
    printf("Points on curves generated in %lds \n", nsec);
    
    //Surfaces
    oofem::IntArray surfaceCounter(3);
    oofem::FloatArray surfaceNormal(3);
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
    printf("Points on surfaces generated in %lds \n", nsec);

    //Regions
    sc = :: clock();
    for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
      ( this->giveRegion(i + 1) )->generatePoints();
      printf( "numberOfVertices = %d\n", this->giveNumberOfVertices() );
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    printf("Points in regions generated in %lds \n", nsec);
    
    nsec = ( ec - start ) / CLOCKS_PER_SEC;
    printf("All points generated in %lds \n", nsec);
        
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
      printf("Points on curves generated in %lds \n", nsec);
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
    printf("Points on surfaces generated in %lds \n", nsec);
    
    //Regions
    sc = :: clock();
    for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
      ( this->giveRegion(i + 1) )->generatePoints();
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    printf("Points in regions generated in %lds \n", nsec);
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
    printf("Points on surfaces generated in %lds \n", nsec);
        
    //Regions
    sc = :: clock();
    for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
      ( this->giveRegion(i + 1) )->generatePoints();
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    printf("Points in regions generated in %lds \n", nsec);
    
  }
  else if(periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1){
    //Region with periodicity in all directions.

    sc = :: clock();
    this->generateControlPoints();
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    printf("Control points generated in %lds \n", nsec);
        
    sc = :: clock();
    //Inclusions
    //@todo: Inclusions need to have a random and regular point generation
    for ( int i = 0; i < this->giveNumberOfInclusions(); i++ ) {
        ( this->giveInclusion(i + 1) )->generatePeriodicPoints();
        printf( "numberOfVertices = %d\n", this->giveNumberOfVertices() );
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    printf("Points for inclusions generated in %lds \n", nsec);

    //Regions
    sc = :: clock();
    for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
        ( this->giveRegion(i + 1) )->generatePoints();
    }
    ec = :: clock();
    nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    printf("Points in regions generated in %lds \n", nsec);    
  }
  else{//should not be here
    printf("Something wrong with mixed point generation and periodic flag in grid->generateRandomPoints\n");
  }

  return 1;
  
}

void Grid :: generateControlPoints()
{
    Vertex *vertex;

    int vertexNumber = generator::size1(this->vertexList);
    
    oofem::FloatArray boundaries;
    defineBoundaries(boundaries);

    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);
    
    //Generate control vertices
    oofem::FloatArray coords(3), mirroredCoords(3);
    
    int controlVertexNumber = generator::size1(this->controlVertexList);
    for ( int i = 0; i < controlVertexNumber; i++ ) {
      this->giveControlVertex(i+1)->giveCoordinates(coords);



      auto *v = new Vertex(vertexNumber+1, this);
      v->setCoordinates(coords);
      this->setVertex(vertexNumber+1, v);
      this->giveGridLocalizer()->insertSequentialNode(vertexNumber+1, coords);
      vertexNumber++;



      if(this->randomFlag == 0){ //Here it is assumed that for Grassl's meshing approach control nodes are on the surface, and for Bolander's approach they are not.
	//Do now the mirroring.
	//It works only for rectangular specimen in a cartesian coordinate system.
	//x-direction


      for ( int k = 1; k <= 2; k++ ) {
	  mirroredCoords.at(1) = coords.at(1) - 2. * ( coords.at(1) - boundaries.at(k) );
	  mirroredCoords.at(2) = coords.at(2);
	  mirroredCoords.at(3) = coords.at(3);


	  auto *v = new Vertex(vertexNumber+1, this);
	  v->setCoordinates(mirroredCoords);
	  this->setVertex(vertexNumber+1, v);
	  this->giveGridLocalizer()->insertSequentialNode(vertexNumber+1, mirroredCoords);
	  vertexNumber++;
	  
	}
	
	//y-direction
	for ( int k = 3; k <= 4; k++ ) {
	  mirroredCoords.at(1) = coords.at(1);
	  mirroredCoords.at(2) = coords.at(2) - 2. * ( coords.at(2) - boundaries.at(k) );
	  mirroredCoords.at(3) = coords.at(3);

	  auto *v = new Vertex(vertexNumber+1, this);
	  v->setCoordinates(mirroredCoords);
	  this->setVertex(vertexNumber+1, v);
	  this->giveGridLocalizer()->insertSequentialNode(vertexNumber+1, mirroredCoords);
	  vertexNumber++;
	  
	}
	
	//z-direction
	for ( int k = 5; k <= 6; k++ ) {
	  mirroredCoords.at(1) = coords.at(1);
	  mirroredCoords.at(2) = coords.at(2);
	  mirroredCoords.at(3) = coords.at(3) - 2. * ( coords.at(3) - boundaries.at(k) );

	  auto *v = new Vertex(vertexNumber+1, this);
	  v->setCoordinates(mirroredCoords);
	  this->setVertex(vertexNumber+1, v);
	  this->giveGridLocalizer()->insertSequentialNode(vertexNumber+1, mirroredCoords);
	  vertexNumber++;
	}
      }
    }
      
  return;
}




void Grid :: generateInputPoints()
{
    Vertex *vertex;

    int vertexNumber = generator::size1(this->vertexList);
    
    //Generate control vertices
    oofem::FloatArray coords(3);
    
    int inputVertexNumber = generator::size1(this->inputVertexList);
    for ( int i = 0; i < inputVertexNumber; i++ ) {
      this->giveInputVertex(i+1)->giveCoordinates(coords);

      auto *v = new Vertex(vertexNumber+1, this);
      v->setCoordinates(coords);
      this->setVertex(vertexNumber+1, v);
      this->giveGridLocalizer()->insertSequentialNode(vertexNumber+1, coords);
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



int Grid :: instanciateYourself(GeneratorDataReader *dr)
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
    auto &ir = dr->giveInputRecord(GeneratorDataReader::GIR_domainRec, 1); 

    IR_GIVE_FIELD(ir, name, _IFT_Grid_type);
    ir.finish();

    // read
    ir = dr->giveInputRecord(GeneratorDataReader :: GIR_controlRec, 1);
    IR_GIVE_FIELD(ir, diameter, _IFT_Grid_diam); // Macro
    TOL = 1.e-6 * diameter;

    //    IR_GIVE_FIELD(ir, density, IFT_dens, "dens"); // Macro

    maxIter = 1.e6;
    IR_GIVE_OPTIONAL_FIELD(ir, maxIter, _IFT_Grid_maxiter); // Macro

    //This will need to be removed later again, once the aggregate generator works
    aggregateFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, aggregateFlag, _IFT_Grid_aggflag); // Macro

    if ( aggregateFlag == 1 ) {
        IR_GIVE_OPTIONAL_FIELD(ir, targetAggregates, _IFT_Grid_aggflag); // Macro
    }

    // We need to write the generator so that it can be used for both periodic and normal meshes.
    // Consequently, we should introduce two flags. One regular flag and one periodicity flag.
    // The standard case should be an irregular generator without periodicity, which would correspond to
    // regularFlag = periodicityFlag = 0

    //when regularFlag=1 regular generator is used.
    this->regularFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->regularFlag, _IFT_Grid_regflag); // Macro
    if ( regularFlag < 0 || regularFlag > 1 ) {
        printf("Error: Unknown regular flag. Should be 0 or 1.\n");
        exit(0);
    }

    //Only needed for the regular generator, but then it should be mandatory
    if ( this->regularFlag == 1 ) {
        xyzEdges.resize(0);
        IR_GIVE_FIELD(ir, xyzEdges, _IFT_Grid_xyzedges); // Macro

        regType = 0;
        IR_GIVE_FIELD(ir, regType, _IFT_Grid_regtype); // Macro
    }

    
    //when periodicityFlag=1 periodic cell generator is used.
    periodicityFlag.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, periodicityFlag, _IFT_Grid_perflag); // Macro
    if ( periodicityFlag.giveSize() != 3 ) {
        printf("Error: Unknown periodicity flag. Should have three components.\n");
        exit(0);
    }

    //two types of random generator
    //0 - Bolander's way of generating meshes. No nodes on boundary. This approach creates problems for coupled analyses, but works well for structural simulations.
    //1 - Grassl's way of generating meshes. Nodes on boundaries. Better for coupled analyses.
    //2 - Grassl's way of generating meshes, but periodic nodes on boundaries.
    randomFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, randomFlag, _IFT_Grid_ranflag); // Macro
    if ( randomFlag < 0 || randomFlag > 2 ) {
        printf("Error: Unknown random flag. Should be 0, 1 or 2.\n");
        exit(0);
    }
    

    //This is useful to be able to regenerate the same mesh.
    //Simply set it to a negative value and this is used to generate the random numbers.
    //If it is not set, then time is used to generate the interger, which is used to generate the random numbers.
    randomInteger = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, randomInteger, _IFT_Grid_ranint); // Macro
    if ( randomInteger >= 0 ) {
        randomInteger = -time(NULL);
        printf("random integer is determined as %d\n", randomInteger);
    }

    ir.finish();

    // read domain description
    int nvertex, ncontrolvertex=0, nelem, ncurve, nsurface, nregion, ninclusion,nrefinement;
    ir = dr->giveInputRecord(GeneratorDataReader :: GIR_domainCompRec, 1);
    IR_GIVE_FIELD(ir, nvertex, _IFT_Grid_nvertex); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, ncontrolvertex, _IFT_Grid_ncontrolvertex); // Macro
    IR_GIVE_FIELD(ir, ncurve, _IFT_Grid_ncurve); // Macro
    IR_GIVE_FIELD(ir, nsurface, _IFT_Grid_nsurface); // Macro
    IR_GIVE_FIELD(ir, nregion, _IFT_Grid_nregion); // Macro


    IR_GIVE_FIELD(ir, ninclusion, _IFT_Grid_ninclusion); // Macro
    IR_GIVE_FIELD(ir, nrefinement, _IFT_Grid_nrefinement); // Macro
    ir.finish();

    //read vertices
    //    inputVertexList->growTo(nvertex);
    inputVertexList.resize(nvertex, nullptr);
    
    for ( i = 0; i < nvertex; i++ ) {


    auto &ir = dr->giveInputRecord(GeneratorDataReader::GIR_vertexRec, i + 1);

    std::string name;
    int num = 0;
    IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

    if (num < 1 || num > nvertex) {
        std::cerr << "instanciateYourself: Invalid vertex number (num=" << num << ")\n";
        std::exit(EXIT_FAILURE);
    }

    if (generator::includes1(inputVertexList, num)) {
        std::cerr << "instanciateYourself: Vertex entry already exists (num=" << num << ")\n";
        std::exit(EXIT_FAILURE);
    }

    Vertex *vertex = new Vertex(num, this);
    vertex->initializeFrom(ir);
    setInputVertex(num, vertex); // 1-based -> vector[num-1]

    ir.finish();

    }


    //read control vertices
    //    controlVertexList->growTo(ncontrolvertex);
    controlVertexList.resize(ncontrolvertex, nullptr);

    for ( i = 0; i < ncontrolvertex; i++ ) {

    auto &ir = dr->giveInputRecord(GeneratorDataReader::GIR_vertexRec, i + 1);

    std::string name;
    int num = 0;
    IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

    // 1) validate number BEFORE touching containers
    if (num < 1 || num > ncontrolvertex) {
        std::cerr << "instanciateYourself: Invalid vertex number (num=" << num << ")\n";
        std::exit(EXIT_FAILURE);
    }

    // 2) duplicate check
    if (generator::includes1(inputVertexList, num)) {
        std::cerr << "instanciateYourself: Vertex entry already exists (num=" << num << ")\n";
        std::exit(EXIT_FAILURE);
    }

    // 3) create, init, and insert
    Vertex *vertex = new Vertex(num, this);
    vertex->initializeFrom(ir);
    setInputVertex(num, vertex); // 1-based -> vector[num-1]

    ir.finish();

    }
    
    // read curves
    //    curveList->growTo(ncurve);
    curveList.resize(ncurve, nullptr);
    for ( i = 0; i < ncurve; i++ ) {

    auto &ir = dr->giveInputRecord(GeneratorDataReader::GIR_curveRec, i + 1);

    std::string name;
    int num = 0;
    IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

    if (num < 1 || num > ncurve) {
        std::cerr << "instanciateYourself: Invalid curve number (num=" << num << ")\n";
        std::exit(EXIT_FAILURE);
    }

    if (generator::includes1(curveList, num)) {
        std::cerr << "instanciateYourself: Curve entry already exists (num=" << num << ")\n";
        std::exit(EXIT_FAILURE);
    }

    Curve *curve = new Curve(num, this);
    curve->initializeFrom(ir);
    setCurve(num, curve);
    ir.finish();
    }

    surfaceList.resize(nsurface, nullptr);
    for ( i = 0; i < nsurface; i++ ) {
      auto &ir = dr->giveInputRecord(GeneratorDataReader::GIR_surfaceRec, i + 1);
      std::string name;
      int num = 0;
      IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

      if (num < 1 || num > nsurface) {
        std::cerr << "instanciateYourself: Invalid curve number (num=" << num << ")\n";
        std::exit(EXIT_FAILURE);
    }

    if (generator::includes1(surfaceList, num)) {
        std::cerr << "instanciateYourself: Curve entry already exists (num=" << num << ")\n";
        std::exit(EXIT_FAILURE);
    }

    Surface *surface = new Surface(num, this);
    surface->initializeFrom(ir);
    setSurface(num, surface);
    ir.finish();
    }


regionList.resize(nregion, nullptr);

for (int i = 0; i < nregion; ++i) {
    auto &ir = dr->giveInputRecord(GeneratorDataReader::GIR_regionRec, i + 1);

    std::string name;
    int num = 0;
    IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

    if (num < 1 || num > nregion) {
        std::cerr << "instanciateYourself: Invalid region number (num=" << num << ")\n";
        std::exit(EXIT_FAILURE);
    }
    if (generator::includes1(regionList, num)) {
        std::cerr << "instanciateYourself: Region entry already exists (num=" << num << ")\n";
        std::exit(EXIT_FAILURE);
    }

    Region* r = nullptr;
    if (name == "prism") {
        r = new Prism(num, this);
    } else if (name == "cylinder") {
        r = new Cylinder(num, this);
    } else if (name == "sphere") {
        r = new Sphere(num, this);
    } else {
        std::cerr << "instanciateYourself: Unknown region type '" << name << "'\n";
        std::exit(EXIT_FAILURE);
    }

    r->initializeFrom(ir);
    setRegion(num, r);
    ir.finish();
}


inclusionList.resize(ninclusion, nullptr);

for (int i = 0; i < ninclusion; ++i) {
    auto &ir = dr->giveInputRecord(GeneratorDataReader::GIR_inclusionRec, i + 1);

    std::string name;
    int num = 0;
    IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

    // 1) validate number
    if (num < 1 || num > ninclusion) {
        std::cerr << "instanciateYourself: Invalid inclusion number (num=" << num << ")\n";
        std::exit(EXIT_FAILURE);
    }

    // 2) duplicate check
    if (generator::includes1(inclusionList, num)) {
        std::cerr << "instanciateYourself: Inclusion entry already exists (num=" << num << ")\n";
        std::exit(EXIT_FAILURE);
    }

    // 3) create the right type
    Inclusion *inc = nullptr;
    if (name == "intersphere") {
        inc = new InterfaceSphere(num, this);
    } else if (name == "interfacecylinder") {
        inc = new InterfaceCylinder(num, this);
    } else {
        std::cerr << "instanciateYourself: Unknown inclusion type '" << name << "'\n";
        std::exit(EXIT_FAILURE);
    }

    // 4) init + store
    inc->initializeFrom(ir);
    setInclusion(num, inc);

    ir.finish();
}
    
 refinementList.resize(nrefinement, nullptr);

 for (int i = 0; i < nrefinement; ++i) {
   auto &ir = dr->giveInputRecord(GeneratorDataReader::GIR_refinementRec, i + 1);
   
   std::string name;
   int num = 0;
   IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);
   
   // 1) validate number
   if (num < 1 || num > nrefinement) {
     std::cerr << "instanciateYourself: Invalid refinement number (num=" << num << ")\n";
     std::exit(EXIT_FAILURE);
   }

   // 2) duplicate check
    if (generator::includes1(refinementList, num)) {
        std::cerr << "instanciateYourself: Refinement entry already exists (num=" << num << ")\n";
        std::exit(EXIT_FAILURE);
    }

    // 3) construct by keyword
    Refinement *ref = nullptr;
    if (name == "refineprism") {
        ref = new RefinePrism(num, this);
    } else {
        std::cerr << "instanciateYourself: Unknown refinement type '" << name << "'\n";
        std::exit(EXIT_FAILURE);
    }

    // 4) init + store
    ref->initializeFrom(ir);
    setRefinement(num, ref);

    ir.finish();
}



    
  
    this->giveGridLocalizer()->init(true);

    return 1;
}




void Grid :: defineBoundaries(oofem::FloatArray &boundaries)
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

      

      
      oofem::FloatArray boundaries;
        defineBoundaries(boundaries);
        FILE *aggregateStream;
        if ( ( aggregateStream = fopen("aggregate.dat", "w") ) == NULL ) {
            printf("Can't open output file aggregate.dat");
            exit(1);
        }

        int outsideFlag;
        oofem::FloatArray coords;
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
