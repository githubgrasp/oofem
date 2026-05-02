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

#include "domain.h"
#include "sphere.h"
#include "interfacesphere.h"
#include "interfacecylinder.h"
#include "refineprism.h"
#include "interfaceplane.h"
#include "interfacesurface.h"
#include "cylinder.h"
#include "notch.h"
#include <iostream>
#include "generatorerror.h"

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

#include <fstream>
#include <iomanip>
#include <sstream>



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



Grid::Grid(int i)
// Constructor. Creates a new domain.
{
    gridLocalizer      = NULL;

    TOL = 0.;
}

Grid::~Grid()
// Destructor.
{
    for (auto v : vertexList) {
        delete v;
    }
    for (auto v : inputVertexList) {
        delete v;
    }
    for (auto v : controlVertexList) {
        delete v;
    }
    for (auto c : curveList) {
        delete c;
    }
    for (auto s : surfaceList) {
        delete s;
    }
    for (auto r : regionList) {
        delete r;
    }
    for (auto i : inclusionList) {
        delete i;
    }
    for (auto r : refinementList) {
        delete r;
    }
    for (auto n : notchList) {
        delete n;
    }

    delete gridLocalizer;
}



GridLocalizer *Grid::giveGridLocalizer()
{
    if ( gridLocalizer == NULL ) {
        gridLocalizer = new OctreeGridLocalizer(1, this);
    }

    return gridLocalizer;
}



Vertex *Grid::giveVertex(int n)
{
    if ( n >= 1 && n <= static_cast < int > ( generator::size1(vertexList) ) && vertexList [ n - 1 ] != nullptr ) {
        return vertexList [ n - 1 ];  // 1-based to 0-based
    }
    generator::errorf("Grid::giveVertex: undefined vertex (%d)", n);
}


Vertex *Grid::giveInputVertex(int n)
{
    if ( n >= 1 && n <= static_cast < int > ( generator::size1(inputVertexList) ) && inputVertexList [ n - 1 ] != nullptr ) {
        return inputVertexList [ n - 1 ];  // 1-based to 0-based
    }
    generator::errorf("Grid::giveInputVertex: undefined inputVertex (%d)", n);
}



Vertex *Grid::giveControlVertex(int n)
{
    if ( n >= 1 && n <= static_cast < int > ( generator::size1(controlVertexList) ) && controlVertexList [ n - 1 ] != nullptr ) {
        return controlVertexList [ n - 1 ];  // 1-based to 0-based
    }
    generator::errorf("Grid::giveControlVertex: undefined controlVertex (%d)", n);
}



Curve *Grid::giveCurve(int n)
{
    if ( n >= 1 && n <= static_cast < int > ( generator::size1(curveList) ) && curveList [ n - 1 ] != nullptr ) {
        return curveList [ n - 1 ];  // 1-based to 0-based
    }
    generator::errorf("Grid::giveCurve: undefined curve (%d)", n);
}

Surface *Grid::giveSurface(int n)
{
    if ( n >= 1 && n <= static_cast < int > ( generator::size1(surfaceList) ) && surfaceList [ n - 1 ] != nullptr ) {
        return surfaceList [ n - 1 ];  // 1-based to 0-based
    }
    generator::errorf("Grid::giveSurface: undefined surface (%d)", n);
}

double Grid::giveDiameter(oofem::FloatArray &coords) {
    double diam = this->diameter;
    for (int n = 1; n <= this->giveNumberOfRefinements(); n++) {
        diam = this->giveRefinement(n)->giveDiameter(coords);
    }

    return diam;
}


Region *Grid::giveRegion(int n)
{
    if ( n >= 1 && n <= static_cast < int > ( generator::size1(regionList) ) && regionList [ n - 1 ] != nullptr ) {
        return regionList [ n - 1 ];  // 1-based to 0-based
    }
    generator::errorf("Grid::giveRegion: undefined region (%d)", n);
}

Inclusion *Grid::giveInclusion(int n)
{
    if ( n >= 1 && n <= static_cast < int > ( generator::size1(inclusionList) ) && inclusionList [ n - 1 ] != nullptr ) {
        return inclusionList [ n - 1 ];  // 1-based to 0-based
    }
    generator::errorf("Grid::giveInclusion: undefined inclusion (%d)", n);
}

Refinement *Grid::giveRefinement(int n)
{
    if ( n >= 1 && n <= static_cast < int > ( generator::size1(refinementList) ) && refinementList [ n - 1 ] != nullptr ) {
        return refinementList [ n - 1 ];  // 1-based to 0-based
    }
    generator::errorf("Grid::giveRefinement: undefined refinement (%d)", n);
}


Notch *Grid::giveNotch(int n)
{
    if ( n >= 1 && n <= static_cast < int > ( generator::size1(notchList) ) && notchList [ n - 1 ] != nullptr ) {
        return notchList [ n - 1 ];
    }
    generator::errorf("Grid::giveNotch: undefined notch (%d)", n);
}


int Grid::giveNumberOfNotches() const { return generator::size1(notchList); }


bool Grid::insideAnyNotch(const oofem::FloatArray &coord) const
{
    for ( const auto *n : notchList ) {
        if ( n != nullptr && n->containsStrictly(coord, this->TOL) ) {
            return true;
        }
    }
    return false;
}


void Grid::resizeVertices(int _newSize) { vertexList.resize(_newSize, nullptr); }
void Grid::resizeInputVertices(int _newSize) { inputVertexList.resize(_newSize, nullptr); }
void Grid::resizeControlVertices(int _newSize) { controlVertexList.resize(_newSize, nullptr); }
void Grid::resizeCurves(int _newSize) { curveList.resize(_newSize, nullptr); }
void Grid::resizeSurfaces(int _newSize) { surfaceList.resize(_newSize, nullptr); }
void Grid::resizeRegions(int _newSize) { regionList.resize(_newSize, nullptr); }
void Grid::resizeInclusions(int _newSize) { inclusionList.resize(_newSize, nullptr); }
void Grid::resizeRefinements(int _newSize) { refinementList.resize(_newSize, nullptr); }


void Grid::setVertex(int i, Vertex *obj) {
    generator::put1_replace(vertexList, i, obj);
    if ( i > vertexCount ) {
        vertexCount = i;
    }
}

void Grid::setInputVertex(int i, Vertex *obj) {
    generator::put1_replace(inputVertexList, i, obj);
}

void Grid::setControlVertex(int i, Vertex *obj) {
    generator::put1_replace(controlVertexList, i, obj);
}

void Grid::setCurve(int i, Curve *obj) {
    generator::put1_replace(curveList, i, obj);
}

void Grid::setSurface(int i, Surface *obj) {
    generator::put1_replace(surfaceList, i, obj);
}

void Grid::setRegion(int i, Region *obj) {
    generator::put1_replace(regionList, i, obj);
}

void Grid::setInclusion(int i, Inclusion *obj) {
    generator::put1_replace(inclusionList, i, obj);
}

void Grid::setRefinement(int i, Refinement *obj) {
    generator::put1_replace(refinementList, i, obj);
}

int Grid::generatePoints()
{
    //@todo: Is the regular generation the same for normal and periodic?
    if ( this->regularFlag == 1 ) { //regular
        this->generateRegularPoints();
    } else { //random (irregular or periodic)
        this->generateRandomPoints();
    }
    return 1;
}

int Grid::generateRandomPoints()
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

    start = ::clock();

    if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0 && randomFlag != 2 ) {//Region without periodicity
        if ( randomFlag == 1 ) { //Input vertices for Grassl's mesh approach
            this->generateInputPoints();
        }


        //generate inclusions
        sc = ::clock();
        //Inclusions
        //@todo: Inclusions need to have a random and regular point generation
        for ( int i = 0; i < this->giveNumberOfInclusions(); i++ ) {
            ( this->giveInclusion(i + 1) )->generatePoints();
            printf("numberOfVertices = %d\n", this->giveNumberOfVertices() );
        }
        ec = ::clock();
        nsec = ( ec - sc ) / CLOCKS_PER_SEC;
        printf("Points for inclusions generated in %lds \n", nsec);


        //Control vertices
        this->generateControlPoints(); //This should now include also periodic shifts

        if ( randomFlag == 1 ) { //Input vertices for Grassl's mesh approach
            //Curves
            sc = ::clock();
            for (int i = 0; i < this->giveNumberOfCurves(); i++ ) {
                ( this->giveCurve(i + 1) )->generatePoints();
                printf("numberOfVertices = %d\n", this->giveNumberOfVertices() );
            }
            ec = ::clock();
            nsec = ( ec - sc ) / CLOCKS_PER_SEC;
            printf("Points on curves generated in %lds \n", nsec);


            //Surfaces
            sc = ::clock();
            for (int i = 0; i < this->giveNumberOfSurfaces(); i++ ) {
                ( this->giveSurface(i + 1) )->generatePoints();
                printf("numberOfVertices = %d\n", this->giveNumberOfVertices() );
            }
            ec = ::clock();
            nsec = ( ec - sc ) / CLOCKS_PER_SEC;
            printf("Points on surfaces generated in %lds \n", nsec);
        }

        //Regions
        sc = ::clock();
        for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
            ( this->giveRegion(i + 1) )->generatePoints();
            printf("numberOfVertices = %d\n", this->giveNumberOfVertices() );
        }
        ec = ::clock();
        nsec = ( ec - sc ) / CLOCKS_PER_SEC;
        printf("Points in regions generated in %lds \n", nsec);

        nsec = ( ec - start ) / CLOCKS_PER_SEC;
        printf("All points generated in %lds \n", nsec);
    } else if ( this->randomFlag == 2 ) {   //Generate periodic surfaces for Adam
        /*Each direction should be generated only once and then shifted by the specimen dimension*/
        oofem::IntArray lineCounter(3);
        oofem::FloatArray lineNormal(3);
        lineCounter.zero();

        this->generateInputPoints();

        this->generateControlPoints(); //This should now include also periodic shifts

        sc = ::clock();
        for (int i = 0; i < this->giveNumberOfCurves(); i++ ) {
            this->giveCurve(i + 1)->giveNormal(lineNormal);
            if ( lineNormal.at(1) == 1 && lineNormal.at(2) == 0 && lineNormal.at(3) == 0 && lineCounter.at(1) == 0 ) {
                ( this->giveCurve(i + 1) )->generatePoints();
                printf("numberOfVertices = %d\n", this->giveNumberOfVertices() );
                lineCounter.at(1) = 1;
            } else if ( lineNormal.at(1) == 0 && lineNormal.at(2) == 1 && lineNormal.at(3) == 0 && lineCounter.at(2) == 0 ) {
                ( this->giveCurve(i + 1) )->generatePoints();
                printf("numberOfVertices = %d\n", this->giveNumberOfVertices() );
                lineCounter.at(2) = 1;
            } else if ( lineNormal.at(1) == 0 && lineNormal.at(2) == 0 && lineNormal.at(3) == 1 && lineCounter.at(3) == 0 ) {
                ( this->giveCurve(i + 1) )->generatePoints();
                printf("numberOfVertices = %d\n", this->giveNumberOfVertices() );
                lineCounter.at(3) = 1;
            }
        }
        ec = ::clock();
        nsec = ( ec - sc ) / CLOCKS_PER_SEC;
        printf("Points on curves generated in %lds \n", nsec);

        //Surfaces
        oofem::IntArray surfaceCounter(3);
        oofem::FloatArray surfaceNormal(3);
        surfaceCounter.zero();

        sc = ::clock();
        for (int i = 0; i < this->giveNumberOfSurfaces(); i++ ) {
            this->giveSurface(i + 1)->giveNormal(surfaceNormal);
            if ( surfaceNormal.at(1) == 1 && surfaceNormal.at(2) == 0 && surfaceNormal.at(3) == 0 && surfaceCounter.at(1) == 0 ) {
                ( this->giveSurface(i + 1) )->generatePoints();
                printf("numberOfVertices = %d\n", this->giveNumberOfVertices() );
                surfaceCounter.at(1) = 1;
            } else if ( surfaceNormal.at(1) == 0 && surfaceNormal.at(2) == 1 && surfaceNormal.at(3) == 0 && surfaceCounter.at(2) == 0 ) {
                ( this->giveSurface(i + 1) )->generatePoints();
                printf("numberOfVertices = %d\n", this->giveNumberOfVertices() );
                surfaceCounter.at(2) = 1;
            } else if ( surfaceNormal.at(1) == 0 && surfaceNormal.at(2) == 0 && surfaceNormal.at(3) == 1 && surfaceCounter.at(3) == 0 ) {
                ( this->giveSurface(i + 1) )->generatePoints();
                printf("numberOfVertices = %d\n", this->giveNumberOfVertices() );
                surfaceCounter.at(3) = 1;
            }
        }
        ec = ::clock();
        nsec = ( ec - sc ) / CLOCKS_PER_SEC;
        printf("Points on surfaces generated in %lds \n", nsec);

        //Regions
        sc = ::clock();
        for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
            ( this->giveRegion(i + 1) )->generatePoints();
            printf("numberOfVertices = %d\n", this->giveNumberOfVertices() );
        }
        ec = ::clock();
        nsec = ( ec - sc ) / CLOCKS_PER_SEC;
        printf("Points in regions generated in %lds \n", nsec);

        nsec = ( ec - start ) / CLOCKS_PER_SEC;
        printf("All points generated in %lds \n", nsec);
    } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0 ) {
        //This is the case of a region with periodicity in x-direction, as is the case for a periodic beam
        sc = ::clock();

        //Curves
        for (int i = 0; i < this->giveNumberOfCurves(); i++ ) {
            this->giveCurve(i + 1)->giveNormal(lineNormal);
            if ( lineNormal.at(1) == 1 && lineNormal.at(2) == 0 && lineNormal.at(3) == 0 ) {// line in x-direction
                ( this->giveCurve(i + 1) )->generatePoints();
                printf("numberOfVertices = %d\n", this->giveNumberOfVertices() );
            }
            ec = ::clock();
            nsec = ( ec - sc ) / CLOCKS_PER_SEC;
            printf("Points on curves generated in %lds \n", nsec);
        }

        //Surfaces
        sc = ::clock();
        for (int i = 0; i < this->giveNumberOfSurfaces(); i++ ) {
            this->giveSurface(i + 1)->giveNormal(surfaceNormal);
            if ( !( surfaceNormal.at(1) == 1 && surfaceNormal.at(2) == 0 && surfaceNormal.at(3) == 0 ) ) {// all surfaces without normal in x-direction
                ( this->giveSurface(i + 1) )->generatePoints();
                printf("numberOfVertices = %d\n", this->giveNumberOfVertices() );
            }
        }
        ec = ::clock();
        nsec = ( ec - sc ) / CLOCKS_PER_SEC;
        printf("Points on surfaces generated in %lds \n", nsec);

        //Regions
        sc = ::clock();
        for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
            ( this->giveRegion(i + 1) )->generatePoints();
        }
        ec = ::clock();
        nsec = ( ec - sc ) / CLOCKS_PER_SEC;
        printf("Points in regions generated in %lds \n", nsec);
    } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
        //Region with periodicity in x and y directions. This could be the case of a slab.

        this->generateControlPoints(); //This should now include also periodic shifts

        //Surfaces
        sc = ::clock();
        for (int i = 0; i < this->giveNumberOfSurfaces(); i++ ) {
            this->giveSurface(i + 1)->giveNormal(surfaceNormal);
            if ( surfaceNormal.at(1) == 0 && surfaceNormal.at(2) == 0 && surfaceNormal.at(3) == 1 ) {// all surfaces without normal in x-direction
                ( this->giveSurface(i + 1) )->generatePoints();
                printf("numberOfVertices = %d\n", this->giveNumberOfVertices() );
            }
        }
        ec = ::clock();
        nsec = ( ec - sc ) / CLOCKS_PER_SEC;
        printf("Points on surfaces generated in %lds \n", nsec);

        //Regions
        sc = ::clock();
        for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
            ( this->giveRegion(i + 1) )->generatePoints();
        }
        ec = ::clock();
        nsec = ( ec - sc ) / CLOCKS_PER_SEC;
        printf("Points in regions generated in %lds \n", nsec);
    } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
        //Region with periodicity in all directions.

        sc = ::clock();
        this->generateControlPoints();
        ec = ::clock();
        nsec = ( ec - sc ) / CLOCKS_PER_SEC;
        printf("Control points generated in %lds \n", nsec);

        sc = ::clock();
        //Inclusions
        //@todo: Inclusions need to have a random and regular point generation
        for ( int i = 0; i < this->giveNumberOfInclusions(); i++ ) {
            ( this->giveInclusion(i + 1) )->generatePeriodicPoints();
            printf("numberOfVertices = %d\n", this->giveNumberOfVertices() );
        }
        ec = ::clock();
        nsec = ( ec - sc ) / CLOCKS_PER_SEC;
        printf("Points for inclusions generated in %lds \n", nsec);

        //Regions
        sc = ::clock();
        for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
            ( this->giveRegion(i + 1) )->generatePoints();
        }
        ec = ::clock();
        nsec = ( ec - sc ) / CLOCKS_PER_SEC;
        printf("Points in regions generated in %lds \n", nsec);
    } else {  //should not be here
        printf("Something wrong with mixed point generation and periodic flag in grid->generateRandomPoints\n");
    }

    return 1;
}

// grid.C
bool Grid::addVertex(const oofem::FloatArray &coords) {
    // Suppress placement strictly inside any #@notch box. Boundary points
    // (within TOL of a face) are kept so the dual mesh can partition cleanly
    // along the notch surface.
    if ( insideAnyNotch(coords) ) {
        return false;
    }

    if ( vertexList.capacity() < static_cast < size_t > ( vertexCount + 1 ) ) {
        vertexList.reserve(vertexCount + 10000); // choose a chunk size that fits your cases
    }

    const int num = vertexCount + 1;         // true next id
    auto * v = new Vertex(num, this);
    v->setCoordinates(coords);               // prefer const-ref setter

    setVertex(num, v);                       // stores (may resize/replace)
    this->vertexCount = num;                       // <-- bump the real count

    if ( auto *loc = giveGridLocalizer() ) {
        loc->insertSequentialNode(num, coords);
    }
    return true;
}


void Grid::generateControlPoints()
{
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
        this->giveControlVertex(i + 1)->giveCoordinates(coords);
        this->addVertex(coords);

        if ( this->randomFlag == 0 ) { //Here it is assumed that for Grassl's meshing approach control nodes are on the surface, and for Bolander's approach they are not.
            //Do now the mirroring.
            //It works only for rectangular specimen in a cartesian coordinate system.
            //x-direction


            for ( int k = 1; k <= 2; k++ ) {
                mirroredCoords.at(1) = coords.at(1) - 2. * ( coords.at(1) - boundaries.at(k) );
                mirroredCoords.at(2) = coords.at(2);
                mirroredCoords.at(3) = coords.at(3);
                this->addVertex(mirroredCoords);
            }

            //y-direction
            for ( int k = 3; k <= 4; k++ ) {
                mirroredCoords.at(1) = coords.at(1);
                mirroredCoords.at(2) = coords.at(2) - 2. * ( coords.at(2) - boundaries.at(k) );
                mirroredCoords.at(3) = coords.at(3);
                this->addVertex(mirroredCoords);
            }

            //z-direction
            for ( int k = 5; k <= 6; k++ ) {
                mirroredCoords.at(1) = coords.at(1);
                mirroredCoords.at(2) = coords.at(2);
                mirroredCoords.at(3) = coords.at(3) - 2. * ( coords.at(3) - boundaries.at(k) );
                this->addVertex(mirroredCoords);
            }
        }
    }

    return;
}




void Grid::generateInputPoints()
{
    printf("At the start of generateInputPoints = %d\n", this->giveNumberOfVertices() );
    //Generate control vertices
    oofem::FloatArray coords(3);

    int inputVertexNumber = generator::size1(this->inputVertexList);
    for ( int i = 0; i < inputVertexNumber; i++ ) {
        this->giveInputVertex(i + 1)->giveCoordinates(coords);
        this->addVertex(coords);
    }

    // Seed notch box boundaries (corners + edges + faces) at edgeRefine /
    // surfaceRefine spacing so the dual mesh partitions cleanly along the
    // notch surface. Done after input vertices so the spatial localiser is
    // already built around the bounding box, and before regions/inclusions
    // run so their distance checks see these as existing neighbours.
    for ( int i = 0; i < this->giveNumberOfNotches(); ++i ) {
        this->giveNotch(i + 1)->generateBoundaryPoints();
    }

    printf("At the end of generateInputPoints = %d\n", this->giveNumberOfVertices() );
    return;
}



int Grid::generateRegularPoints()
{
    //Regular approach for both normal and periodic
    clock_t start;
    clock_t sc;

    start = ::clock();
    //The old method is reviced, the only generation needed is the region one.
    //The user has to set in the input a larger region, than the one in the converter, in order to avoid boundary discrepancies.
    //Regions
    for ( int i = 0; i < this->giveNumberOfRegions(); i++ ) {
        if ( regType == 1 ) { //BCC
            ( this->giveRegion(i + 1) )->generateRegularPoints1();
        } else if ( regType == 2 ) {     //FCC
            ( this->giveRegion(i + 1) )->generateRegularPoints2();
        } else {
            generator::error("Grid::generateRegularPoints: regular grid type was not determined");
        }
    }
    sc = ::clock();
    double elapsed_sec = double( sc - start ) / CLOCKS_PER_SEC;
    std::cout << "Generation took " << elapsed_sec << " seconds\n";

    return 1;
}


int Grid::readControlRecords(const std::string &controlFile)
{
    std::ifstream in(controlFile);
    if ( !in ) {
        generator::errorf("Grid::readControlRecords: cannot open control file '%s'", controlFile.c_str() );
    }

    // Match legacy defaults from instanciateYourself().
    maxIter       = 1000000;
    regularFlag   = 0;
    randomFlag    = 0;
    vtkFlag       = 0;
    randomInteger = 0;
    periodicityFlag.resize(3);
    periodicityFlag.zero();

    std::string line;
    while ( std::getline(in, line) ) {
        std::istringstream iss(line);
        std::string tag;
        if ( !( iss >> tag ) ) {
            continue;
        }
        if ( tag == "#@output" ) {
            iss >> outputFileName;
        } else if ( tag == "#@domain" ) {
            iss >> dimension;
        } else if ( tag == "#@diam" ) {
            iss >> diameter;
            TOL = 1.e-6 * diameter;
        } else if ( tag == "#@maxiter" ) {
            iss >> maxIter;
        } else if ( tag == "#@ranint" ) {
            iss >> randomInteger;
            if ( randomInteger >= 0 ) {
                randomInteger = -time(NULL);
                printf("random integer is determined as %d\n", randomInteger);
            }
        } else if ( tag == "#@perflag" ) {
            int n;
            iss >> n;
            periodicityFlag.resize(n);
            for ( int i = 1; i <= n; ++i ) {
                iss >> periodicityFlag.at(i);
            }
            if ( periodicityFlag.giveSize() != 3 ) {
                generator::error("#@perflag must have three components");
            }
        } else if ( tag == "#@ranflag" ) {
            iss >> randomFlag;
            if ( randomFlag < 0 || randomFlag > 2 ) {
                generator::error("#@ranflag must be 0, 1, or 2");
            }
        } else if ( tag == "#@vtk" ) {
            // Opt in to writing points.vtk alongside the nodes.dat output.
            vtkFlag = 1;
        } else if ( tag == "#@vertex" ) {
            int num;
            iss >> num;
            if ( num < 1 ) {
                generator::errorf("#@vertex: invalid number %d", num);
            }
            auto *v = new Vertex(num, this);
            v->initializeFromTokens(iss);
            if ( ( int ) inputVertexList.size() < num ) {
                inputVertexList.resize(num, nullptr);
            }
            setInputVertex(num, v);
        } else if ( tag == "#@controlvertex" ) {
            int num;
            iss >> num;
            if ( num < 1 ) {
                generator::errorf("#@controlvertex: invalid number %d", num);
            }
            auto *v = new Vertex(num, this);
            v->initializeFromTokens(iss);
            if ( ( int ) controlVertexList.size() < num ) {
                controlVertexList.resize(num, nullptr);
            }
            setControlVertex(num, v);
        } else if ( tag == "#@curve" ) {
            int num;
            iss >> num;
            auto *c = new Curve(num, this);
            c->initializeFromTokens(iss);
            if ( ( int ) curveList.size() < num ) {
                curveList.resize(num, nullptr);
            }
            setCurve(num, c);
        } else if ( tag == "#@surface" ) {
            int num;
            iss >> num;
            auto *s = new Surface(num, this);
            s->initializeFromTokens(iss);
            if ( ( int ) surfaceList.size() < num ) {
                surfaceList.resize(num, nullptr);
            }
            setSurface(num, s);
        } else if ( tag == "#@prism" ) {
            int num;
            iss >> num;
            auto *r = new Prism(num, this);
            r->initializeFromTokens(iss);
            if ( ( int ) regionList.size() < num ) {
                regionList.resize(num, nullptr);
            }
            setRegion(num, r);
        } else if ( tag == "#@cylinder" ) {
            int num;
            iss >> num;
            auto *r = new Cylinder(num, this);
            r->initializeFromTokens(iss);
            if ( ( int ) regionList.size() < num ) {
                regionList.resize(num, nullptr);
            }
            setRegion(num, r);
        } else if ( tag == "#@sphere" ) {
            int num;
            iss >> num;
            auto *r = new Sphere(num, this);
            r->initializeFromTokens(iss);
            if ( ( int ) regionList.size() < num ) {
                regionList.resize(num, nullptr);
            }
            setRegion(num, r);
        } else if ( tag == "#@intersphere" ) {
            int num;
            iss >> num;
            auto *inc = new InterfaceSphere(num, this);
            inc->initializeFromTokens(iss);
            if ( ( int ) inclusionList.size() < num ) {
                inclusionList.resize(num, nullptr);
            }
            setInclusion(num, inc);
        } else if ( tag == "#@interfacecylinder" ) {
            int num;
            iss >> num;
            auto *inc = new InterfaceCylinder(num, this);
            inc->initializeFromTokens(iss);
            if ( ( int ) inclusionList.size() < num ) {
                inclusionList.resize(num, nullptr);
            }
            setInclusion(num, inc);
        } else if ( tag == "#@interfaceplane" ) {
            int num;
            iss >> num;
            auto *inc = new InterfacePlane(num, this);
            inc->initializeFromTokens(iss);
            if ( ( int ) inclusionList.size() < num ) {
                inclusionList.resize(num, nullptr);
            }
            setInclusion(num, inc);
        } else if ( tag == "#@interfacesurface" ) {
            int num;
            iss >> num;
            auto *inc = new InterfaceSurface(num, this);
            inc->initializeFromTokens(iss);
            if ( ( int ) inclusionList.size() < num ) {
                inclusionList.resize(num, nullptr);
            }
            setInclusion(num, inc);
        } else if ( tag == "#@refineprism" ) {
            int num;
            iss >> num;
            auto *ref = new RefinePrism(num, this);
            ref->initializeFromTokens(iss);
            if ( ( int ) refinementList.size() < num ) {
                refinementList.resize(num, nullptr);
            }
            setRefinement(num, ref);
        } else if ( tag == "#@notch" ) {
            int num;
            iss >> num;
            auto *nh = new Notch(num, this);
            nh->initializeFromTokens(iss);
            if ( ( int ) notchList.size() < num ) {
                notchList.resize(num, nullptr);
            }
            notchList [ num - 1 ] = nh;
        } else if ( tag == "#@inclusionfile" ) {
            std::string path;
            iss >> path;
            double itz = 0.0;
            double refine = 1.0;
            std::string sub;
            while ( iss >> sub ) {
                if ( sub == "itz" ) {
                    iss >> itz;
                } else if ( sub == "refine" ) {
                    iss >> refine;
                } else {
                    generator::errorf("Grid::readControlRecords: '#@inclusionfile' unknown sub-keyword '%s'", sub.c_str());
                }
            }
            this->readInclusionFile(path, itz, refine);
        } else if ( tag.rfind("#@", 0) == 0 ) {
            generator::errorf("Grid::readControlRecords: unknown directive '%s'", tag.c_str() );
        }
        // Non-directive lines (blank, `#` comments, anything not starting
        // with `#@`) are skipped silently.
    }

    if ( outputFileName.empty() ) {
        generator::error("Grid::readControlRecords: no #@output directive found");
    }
    if ( diameter <= 0. ) {
        generator::error("Grid::readControlRecords: missing or invalid #@diam");
    }

    this->giveGridLocalizer()->init(true);

    return 1;
}


void Grid::readInclusionFile(const std::string &path, double itz, double refine)
{
    std::ifstream in(path);
    if ( !in ) {
        generator::errorf("Grid::readInclusionFile: cannot open '%s'", path.c_str());
    }

    int sphereCount = 0;
    int ellipsoidCount = 0;
    int fibreCount = 0;
    std::string line;
    while ( std::getline(in, line) ) {
        if ( line.empty() || line[0] == '#' ) {
            continue; // skip metadata comments
        }
        std::istringstream iss(line);
        std::string keyword;
        if ( !( iss >> keyword ) ) {
            continue;
        }
        if ( keyword == "sphere" ) {
            int packingId;
            iss >> packingId;
            // Append the directive's refine/itz so the existing
            // InterfaceSphere parser sees one well-formed token stream.
            std::string remainder;
            std::getline(iss, remainder);
            std::ostringstream merged;
            merged << remainder << " refine " << refine << " itz " << itz;
            std::istringstream mergedStream(merged.str());

            // Renumber to avoid collisions with any inclusions added by
            // other directives (`#@intersphere`, etc.).
            const int newNumber = static_cast<int>(inclusionList.size()) + 1;
            auto *inc = new InterfaceSphere(newNumber, this);
            inc->initializeFromTokens(mergedStream);
            inclusionList.resize(newNumber, nullptr);
            setInclusion(newNumber, inc);
            ++sphereCount;
        } else if ( keyword == "ellipsoid" ) {
            ++ellipsoidCount;
        } else if ( keyword == "fibre" ) {
            ++fibreCount; // intentionally not loaded — converter handles fibres
        }
        // Unknown keywords are ignored to stay forward-compatible.
    }

    if ( ellipsoidCount > 0 ) {
        std::fprintf(stderr,
                     "Warning: Grid::readInclusionFile: %d ellipsoid line(s) in '%s' ignored "
                     "(generator only handles spheres from packing files)\n",
                     ellipsoidCount, path.c_str());
    }
    std::printf("readInclusionFile('%s'): %d sphere(s) loaded, %d fibre(s) skipped\n",
                path.c_str(), sphereCount, fibreCount);
}



void Grid::exportVTK(const std::string &path)
{
    std::ofstream out(path);
    if ( !out ) {
        std::cerr << "Failed to open " << path << " for writing\n";
        return;
    }

    const int n = this->giveNumberOfVertices();
    out << "# vtk DataFile Version 3.0\n";
    out << "grid points\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << n << " float\n";

    for (int i = 1; i <= n; ++i) {
        const oofem::FloatArray *c = this->giveVertex(i)->giveCoordinates();
        double x = c->giveSize() >= 1 ? c->at(1) : 0.0;
        double y = c->giveSize() >= 2 ? c->at(2) : 0.0;
        double z = c->giveSize() >= 3 ? c->at(3) : 0.0; // pad to 3D
        out << std::setprecision(9) << x << " " << y << " " << z << "\n";
    }

    // Add a vertex cell for each point so they render as glyphs
    out << "VERTICES " << n << " " << 2 * n << "\n";
    for (int i = 0; i < n; ++i) {
        out << "1 " << i << "\n";
    }

    // Optional: save the original 1-based IDs as point data
    out << "POINT_DATA " << n << "\n";
    out << "SCALARS id int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int i = 1; i <= n; ++i) {
        out << i << "\n";
    }
}


void Grid::defineBoundaries(oofem::FloatArray &boundaries)
//Determine the boundaries of the domain
{
    if ( this->giveNumberOfRegions() > 0 ) {
        if ( this->giveNumberOfRegions() != 1 ) {
            generator::error("Grid::defineBoundaries: cannot define boundaries for multiple regions yet");
        }
        this->giveRegion(1)->defineBoundaries(boundaries);
    } else if ( this->giveNumberOfSurfaces() > 0 ) {
        if ( this->giveNumberOfSurfaces() != 1 ) {
            generator::error("Grid::defineBoundaries: cannot define boundaries for multiple surfaces yet");
        }
        this->giveSurface(1)->defineBoundaries(boundaries);
        return;
    }
    return;
}

void Grid::giveOutput(FILE *outputStream)
{
    int dimension = 3;
    fprintf(outputStream, "%d\n", dimension);

    fprintf(outputStream, "%d\n", this->giveNumberOfVertices() );

    for ( int i = 0; i < this->giveNumberOfVertices(); i++ ) {
        fprintf(outputStream, "%.16e %.16e %.16e\n", ( this->giveVertex(i + 1) )->giveCoordinate(1), ( this->giveVertex(i + 1) )->giveCoordinate(2), ( this->giveVertex(i + 1) )->giveCoordinate(3) );
    }
}


double Grid::ran1(int *idum)
{
    // random number generator from "Numerical recipes book"
    int j;
    long k;
    static long iy = 0;
    static long iv[ NTAB ];
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
