#include "fibre.h"
#include "vertex.h"

#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include <iostream>
 #include "line.h"
 #include <vector>
#endif


Fibre::Fibre(int n, Grid *aGrid) : GridComponent(n, aGrid)
{
    this->number = n;
    this->tolerance = 1.0e-10;
}




void
Fibre::initializeFromCoords(const oofem::FloatArray &endpts, double diam)
{
    endpoints = endpts;
    diameter  = diam;
    coordP.resize(3);
    coordQ.resize(3);
    coordP.at(1) = endpts.at(1);
    coordP.at(2) = endpts.at(2);
    coordP.at(3) = endpts.at(3);
    coordQ.at(1) = endpts.at(4);
    coordQ.at(2) = endpts.at(5);
    coordQ.at(3) = endpts.at(6);

    pointP = grid->createInterNode(coordP);
    pointQ = grid->createInterNode(coordQ);

    length          = computeDistance(coordP, coordQ);
    directionVector = ( 1. / length ) * ( coordQ - coordP );
}


void
Fibre::findDelaunayVertex()
{
    const int n = ( int ) intersectionPoints.size();
    if ( n <= 0 ) {
        printf("error in findDelaunayVertex in discretization of fibres\n");
        return;
    }

    oofem::FloatArray coordsOne, coordsTwo;
    Vertex *interP = grid->giveInterNode(intersectionPoints.at(n - 1));
    interP->giveCoordinates(coordsOne);

    oofem::IntArray nodeCandidates =
        grid->findDelaunayNodesWithinBox(coordsOne, 2. * grid->giveDiameter());

    const int numberOfDelaunayCandidates = nodeCandidates.giveSize();
    if ( numberOfDelaunayCandidates < 1 ) {
        printf("error in discretization of fibre: preselection in localization of first Delaunay vertex must be wider\n");
        exit(1);
    }

    int    bestIndex   = 1;
    double bestDistance = 2. * grid->giveDiameter();
    for ( int i = 1; i <= numberOfDelaunayCandidates; i++ ) {
        grid->giveDelaunayVertex(nodeCandidates.at(i))->giveCoordinates(coordsTwo);
        const double d = computeDistance(coordsOne, coordsTwo);
        if ( d < bestDistance ) {
            bestIndex    = nodeCandidates.at(i);
            bestDistance = d;
        }
    }
    delaunayVertices.push_back(bestIndex);
}


void
Fibre::findIntersect()
{
    const int numberOfIntersections = ( int ) intersectionPoints.size();

    if ( numberOfIntersections == 0 ) {
        // First call: seed the walk with the start endpoint.
        intersectionPoints.push_back(pointP->giveNumber());
        return;
    }

    const int numberOfDelaunayPicked = ( int ) delaunayVertices.size();
    const int currentDelaunayNumber  = delaunayVertices.at(numberOfDelaunayPicked - 1);

    oofem::IntArray neighbourLines;
    grid->giveDelaunayVertex(delaunayVertices.at(numberOfDelaunayPicked - 1))->giveLocalLines(neighbourLines);

    const int numberOfNeighbours = neighbourLines.giveSize();
    oofem::IntArray neighbourNumbers(numberOfNeighbours);
    oofem::IntArray nodes(2);

    for ( int i = 0; i < numberOfNeighbours; i++ ) {
        grid->giveDelaunayLine(neighbourLines(i))->giveLocalVertices(nodes);
        if ( nodes.at(1) == currentDelaunayNumber ) {
            neighbourNumbers.at(i + 1) = nodes.at(2);
        } else if ( nodes.at(2) == currentDelaunayNumber ) {
            neighbourNumbers.at(i + 1) = nodes.at(1);
        } else {
            printf("error in the discretization of fibre, bad use of data\n");
        }
    }

    oofem::FloatArray coordI, coordJ, coordM, coordU, diffIJ, diffMU, diffQU, coordS;
    grid->giveDelaunayVertex(delaunayVertices.at(numberOfDelaunayPicked - 1))->giveCoordinates(coordI);
    grid->giveInterNode(intersectionPoints.at(numberOfIntersections - 1))->giveCoordinates(coordU);

    double rhoMin = 2.;        // any computed rho ≤ 1 will replace this
    int    bestNeighbour = 0;

    for ( int i = 0; i < numberOfNeighbours; i++ ) {
        grid->giveDelaunayVertex(neighbourNumbers.at(i + 1))->giveCoordinates(coordJ);
        coordM = 0.5 * ( coordI + coordJ );
        diffMU = coordM - coordU;
        diffIJ = coordI - coordJ;
        diffQU = coordQ - coordU;
        const double scalar1 = diffIJ.dotProduct(diffMU, 3);
        const double scalar2 = diffIJ.dotProduct(diffQU, 3);
        if ( scalar2 == 0 ) continue;        // parallel — skip

        const double rho = scalar1 / scalar2;
        if ( rho < rhoMin && rho > tolerance ) {
            rhoMin        = rho;
            bestNeighbour = neighbourNumbers.at(i + 1);
        }
    }

    if ( rhoMin < 1. ) {
        coordS = coordU + rhoMin * diffQU;
        intersectionPoints.push_back(grid->createInterNode(coordS)->giveNumber());
        delaunayVertices.push_back(grid->giveDelaunayVertex(bestNeighbour)->giveNumber());
    } else if ( rhoMin >= 1. ) {
        // Reached the end of the fibre — terminate with the end endpoint.
        intersectionPoints.push_back(pointQ->giveNumber());
    } else {
        printf("error in the discretization of fibre\n");
    }
}


void
Fibre::placeReinforcementNode()
{
    oofem::FloatArray coord1, coord2, coordR;
    const int n = ( int ) intersectionPoints.size();
    grid->giveInterNode(intersectionPoints.at(n - 2))->giveCoordinates(coord1);
    grid->giveInterNode(intersectionPoints.at(n - 1))->giveCoordinates(coord2);
    coordR = 0.5 * ( coord1 + coord2 );

    reinforcementPoints.push_back(grid->createReinfNode(coordR)->giveNumber());

    const double distFromP = computeDistance(coordR, coordP);
    endLengths.push_back(std::min(distFromP, length - distFromP));
}


void
Fibre::discretize()
{
    findIntersect();        // seed the walk with start endpoint
    findDelaunayVertex();

    int iter = 0;
    while ( true ) {
        if ( ++iter > 1000 ) {
            printf("Error in discretization: max iteration reached. Adjust `tolerance` to avoid instabilities.\n");
            exit(1);
        }
        const int n = ( int ) intersectionPoints.size();
        if ( intersectionPoints.at(n - 1) == pointQ->giveNumber() ) break;
        findIntersect();
        placeReinforcementNode();
    }
}


Fibre *
Fibre::ofType()
{
    return new Fibre(number, grid);
}


Vertex *
Fibre::findNearestReinforcementNode(const oofem::FloatArray &coord, Grid *grid, double tol)
{
    Vertex *bestNode = nullptr;
    oofem::FloatArray coordN;

    // First pass: any node within `tol`.
    for ( int i = 0; i < grid->giveNumberOfReinforcementNode(); i++ ) {
        grid->giveReinforcementNode(i + 1)->giveCoordinates(coordN);
        if ( computeDistance(coord, coordN) <= tol ) {
            return grid->giveReinforcementNode(i + 1);
        }
    }

    // Fallback: nearest by distance, with a warning.
    double bestDistance = std::numeric_limits< double >::infinity();
    for ( int i = 0; i < grid->giveNumberOfReinforcementNode(); i++ ) {
        grid->giveReinforcementNode(i + 1)->giveCoordinates(coordN);
        const double d = computeDistance(coord, coordN);
        if ( d <= bestDistance ) {
            bestDistance = d;
            bestNode     = grid->giveReinforcementNode(i + 1);
        }
    }
    printf("findNearestReinforcementNode: no node within tol=%e; nearest is at %e\n", tol, bestDistance);
    return bestNode;
}


double
Fibre::computeDistance(const oofem::FloatArray &a, const oofem::FloatArray &b)
{
    return sqrt(pow(a.at(1) - b.at(1), 2.) +
                pow(a.at(2) - b.at(2), 2.) +
                pow(a.at(3) - b.at(3), 2.));
}
