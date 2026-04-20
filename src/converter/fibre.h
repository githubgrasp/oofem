#ifndef FIBRE_H_INCLUDED
#define FIBRE_H_INCLUDED

#include <vector>
#include "grid.h"
#include "vertex.h"
#include "gridcomponent.h"
#include "floatarray.h"
#include "intarray.h"



#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

// Notes:
//   * Only straight fibres are supported.
//   * All fibres must lie entirely within the 27-cell periodic structure.

class Fibre : public GridComponent
{
protected:
    int    number;
    double diameter;
    double length;
    double tolerance;                       // guard against floating-point instabilities in the discretisation walk

    oofem::FloatArray endpoints;            // 6 doubles: (xP, yP, zP, xQ, yQ, zQ)
    oofem::FloatArray coordP, coordQ;       // endpoint coordinates as 3-vectors
    oofem::FloatArray directionVector;      // unit vector from P to Q

    // Vertices created in the grid for the fibre endpoints.
    Vertex *pointP;
    Vertex *pointQ;

    // Discretisation outputs: indices into the grid's lists.
    std::vector< int >    reinforcementPoints;
    std::vector< int >    delaunayVertices;       // matrix Delaunay vertex paired with each reinforcement node
    std::vector< int >    intersectionPoints;     // points where the fibre crosses Voronoi facets
    std::vector< double > endLengths;             // distance from each reinforcement node to the nearer fibre endpoint

    // Discretisation helpers.
    void findDelaunayVertex();
    void findIntersect();
    void placeReinforcementNode();

public:
    Fibre(int n, Grid *aGrid);
    ~Fibre() override = default;

    Fibre *ofType();

    /// Initialise from already-parsed endpoints + diameter (qhull-template path).
    void initializeFromCoords(const oofem::FloatArray &endpoints, double diam);

    /// Walk the fibre, place reinforcement nodes, and record the matched
    /// matrix Delaunay vertices and Voronoi-cell intersection points.
    void discretize();

    int               giveNumber() { return number; }
    double            giveDiameter() { return diameter; }
    oofem::FloatArray giveDirectionVector() { return directionVector; }

    int    giveNumberOfReinforcementNodes() { return reinforcementPoints.size(); }
    int    giveNumberOfDelaunayNodes() { return delaunayVertices.size(); }
    int    giveNumberOfIntersectionPoints() { return intersectionPoints.size(); }
    int    giveNumberOfEndLengths() { return endLengths.size(); }
    int    giveReinforcementNodeNumber(int i) { return reinforcementPoints.at(i - 1); }
    int    giveDelaunayNodeNumber(int i) { return delaunayVertices.at(i - 1); }
    int    giveIntersectionPointNumber(int i) { return intersectionPoints.at(i - 1); }
    double giveEndLength(int i) { return endLengths.at(i - 1); }

    /// Find the reinforcement node nearest to `coord`, falling back to the
    /// closest one when no node lies within `tol`. Returns a pointer to an
    /// existing grid vertex (no allocation).
    static Vertex *findNearestReinforcementNode(const oofem::FloatArray &coord, Grid *grid, double tol = 1.0e-14);

    static double computeDistance(const oofem::FloatArray &a, const oofem::FloatArray &b);
};

#endif // FIBRE_H_INCLUDED
