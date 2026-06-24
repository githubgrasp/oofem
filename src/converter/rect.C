#include "rect.h"
#include "vertex.h"
#include "line.h"
#include "octreegridlocalizer.h"
#include "gridlocalizer.h"

#include <cmath>
#include <list>
#include <vector>


Rect::Rect(int n, Grid *aGrid) : Region(n, aGrid)
{
    this->number = n;
}


void Rect::defineBoundaries(oofem::FloatArray &boundaries)
{
    boundaries.resize(6);
    boundaries.at(1) = box.at(1);   // xmin
    boundaries.at(2) = box.at(3);   // xmax
    boundaries.at(3) = box.at(2);   // ymin
    boundaries.at(4) = box.at(4);   // ymax
    boundaries.at(5) = 0.;          // zmin
    boundaries.at(6) = 0.;          // zmax
}


bool Rect::contains(double x, double y, double tol) const
{
    return x > box.at(1) - tol && x < box.at(3) + tol &&
           y > box.at(2) - tol && y < box.at(4) + tol;
}


void Rect::findOutsiders(oofem::FloatArray &boundaries)
{
    const double tol = grid->giveTol();
    const double xmin = boundaries.at(1), xmax = boundaries.at(2);
    const double ymin = boundaries.at(3), ymax = boundaries.at(4);
    const double dx = xmax - xmin, dy = ymax - ymin;

    oofem::IntArray periodicityFlag;
    grid->givePeriodicityFlag(periodicityFlag);
    const bool perX = ( periodicityFlag.giveSize() >= 1 && periodicityFlag.at(1) == 1 );
    const bool perY = ( periodicityFlag.giveSize() >= 2 && periodicityFlag.at(2) == 1 );
    const bool periodic = perX || perY;

    auto classify = [&]( double x, double y ) -> int {
        // A vertex inside a `delete` (hole) inclusion is treated as outside the
        // domain — carving a void through the same outside-node machinery used
        // for boundary mirror / periodic image nodes.
        if ( grid->isInsideDeleteHole(x, y) ) return 1;
        const bool sx = ( x + tol < xmin ) || ( x - tol > xmax );
        const bool sy = ( y + tol < ymin ) || ( y - tol > ymax );
        if ( sx || sy ) return 1;
        const bool bx = std::fabs(x - xmin) < tol || std::fabs(x - xmax) < tol;
        const bool by = std::fabs(y - ymin) < tol || std::fabs(y - ymax) < tol;
        if ( bx || by ) return 2;
        return 0;
    };

    auto axisShift = [&]( double v, double vmin, double vmax ) -> int {
        if ( v + tol < vmin ) return -1;
        if ( v - tol > vmax ) return +1;
        return 0;
    };

    // Compass-direction location code (1..8) used by Lattice2dBoundary.
    auto encodeLocation = []( int sx, int sy ) -> int {
        if ( sx == +1 && sy ==  0 ) return 1;
        if ( sx == +1 && sy == +1 ) return 2;
        if ( sx ==  0 && sy == +1 ) return 3;
        if ( sx == -1 && sy == +1 ) return 4;
        if ( sx == -1 && sy ==  0 ) return 5;
        if ( sx == -1 && sy == -1 ) return 6;
        if ( sx ==  0 && sy == -1 ) return 7;
        if ( sx == +1 && sy == -1 ) return 8;
        return 0;
    };

    auto findPartner = [&]( double x, double y, int nodeType ) -> int {
        oofem::FloatArray probe(3);
        probe.at(1) = x; probe.at(2) = y; probe.at(3) = 0.;
        std::list< int > nodeSet;
        if ( nodeType == 0 && grid->delaunayLocalizer ) {
            grid->delaunayLocalizer->giveAllNodesWithinBox(nodeSet, probe, tol, 0);
        } else if ( nodeType == 1 && grid->voronoiLocalizer ) {
            grid->voronoiLocalizer->giveAllNodesWithinBox(nodeSet, probe, tol, 1);
        }
        return nodeSet.empty() ? 0 : *nodeSet.begin();
    };

    auto markVertex = [&]( int id, int nodeType, oofem::FloatArray &c ) {
        Vertex *v = ( nodeType == 0 ) ? grid->giveDelaunayVertex(id) : grid->giveVoronoiVertex(id);
        v->giveCoordinates(c);
        const int flag = classify(c.at(1), c.at(2));
        v->setOutsideFlag(flag);
        if ( flag == 1 && periodic ) {
            const int sx = perX ? axisShift(c.at(1), xmin, xmax) : 0;
            const int sy = perY ? axisShift(c.at(2), ymin, ymax) : 0;
            if ( sx == 0 && sy == 0 ) return;   // truly outside on a non-periodic axis
            const double px = c.at(1) - sx * dx;
            const double py = c.at(2) - sy * dy;
            const int partner = findPartner(px, py, nodeType);
            if ( partner != 0 ) {
                v->setPeriodicNode(partner);
                v->setLocation(encodeLocation(sx, sy));
            }
        }
    };

    oofem::FloatArray c(3);
    for ( int i = 1; i <= grid->giveNumberOfDelaunayVertices(); ++i ) {
        markVertex(i, 0, c);
    }
    for ( int i = 1; i <= grid->giveNumberOfVoronoiVertices(); ++i ) {
        markVertex(i, 1, c);
    }

    auto lineFlag = []( int f1, int f2 ) -> int {
        // 0 inside, 1 fully outside, 2 crossing, 3 on boundary
        if ( f1 == 0 && f2 == 0 ) return 0;
        if ( f1 == 1 && f2 == 1 ) return 1;
        if ( f1 == 2 && f2 == 2 ) return 3;
        if ( f1 == 1 || f2 == 1 ) return 2;
        return 0;   // mixed inside / on-boundary stays inside
    };

    oofem::IntArray ep;
    for ( int i = 1; i <= grid->giveNumberOfDelaunayLines(); ++i ) {
        grid->giveDelaunayLine(i)->giveLocalVertices(ep);
        if ( ep.giveSize() != 2 ) continue;
        const int f1 = grid->giveDelaunayVertex(ep.at(1))->giveOutsideFlag();
        const int f2 = grid->giveDelaunayVertex(ep.at(2))->giveOutsideFlag();
        grid->giveDelaunayLine(i)->setOutsideFlag(lineFlag(f1, f2));
    }
    for ( int i = 1; i <= grid->giveNumberOfVoronoiLines(); ++i ) {
        grid->giveVoronoiLine(i)->giveLocalVertices(ep);
        if ( ep.giveSize() != 2 ) continue;
        if ( ep.at(1) == 0 || ep.at(2) == 0 ) {
            grid->giveVoronoiLine(i)->setOutsideFlag(1);
            continue;
        }
        const int f1 = grid->giveVoronoiVertex(ep.at(1))->giveOutsideFlag();
        const int f2 = grid->giveVoronoiVertex(ep.at(2))->giveOutsideFlag();
        grid->giveVoronoiLine(i)->setOutsideFlag(lineFlag(f1, f2));
    }

    // The region owns its Voronoi cross-section adjustment (mirrors
    // Cylinder/Disk findOutsiders): clamp boundary-crossing Voronoi vertices
    // onto the rectangle faces. Notch voids are handled separately.
    this->projectVoronoiToBoundary();
}


void Rect::projectVoronoiToBoundary()
{
    // Clamp the outside endpoint of every crossing Voronoi edge (one endpoint
    // inside the rectangle, one outside) onto the nearest rectangle face. Only
    // crossing-edge endpoints are moved — far-outside vertices are left where
    // qhull put them (they are never emitted), avoiding coincident corner nodes.
    const double xmin = box.at(1), ymin = box.at(2);
    const double xmax = box.at(3), ymax = box.at(4);
    const double tol = grid->giveTol();

    auto outsideRect = [&]( const oofem::FloatArray &c ) {
        return c.at(1) < xmin - tol || c.at(1) > xmax + tol ||
               c.at(2) < ymin - tol || c.at(2) > ymax + tol;
    };
    auto clampToRect = [&]( oofem::FloatArray &c ) {
        if ( c.at(1) < xmin ) c.at(1) = xmin;
        if ( c.at(1) > xmax ) c.at(1) = xmax;
        if ( c.at(2) < ymin ) c.at(2) = ymin;
        if ( c.at(2) > ymax ) c.at(2) = ymax;
    };

    const int nVorV = grid->giveNumberOfVoronoiVertices();
    const int nVorL = grid->giveNumberOfVoronoiLines();

    std::vector< int > vertOutside( nVorV + 1, 0 );
    oofem::FloatArray c(3);
    for ( int i = 1; i <= nVorV; ++i ) {
        grid->giveVoronoiVertex(i)->giveCoordinates(c);
        if ( outsideRect(c) ) vertOutside[i] = 1;
    }

    std::vector< char > toProject( nVorV + 1, 0 );
    oofem::IntArray ep;
    for ( int i = 1; i <= nVorL; ++i ) {
        grid->giveVoronoiLine(i)->giveLocalVertices(ep);
        if ( ep.giveSize() != 2 ) continue;
        const int a = ep.at(1), b = ep.at(2);
        if ( a == 0 || b == 0 ) continue;   // qhull at-infinity marker
        if ( vertOutside[a] == 0 && vertOutside[b] == 1 ) toProject[b] = 1;
        else if ( vertOutside[a] == 1 && vertOutside[b] == 0 ) toProject[a] = 1;
    }

    for ( int i = 1; i <= nVorV; ++i ) {
        if ( !toProject[i] ) continue;
        grid->giveVoronoiVertex(i)->giveCoordinates(c);
        clampToRect(c);
        grid->giveVoronoiVertex(i)->setCoordinates(c);
    }
}
