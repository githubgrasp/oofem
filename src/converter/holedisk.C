#include "holedisk.h"
#include "grid.h"
#include "line.h"
#include "vertex.h"
#include "intarray.h"

#include <cmath>
#include <map>


HoleDisk::HoleDisk(int n, double cx, double cy, double radius)
    : number(n), cx(cx), cy(cy), radius(radius)
{ }


bool HoleDisk::isStrictlyInside(double x, double y) const
{
    const double dx = x - cx, dy = y - cy;
    return std::sqrt(dx * dx + dy * dy) < radius;
}


bool HoleDisk::onRim(double x, double y, double tol) const
{
    const double dx = x - cx, dy = y - cy;
    return std::fabs(std::sqrt(dx * dx + dy * dy) - radius) < tol;
}


void HoleDisk::projectVoronoiCrossSection(Grid *grid, double tol) const
{
    // In 2D the cross-section of a Delaunay edge is a single dual Voronoi edge
    // (2 vertices); for a rim edge one endpoint is in the matrix and the other
    // inside the void — the latter is projected radially onto the circle and
    // kept (flag 2) as the boundary transport vertex.
    const double rimTol = 1.e3 * tol;

    oofem::IntArray ep, csVerts;
    oofem::FloatArray ca, cb, coords;

    auto onRimLocal = [&]( const oofem::FloatArray &c ) -> bool {
        const double d = std::sqrt(( c.at(1) - cx ) * ( c.at(1) - cx ) +
                                   ( c.at(2) - cy ) * ( c.at(2) - cy ) );
        return std::fabs(d - radius) < rimTol;
    };

    for ( int i = 1; i <= grid->giveNumberOfDelaunayLines(); ++i ) {
        grid->giveDelaunayLine(i)->giveLocalVertices(ep);
        if ( ep.giveSize() != 2 ) continue;
        grid->giveDelaunayVertex(ep.at(1) )->giveCoordinates(ca);
        grid->giveDelaunayVertex(ep.at(2) )->giveCoordinates(cb);
        if ( !onRimLocal(ca) || !onRimLocal(cb) ) continue;   // a rim mechanical edge of this hole

        grid->giveDelaunayLine(i)->giveCrossSectionVertices(csVerts);
        for ( int k = 1; k <= csVerts.giveSize(); ++k ) {
            if ( csVerts.at(k) == 0 ) continue;
            Vertex *vv = grid->giveVoronoiVertex(csVerts.at(k) );
            vv->giveCoordinates(coords);
            const double d = std::sqrt(( coords.at(1) - cx ) * ( coords.at(1) - cx ) +
                                       ( coords.at(2) - cy ) * ( coords.at(2) - cy ) );
            if ( d < radius - tol && d > 0. ) {            // strictly inside the void → project onto rim
                coords.at(1) = cx + radius / d * ( coords.at(1) - cx );
                coords.at(2) = cy + radius / d * ( coords.at(2) - cy );
                vv->setCoordinates(coords);
                vv->setOutsideFlag(2);
            }
        }
    }
}


bool HoleDisk::boundaryNodeCoord(const oofem::FloatArray &dA, const oofem::FloatArray &dB,
                                 double tol, oofem::FloatArray &m) const
{
    const double da = std::sqrt(( dA.at(1) - cx ) * ( dA.at(1) - cx ) + ( dA.at(2) - cy ) * ( dA.at(2) - cy ) );
    const double db = std::sqrt(( dB.at(1) - cx ) * ( dB.at(1) - cx ) + ( dB.at(2) - cy ) * ( dB.at(2) - cy ) );
    if ( std::fabs(da - radius) >= tol || std::fabs(db - radius) >= tol ) {
        return false;
    }
    m.resize(3);
    m.at(1) = 0.5 * ( dA.at(1) + dB.at(1) );
    m.at(2) = 0.5 * ( dA.at(2) + dB.at(2) );
    m.at(3) = 0.;
    const double ex = m.at(1) - cx, ey = m.at(2) - cy;
    const double d = std::sqrt(ex * ex + ey * ey);
    if ( d > 0. ) {
        m.at(1) = cx + radius / d * ex;
        m.at(2) = cy + radius / d * ey;
    }
    return true;
}


void HoleDisk::computeRimCoupling(Grid *grid, std::vector< RimCouplingEntry > &entries) const
{
    entries.clear();

    const double rimTol = 1.e3 * grid->giveTol();

    auto rad = [&]( const oofem::FloatArray &c ) -> double {
        return std::sqrt(( c.at(1) - cx ) * ( c.at(1) - cx ) + ( c.at(2) - cy ) * ( c.at(2) - cy ) );
    };

    // Accumulate the tributary length per rim vertex over its incident rim
    // Delaunay edges (the boundary mechanical elements). No Delaunay edge
    // crosses the empty hole, so the only edges with both endpoints on the
    // circle are the boundary-polygon edges — each rim node sees its two
    // neighbours and gets half of each edge.
    std::map< int, int > vertexToEntry;
    oofem::IntArray ep;
    oofem::FloatArray ca, cb;

    for ( int i = 1; i <= grid->giveNumberOfDelaunayLines(); ++i ) {
        if ( grid->giveDelaunayLine(i)->giveOutsideFlag() == 1 ) continue;
        grid->giveDelaunayLine(i)->giveLocalVertices(ep);
        if ( ep.giveSize() != 2 ) continue;
        grid->giveDelaunayVertex(ep.at(1) )->giveCoordinates(ca);
        grid->giveDelaunayVertex(ep.at(2) )->giveCoordinates(cb);
        if ( std::fabs(rad(ca) - radius) >= rimTol || std::fabs(rad(cb) - radius) >= rimTol ) continue;

        const double dx = cb.at(1) - ca.at(1);
        const double dy = cb.at(2) - ca.at(2);
        const double halfLen = 0.5 * std::sqrt(dx * dx + dy * dy);

        for ( int k = 1; k <= 2; ++k ) {
            const int gv = ep.at(k);
            auto it = vertexToEntry.find(gv);
            int idx;
            if ( it == vertexToEntry.end() ) {
                oofem::FloatArray cc;
                grid->giveDelaunayVertex(gv)->giveCoordinates(cc);
                const double d = rad(cc);
                RimCouplingEntry e;
                e.delaunayVertex = gv;
                e.dirX = ( cc.at(1) - cx ) / d;   // outward radial unit vector
                e.dirY = ( cc.at(2) - cy ) / d;
                e.tributary = 0.;
                idx = (int) entries.size();
                entries.push_back(e);
                vertexToEntry[ gv ] = idx;
            } else {
                idx = it->second;
            }
            entries[ idx ].tributary += halfLen;
            entries[ idx ].neighbourVertices.push_back(ep.at(3 - k) );  // other endpoint of this rim edge
        }
    }
}
