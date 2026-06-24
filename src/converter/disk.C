#include "disk.h"
#include "vertex.h"
#include "line.h"
#include "convertererror.h"

#include <cmath>
#include <map>


Disk::Disk(int n, Grid *aGrid) : Region(n, aGrid)
{
    this->number = n;
}


void Disk::defineBoundaries(oofem::FloatArray &boundaries)
{
    boundaries.resize(6);
    boundaries.at(1) = centre.at(1) - radius;   // xmin
    boundaries.at(2) = centre.at(1) + radius;   // xmax
    boundaries.at(3) = centre.at(2) - radius;   // ymin
    boundaries.at(4) = centre.at(2) + radius;   // ymax
    boundaries.at(5) = 0.;                       // zmin
    boundaries.at(6) = 0.;                       // zmax
}


bool Disk::contains(double x, double y, double tol) const
{
    const double dx = x - centre.at(1), dy = y - centre.at(2);
    const double d = std::sqrt(dx * dx + dy * dy);
    return d < radius + tol && ( innerRadius <= 0. || d > innerRadius - tol );
}


void Disk::findOutsiders(oofem::FloatArray &boundaries)
{
    const double tol = grid->giveTol();
    const double cx = centre.at(1), cy = centre.at(2);

    auto classify = [&]( double x, double y ) -> int {
        const double d = std::sqrt(( x - cx ) * ( x - cx ) + ( y - cy ) * ( y - cy ) );
        if ( grid->isInsideDeleteHole(x, y) ) return 1;                            // inside a #@holedisk void
        if ( grid->onDeleteHoleRim(x, y, tol) ) return 2;                          // on a #@holedisk rim
        if ( d - tol > radius ) return 1;                                          // beyond outer circle
        if ( innerRadius > 0. && d + tol < innerRadius ) return 1;                 // inside the hole
        if ( std::fabs(d - radius) < tol ) return 2;                               // on outer circle
        if ( innerRadius > 0. && std::fabs(d - innerRadius) < tol ) return 2;      // on inner circle
        return 0;                                                                   // inside the annulus
    };

    oofem::FloatArray c(3);
    for ( int i = 1; i <= grid->giveNumberOfDelaunayVertices(); ++i ) {
        grid->giveDelaunayVertex(i)->giveCoordinates(c);
        grid->giveDelaunayVertex(i)->setOutsideFlag( classify(c.at(1), c.at(2)) );
    }
    for ( int i = 1; i <= grid->giveNumberOfVoronoiVertices(); ++i ) {
        grid->giveVoronoiVertex(i)->giveCoordinates(c);
        const int flag = classify(c.at(1), c.at(2));
        if ( flag == 1 ) {
            // Pull far-outside Voronoi vertices in to 2*radius so the dual cells
            // stay bounded (mirror of Cylinder::findOutsiders).
            const double d = std::sqrt(( c.at(1) - cx ) * ( c.at(1) - cx ) + ( c.at(2) - cy ) * ( c.at(2) - cy ) );
            if ( d > 2. * radius ) {
                c.at(1) = cx + 2. * radius / d * ( c.at(1) - cx );
                c.at(2) = cy + 2. * radius / d * ( c.at(2) - cy );
                grid->giveVoronoiVertex(i)->setCoordinates(c);
            }
        }
        grid->giveVoronoiVertex(i)->setOutsideFlag(flag);
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
        grid->giveDelaunayLine(i)->setOutsideFlag( lineFlag(f1, f2) );
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
        grid->giveVoronoiLine(i)->setOutsideFlag( lineFlag(f1, f2) );
    }

    // The region owns the cross-section adjustment (mirrors
    // Cylinder::findOutsiders): for each boundary Delaunay line, project its
    // Voronoi cross-section vertices radially onto the circle.
    for ( int i = 1; i <= grid->giveNumberOfDelaunayLines(); ++i ) {
        const int lf = grid->giveDelaunayLine(i)->giveOutsideFlag();
        if ( lf == 2 || lf == 3 ) {
            this->modifyVoronoiCrossSection(i);
        }
    }
}


void Disk::computeInnerRimCoupling(std::vector< RimCouplingEntry > &entries)
{
    entries.clear();

    if ( innerRadius <= 0. ) {
        converter::error("Disk::computeInnerRimCoupling requires innerradius > 0 (an annulus)");
    }
    const double cx = centre.at(1), cy = centre.at(2);
    const double rimTol = 1.e3 * this->grid->giveTol();

    auto radius = [&]( const oofem::FloatArray &c ) -> double {
        return std::sqrt(( c.at(1) - cx ) * ( c.at(1) - cx ) + ( c.at(2) - cy ) * ( c.at(2) - cy ) );
    };

    // Accumulate the tributary length per inner-rim vertex over its incident
    // inner-rim Delaunay edges (the boundary mechanical elements). No Delaunay
    // edge crosses the empty hole, so the only edges with both endpoints on the
    // inner circle are the boundary-polygon edges — each rim node sees its two
    // neighbours and gets half of each edge.
    std::map< int, int > vertexToEntry;
    oofem::IntArray ep;
    oofem::FloatArray ca, cb;

    for ( int i = 1; i <= this->grid->giveNumberOfDelaunayLines(); ++i ) {
        if ( this->grid->giveDelaunayLine(i)->giveOutsideFlag() == 1 ) continue;
        this->grid->giveDelaunayLine(i)->giveLocalVertices(ep);
        if ( ep.giveSize() != 2 ) continue;
        this->grid->giveDelaunayVertex(ep.at(1) )->giveCoordinates(ca);
        this->grid->giveDelaunayVertex(ep.at(2) )->giveCoordinates(cb);
        if ( std::fabs(radius(ca) - innerRadius) >= rimTol ||
             std::fabs(radius(cb) - innerRadius) >= rimTol ) continue;

        const double dx = cb.at(1) - ca.at(1);
        const double dy = cb.at(2) - ca.at(2);
        const double halfLen = 0.5 * std::sqrt(dx * dx + dy * dy);

        for ( int k = 1; k <= 2; ++k ) {
            const int gv = ep.at(k);
            auto it = vertexToEntry.find(gv);
            int idx;
            if ( it == vertexToEntry.end() ) {
                oofem::FloatArray cc;
                this->grid->giveDelaunayVertex(gv)->giveCoordinates(cc);
                const double d = radius(cc);
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
        }
    }
}


int Disk::modifyVoronoiCrossSection(int elementNumber)
{
    // 2D analog of Cylinder::modifyVoronoiCrossSection: for the cross-section
    // (dual Voronoi edges) of a boundary-crossing mechanical element, drop the
    // fully-outside dual edges and shift the dual vertices of crossing edges
    // radially onto the circle, so transport nodes lie on the boundary.
    oofem::IntArray crossSectionElements;
    oofem::IntArray crossSectionVertices;
    oofem::IntArray nodes(2);
    oofem::FloatArray coords(3);

    const double newTol = 1.e3 * this->grid->giveTol();
    const double cx = centre.at(1), cy = centre.at(2);

    this->grid->giveDelaunayLine(elementNumber)->giveCrossSectionElements(crossSectionElements);
    this->grid->giveDelaunayLine(elementNumber)->giveCrossSectionVertices(crossSectionVertices);

    const int elementSize = crossSectionElements.giveSize();
    const int vertexSize = crossSectionVertices.giveSize();
    int elementCounter = 0;

    for ( int i = 0; i < elementSize; i++ ) {
        this->grid->giveVoronoiLine(crossSectionElements.at(i + 1) )->giveLocalVertices(nodes);
        const int lflag = this->grid->giveVoronoiLine(crossSectionElements.at(i + 1) )->giveOutsideFlag();
        if ( lflag == 1 || lflag == 3 ) {            // fully outside or on the circle
            crossSectionElements.at(i + 1) = 0;
            elementCounter++;
        } else if ( lflag == 2 ) {                   // crossing the boundary
            for ( int k = 0; k < 2; k++ ) {
                if ( this->grid->giveVoronoiVertex(nodes.at(k + 1) )->giveOutsideFlag() == 1 ) {
                    this->grid->giveVoronoiVertex(nodes.at(k + 1) )->giveCoordinates(coords);
                    const double d = std::sqrt(( coords.at(1) - cx ) * ( coords.at(1) - cx ) +
                                               ( coords.at(2) - cy ) * ( coords.at(2) - cy ) );
                    // Project onto whichever circle the vertex is beyond: the
                    // outer circle if past it, the inner (hole) circle if in the hole.
                    double targetR = 0.;
                    if ( d - newTol > radius ) {
                        targetR = this->radius;
                    } else if ( innerRadius > 0. && d + newTol < innerRadius ) {
                        targetR = innerRadius;
                    }
                    if ( targetR > 0. && d > 0. ) {
                        coords.at(1) = cx + targetR / d * ( coords.at(1) - cx );
                        coords.at(2) = cy + targetR / d * ( coords.at(2) - cy );
                    }
                    this->grid->giveVoronoiVertex(nodes.at(k + 1) )->setCoordinates(coords);
                    this->grid->giveVoronoiVertex(nodes.at(k + 1) )->setOutsideFlag(2);
                }
            }
            this->grid->giveVoronoiLine(crossSectionElements.at(i + 1) )->setOutsideFlag(0);
        }
    }

    // Compact the cross-section element list (drop the zeroed entries).
    oofem::IntArray modifiedCrossSectionElements(elementSize - elementCounter);
    int help = 0;
    for ( int i = 0; i < elementSize; i++ ) {
        if ( crossSectionElements.at(i + 1) == 0 ) {
            help++;
        } else {
            modifiedCrossSectionElements.at(i - help + 1) = crossSectionElements.at(i + 1);
        }
    }
    this->grid->giveDelaunayLine(elementNumber)->setCrossSectionElements(modifiedCrossSectionElements);

    // Compact the cross-section vertex list (drop fully-outside dual vertices).
    int nodeCounter = 0;
    for ( int i = 0; i < vertexSize; i++ ) {
        if ( this->grid->giveVoronoiVertex(crossSectionVertices.at(i + 1) )->giveOutsideFlag() == 1 ) {
            crossSectionVertices.at(i + 1) = 0;
            nodeCounter++;
        }
    }
    oofem::IntArray modifiedCrossSectionVertices(vertexSize - nodeCounter);
    help = 0;
    for ( int i = 0; i < vertexSize; i++ ) {
        if ( crossSectionVertices.at(i + 1) == 0 ) {
            help++;
        } else {
            modifiedCrossSectionVertices.at(i - help + 1) = crossSectionVertices.at(i + 1);
        }
    }
    this->grid->giveDelaunayLine(elementNumber)->setCrossSectionVertices(modifiedCrossSectionVertices);

    // In 2D the cross-section of a Delaunay edge is a single dual Voronoi edge
    // (2 vertices); it is valid as long as both survive. (The 3D Cylinder uses
    // < 3 because there the cross-section is a polygon.)
    int flag = 1;
    if ( modifiedCrossSectionVertices.giveSize() < 2 ) {
        flag = 0;
        this->grid->giveDelaunayLine(elementNumber)->setOutsideFlag(1);
    }
    return flag;
}
