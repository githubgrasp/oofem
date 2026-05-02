#include "notch.h"
#include "generatorerror.h"

#include <algorithm>
#include <cmath>


Notch::Notch(int n, Grid *aGrid) : GridComponent(n, aGrid) {}

Notch::~Notch() {}


void Notch::initializeFromTokens(std::istringstream &iss)
{
    bool gotBox = false;
    std::string tok;
    while ( iss >> tok ) {
        if ( tok == "box" ) {
            int n;
            iss >> n;
            box.resize(n);
            for ( int i = 1; i <= n; ++i ) {
                iss >> box.at(i);
            }
            gotBox = true;
        } else if ( tok == "edgerefine" ) {
            iss >> edgeRefine;
        } else if ( tok == "surfacerefine" ) {
            iss >> surfaceRefine;
        } else {
            generator::errorf("Notch::initializeFromTokens: unknown keyword '%s'", tok.c_str() );
        }
    }
    if ( !gotBox ) {
        generator::error("Notch::initializeFromTokens: missing 'box' keyword");
    }
    if ( box.giveSize() != 4 && box.giveSize() != 6 ) {
        generator::errorf("Notch::initializeFromTokens: 'box' must have 4 (2D) or 6 (3D) entries, got %d",
                          box.giveSize() );
    }
}


bool Notch::containsStrictly(const oofem::FloatArray &coord, double tol) const
{
    if ( box.giveSize() == 6 ) {
        return coord.at(1) > box.at(1) + tol && coord.at(1) < box.at(4) - tol &&
               coord.at(2) > box.at(2) + tol && coord.at(2) < box.at(5) - tol &&
               coord.at(3) > box.at(3) + tol && coord.at(3) < box.at(6) - tol;
    }
    // 2D box: xmin ymin xmax ymax
    return coord.at(1) > box.at(1) + tol && coord.at(1) < box.at(3) - tol &&
           coord.at(2) > box.at(2) + tol && coord.at(2) < box.at(4) - tol;
}


void Notch::generateBoundaryPoints()
{
    auto nIntervals = [](double length, double spacing) {
        if ( spacing <= 0. || length <= 0. ) return 1;
        return std::max(1, static_cast< int >( std::ceil(length / spacing) ));
    };

    const double diam = grid->diameter;
    const double sEdge = edgeRefine * diam;
    const double sFace = surfaceRefine * diam;

    if ( box.giveSize() == 6 ) {
        // 3D box: 8 corners + 12 edges + 6 faces.
        const double xmin = box.at(1), ymin = box.at(2), zmin = box.at(3);
        const double xmax = box.at(4), ymax = box.at(5), zmax = box.at(6);
        const double dx = xmax - xmin, dy = ymax - ymin, dz = zmax - zmin;
        oofem::FloatArray p(3);

        // Corners
        for ( int ix = 0; ix < 2; ++ix ) {
            p.at(1) = ix == 0 ? xmin : xmax;
            for ( int iy = 0; iy < 2; ++iy ) {
                p.at(2) = iy == 0 ? ymin : ymax;
                for ( int iz = 0; iz < 2; ++iz ) {
                    p.at(3) = iz == 0 ? zmin : zmax;
                    grid->addVertex(p);
                }
            }
        }

        // Edges parallel to x — 4 of them, one per (y,z) corner pair
        const int nx = nIntervals(dx, sEdge);
        for ( int iy = 0; iy < 2; ++iy ) {
            for ( int iz = 0; iz < 2; ++iz ) {
                p.at(2) = iy == 0 ? ymin : ymax;
                p.at(3) = iz == 0 ? zmin : zmax;
                for ( int k = 1; k < nx; ++k ) {
                    p.at(1) = xmin + ( double ) k * dx / nx;
                    grid->addVertex(p);
                }
            }
        }
        // Edges parallel to y
        const int ny = nIntervals(dy, sEdge);
        for ( int ix = 0; ix < 2; ++ix ) {
            for ( int iz = 0; iz < 2; ++iz ) {
                p.at(1) = ix == 0 ? xmin : xmax;
                p.at(3) = iz == 0 ? zmin : zmax;
                for ( int k = 1; k < ny; ++k ) {
                    p.at(2) = ymin + ( double ) k * dy / ny;
                    grid->addVertex(p);
                }
            }
        }
        // Edges parallel to z
        const int nz = nIntervals(dz, sEdge);
        for ( int ix = 0; ix < 2; ++ix ) {
            for ( int iy = 0; iy < 2; ++iy ) {
                p.at(1) = ix == 0 ? xmin : xmax;
                p.at(2) = iy == 0 ? ymin : ymax;
                for ( int k = 1; k < nz; ++k ) {
                    p.at(3) = zmin + ( double ) k * dz / nz;
                    grid->addVertex(p);
                }
            }
        }

        // Faces — interior grid (excluding edges, which are already done).
        const int nxF = nIntervals(dx, sFace);
        const int nyF = nIntervals(dy, sFace);
        const int nzF = nIntervals(dz, sFace);

        // Faces perpendicular to x (x = xmin and x = xmax)
        for ( int ix = 0; ix < 2; ++ix ) {
            p.at(1) = ix == 0 ? xmin : xmax;
            for ( int j = 1; j < nyF; ++j ) {
                p.at(2) = ymin + ( double ) j * dy / nyF;
                for ( int k = 1; k < nzF; ++k ) {
                    p.at(3) = zmin + ( double ) k * dz / nzF;
                    grid->addVertex(p);
                }
            }
        }
        // Faces perpendicular to y
        for ( int iy = 0; iy < 2; ++iy ) {
            p.at(2) = iy == 0 ? ymin : ymax;
            for ( int i = 1; i < nxF; ++i ) {
                p.at(1) = xmin + ( double ) i * dx / nxF;
                for ( int k = 1; k < nzF; ++k ) {
                    p.at(3) = zmin + ( double ) k * dz / nzF;
                    grid->addVertex(p);
                }
            }
        }
        // Faces perpendicular to z
        for ( int iz = 0; iz < 2; ++iz ) {
            p.at(3) = iz == 0 ? zmin : zmax;
            for ( int i = 1; i < nxF; ++i ) {
                p.at(1) = xmin + ( double ) i * dx / nxF;
                for ( int j = 1; j < nyF; ++j ) {
                    p.at(2) = ymin + ( double ) j * dy / nyF;
                    grid->addVertex(p);
                }
            }
        }
        return;
    }

    // 2D box: 4 corners + 4 edges. (Faces don't apply.)
    const double xmin = box.at(1), ymin = box.at(2);
    const double xmax = box.at(3), ymax = box.at(4);
    const double dx = xmax - xmin, dy = ymax - ymin;
    oofem::FloatArray p(2);

    // Corners
    for ( int ix = 0; ix < 2; ++ix ) {
        p.at(1) = ix == 0 ? xmin : xmax;
        for ( int iy = 0; iy < 2; ++iy ) {
            p.at(2) = iy == 0 ? ymin : ymax;
            grid->addVertex(p);
        }
    }

    // Edges parallel to x (top and bottom)
    const int nx = nIntervals(dx, sEdge);
    for ( int iy = 0; iy < 2; ++iy ) {
        p.at(2) = iy == 0 ? ymin : ymax;
        for ( int k = 1; k < nx; ++k ) {
            p.at(1) = xmin + ( double ) k * dx / nx;
            grid->addVertex(p);
        }
    }
    // Edges parallel to y (left and right)
    const int ny = nIntervals(dy, sEdge);
    for ( int ix = 0; ix < 2; ++ix ) {
        p.at(1) = ix == 0 ? xmin : xmax;
        for ( int k = 1; k < ny; ++k ) {
            p.at(2) = ymin + ( double ) k * dy / ny;
            grid->addVertex(p);
        }
    }
}
