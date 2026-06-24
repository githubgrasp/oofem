#ifndef converter_disk_h
#define converter_disk_h

#include "grid.h"
#include "region.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <vector>
#endif


/// 2D circular domain region — the converter counterpart of the generator's
/// `#@disk` and the 2D analog of the 3D `Cylinder`. `#@disk <id> centre 2 cx cy
/// radius r`. Vertices outside the circle are marked outside the domain;
/// `modifyVoronoiCrossSection` projects the Voronoi cross-section vertices of
/// boundary-crossing mechanical elements radially onto the circle, so transport
/// nodes lie on the boundary and no transport element runs tangentially along
/// it (the requirement of the coupled transport-mechanical lattice).
class Disk : public Region
{
protected:
    /// Circle centre `[cx cy]`.
    oofem::FloatArray centre;
    /// Outer circle radius.
    double radius = 0.;
    /// Inner (hole) radius; 0 = full disc, > 0 = annulus.
    double innerRadius = 0.;

public:
    Disk(int n, Grid *aGrid);
    ~Disk() override = default;

    /// Set the circle centre `[cx cy]`.
    void setCentre(const oofem::FloatArray &c) { centre = c; }
    /// Returns the circle centre `[cx cy]`.
    const oofem::FloatArray &giveCentre() const { return centre; }
    /// Set the outer circle radius.
    void setRadius(double r) { radius = r; }
    /// Returns the outer circle radius.
    double giveRadius() const { return radius; }
    /// Set the inner (hole) radius; 0 = full disc.
    void setInnerRadius(double r) { innerRadius = r; }
    /// Returns the inner (hole) radius (0 = full disc).
    double giveInnerRadius() const { return innerRadius; }

    /// Fill `boundaries` with the disc bounding box `[xmin xmax ymin ymax 0 0]`.
    void defineBoundaries(oofem::FloatArray &boundaries) override;

    /// Radial analog of `Rect::findOutsiders`: sets `outsideFlag` on every
    /// Delaunay / Voronoi vertex and line by distance to the centre
    /// (0=inside, 1=outside, 2=on-circle for vertices; lines via endpoints).
    /// Far-outside Voronoi vertices are pulled in to 2*radius (as `Cylinder`).
    void findOutsiders(oofem::FloatArray &boundaries) override;

    /// Project the Voronoi cross-section vertices of a boundary-crossing
    /// mechanical element radially onto the circle (2D analog of
    /// `Cylinder::modifyVoronoiCrossSection`).
    int modifyVoronoiCrossSection(int elementNumber) override;

    /// True iff `(x, y)` lies inside the circle (with `tol` margin).
    bool contains(double x, double y, double tol) const;

    /// One inner-rim mechanical node receiving the hydro-mechanical coupling load.
    struct RimCouplingEntry {
        int delaunayVertex = 0;      ///< global Delaunay-vertex id (mechanical node)
        double dirX = 0., dirY = 0.; ///< outward radial unit direction at the node
        double tributary = 0.;       ///< tributary boundary length (sum of incident half-edges)
    };

    /// Identify the mechanical nodes on this disk's inner circle (the annulus
    /// hole rim) and fill `entries` with each node's outward radial direction
    /// and tributary boundary length (half of every incident inner-rim Delaunay
    /// edge). Pure geometry on the region's own circle — the grid turns the
    /// result into coupling boundary-condition records. Requires `innerRadius > 0`.
    void computeInnerRimCoupling(std::vector< RimCouplingEntry > &entries);

    const char *giveClassName() const override { return "Disk"; }
};


#endif // converter_disk_h
