#ifndef converter_disk_h
#define converter_disk_h

#include "grid.h"
#include "region.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif


/// 2D circular domain region — the converter counterpart of the generator's
/// `#@disk` and the 2D analog of the 3D `Cylinder`. `#@disk <id> centre 2 cx cy
/// radius r`. Vertices outside the circle are marked outside the domain;
/// `modifyVoronoiCrossSection` projects the Voronoi cross-section vertices of
/// boundary-crossing mechanical elements radially onto the circle, so transport
/// nodes lie on the boundary and no transport element runs tangentially along
/// it (the requirement of the coupled transport-mechanical lattice). Holes
/// (annulus, multiple voids) are separate `#@holedisk` inclusions, not a
/// property of the disc.
class Disk : public Region
{
protected:
    /// Circle centre `[cx cy]`.
    oofem::FloatArray centre;
    /// Outer circle radius.
    double radius = 0.;

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

    const char *giveClassName() const override { return "Disk"; }
};


#endif // converter_disk_h
