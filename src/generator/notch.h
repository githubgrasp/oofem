#ifndef notch_h
#define notch_h

#include "gridcomponent.h"
#include "grid.h"

#include "floatarray.h"

#ifndef __MAKEDEPEND
 #include <sstream>
 #include <stdio.h>
#endif


/// An axis-aligned rectangular/cuboid notch declared via `#@notch`. The grid
/// suppresses point placement strictly inside the box (boundary points are
/// retained, so the dual mesh can partition cleanly along the notch surface).
/// The companion converter-side `#@notch ... delete` directive then drops
/// elements whose midpoint falls inside the same box. Optional `edgerefine`
/// and `surfacerefine` factors are stored for the (forthcoming) explicit
/// face / edge seeding pass.
class Notch : public GridComponent
{
protected:
    /// `[xmin ymin zmin xmax ymax zmax]`. For 2D grids the z entries are
    /// ignored.
    oofem::FloatArray box;

    /// Spacing factor for points seeded on the notch box edges
    /// (relative to `grid->diameter`). Default 1.
    double edgeRefine = 1.;

    /// Spacing factor for points seeded on the notch box faces
    /// (relative to `grid->diameter`). Default 1.
    double surfaceRefine = 1.;

public:
    /**
     * Constructor. Creates a notch entry belonging to `aGrid`.
     * @param n notch number in the grid
     * @param aGrid grid to which the notch belongs
     */
    Notch(int n, Grid *aGrid);

    /// Destructor.
    ~Notch();

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@notch <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);

    /// True iff `coord` lies strictly inside the notch box (more than
    /// `tol` away from every face). Boundary points return false.
    bool containsStrictly(const oofem::FloatArray &coord, double tol) const;

    /// Returns the receiver's box as `[xmin ymin zmin xmax ymax zmax]`.
    const oofem::FloatArray &giveBox() const { return box; }

    /// Returns the edge spacing factor (relative to `grid->diameter`).
    double giveEdgeRefine() const { return edgeRefine; }

    /// Returns the face spacing factor (relative to `grid->diameter`).
    double giveSurfaceRefine() const { return surfaceRefine; }

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Notch"; }
};

#endif // notch_h
