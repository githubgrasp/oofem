#ifndef holedisk_h
#define holedisk_h

#include "inclusion.h"
#include "grid.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <sstream>
#endif


/// 2D circular hole (void) inclusion. Reads
/// `#@holedisk <id> centre 2 cx cy radius r [refine <r>]`.
/// Unlike `InterfaceDisk` (a solid inclusion with an ITZ sister ring), a hole
/// is meshed with mechanical nodes placed ON the hole/matrix boundary (a single
/// rim ring) plus sacrificial points filling the hole interior. The sacrificial
/// points only control the Delaunay triangulation so the rim edges come out
/// clean (no chords spanning the void); the converter then drops everything
/// strictly inside the hole and shifts the dual cross-section vertices onto the
/// circle, leaving mechanical AND transport nodes on the boundary.
/// Matrix fill points that would land inside the hole are rejected by
/// `Grid::addVertex` (see `Grid::insideAnyHole`).
class HoleDisk : public Inclusion
{
protected:
    /// `[cx cy]`.
    oofem::FloatArray centre;

public:
    /**
     * Constructor. Creates a 2D circular hole.
     * @param n inclusion number in the grid
     * @param aGrid grid to which the hole belongs
     */
    HoleDisk(int n, Grid *aGrid);

    /// Destructor.
    ~HoleDisk();

    /// Returns the hole radius.
    double giveRadius() { return this->radius; }

    /// Returns the hole centre `[cx cy]`.
    const oofem::FloatArray &giveCentre() const { return this->centre; }

    /// True iff `(coord)` lies strictly inside the hole (more than `tol`
    /// inside the circle). Boundary points within `tol` of the rim are kept.
    bool containsStrictly(const oofem::FloatArray &coord, double tol) const;

    /// Place a ring of mechanical nodes on the circle plus sacrificial points
    /// filling the interior (the latter are dropped by the converter).
    int generatePoints() override;

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "HoleDisk"; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@holedisk <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);
};

#endif // holedisk_h
