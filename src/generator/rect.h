#ifndef rect_h
#define rect_h


#include "grid.h"
#include "region.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <sstream>
#endif


/// Axis-aligned 2D rectangular region — the 2D analog of `Prism`. Reads
/// `#@rect <id> box 4 xmin ymin xmax ymax [refine <r>] [edgerefine <re>]
/// [regionrefine <rr>]`. Coordinates are stored internally as 3D with
/// `z = 0` so the existing localiser / vertex / output infrastructure
/// works unchanged; only `Grid::giveOutput` strips the z column when
/// `#@domain 2` is set.
class Rect : public Region
{
protected:
    /// `[xmin ymin xmax ymax]`.
    oofem::FloatArray box;

    /// Spacing factor for points seeded on the rectangle edges
    /// (relative to `grid->diameter`). Default 1.
    double edgeRefine = 1.;

    /// Spacing factor for the random interior fill
    /// (relative to `grid->diameter`). Default 1.
    double regionRefine = 1.;

public:

    /**
     * Constructor. Creates an axis-aligned rectangle region.
     * @param n region number in the grid
     * @param aGrid grid to which the region belongs
     */
    Rect(int n, Grid *aGrid);

    /// Destructor.
    ~Rect();

    /// Fill `boundaries` with `[xmin xmax ymin ymax 0 0]`. The z slot is
    /// kept zero so the (3D-internal) localiser sees a degenerate z axis.
    void defineBoundaries(oofem::FloatArray &boundaries) override;

    /// Place edge vertices along the rectangle's faces (uniform spacing
    /// `edgeRefine * grid->diameter`). Called separately and earlier than
    /// `generatePoints()` so that the boundary discretisation is laid down
    /// before inclusions and interior fill.
    int generateBoundaryPoints() override;

    /// Random interior fill. Edge vertices are assumed to already exist
    /// (`generateBoundaryPoints()` is invoked beforehand by the pipeline);
    /// any duplicate edge placements here would be rejected by the
    /// neighbour-distance check anyway.
    int generatePoints() override;

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Rect"; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@rect <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);
};


#endif // rect_h
