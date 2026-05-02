#ifndef converter_rect_h
#define converter_rect_h

#include "grid.h"
#include "region.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif


/// 2D analog of `Prism` in the converter — an axis-aligned rectangle. For
/// Stage 6 MVP the rectangle is purely a bounding-box descriptor (used by
/// `defineBoundaries` and the inside-rectangle test). Boundary clipping of
/// Voronoi cross-sections (the equivalent of `Prism::modifyVoronoiCrossSection`
/// in 3D) is not implemented yet — Delaunay edges with any endpoint outside
/// the rectangle are dropped during emission.
class Rect : public Region
{
protected:
    /// `[xmin ymin xmax ymax]`.
    oofem::FloatArray box;

public:

    /**
     * Constructor.
     * @param n region number in grid
     * @param aGrid grid to which the region belongs
     */
    Rect(int n, Grid *aGrid);

    /// Destructor.
    ~Rect() override = default;

    /// Set the bounding box `[xmin ymin xmax ymax]`.
    void setBox(const oofem::FloatArray &b) { box = b; }

    /// Returns the bounding box `[xmin ymin xmax ymax]`.
    const oofem::FloatArray &giveBox() const { return box; }

    /// Fill `boundaries` with `[xmin xmax ymin ymax 0 0]` — the same 6-slot
    /// convention as the 3D regions, with zero z extent.
    void defineBoundaries(oofem::FloatArray &boundaries) override;

    /// No-op for the MVP — boundary handling defers to the inside-rect
    /// midpoint test in the emitter.
    void findOutsiders(oofem::FloatArray &boundaries) override {}

    /// True iff `(x, y)` lies inside the rectangle (with `tol` margin).
    bool contains(double x, double y, double tol) const;

    /// Returns class name of the receiver.
    const char *giveClassName() const override { return "Rect"; }
};


#endif // converter_rect_h
