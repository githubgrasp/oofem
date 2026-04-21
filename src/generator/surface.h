#ifndef surface_h
#define surface_h


#include "grid.h"
#include "gridcomponent.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <sstream>
#endif




//class FloatArray;
//class IntArray;

class Surface : public GridComponent
{
protected:
    /// Array storing nodal coordinates.
    oofem::IntArray curves;
    int boundaryFlag;
    oofem::FloatArray boundaryShift;
    oofem::FloatArray normal;

public:

    /**
     * Constructor. Creates a surface belonging to `aGrid`.
     * @param n surface number in the grid
     * @param aGrid grid to which the surface belongs
     */
    Surface(int n, Grid *aGrid);
    /// Destructor.
    ~Surface();

    /// Returns the `i`-th (1-based) bounding-curve id of the surface.
    int giveLocalCurve(int i);
    /// Returns a pointer to the surface's bounding-curve id array.
    oofem::IntArray *giveLocalCurves() { return & curves; }

    /// Returns the number of bounding curves of the surface.
    int giveNumberOfLocalCurves();

    /// Fill `boundaries` with the surface bounding box in
    /// `[xmin xmax ymin ymax zmin zmax]` order.
    void defineBoundaries(oofem::FloatArray &boundaries);

    /// Copies the surface's normal vector into `answer`.
    void giveNormal(oofem::FloatArray &answer) { answer = this->normal; }

    /// Generate random points on the surface, respecting periodic shifts
    /// and mirroring when the grid is periodic. Returns 1 on success.
    int generatePoints();

    /// Apply periodic mirror/shift to a candidate point. Used by the
    /// surface point generator to reflect a proposed point across
    /// periodic boundaries before final acceptance.
    void mirrorShift(oofem::FloatArray &random, oofem::FloatArray &normal, oofem::FloatArray &specimenDimension, oofem::FloatArray &boundaries, oofem::IntArray &periodicityFlag);

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Surface"; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@surface <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);

};


#endif // node_h
