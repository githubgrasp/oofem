#ifndef curve_h
#define curve_h


#include "gridcomponent.h"
#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <sstream>
#endif


class Curve : public GridComponent
{
protected:
    /// Array storing nodal coordinates.
    oofem::IntArray vertices;
    oofem::FloatArray normal;

public:

    /**
     * Constructor. Creates a curve belonging to `aGrid`.
     * @param n curve number in the grid
     * @param aGrid grid to which the curve belongs
     */
    Curve(int n, Grid *aGrid);

    /// Destructor.
    virtual ~Curve();

    /// Returns the `i`-th (1-based) vertex id of the curve.
    int giveLocalVertex(int i);
    /// Returns a pointer to the curve's vertex-id array.
    oofem::IntArray *giveLocalVertices() { return & vertices; }

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Curve"; }

    /// Generate random points along the curve using the current grid's
    /// diameter / refinement. Returns 1 on success.
    int generatePoints();

    /// Generate periodic-image points for the curve when the grid is
    /// periodic. Returns 1 on success.
    int generatePeriodicPoints();

    /// Apply periodic mirror/shift to a candidate point. Used by the
    /// periodic point generator to reflect a proposed point across
    /// periodic boundaries before final acceptance.
    void mirrorShift(oofem::FloatArray &random, oofem::FloatArray &normal, oofem::FloatArray &specimenDimension, oofem::FloatArray &boundaries, oofem::IntArray &periodicityFlag);

    /// Copies the curve's normal vector into `answer`.
    void giveNormal(oofem::FloatArray &answer) { answer = this->normal; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@curve <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);

};


#endif // node_h
