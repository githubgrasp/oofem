#ifndef prism_h
#define prism_h


#include "grid.h"

#include "region.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <sstream>
#endif

class Prism : public Region
{
protected:
    /// Array storing nodal coordinates.
    // per-stage spacing ratios applied to grid->diameter; default 1.
    double edgeRefine;
    double surfaceRefine;
    double regionRefine;
    oofem::FloatArray box;

public:

    /**
     * Constructor. Creates an axis-aligned box region.
     * @param n region number in the grid
     * @param aGrid grid to which the region belongs
     */
    Prism(int n, Grid *aGrid);
    /// Destructor.
    ~Prism();

    /// Returns the `i`-th (1-based) bounding-surface id.
    int giveLocalSurface(int i);
    /// Copies the bounding-surface id array into `surf`.
    void giveLocalSurfaces(oofem::IntArray &surf) { surf = this->surfaces; }

    void defineBoundaries(oofem::FloatArray &boundaries) override;

    /// Regular point-grid seeding, variant 1 (BCC-style).
    int generateRegularPoints1();
    /// Regular point-grid seeding, variant 2 (FCC-style).
    int generateRegularPoints2();

    /// Generate random points in the box treating it as fully periodic.
    int generatePeriodicPoints();

    /// Apply periodic mirror/shift to a candidate point on the box
    /// boundary. Used by the periodic point generator to reflect a
    /// proposed point across periodic faces before final acceptance.
    void mirrorShiftSurface(oofem::FloatArray &random, oofem::FloatArray &normal, oofem::FloatArray &specimenDimension, oofem::FloatArray &boundaries, oofem::IntArray &periodicityFlag);

    /// Generate random points in the box (non-periodic fallback).
    int generatePoints();

    /// Generate random points when some axes are periodic and others
    /// are not.
    int generateMixedPoints();

    /// Returns the extent of the prism along x.
    double giveXLength() { return this->xlength; }
    /// Returns the extent of the prism along y.
    double giveYLength() { return this->ylength; }
    /// Returns the extent of the prism along z.
    double giveZLength() { return this->zlength; }


    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Prism"; }

    /// Returns the region number (1-based).
    int giveNumber() { return this->number; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@prism <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);
};


#endif // node_h
