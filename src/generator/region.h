#ifndef region_h
#define region_h


#include "grid.h"
#include "gridcomponent.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif


//class FloatArray;
//class IntArray;

class Region : public GridComponent
{
protected:
    /// Array storing nodal coordinates.
    oofem::IntArray surfaces;
    double xlength, ylength, zlength;

public:

    /**
     * Constructor. Creates a region belonging to `aGrid`.
     * @param n region number in the grid
     * @param aGrid grid to which the region belongs
     */
    Region(int n, Grid *aGrid);
    /// Destructor.
    virtual ~Region();

    /// Returns the `i`-th (1-based) bounding-surface id of the region.
    int giveLocalSurface(int i);
    /// Copies the region's bounding-surface id array into `surf`.
    void giveLocalSurfaces(oofem::IntArray &surf) { surf = this->surfaces; }

    /// Fill `boundaries` with the region bounding box in
    /// `[xmin xmax ymin ymax zmin zmax]` order.
    virtual void defineBoundaries(oofem::FloatArray &boundaries);

    /// Regular point-grid seeding, variant 1 (BCC-style).
    int generateRegularPoints1();
    /// Regular point-grid seeding, variant 2 (FCC-style).
    int generateRegularPoints2();

    /// Generate random points in the region treating it as fully
    /// periodic on all active axes.
    int generatePeriodicPoints();

    /// Generate random points in the region (non-periodic fallback).
    /// Subclasses override to provide shape-specific acceptance tests.
    virtual int generatePoints();

    /// Generate random points in the region when some axes are periodic
    /// and others are not.
    int generateMixedPoints();

    /// Returns the region's extent along x.
    double giveXLength() { return this->xlength; }

    /// Returns the region's extent along y.
    double giveYLength() { return this->ylength; }

    /// Returns the region's extent along z.
    double giveZLength() { return this->zlength; }

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Region"; }

    /// Returns the region number (1-based).
    int giveNumber() { return this->number; }

};


#endif // node_h
