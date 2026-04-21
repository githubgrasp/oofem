#ifndef sphere_h
#define sphere_h

#include "grid.h"
#include "region.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <sstream>
#endif

//class FloatArray;
//class IntArray;

class Sphere : public Region
{
protected:
    /// Array storing nodal coordinates.
    oofem::FloatArray centre;


public:

    /**
     * Constructor. Creates a spherical region.
     * @param n region number in the grid
     * @param aGrid grid to which the region belongs
     */
    Sphere(int n, Grid *aGrid);
    /// Destructor.
    ~Sphere();

    /// Returns the sphere radius.
    double giveRadius() { return this->radius; }

    /// Generate random points inside the sphere.
    virtual int generatePoints();

    /// Generate periodic-image points when the grid is periodic.
    int generatePeriodicPoints();

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Sphere"; }

    /// Returns the region number (1-based).
    int giveNumber() { return this->number; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@sphere <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);

};

#endif // sphere_h
