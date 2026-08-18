#ifndef cylinder_h
#define cylinder_h

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

class Cylinder : public Region
{
protected:
    /// Array storing nodal coordinates.
    oofem::FloatArray line; //Storing start and end point of axis of cylinder
    double diameter;

public:

    /**
     * Constructor. Creates a cylindrical region.
     * @param n region number in the grid
     * @param aGrid grid to which the region belongs
     */
    Cylinder(int n, Grid *aGrid);
    /// Destructor.
    ~Cylinder();

    /// Returns the cylinder diameter.
    double giveDiameter() { return this->diameter; }

    /// Generate random points inside the cylinder.
    virtual int generatePoints();


    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Cylinder"; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@cylinder <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);

};


#endif // cylinder_h
