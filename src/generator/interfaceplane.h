#ifndef interfaceplane_h
#define interfacecplane_h


#include "grid.h"
#include "inclusion.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <sstream>
#endif


//class FloatArray;
//class IntArray;

class InterfacePlane : public Inclusion
{
protected:
    /// Array storing nodal coordinates.
    oofem::FloatArray line; //Storing start and end point of axis of plane
    double diameter;
    double itzThickness;

public:

    /**
     * Constructor. Creates a planar interface inclusion.
     * @param n inclusion number in the grid
     * @param aGrid grid to which the inclusion belongs
     */
    InterfacePlane(int n, Grid *aGrid);
    /// Destructor.
    ~InterfacePlane();

    /// Returns the inclusion diameter.
    double giveDiameter() { return this->diameter; }
    /// Returns the thickness of the ITZ halo around the inclusion.
    double giveITZThickness() { return this->itzThickness; }

    /// Place points on the interface plane plus an ITZ halo of thickness
    /// `itzThickness`. Returns 1 on success.
    int generatePoints();


    /// Returns class name of the receiver.
    const char *giveClassName() const { return "InterfacePlane"; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@interfaceplane <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);
};


#endif // node_h
