#ifndef interfacecylinder_h
#define interfacecylinder_h


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

class InterfaceCylinder : public Inclusion
{
protected:
    /// Array storing nodal coordinates.
    oofem::FloatArray line; //Storing start and end point of axis of cylinder
    double diameter;
    double itzThickness;
    int nInterval;

public:

    /**
     * Constructor. Creates a cylindrical inclusion.
     * @param n inclusion number in the grid
     * @param aGrid grid to which the inclusion belongs
     */
    InterfaceCylinder(int n, Grid *aGrid);
    /// Destructor.
    ~InterfaceCylinder();

    /// Returns the inclusion diameter.
    double giveDiameter() { return this->diameter; }
    /// Returns the thickness of the ITZ halo around the inclusion.
    double giveITZThickness() { return this->itzThickness; }

    /// Place points on the cylinder surface plus an ITZ halo of thickness
    /// `itzThickness`. Returns 1 on success.
    int generatePoints();


    /// Returns class name of the receiver.
    const char *giveClassName() const { return "InterfaceCylinder"; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@interfacecylinder <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);

};


#endif // node_h
