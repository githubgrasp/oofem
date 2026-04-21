#ifndef interfacesphere_h
#define interfacesphere_h

#include "inclusion.h"

#include "grid.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <sstream>
#endif

class InterfaceSphere : public Inclusion
{
protected:
    /// Array storing nodal coordinates.
    oofem::FloatArray centre;
    double diameter;
    double itzThickness;

public:

    /**
     * Constructor. Creates a spherical inclusion.
     * @param n inclusion number in the grid
     * @param aGrid grid to which the inclusion belongs
     */
    InterfaceSphere(int n, Grid *aGrid);
    /// Destructor.
    ~InterfaceSphere();

    /// Returns the inclusion diameter.
    double giveDiameter() { return this->diameter; }
    /// Returns the thickness of the ITZ halo around the inclusion.
    double giveITZThickness() { return this->itzThickness; }
    /// Place points on the sphere surface plus an ITZ halo of thickness
    /// `itzThickness`. Returns 1 on success.
    int generatePoints();


    /// Returns class name of the receiver.
    const char *giveClassName() const { return "InterfaceSphere"; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@intersphere <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);

};

#endif // interfasesphere_h
