#ifndef interfacesurface_h
#define interfacesurface_h

#include "grid.h"
#include "inclusion.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <sstream>
#endif


class
InterfaceSurface : public Inclusion
{
protected:
    /// Array storing nodal coordinates.
    oofem::IntArray curves;

public:

    /**
     * Constructor. Creates a curve-bounded interface surface.
     * @param n inclusion number in the grid
     * @param aGrid grid to which the inclusion belongs
     */
    InterfaceSurface(int n, Grid *aGrid);
    /// Destructor.
    ~InterfaceSurface();

    /// Place points on the interface surface. Returns 1 on success.
    int generatePoints();


    /// Returns class name of the receiver.
    const char *giveClassName() const { return "InterfaceSurface"; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@interfacesurface <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);

};

#endif // node_h
