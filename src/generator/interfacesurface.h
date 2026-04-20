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
    int number;
    double refinement;

public:

    /**
     * Constructor. Creates a interface surfacec belonging to grid.
     * @param n node number in grid aGrid
     * @param aGrid grid to which node belongs
     */

    InterfaceSurface(int n, Grid *aGrid);                      // constructor

    /// Destructor.
    ~InterfaceSurface();                                           // destructor

    int generatePoints();

    InterfaceSurface *ofType();

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "InterfaceSurface"; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@interfacesurface <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);

    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();
};

#endif // node_h
