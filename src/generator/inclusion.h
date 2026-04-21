#ifndef inclusion_h
#define inclusion_h


#include "grid.h"
#include "gridcomponent.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

//class FloatArray;
//class IntArray;

class Inclusion : public GridComponent
{
protected:
    /// Array storing nodal coordinates.



public:

    /**
     * Constructor. Creates an inclusion belonging to `aGrid`.
     * @param n inclusion number in the grid
     * @param aGrid grid to which the inclusion belongs
     */
    Inclusion(int n, Grid *aGrid);
    /// Destructor.
    virtual ~Inclusion();

    /// Generate points on the inclusion surface (and ITZ halo, if any).
    /// Pure virtual — each concrete inclusion implements its own shape.
    virtual int generatePoints() = 0;

    /// Generate periodic-image points for the inclusion when the grid
    /// is periodic. Default implementation is a no-op (returns 0).
    virtual int generatePeriodicPoints() { return 0; }
};


#endif // node_h
