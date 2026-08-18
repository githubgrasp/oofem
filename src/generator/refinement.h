#ifndef refinement_h
#define refinement_h


#include "grid.h"
#include "gridcomponent.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

class Refinement : public GridComponent
{
protected:
    /// Array storing nodal coordinates.



public:

    /**
     * Constructor. Creates a refinement entry belonging to `aGrid`.
     * @param n refinement number in the grid
     * @param aGrid grid to which the refinement belongs
     */
    Refinement(int n, Grid *aGrid);
    /// Destructor.
    virtual ~Refinement();

    /// Returns the target diameter at point `coord`. Default implementation
    /// returns 0 (i.e. no override); subclasses override to return a local
    /// target spacing inside their region of influence.
    virtual double giveDiameter(oofem::FloatArray &coord) { return 0.; }
};


#endif // node_h
