#ifndef refineprism_h
#define refineprism_h

#include "refinement.h"

#include "grid.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <sstream>
#endif


class RefinePrism : public Refinement
{
protected:
    /// Array storing nodal coordinates.
    oofem::FloatArray box;


public:

    /**
     * Constructor. Creates a box-shaped local refinement region.
     * @param n refinement number in the grid
     * @param aGrid grid to which the refinement belongs
     */
    RefinePrism(int n, Grid *aGrid);
    /// Destructor.
    ~RefinePrism();

    virtual double giveDiameter(oofem::FloatArray &coord) override;


    /// Returns class name of the receiver.
    const char *giveClassName() const { return "RefinePrism"; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@refineprism <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);
};

#endif // refineprism_h
