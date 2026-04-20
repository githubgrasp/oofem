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
    int number;
    double refinement;
    oofem::FloatArray box;


public:

    /**
     * Constructor. Creates a node belonging to domain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    RefinePrism(int n, Grid *aGrid);                      // constructor
    /// Destructor.
    ~RefinePrism();                                           // destructor

    virtual double giveDiameter(oofem::FloatArray &coord) override;

    RefinePrism * ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "RefinePrism"; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@refineprism <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);
    //virtual IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();
};

#endif // refineprism_h
