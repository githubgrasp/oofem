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
    int number;
    double refinement;
    double itzThickness;

public:

    /**
     * Constructor. Creates a node belonging to domain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    InterfaceSphere(int n, Grid *aGrid);                      // constructor
    /// Destructor.
    ~InterfaceSphere();                                           // destructor

    double giveDiameter() { return this->diameter; }
    double giveITZThickness() { return this->itzThickness; }
    int generatePoints();

    InterfaceSphere *ofType();

    const char *giveClassName() const { return "InterfaceSphere"; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@intersphere <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);

    void         printYourself();
};

#endif // interfasesphere_h
