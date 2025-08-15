#ifndef interfacesphere_h
#define interfacesphere_h

#include "inclusion.h"

#include "grid.h"

#include "floatarray.h"
#include "intarray.h"


#include "inputrecord.h"

#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "oofemtxtinputrecord.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

#define _IFT_InterfaceSphere_centre "centre"
#define _IFT_InterfaceSphere_refine "refine"
#define _IFT_InterfaceSphere_radius "radius"
#define _IFT_InterfaceSphere_itz "itz"

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

    void initializeFrom(GeneratorInputRecord &ir);

    void         printYourself();
};

#endif // interfasesphere_h
