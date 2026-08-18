#ifndef interfacecylinder_h
#define interfacecylinder_h


#include "grid.h"
#include "inclusion.h"

#include "floatarray.h"
#include "intarray.h"



#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif



class InterfaceCylinder : public Inclusion
{
protected:
    /// Array storing nodal coordinates.
    oofem::FloatArray line;
    double radius;
    int number;
    double refinement;
    double itzThickness;

public:

    /**
     * Constructor. Creates a node belonging to domain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    InterfaceCylinder(int n, Grid *aGrid);                      // constructor
    /// Destructor.
    ~InterfaceCylinder() override = default;                                           // destructor


    double giveRadius() { return this->radius; }
    double giveITZThickness() { return this->itzThickness; }
    void giveLine(oofem::FloatArray &lin) { lin = line; }
    void setLine(const oofem::FloatArray &l) { line = l; }
    void setRadius(double r) { radius = r; }
    void setITZThickness(double t) { itzThickness = t; }



    InterfaceCylinder *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const override { return "InterfaceCylinder"; }


    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();
};


#endif // interfacecylinder_h
