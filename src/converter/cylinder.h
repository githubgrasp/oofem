#ifndef cylinder_h
#define cylinder_h

#include "grid.h"
#include "region.h"

#include "floatarray.h"
#include "intarray.h"



#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif


class Cylinder : public Region
{
protected:
    /// Array storing nodal coordinates.
    oofem::FloatArray line;
    double radius;
    int number;
    double refinement;

public:

    /**
     * Constructor. Creates a node belonging to domain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    Cylinder(int n, Grid *aGrid);                      // constructor
    /// Destructor.
    ~Cylinder() override = default;                                           // destructor


    double giveRadius() { return this->radius; }
    void giveLine(oofem::FloatArray &lin) { lin = line; }
    void setLine(const oofem::FloatArray &l) { line = l; }
    void setRadius(double r) { radius = r; }


    Cylinder *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const override { return "Cylinder"; }


    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();

    void defineBoundaries(oofem::FloatArray &boundaries) override;

    void findOutsiders(oofem::FloatArray &boundaries) override;

    int modifyVoronoiCrossSection(int elementNumber) override;

     int areaCheck(int elementNumber) override;
};


#endif // cylinder_h
