#ifndef boundarysphere_h
#define boundarysphere_h


#include "grid.h"
#include "region.h"

#include "floatarray.h"
#include "intarray.h"



#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif


class BoundarySphere : public Region
{
protected:
    /// Array storing nodal coordinates.
    oofem::FloatArray centre;
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
    BoundarySphere(int n, Grid *aGrid);                      // constructor
    /// Destructor.
    ~BoundarySphere() override = default;                                           // destructor


    double giveRadius() { return this->radius; }
    double giveITZThickness() { return this->itzThickness; }
    void giveCentre(oofem::FloatArray &cent) { cent = centre; }


    BoundarySphere *ofType();

    const char *giveClassName() const override { return "BoundarySphere"; }


    void printYourself();

    virtual void defineBoundaries(oofem::FloatArray &boundaries) override;

    virtual void findOutsiders(oofem::FloatArray &boundaries) override;

    virtual int modifyVoronoiCrossSection(int elementNumber) override;
};


#endif // boundarysphere_h
