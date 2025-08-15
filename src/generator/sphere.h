#ifndef sphere_h
#define sphere_h

#include "grid.h"
#include "region.h"

#include "floatarray.h"
#include "intarray.h"

#include "inputrecord.h"

#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "oofemtxtinputrecord.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

#define _IFT_Sphere_centre "centre"
#define _IFT_Sphere_refine "refine"
#define _IFT_Sphere_radius "radius"

//class FloatArray;
//class IntArray;

class Sphere : public Region
{
protected:
    /// Array storing nodal coordinates.
    oofem::FloatArray centre;
    double radius;
    int number;
    double refinement;


public:

    /**
     * Constructor. Creates a node belonging to domain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    Sphere(int n, Grid *aGrid);                      // constructor
    /// Destructor.
    ~Sphere();                                           // destructor

    double giveRadius() { return this->radius; }

    virtual int generatePoints();

    int generatePeriodicPoints();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Sphere"; }

    ///Returns the number of region
    int giveNumber() { return this->number; }


    void initializeFrom(GeneratorInputRecord &ir);

    void printYourself();
};

#endif // sphere_h
