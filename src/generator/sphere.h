#ifndef sphere_h
#define sphere_h

#include "grid.h"
#include "region.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <sstream>
#endif

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


    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@sphere <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);

    void printYourself();
};

#endif // sphere_h
