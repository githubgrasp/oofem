#ifndef sphere_h
#define sphere_h


#include "grid.h"
#include "region.h"

#include "flotarry.h"
#include "intarray.h"

#include "inputrecord.h"

#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "oofemtxtinputrecord.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

//class FloatArray;
//class IntArray;

class Sphere : public Region
{

protected:
    /// Array storing nodal coordinates.
  FloatArray centre;
  double radius;
  int number;
  double refinement;


public:

    /**
     * Constructor. Creates a node belonging to domain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    Sphere(int n, Grid* aGrid);                      // constructor
    /// Destructor.
    ~Sphere();                                           // destructor

    double giveRadius(){return this->radius;}
    
    virtual int generatePoints();

    int generatePeriodicPoints();

    Sphere *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Sphere"; }

    ///Returns the number of region
    int giveNumber() { return this->number; }

    
    IRResultType initializeFrom(InputRecord *ir);
    //virtual IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();

};

#endif // sphere_h






