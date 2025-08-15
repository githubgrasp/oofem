#ifndef cylinder_h
#define cylinder_h

#include "grid.h"
#include "region.h"

#include "floatarry.h"
#include "intarray.h"

#include "converterinputrecord.h"

#include "converterdatareader.h"
#include "convertertxtdatareader.h"
#include "convertertxtinputrecord.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

//class oofem::FloatArray;
//class IntArray;

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
    Cylinder(int n, Grid* aGrid);                      // constructor
    /// Destructor.
    ~Cylinder();                                           // destructor


    double giveRadius(){return this->radius;}
    void giveLine(oofem::FloatArray& lin){lin = line;}
    
    
    Cylinder *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Cylinder"; }

    IRResultType initializeFrom(InputRecord *ir);
    //virtual IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();

    virtual void defineBoundaries(oofem::FloatArray &boundaries);

    virtual void findOutsiders(oofem::FloatArray &boundaries);

    virtual int modifyVoronoiCrossSection(int elementNumber);
    
    virtual int areaCheck(int elementNumber);
    
};


#endif // cylinder_h






