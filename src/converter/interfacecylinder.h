#ifndef interfacecylinder_h
#define interfacecylinder_h


#include "grid.h"
#include "inclusion.h"

#include "floatarray.h"
#include "intarray.h"

#include "converterinputrecord.h"

#include "converterdatareader.h"
#include "convertertxtdatareader.h"
#include "convertertxtinputrecord.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

//class FloatArray;
//class oofem::IntArray;

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
    InterfaceCylinder(int n, Grid* aGrid);                      // constructor
    /// Destructor.
    ~InterfaceCylinder();                                           // destructor


    double giveRadius(){return this->radius;}
    double giveITZThickness(){return this->itzThickness;}
    void giveLine(oofem::FloatArray& lin){lin = line;}

    

    InterfaceCylinder *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "InterfaceCylinder"; }

    IRResultType initializeFrom(InputRecord *ir);
    //virtual oofem::IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();

};


#endif // node_h






