#ifndef interfacesphere_h
#define interfacesphere_h


#include "grid.h"
#include "inclusion.h"

#include "floatarray.h"
#include "intarray.h"

#include "inputrecord.h"

#include "converterdatareader.h"
#include "convertertxtdatareader.h"
#include "convertertxtinputrecord.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

//class FloatArray;
//class IntArray;

class InterfaceSphere : public Inclusion
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
    InterfaceSphere(int n, Grid* aGrid);                      // constructor
    /// Destructor.
    ~InterfaceSphere();                                           // destructor


    double giveRadius(){return this->radius;}
    double giveITZThickness(){return this->itzThickness;}
    void giveCentre(FloatArray& cent){cent = centre;}

    

    InterfaceSphere *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "InterfaceSphere"; }

    IRResultType initializeFrom(InputRecord *ir);
    //virtual IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();

};


#endif // node_h






