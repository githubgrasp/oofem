#ifndef interfacesphere_h
#define interfacesphere_h

#include "inclusion.h"

#include "grid.h"

#include "floatarray.h"
#include "intarray.h"
using oofem::FloatArray;
using oofem::IntArray;


#include "inputrecord.h"

#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "oofemtxtinputrecord.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

class InterfaceSphere : public Inclusion
{

protected:
    /// Array storing nodal coordinates.
  FloatArray centre;
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
    InterfaceSphere(int n, Grid* aGrid);                      // constructor
    /// Destructor.
    ~InterfaceSphere();                                           // destructor

    double giveDiameter(){return this->diameter;}
    double giveITZThickness(){return this->itzThickness;}
    int generatePoints();

    InterfaceSphere *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "InterfaceSphere"; }

  void initializeFrom(oofem::InputRecord *ir);
    //virtual IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();

};

#endif // interfasesphere_h






