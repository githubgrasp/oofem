#ifndef interfaceplane_h
#define interfacecplane_h


#include "grid.h"
#include "inclusion.h"

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

class InterfacePlane : public Inclusion
{

protected:
    /// Array storing nodal coordinates.
  FloatArray line; //Storing start and end point of axis of plane
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
    InterfacePlane(int n, Grid* aGrid);                      // constructor
    /// Destructor.
    ~InterfacePlane();                                           // destructor


    double giveDiameter(){return this->diameter;}
    double giveITZThickness(){return this->itzThickness;}
    
    int generatePoints();

    InterfacePlane *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "InterfacePlane"; }

    IRResultType initializeFrom(InputRecord *ir);
    //virtual IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();

};


#endif // node_h






