#ifndef interfacecylinder_h
#define interfacecylinder_h


#include "grid.h"
#include "inclusion.h"

#include "floatarray.h"
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

class InterfaceCylinder : public Inclusion
{

protected:
    /// Array storing nodal coordinates.
  FloatArray line; //Storing start and end point of axis of cylinder
  double diameter;
  int number;
  double refinement;
  double itzThickness;
  int nInterval;

public:

    /**
     * Constructor. Creates a node belonging to domain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    InterfaceCylinder(int n, Grid* aGrid);                      // constructor
    /// Destructor.
    ~InterfaceCylinder();                                           // destructor


    double giveDiameter(){return this->diameter;}
    double giveITZThickness(){return this->itzThickness;}
    
    int generatePoints();

    InterfaceCylinder *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "InterfaceCylinder"; }

  void initializeFrom(oofem::InputRecord *ir);
    //virtual IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();

};


#endif // node_h






