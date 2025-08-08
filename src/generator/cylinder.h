#ifndef cylinder_h
#define cylinder_h

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


#define _IFT_Cylinder_line "line"
#define _IFT_Cylinder_radius "radius"
#define _IFT_Cylinder_refine "refine"

//class FloatArray;
//class IntArray;

class Cylinder : public Region
{

protected:
    /// Array storing nodal coordinates.
  oofem::FloatArray line; //Storing start and end point of axis of cylinder
  double diameter;
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


    double giveDiameter(){return this->diameter;}
    
    virtual int generatePoints();

    Cylinder *ofType();

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Cylinder"; }

  void initializeFrom(GeneratorInputRecord &ir);

    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();

};


#endif // cylinder_h






