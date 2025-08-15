#ifndef boundarysphere_h
#define boundarysphere_h


#include "grid.h"
#include "region.h"

#include "floatarray.h"
#include "intarray.h"

#include "converterinputrecord.h"

#include "converterdatareader.h"
#include "convertertxtdatareader.h"
#include "convertertxtinputrecord.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

class BoundarySphere : public Region
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
    BoundarySphere(int n, Grid* aGrid);                      // constructor
    /// Destructor.
    ~BoundarySphere();                                           // destructor


    double giveRadius(){return this->radius;}
    double giveITZThickness(){return this->itzThickness;}
    void giveCentre(oofem::FloatArray& cent){cent = centre;}

    

    BoundarySphere *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "BoundarySphere"; }

    IRResultType initializeFrom(InputRecord *ir);
    //virtual IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();
    
    virtual void defineBoundaries(oofem::FloatArray &boundaries);

    virtual void findOutsiders(oofem::FloatArray &boundaries);

    virtual int modifyVoronoiCrossSection(int elementNumber);

    
};


#endif // node_h






