#ifndef region_h
#define region_h


#include "grid.h"
#include "gridcomponent.h"

#include "floatarray.h"
#include "intarray.h"

#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "oofemtxtinputrecord.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

//class FloatArray;
//class IntArray;

class Region : public GridComponent
{
protected:
    /// Array storing nodal coordinates.
  oofem::IntArray surfaces;
    int number;
    double refinement;
    double xlength, ylength, zlength;

public:

    /**
     * Constructor. Creates a node belonging to grid.
     * @param n node number in grid aGrid
     * @param aGrid grid to which node belongs
     */
    Region(int n, Grid *aGrid);                      // constructor
    /// Destructor.
    virtual ~Region();                                           // destructor

    /// Returns i-th vertex of curve.
    int      giveLocalSurface(int i);
    /// Returns pointer to curve vertex array.

  void giveLocalSurfaces(oofem::IntArray &surf) { surf = this->surfaces; }

    /// Define boundaries
  void defineBoundaries(oofem::FloatArray &boundaries);
    
    //generate regular points
    int generateRegularPoints1();
    int generateRegularPoints2();

    //generate random points in periodic cell
    int generatePeriodicPoints();

    //generate random points
    virtual int generatePoints();

    //generate mixed points.
    //This means that some direction are periodic and others are not
    int generateMixedPoints();
    
    ///Returns the x length
    double giveXLength() { return this->xlength; }

    ///Returns the y length
    double giveYLength() { return this->ylength; }

    ///Returns the z length
    double giveZLength() { return this->zlength; }

    Region *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Region"; }

    int giveNumber() { return this->number; }

  void initializeFrom(oofem::InputRecord *ir);

    void         printYourself();
};


#endif // node_h
