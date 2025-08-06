#ifndef prism_h
#define prism_h


#include "grid.h"

#include "region.h"

#include "floatarray.h"
#include "intarray.h"

using oofem::FloatArray;
using oofem::IntArray;

#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "oofemtxtinputrecord.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif


class Prism : public Region
{
protected:
    /// Array storing nodal coordinates.
    int number;
    double refinement;
  FloatArray box;
    double xlength, ylength, zlength;

public:

    /**
     * Constructor. Creates a node belonging to grid.
     * @param n node number in grid aGrid
     * @param aGrid grid to which node belongs
     */
    Prism(int n, Grid *aGrid);                      // constructor
    /// Destructor.
    ~Prism();                                           // destructor

    /// Returns i-th vertex of curve.
    int      giveLocalSurface(int i);
    /// Returns pointer to curve vertex array.

    void giveLocalSurfaces(IntArray &surf) { surf = this->surfaces; }

    /// Define boundaries
  void defineBoundaries(FloatArray &boundaries);
    
    //generate regular points
    int generateRegularPoints1();
    int generateRegularPoints2();

    //generate random points in periodic cell
    int generatePeriodicPoints();

  void mirrorShiftSurface(oofem::FloatArray& random, oofem::FloatArray& normal,FloatArray& specimenDimension,FloatArray& boundaries, int& vertexNumber, IntArray& periodicityFlag);
    
    //generate random points
    int generatePoints();

    //generate mixed points.
    //This means that some direction are periodic and others are not
    int generateMixedPoints();
    
    ///Returns the x length
    double giveXLength() { return this->xlength; }

    ///Returns the y length
    double giveYLength() { return this->ylength; }

    ///Returns the z length
    double giveZLength() { return this->zlength; }

    Prism *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Prism"; }

    ///Returns the number of region
    int giveNumber() { return this->number; }

  void initializeFrom(oofem::InputRecord *ir);
    //virtual IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();
};


#endif // node_h
