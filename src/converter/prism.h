#ifndef prism_h
#define prism_h


#include "grid.h"
#include "region.h"

#include "floatarray.h"
#include "intarray.h"

#include "converterdatareader.h"
#include "convertertxtdatareader.h"
#include "convertertxtinputrecord.h"
#include <list>

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

class GridLocalizer;

//class oofem::FloatArray;
//class oofem::IntArray;

class Prism : public Region
{
protected:
    /// Array storing nodal coordinates.
  oofem::FloatArray box;

    int number;
    double refinement;

    typedef std :: list< int > nodeContainerType;
    
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
    void giveLocalSurfaces(oofem::IntArray &surf) { surf = this->surfaces; }


    Prism *ofType();

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Prism"; }

    IRResultType initializeFrom(InputRecord *ir);
    //virtual oofem::IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();

    virtual int giveSwitches(oofem::IntArray& switches, const oofem::FloatArray& coords);
    virtual int giveCornerSwitches(oofem::IntArray &switches, const oofem::FloatArray &coords);
    virtual int giveEdgeSwitches(oofem::IntArray &switches, const oofem::FloatArray &coords);
    
    /** This function deletes cross-section elements which lie outside the specimen and
     *  shifts not the boundary for elements which cross the boundary*/
    virtual int modifyVoronoiCrossSection(int elementNumber);
    
    virtual void defineBoundaries(oofem::FloatArray &boundaries);

    virtual void findOutsiders(oofem::FloatArray &boundaries);
};


#endif // prism_h
