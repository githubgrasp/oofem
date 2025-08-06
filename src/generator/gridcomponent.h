#ifndef gridcomponent_h
#define gridcomponent_h

#include "grid.h"

#include "floatarray.h"
#include "intarray.h"

#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "oofemtxtinputrecord.h"

#include "inputrecord.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

//class FloatArray;
//class IntArray;

class GridComponent
{

protected:
    /// Array storing nodal coordinates.
  oofem::FloatArray coordinates;
    int number;
    double refinement;
    double radius;
    
    Grid* grid;
public:

    /**
     * Constructor. Creates a node belonging to grid.
     * @param n node number in grid aGrid
     * @param aGrid grid to which node belongs
     */
    GridComponent(int n, Grid* aGrid);                      // constructor
    /// Destructor.
    ~GridComponent();                                           // destructor

    GridComponent *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "GridComponent"; }

  initializeFrom(InputRecord *ir);
    //virtual IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();

};


#endif // node_h






