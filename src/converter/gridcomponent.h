#ifndef gridcomponent_h
#define gridcomponent_h

#include "grid.h"

#include "floatarray.h"
#include "intarray.h"

#include "converterdatareader.h"
#include "convertertxtdatareader.h"
#include "convertertxtinputrecord.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

class GridComponent
{
protected:
    /// Array storing nodal coordinates.
    int number;

    Grid *grid;
public:

    /**
     * Constructor. Creates a node belonging to grid.
     * @param n node number in grid aGrid
     * @param aGrid grid to which node beloy
     */
    GridComponent(int n, Grid *aGrid);                      // constructor
    /// Destructor.
    virtual ~GridComponent() = default;                                           // destructor

    /// Returns class name of the receiver.
    virtual const char *giveClassName() const { return "GridComponent"; }
};


#endif // gridcomponent_h
