#ifndef inclusion_h
#define inclusion_h


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

class Inclusion : public GridComponent
{
protected:
    /// Array storing nodal coordinates.

    int number;


public:

    /**
     * Constructor. Creates a node belonging to domain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    Inclusion(int n, Grid *aGrid);                      // constructor
    /// Destructor.
    virtual ~Inclusion();                                           // destructor

    virtual int generatePoints() = 0;

    virtual int generatePeriodicPoints() { return 0; }
};


#endif // node_h
