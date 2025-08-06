#ifndef refinement_h
#define refinement_h


#include "grid.h"
#include "gridcomponent.h"

#include "flotarry.h"
#include "intarray.h"

#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "oofemtxtinputrecord.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

class Refinement : public GridComponent
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
    Refinement(int n, Grid* aGrid);                      // constructor
    /// Destructor.
    ~Refinement();                                           // destructor

    virtual double giveDiameter(FloatArray &coord) {return 0.;}

};


#endif // node_h






