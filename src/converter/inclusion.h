#ifndef inclusion_h
#define inclusion_h


#include "grid.h"
#include "gridcomponent.h"

#include "flotarry.h"
#include "intarray.h"

#include "converterdatareader.h"
#include "convertertxtdatareader.h"
#include "convertertxtinputrecord.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif


class Inclusion : public GridComponent
{

protected:

    int number;


public:

    /**
     * Constructor. Creates a node belonging to domain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    Inclusion(int n, Grid* aGrid);                      // constructor
    /// Destructor.
    ~Inclusion();                                           // destructor

    virtual const char *giveClassName() const = 0;
};


#endif // node_h






