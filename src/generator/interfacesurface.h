#ifndef interfacesurface_h
#define interfacesurface_h

#include "grid.h"
#include "inclusion.h"

#include "floatarray.h"
#include "intarray.h"

#include "inputrecord.h"

#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "oofemtxtinputrecord.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif


#define _IFT_InterfaceSurface_curves "curves"
#define _IFT_InterfaceSurface_refine "refine"
 

class
InterfaceSurface : public Inclusion
{

protected:
    /// Array storing nodal coordinates.
  oofem::IntArray curves;
    int number;
    double refinement;

public:

    /**
     * Constructor. Creates a interface surfacec belonging to grid.
     * @param n node number in grid aGrid
     * @param aGrid grid to which node belongs
     */

    InterfaceSurface(int n, Grid* aGrid);                      // constructor

    /// Destructor.
    ~InterfaceSurface();                                           // destructor
    
    int generatePoints();

    InterfaceSurface *ofType();

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "InterfaceSurface"; }

    void initializeFrom(GeneratorInputRecord &ir);

    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();

};

#endif // node_h






