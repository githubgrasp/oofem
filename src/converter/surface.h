#ifndef surface_h
#define surface_h


#include "grid.h"
#include "gridcomponent.h"

#include "floatarray.h"
#include "intarray.h"

#include "converterdatareader.h"
#include "convertertxtdatareader.h"
#include "convertertxtinputrecord.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

#define _IFT_Surface_curves "curves"
#define _IFT_Surface_refine "refine"

class Surface : public GridComponent
{
protected:
    /// Array storing nodal coordinates.
    oofem::IntArray curves;
    int number;
    double refinement;

public:

    /**
     * Constructor. Creates a node belonging to grid.
     * @param n node number in grid aGrid
     * @param aGrid grid to which node belongs
     */
    Surface(int n, Grid *aGrid);

    /// Destructor.
    ~Surface() override = default;

    /// Returns i-th vertex of curve.
    int      giveLocalCurve(int i);

    /// Returns pointer to curve vertex array.
    oofem::IntArray *giveLocalCurves() { return & curves; }

    int giveNumberOfLocalCurves();




    Surface *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Surface"; }

    void initializeFrom(ConverterInputRecord &ir);

    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();
};


#endif // surface_h
