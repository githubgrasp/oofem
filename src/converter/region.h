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

class Region : public GridComponent
{
protected:
    /// Array storing nodal coordinates.
    oofem::IntArray surfaces;
    int number;
    double refinement;

public:

    /**
     * Constructor. Creates a region belonging to grid.
     * @param n region number in grid aGrid
     * @param aGrid grid to which node belongs
     */
    Region(int n, Grid *aGrid);                      // constructor
    /// Destructor.
    virtual ~Region() = default;                                           // destructor

    virtual void initializeFrom(ConverterInputRecord &ir) = 0;

    virtual void defineBoundaries(oofem::FloatArray &boundaries) = 0;

    virtual void findOutsiders(oofem::FloatArray &boundaries) = 0;

    virtual int giveSwitches(oofem::IntArray &switches, oofem::FloatArray &coords){return 0;}

    virtual int modifyVoronoiCrossSection(int elementNumber){return 1;}

    virtual int areaCheck(int elementNumber){return 0;};
};


#endif // region_h
