#ifndef tetra_h
#define tetra_h


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

class Tetra : public GridComponent
{
protected:
    /// Array storing nodal coordinates.
    oofem::IntArray vertices;
    oofem::IntArray crossSectionVertices;
    oofem::IntArray crossSectionElements;

    int periodicElement;

    int outsideFlag;
    oofem::IntArray boundaryFlags;
    oofem::IntArray infinityFlags;

    int number;
    oofem::IntArray localVertexFlag;
    oofem::IntArray boundaryElements;
    int mechanicalLineFlag;
    int globalPossition;

    double radius;

    int edgeFlag;

public:

    /**
     * Constructor. Creates a tetra belonging to grid.
     * @param n element number in grid aGrid
     * @param aGrid grid to which node belongs
     */

    Tetra(int n, Grid *aGrid);                   // constructor
    /// Destructor.
    ~Tetra() override = default;                                           // destructor

    /// Returns i-th vertex of curve.
    int      giveLocalVertex(int i);
    /// Returns pointer to curve vertex array.

    void giveCoordinates(oofem::FloatArray &coords);

    void giveLocalVertices(oofem::IntArray &nodes) { nodes = this->vertices; }
    void setPeriodicElement(int element) { this->periodicElement = element; }
    int givePeriodicElement() { return this->periodicElement; }

    void setLocalVertices(oofem::IntArray &nodes) { this->vertices = nodes; }
    void setLocalVertexFlag(oofem::IntArray &flag) { this->localVertexFlag = flag; }

    Tetra *ofType();

    int giveOutsideFlag() { return this->outsideFlag; }
    void setOutsideFlag(int flag) { this->outsideFlag = flag; }

    void giveBoundaryFlags(oofem::IntArray &flags) { flags = this->boundaryFlags; }
    void setBoundaryFlags(oofem::IntArray &flags) { this->boundaryFlags = flags; }

    void giveInfinityFlags(oofem::IntArray &flags) { flags = this->infinityFlags; }
    void setInfinityFlags(oofem::IntArray &flags) { this->infinityFlags = flags; }


    /// Returns class name of the receiver.
    const char *giveClassName() const override { return "Line"; }

    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();
};


#endif //tetra_h
