#ifndef curve_h
#define curve_h


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

//class oofem::FloatArray;
//class oofem::IntArray;

class Curve : public GridComponent
{

protected:
    /// Array storing nodal coordinates.
    oofem::IntArray vertices;
    int number;
    double refinement;

public:

    /**
     * Constructor. Creates a node belonging to grid.
     * @param n node number in grid aGrid
     * @param aGrid grid to which node belongs
     */
 Curve(int n, Grid* aGrid);                      // constructor
    /// Destructor.
    ~Curve();                                           // destructor

    /// Returns i-th vertex of curve.
    int      giveLocalVertex(int i);
    /// Returns pointer to curve vertex array.
    oofem::IntArray *giveLocalVertices() { return & vertices; }

    Vertex *giveGlobalVertex(int i, AList <Vertex> *vertexList);

    Curve *ofType();

    const char *giveClassName() const { return "Curve"; }
    
    IRResultType initializeFrom(InputRecord *ir);

    void         printYourself();

};


#endif // node_h






