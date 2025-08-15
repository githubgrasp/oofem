#ifndef vertex_h
#define vertex_h

#include "grid.h"

#include "floatarray.h"
#include "intarray.h"

#include "generatordatareader.h"
#include "generatortxtdatareader.h"
#include "generatortxtinputrecord.h"
#include "gridcomponent.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

#define _IFT_Vertex_coords "coords"
#define _IFT_Vertex_refine "refine"
#define _IFT_Vertex_radius "radius"


//class FloatArray;
//class IntArray;

class Vertex : public GridComponent
{
protected:
    /// Array storing nodal coordinates.
    oofem::FloatArray coordinates;
    int number;
    double refinement;
    double radius;
    int randomswitch;
public:




    /**
     * Constructor. Creates a node belonging to domain.
     * @param n node number in grid aGrid
     * @param aGrid grid to which node belongs
     */
    Vertex(int n, Grid *aGrid);                      // constructor
    /// Destructor.
    virtual ~Vertex();                                           // destructor

    // coordinates
    bool        hasCoordinates() { return true; }
    /// Returns i-th coordinate of node.
    double      giveCoordinate(int i);


    double giveRefinement() { return this->refinement; }

    double giveRadius() { return this->radius; }

    /// Returns pointer to node coordinate array
    void giveCoordinates(oofem::FloatArray &coord) { coord = this->coordinates; }
    oofem::FloatArray *giveCoordinates() { return & coordinates; }

    void setCoordinates(const oofem::FloatArray &coords) { this->coordinates = coords; }

    /// Sets i-th componet. The component will be futher managed and maintained by grid object.


    Vertex *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Vertex"; }

    void initializeFrom(GeneratorInputRecord &ir);
    //virtual IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();
};


#endif // node_h
