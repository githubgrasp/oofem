#ifndef vertex_h
#define vertex_h

#include "grid.h"

#include "floatarray.h"
#include "intarray.h"

#include "gridcomponent.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <sstream>
#endif


//class FloatArray;
//class IntArray;

class Vertex : public GridComponent
{
public:




    /**
     * Constructor. Creates a vertex belonging to `aGrid`.
     * @param n vertex number in the grid
     * @param aGrid grid to which the vertex belongs
     */
    Vertex(int n, Grid *aGrid);
    /// Destructor.
    virtual ~Vertex();

    /// Always true — vertices carry coordinates by construction.
    bool hasCoordinates() { return true; }
    /// Returns the `i`-th (1-based) coordinate of the vertex.
    double giveCoordinate(int i);

    /// Local refinement factor (multiplier on `Grid::diameter`).
    double giveRefinement() { return this->refinement; }

    /// Exclusion radius around the vertex during random placement.
    double giveRadius() { return this->radius; }

    /// Copies the vertex coordinates into `coord`.
    void giveCoordinates(oofem::FloatArray &coord) { coord = this->coordinates; }
    /// Returns a pointer to the internal coordinate array.
    oofem::FloatArray *giveCoordinates() { return & coordinates; }

    /// Overwrites the vertex coordinates.
    void setCoordinates(const oofem::FloatArray &coords) { this->coordinates = coords; }


    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Vertex"; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@vertex <num>` (or `#@controlvertex <num>`) prefix.
    void initializeFromTokens(std::istringstream &iss);
};


#endif // node_h
