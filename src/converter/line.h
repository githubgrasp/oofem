#ifndef line_h
#define line_h


#include "grid.h"
#include "gridcomponent.h"
#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

class Line : public GridComponent
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
    int globalPosition;

    double radius;

    int edgeFlag;

    // material assigned to the line by the writer.
    int materialType;

    // For fibre/matrix coupling links: portion of fibre length attributed to this link.
    double associatedLength;

    // Fibre-element parameters (also used by coupling links).
    double diameter;
    oofem::FloatArray directionVector;        // unit vector along the fibre, not necessarily along the element

    double endLength;                         // distance to the nearer fibre endpoint

public:

    /**
     * Constructor. Creates a node belonging to grid.
     * @param n node number in grid aGrid
     * @param aGrid grid to which node belongs
     */
    Line(int n, Grid *aGrid);                   // constructor
    /// Destructor.
    ~Line() override = default;                                           // destructor

    /// Returns i-th vertex of curve.
    int      giveLocalVertex(int i);
    /// Returns pointer to curve vertex array.
    void giveLocalVertices(oofem::IntArray &nodes) { nodes = vertices; }
    void setPeriodicElement(int element) { this->periodicElement = element; }
    int givePeriodicElement() { return this->periodicElement; }

    void setLocalVertices(oofem::IntArray &nodes) { this->vertices = nodes; }
    void setLocalVertexFlag(oofem::IntArray &flag) { this->localVertexFlag = flag; }

    int delaunayAreaCheck();

    Line *ofType();

    int giveOutsideFlag() { return this->outsideFlag; }
    void setOutsideFlag(int flag) { this->outsideFlag = flag; }

    void giveBoundaryFlags(oofem::IntArray &flags) { flags = this->boundaryFlags; }
    void setBoundaryFlags(oofem::IntArray &flags) { this->boundaryFlags = flags; }

    void giveInfinityFlags(oofem::IntArray &flags) { flags = this->infinityFlags; }
    void setInfinityFlags(oofem::IntArray &flags) { this->infinityFlags = flags; }

    void setAssociatedLength(double length) { associatedLength = length; }
    double giveAssociatedLength() { return associatedLength; }

    void updateCrossSectionElement(int element);
    void updateCrossSectionElements(oofem::IntArray &elements);

    void setRadius(double &rad) { this->radius = rad; }
    double giveRadius() { return this->radius; }

    void giveCrossSectionElements(oofem::IntArray &elements) { elements = this->crossSectionElements; }
    void setCrossSectionElements(oofem::IntArray &elements) { this->crossSectionElements = elements; }

    void setVertices(oofem::IntArray &_nodes) { this->vertices = _nodes; }
    void setCrossSectionVertices(oofem::IntArray &_nodes) { this->crossSectionVertices = _nodes; }

    void updateCrossSectionVertices(oofem::IntArray &nodes);


    void  giveCrossSectionVertices(oofem::IntArray &answer) const
    { answer = this->crossSectionVertices; }

    void  giveCrossSectionElements(oofem::IntArray &answer) const
    { answer = this->crossSectionElements; }


    /// Returns class name of the receiver.
    const char *giveClassName() const override { return "Line"; }

    /// Print receiver state on stdout — useful for debugging.
    void printYourself();

    void updateMaterial(int typeOfMaterial);
    int  giveMaterial() { return this->materialType; }

    // Fibre-element parameters (also used by coupling links).
    void setDiameter(double diameterFibre) { diameter = diameterFibre; }
    double giveDiameter() { return diameter; }
    void setDirectionVector(const oofem::FloatArray &dir) { directionVector = dir; }
    oofem::FloatArray giveDirectionVector() { return directionVector; }
    void setEndLength(double l) { endLength = l; }
    double giveEndLength() { return endLength; }
};


#endif //line_h
