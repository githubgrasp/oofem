#ifndef vertex_h
#define vertex_h

#include "grid.h"

#include "floatarray.h"
#include "intarray.h"

#include "converterdatareader.h"
#include "convertertxtdatareader.h"
#include "convertertxtinputrecord.h"
#include "gridcomponent.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

#define _IFT_Vertex_coords "coords"
#define _IFT_Vertex_refine "refine"
#define _IFT_Vertex_radius "radius"

class Vertex : public GridComponent
{
protected:
  int printFlag;
  
    int outsideFlag;

    oofem::IntArray helpOutsideFlag;

    int boundaryFlag;
    
    oofem::IntArray cellVertices;

    oofem::IntArray cellElements;

    oofem::FloatArray coordinates;

    oofem::IntArray localLines, localTetras, localLinks;
    // localLines to stock elements linking elements of the same nature
    // localLink to stock elements linking elements of different natures
    // this distinction is set to avoid confusion between elements

    int location;

    int periodicNode;

    int number;

    double refinement;

    double radius;

public:

    /**
     * Constructor. Creates a node belonging to grid.
     * @param n node number in grid aGrid
     * @param aGrid grid to which node belongs
     */
    Vertex(int n, Grid *aGrid);                      // constructor
    
   // Vertex(Grid *aGrid);    // constructor, for temporary points, without numerotation
    
    /// Destructor.
    ~Vertex() override = default;                                           // destructor

    int giveNumber() { return number; }
    void updateNumber(int i) {this->number=i ;}

    // coordinates
    bool        hasCoordinates() { return true; }

    /// Returns i-th coordinate of node.
    double      giveCoordinate(int i);

    int giveLocation() { return this->location; }

    void setLocation(int loc) { this->location = loc; }

    void setPeriodicNode(int node) { this->periodicNode = node; }

    int givePeriodicNode() { return this->periodicNode; }

    double giveRadius() { return this->radius; }

    void setRadius(double rad) { this->radius = rad; }

    /// Returns pointer to node coordinate array.
    void giveCoordinates(oofem::FloatArray &coord) { coord = this->coordinates; }

    oofem::FloatArray *giveCoordinates() { return & coordinates; }

    void giveLocalLines(oofem::IntArray &elem) { elem = this->localLines; }

    void giveLocalTetras(oofem::IntArray &elem) { elem = this->localTetras; }
	
     void giveLocalLinks(oofem::IntArray &elem) { elem = this->localLinks; }

    void setLocalLine(int elem);

    void setLocalTetra(int elem);

    void setLocalLink(int elem);

    void setOutsideFlag(int flag) { this->outsideFlag = flag; }

    int giveOutsideFlag() { return this->outsideFlag; }

    void setHelpOutsideFlag(oofem::IntArray &_flags) { this->helpOutsideFlag = _flags; }

    void giveHelpOutsideFlag(oofem::IntArray &_flags) { _flags = this->helpOutsideFlag; }

    void setBoundaryFlag(int flag) { this->boundaryFlag = flag; }
    
    int giveBoundaryFlag() { return this->boundaryFlag; }
    
    void setCoordinates(oofem::FloatArray &_coords) { this->coordinates = _coords; }

    void setPrintFlag(int flag) { this->printFlag = flag; }

    int givePrintFlag() { return this->printFlag; }
    
    void updateCellVertices(oofem::IntArray &nodes);

    void giveCellVertices(oofem::IntArray &nodes) { nodes = this->cellVertices; }

    void updateCellElements(oofem::IntArray &elements);

    void giveCellElements(oofem::IntArray &elements) { elements = this->cellElements; }
    
    Vertex *ofType();

    const char *giveClassName() const { return "Vertex"; }

    void initializeFrom(ConverterInputRecord &ir);

    void  printYourself();
};

#endif // node_h
