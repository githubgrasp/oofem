#ifndef line_h
#define line_h


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
    int globalPossition;
    
    



    /* int keepingVoronoiMechanicalFlag; */
    /* int keepingVoronoiTransportFlag; */
    /* int keepingTransportFlag; */
    /* int delaunayTransportLineFlag; */
    /* int voronoiTransportBoundaryFlag; */
    /* int voronoiTransportCrossingFlag; */
    /* int voronoiTransportCorrespondingLine; */
    /* int voronoiMechanicalCorrespondingLine; */
    /* int delaunayMechanicalBoundaryFlag; */
    /* int voronoiTransportCorrespondancePassedLineFlag; */
    /* int voronoiMechanicalCorrespondancePassedLineFlag; */
    /* int thisTransportElementHasOutsideMechanicalCrossS; */
    /* oofem::IntArray theCouplingFlagsOfThisTransportElementCS; */
    /* oofem::IntArray theCouplingNumberOfThisTransportElementCS; */
    /* int thisMechanicalElementHasOutsideTransportCrossS; */
    /* oofem::IntArray theCouplingFlagsOfThisMechanicalElementCS; */
    /* oofem::IntArray theCouplingNumberOfThisMechanicalElementCS; */

    double radius;

    int edgeFlag;
    
    // info about material
    int m_typeOfMaterial;
    
    //for the links (case with fibre)
    double associated_length;

    
    // for the use of fibre elenents
    double diameter;
    oofem::FloatArray dir_vector; // direction vector of the fibre (not necessarily those of the elements...)
    
    double L_end;// distance to the nearest fibre endpoint

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

    //void setMechanicalLineFlag(int flag) { this->mechanicalLineFlag = flag; }

    //    void setDelaunayTransportLineFlag(int flag) { this->delaunayTransportLineFlag = flag; }
    //void setPossitionOfMechanicalLineInGlobal(oofem::IntArray &flag) { this->globalPossition = flag;}
    //    int giveLocalVertexFlag(oofem::IntArray &flag) { flag =  this->localVertexFlag; }
    //    int giveMechanicalLineFlag() { return this->mechanicalLineFlag; }
    //    int giveDelaunayTransportLineFlag() { return this->delaunayTransportLineFlag; }
    // int givePossitionOfMechanicalLineInGlobal(oofem::IntArray &flag) { flag = this->globalPossition;}
    Line *ofType();

    //    void setKeepingVoronoiMechanicalFlag(int flag) { this->keepingVoronoiMechanicalFlag = flag; }


    int giveOutsideFlag() { return this->outsideFlag; }
    void setOutsideFlag(int flag) { this->outsideFlag = flag; }

    void giveBoundaryFlags(oofem::IntArray &flags) { flags = this->boundaryFlags; }
    void setBoundaryFlags(oofem::IntArray &flags) { this->boundaryFlags = flags; }

    void giveInfinityFlags(oofem::IntArray &flags) { flags = this->infinityFlags; }
    void setInfinityFlags(oofem::IntArray &flags) { this->infinityFlags = flags; }
    
    void setAssociatedLength(double length){associated_length=length;}
    double giveAssociatedLength(){return associated_length;}


    //    void setVoronoiTransportBoundaryFlag(int flag) { this->voronoiTransportBoundaryFlag = flag; }

    //    void setVoronoiTransportCrossingFlag(int flag) { this->voronoiTransportCrossingFlag = flag; }

    //    void setVoronoiTransportCorrespondingLine(int flag) { this->voronoiTransportCorrespondingLine = flag; }

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
    const char *giveClassName() const { return "Line"; }

    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();
    
    
    // update material
    void updateMaterial(int typeOfMaterial);
    // access to material
    int giveMaterial() { return this->m_typeOfMaterial; }
    
    // for the use of fibre elenents
    void setDiameter(double diameter_fibre){diameter=diameter_fibre;}
    double giveDiameter(){return diameter;}
    void setDirVector(oofem::FloatArray dir_vector_fibre){dir_vector=dir_vector_fibre;}
    oofem::FloatArray giveDirectionVector(){return dir_vector;}
    void setL_end(double L){L_end=L;}
    double giveL_end(){return L_end;}

    
};


#endif //line_h
