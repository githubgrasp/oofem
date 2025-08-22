#ifndef FIBRE_H_INCLUDED
#define FIBRE_H_INCLUDED



#include <vector>
#include "grid.h"
#include "vertex.h"
#include "gridcomponent.h"
#include "floatarray.h"
#include "intarray.h"

#include "converterdatareader.h"
#include "convertertxtdatareader.h"
#include "convertertxtinputrecord.h"

#define  _IFT_Fibre_endpoints "endpoints"
#define  _IFT_Fibre_diameter "diameter"


#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif


// warning 1 : vector indices starts at O, for the other type : at 1
// warning 2 : only staight fibres are here taken into account
// warning 3 : all fibres in input file must be entirely within the 27 cell structure

class Fibre : public GridComponent
{
protected:

    double m_TOL;// used to avoid computation instabilities during discretization

    /// Array storing nodal coordinates.
    oofem::FloatArray m_endpoints;
    oofem::FloatArray coordP, coordQ;// coord endpoints


    int periodicElement;

    double diameter;
    oofem::FloatArray direction_vector;

    double length;


    int number;




    // data from Discretization
    std::vector< int >m_ReinforcementPoints;
    std::vector< int >m_DelaunayVertices; // Delaunay corresponding vertices (to create bond link)
    std::vector< int >m_DelaunayVerticesTwo; // Delaunay corresponding vertices (to create bond link)
    std::vector< int >m_IntersectionPoints; // Intersection points with Voroinoi cells

    std::vector< double >m_listOfL_end; // list of L_end for use of Link elements


    // vertices associated to endpoints
    Vertex *PointP;
    Vertex *PointQ;



    // tool functions
    void findDelaunayVertex(); // to find nearest Delaunay vertex : called at the first iteration of the discretisation process
    void findVoronoiVertex();     // to find nearest Delaunay vertex : called at the first iteration of the discretisation process
    void findIntersect();//find the intersection point with the cell and the next delaunay vertex to consider
    void placeReinforcementNode();




public:


    /**
     * Constructor. Creates a fibre belonging to grid.
     * @param n fibre number in grid aGrid
     * @param aGrid grid to which node belongs
     */

    Fibre(int n, Grid *aGrid);                   // constructor
    /// Destructor.
    ~Fibre() override = default;

    Fibre *ofType();

    void initializeFrom(ConverterInputRecord &ir);
    void discretizeYouself();//


    //accessors
    int giveNumberReinforcementNode(int i);
    int NbOfReinfNodes();
    int giveNumberDelaunayNode(int i);
    int NbOfDelNodes();
    int giveNumberIntersectionPoint(int i);
    int NbOfIntersectionPoints();
    oofem::FloatArray giveDirVector() { return direction_vector; }
    double giveDiameter() { return diameter; }
    double giveL_end(int i);
    int NbOfL_end();
    int giveNumber() { return number; };


    static Vertex *reinforcementLocalizer(oofem::FloatArray coord, Grid *agrid, double TOL = 1.0e-14);
    // find the reinforcement node localized around the specified coordinates (with TOL as tolerance)
    // error if no node found
    // TOL must be chosen enough small, or the returned node will be the first found

    static double computedistance(oofem::FloatArray coordsOne, oofem::FloatArray coordsTwo);
};


#endif // FIBRE_H_INCLUDED
