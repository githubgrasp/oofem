#ifndef curve_h
#define curve_h


#include "gridcomponent.h"
#include "flotarry.h"
#include "intarray.h"


#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "oofemtxtinputrecord.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

//class FloatArray;
//class IntArray;

class Curve : public GridComponent
{
protected:
    /// Array storing nodal coordinates.
    IntArray vertices;
    int number;
    double refinement;
    int randomswitch;
    FloatArray normal;

public:

    /**
     * Constructor. Creates a node belonging to domain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    Curve(int n, Grid *aGrid);                   // constructor
    /// Destructor.
    ~Curve();                                           // destructor

    /// Returns i-th vertex of curve.
    int      giveLocalVertex(int i);
    /// Returns pointer to curve vertex array.
    IntArray *giveLocalVertices() { return & vertices; }

    Vertex *giveGlobalVertex(int i, AList< Vertex > *vertexList);

    Curve *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Curve"; }

    //Random approach to generate point on curves
    int generatePoints();

    //Generate periodicPoints
    int generatePeriodicPoints();

    void mirrorShift(FloatArray& random, FloatArray& normal,FloatArray& specimenDimension,FloatArray& boundaries, int& vertexNumber, IntArray& periodicityFlag);
   

    //Give normal of curve
    void giveNormal(FloatArray &answer){answer = this->normal;}
    
    IRResultType initializeFrom(InputRecord *ir);
    //virtual IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();
};


#endif // node_h
