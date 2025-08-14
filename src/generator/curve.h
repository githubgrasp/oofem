#ifndef curve_h
#define curve_h


#include "gridcomponent.h"
#include "floatarray.h"
#include "intarray.h"


#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "oofemtxtinputrecord.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif


#define _IFT_Curve_vertices "vertices"
#define _IFT_Curve_refine "refine"
#define _IFT_Curve_normal "normal"


class Curve : public GridComponent
{
protected:
    /// Array storing nodal coordinates.
  oofem::IntArray vertices;
    int number;
    double refinement;
    int randomswitch;
  oofem::FloatArray normal;

public:

    /**
     * Constructor. Creates a node belonging to domain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    Curve(int n, Grid *aGrid);                   

    virtual ~Curve();

    /// Returns i-th vertex of curve.
    int      giveLocalVertex(int i);
    /// Returns pointer to curve vertex array.
  oofem::IntArray *giveLocalVertices() { return & vertices; }


  //    Curve *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Curve"; }

    //Random approach to generate point on curves
    int generatePoints();

    //Generate periodicPoints
    int generatePeriodicPoints();
  
  void mirrorShift(oofem::FloatArray& random, oofem::FloatArray& normal,oofem::FloatArray& specimenDimension,oofem::FloatArray& boundaries, int& vertexNumber, oofem::IntArray& periodicityFlag);
   
  //Give normal of curve
  void giveNormal(oofem::FloatArray &answer){answer = this->normal;}
    
  void initializeFrom(GeneratorInputRecord &ir);

  void printYourself();
};


#endif // node_h
