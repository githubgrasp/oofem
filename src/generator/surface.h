#ifndef surface_h
#define surface_h


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


#define _IFT_Surface_curves "curves"
#define _IFT_Surface_refine "refine"
#define _IFT_Surface_normal "normal"

#define _IFT_Surface_boundaryflag "boundaryflag"
#define _IFT_Surface_boundaryshift "boundaryshift"




//class FloatArray;
//class IntArray;

class Surface : public GridComponent
{
protected:
    /// Array storing nodal coordinates.
  oofem::IntArray curves;
    int number;
    double refinement;
    double xedges, yedges, zedges;
    int boundaryFlag;
  oofem::FloatArray boundaryShift;
  oofem::FloatArray normal;
    
public:

    /**
     * Constructor. Creates a node belonging to domain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    Surface(int n, Grid *aGrid);                      // constructor
    /// Destructor.
    ~Surface();                                           // destructor

    /// Returns i-th vertex of curve.
    int      giveLocalCurve(int i);
    /// Returns pointer to curve vertex array.
  oofem::IntArray *giveLocalCurves() { return & curves; }

    int giveNumberOfLocalCurves();

    /// Define boundaries
  void defineBoundaries(oofem::FloatArray &boundaries);
    
  void giveNormal(oofem::FloatArray &answer){answer = this->normal;}
    
    //Random approach to generate points on surfaces including period shift and mirroring
    int generatePoints();
    
    //Shorten code by putting routine to shift and mirror point separately
  void mirrorShift(oofem::FloatArray& random, oofem::FloatArray& normal,oofem::FloatArray& specimenDimension,oofem::FloatArray& boundaries, int& vertexNumber, oofem::IntArray& periodicityFlag);
        
  Surface *ofType();

  // miscellaneous
  /// Returns class name of the receiver.
  const char *giveClassName() const { return "Surface"; }

  void initializeFrom(GeneratorInputRecord &ir);

  /// prints receiver state on stdout. Usefull for debuging.
  void         printYourself();
};


#endif // node_h
