#ifndef ELLIPSOID_H_INCLUDED
#define ELLIPSOID_H_INCLUDED


#include "grid.h"
#include "inclusion.h"

#include "flotarry.h"
#include "intarray.h"
#include "floatmatrix.h"

#include "converterinputrecord.h"

#include "converterdatareader.h"
#include "convertertxtdatareader.h"
#include "convertertxtinputrecord.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif


class Ellipsoid : public Inclusion
{

protected:
    /// Array storing nodal coordinates.
  oofem::FloatArray centre;
  oofem::FloatArray angles;
  oofem::FloatArray radii;
  oofem::FloatMatrix matrixA;
  oofem::FloatMatrix matrixA_wITZ; // matrix for ellipsoid+ITZ (if itzThickness>0)

  int number;
  double refinement;
  double itzThickness; //

public:

    /**
     * Constructor. Creates a node belonging to domain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    Ellipsoid(int n, Grid* aGrid);                      // constructor
    /// Destructor.
    ~Ellipsoid();                                           // destructor


    oofem::FloatArray giveCenter(){return this->centre;}
     oofem::FloatArray giveRadii(){return this->radii;}
      oofem::FloatArray giveAngles(){return this->angles;}
    double giveITZThickness(){return this->itzThickness;}
    void giveCentre(oofem::FloatArray& cent){cent = centre;}
    void giveRadii(oofem::FloatArray& rad){rad = radii;}



    Ellipsoid *ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Ellipsoid"; }

    IRResultType initializeFrom(InputRecord *ir);
    //virtual IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();

    // to know whether a point is inside the ellipsoid
    int isInside(double coord_x,double coord_y,double coord_z);
    // 1 means inside , 0 outisde , 2 in the ITZ (if specified)

};

#endif // ELLIPSOID_H_INCLUDED
