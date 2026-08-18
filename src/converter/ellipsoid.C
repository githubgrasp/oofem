#include "ellipsoid.h"
#include "curve.h"
#include "vertex.h"
#include <iostream>


#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
#endif

using namespace std;


Ellipsoid::Ellipsoid(int n, Grid *aGrid) : Inclusion(n, aGrid) //, coordinates()
{
    this->number = n;
}




Ellipsoid *Ellipsoid::ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    Ellipsoid *ellipsoid;

    ellipsoid = new Ellipsoid(number, grid);

    return ellipsoid;
}

int Ellipsoid::isInside(double coord_x, double coord_y, double coord_z)
{
    double scalar(0);
    oofem::FloatMatrix point(4, 1);
    point.at(1, 1) = coord_x;
    point.at(2, 1) = coord_y;
    point.at(3, 1) = coord_z;
    point.at(4, 1) = 1.;
    oofem::FloatMatrix Tpoint;
    Tpoint.beTranspositionOf(point);
    oofem::FloatMatrix P1;
    P1.beProductOf(matrixA, point);
    oofem::FloatMatrix P2;
    P2.beProductOf(Tpoint, P1);
    scalar = P2.at(1, 1);
    if ( scalar <= 0. ) {
        return 1;
    } else   {
        if ( itzThickness > 0. )          {
            oofem::FloatMatrix P1ItZ;
            P1ItZ.beProductOf(matrixA_wITZ, point);
            oofem::FloatMatrix P2ItZ;
            P2ItZ.beProductOf(Tpoint, P1ItZ);
            double scalarItZ(0);
            scalarItZ = P2ItZ.at(1, 1);
            if ( scalarItZ <= 0. ) {
                return 2;
            } else   {
                return 0;
            }
        }   else {
            return 0;
        }
    }
}
