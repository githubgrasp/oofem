#include "ellipsoid.h"
#include "curve.h"
#include "vertex.h"
#include <iostream>


#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#endif

using namespace std;


Ellipsoid :: Ellipsoid(int n,Grid* aGrid) : Inclusion(n,aGrid) //, coordinates()
{
  this->number = n;
}

Ellipsoid :: ~Ellipsoid()
// Destructor.
{
}


void
Ellipsoid :: initializeFrom(ConverterInputRecord &ir)
// Gets from the source line from the data file all the data of the receiver.
{
    
    IR_GIVE_FIELD(ir, centre, _IFT_Ellipsoid_centre);

    IR_GIVE_FIELD(ir, radii, _IFT_Ellipsoid_radii); // Macro
    IR_GIVE_FIELD(ir, angles, _IFT_Ellipsoid_angles); // Macro
    refinement = 1.;
   
    IR_GIVE_OPTIONAL_FIELD(ir, refinement, _IFT_Ellipsoid_refine); // Macro
    itzThickness = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, itzThickness, _IFT_Ellipsoid_itz); // Macro
     
   // compute the matrix of the cartesian equation of the ellipsoid
    
    oofem::FloatMatrix matrixAloc(4,4);
    matrixAloc.at(1,1)=(1./pow(radii.at(1),2.));
    matrixAloc.at(2,2)=(1./pow(radii.at(2),2.));
    matrixAloc.at(3,3)=(1./pow(radii.at(3),2.));
    matrixAloc.at(4,4)=-1;
    
    
   oofem::FloatMatrix Rx(4,4);
   oofem::FloatMatrix Ry(4,4);
   oofem::FloatMatrix Rz(4,4);
   oofem::FloatMatrix T(4,4);

   T.at(4,4)=1;
   T.at(1,1)=1;
   T.at(2,2)=1;
   T.at(3,3)=1;
   T.at(4,1)=-centre.at(1);
   T.at(4,2)=-centre.at(2);
   T.at(4,3)=-centre.at(3);

   Rx.at(1,1)=1;
   Rx.at(4,4)=1;
   Rx.at(2,2)=cos(angles.at(1));
   Rx.at(2,3)=-sin(angles.at(1));
   Rx.at(3,2)=sin(angles.at(1));
   Rx.at(3,3)=cos(angles.at(1));
    
    Ry.at(2,2)=1;
    Ry.at(4,4)=1;
    Ry.at(1,1)=cos(angles.at(2));
    Ry.at(1,3)=sin(angles.at(2));
    Ry.at(3,1)=-sin(angles.at(2));
    Ry.at(3,3)=cos(angles.at(2));
    
    Rz.at(3,3)=1;
    Rz.at(4,4)=1;
    Rz.at(1,1)=cos(angles.at(3));
    Rz.at(1,2)=-sin(angles.at(3));
    Rz.at(2,1)=sin(angles.at(3));
    Rz.at(2,2)=cos(angles.at(3));

    oofem::FloatMatrix Rprov;
    Rprov.beProductOf(Rz,Ry);
    oofem::FloatMatrix R;
    R.beProductOf(Rprov,Rx);
    oofem::FloatMatrix tR;
    tR.beTranspositionOf(R);
    oofem::FloatMatrix tT;
    tT.beTranspositionOf(T);
    
    
    oofem::FloatMatrix tempM;
    tempM.beProductOf(T,R);
    oofem::FloatMatrix tempM1;
    tempM1.beProductOf(tempM,matrixAloc);
    oofem::FloatMatrix tempM2;
    tempM2.beProductOf(tempM1,tR);
    matrixA.beProductOf(tempM2,tT);
    
    if (itzThickness>0.)

    {
        oofem::FloatMatrix matrixAlocItZ(4,4);
        matrixAlocItZ.at(1,1)=(1./pow(radii.at(1)+itzThickness,2.));
        matrixAlocItZ.at(2,2)=(1./pow(radii.at(2)+itzThickness,2.));
        matrixAlocItZ.at(3,3)=(1./pow(radii.at(3)+itzThickness,2.));
        matrixAlocItZ.at(4,4)=-1;
        
        oofem::FloatMatrix tempM1ItZ;
        tempM1ItZ.beProductOf(tempM,matrixAlocItZ);
        oofem::FloatMatrix tempM2ItZ;
        tempM2ItZ.beProductOf(tempM1ItZ,tR);
        matrixA_wITZ.beProductOf(tempM2ItZ,tT);
                
    }
    
    return;

}


Ellipsoid *Ellipsoid :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    Ellipsoid *ellipsoid;

    ellipsoid = new Ellipsoid(number,grid);

    return ellipsoid;
}

int Ellipsoid::isInside(double coord_x,double coord_y,double coord_z)
{
double scalar(0);
oofem::FloatMatrix point(4,1);
point.at(1,1)=coord_x;
point.at(2,1)=coord_y;
point.at(3,1)=coord_z;
point.at(4,1)=1.;
oofem::FloatMatrix Tpoint;
Tpoint.beTranspositionOf(point);
oofem::FloatMatrix P1;
P1.beProductOf(matrixA,point);
oofem::FloatMatrix P2;
P2.beProductOf(Tpoint,P1);
scalar=P2.at(1,1);
if (scalar<=0.)
{
    return 1;
}
else
{   if (itzThickness>0.)
{
    oofem::FloatMatrix P1ItZ;
    P1ItZ.beProductOf(matrixA_wITZ,point);
    oofem::FloatMatrix P2ItZ;
    P2ItZ.beProductOf(Tpoint,P1ItZ);
    double scalarItZ(0);
    scalarItZ=P2ItZ.at(1,1);
    if (scalarItZ<=0.)
    {
        return 2 ;
    }
    else {
        return 0;
    }
    
}
    
    else{
         return 0;
    }
}

}
