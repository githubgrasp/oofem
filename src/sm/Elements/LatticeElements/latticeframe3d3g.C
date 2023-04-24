/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2019   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "domain.h"
#include "latticeframe3d3g.h"
#include "../sm/Materials/LatticeMaterials/latticematstatus.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "mathfem.h"
#include "latticeframe3d.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"
#include "sm/CrossSections/latticecrosssection.h"
#include "engngm.h"

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "../sm/Materials/structuralmaterial.h"
#endif

namespace oofem {
REGISTER_Element(LatticeFrame3d3g);

LatticeFrame3d3g::LatticeFrame3d3g(int n, Domain *aDomain) : LatticeFrame3d(n, aDomain)
{
    numberOfDofMans = 2;
}

LatticeFrame3d3g::~LatticeFrame3d3g()
{}

 
void
LatticeFrame3d3g::computeBDmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the strain matrix of the receiver.
{
    //Assemble Bmatrix (used to compute strains and rotations)
    answer.resize(6, 12);
    answer.zero();
    TimeStep *tStep = this->domain->giveEngngModel()->giveCurrentStep();

    this->length = computeLength();
    double tol=1.e-16;
    FloatArray u;
    //    this->computeVectorOf(VM_Total, tStep, u);
    this->computeVectorOf(VM_Incremental, tStep, u);
    double l1 = this->length*(1.-this->s)/2;
    double l2 = this->length*(1.+this->s)/2;
    //Normal displacement jump in x-direction
    //First node
    answer.at(1, 1) = -1.;
    answer.at(1, 2) = 0.;
    answer.at(1, 3) = 0.;
    answer.at(1, 4) = 0.;
    if ( fabs(u.at(5)) <= tol ) {
        answer.at(1, 5) =  0.;
    } else {
        answer.at(1, 5) =  l1*(1. - cos(u.at(5))*cos(u.at(6)))/u.at(5);
    }
    if ( fabs(u.at(6)) <= tol ) {
        answer.at(1, 6) =  0.;
    } else {
        answer.at(1, 6) =  l1*(1. - cos(u.at(5))*cos(u.at(6)))/u.at(6);
    }
    //
    //Second node
    answer.at(1, 7) = 1.;
    answer.at(1, 8) = 0.;
    answer.at(1, 9) = 0.;
    answer.at(1, 10) = 0.;
    if ( fabs(u.at(11)) <= tol ) {
        answer.at(1, 11) =  0;
    } else {
        answer.at(1, 11) =  l2*(1.-cos(u.at(11))*cos(u.at(12)))/u.at(11);
    }
    if ( fabs(u.at(12)) <= tol ) {
        answer.at(1, 12) =  0.;
    } else {
        answer.at(1, 12) =  l2*(1.-cos(u.at(11))*cos(u.at(12)))/u.at(12);
    }

    //Shear displacement jump in y-plane
    //first node
    answer.at(2, 1) = 0.;
    answer.at(2, 2) = -1.;
    answer.at(2, 3) =  0.;
    answer.at(2, 4) = 0.;
    answer.at(2, 5) = 0;
    if ( fabs(u.at(6)) <= tol ) {
        answer.at( 2, 6 ) = l1;
    } else {
        answer.at( 2, 6 ) = -(cos( u.at( 4 ) )*sin( u.at( 6 ) )+sin( u.at( 4 ) )*sin( u.at( 5 ) )*cos( u.at( 6 ) ) )* l1/u.at(6);
    }
    //Second node
    answer.at(2, 7) = 0.;
    answer.at(2, 8) = 1.;
    answer.at(2, 9) =  0.;
    answer.at(2, 10) = 0.;
    answer.at(2, 11) = 0;
    if ( fabs(u.at(12)) <= tol ) {
        answer.at(2, 12) =  l2;
    } else {
        answer.at( 2, 12 ) = -(cos( u.at( 10 ) )*sin( u.at( 12 ) )+sin( u.at( 10 ) )*sin( u.at( 11 ) )*cos( u.at( 12 ) ) ) * l2/u.at(12);
    }

    //Shear displacement jump in z-plane
    //first node
    answer.at(3, 1) = 0.;
    answer.at(3, 2) = 0.;
    answer.at(3, 3) = -1.;
    answer.at(3, 4) = 0.;
    if ( fabs(u.at(5)) <= tol ) {
        answer.at(3, 5) =  l1;
    } else {
        answer.at( 3, 5 ) = -(sin( u.at( 4 ) )*sin( u.at( 6 ) )-cos( u.at( 4 ) )*sin( u.at( 5 ) )*cos( u.at( 6 ) )) * l1/u.at(5);
    }
    answer.at(3, 6) = 0.;
    //Second node
    answer.at(3, 7) = 0.;
    answer.at(3, 8) = 0.;
    answer.at(3, 9) =  1.;
    answer.at(3, 10) = 0.;
    if ( fabs(u.at(11)) <= tol ) {
        answer.at(3, 11) =  l2;
    } else {
        answer.at( 3, 11 ) = -(sin( u.at( 10 ) )*sin( u.at( 12 ) )-cos( u.at( 10 ) )*sin( u.at( 11 ) )*cos( u.at( 12 ) )) *l2/u.at(11);
    }
    answer.at(3, 12) = 0.;

    //Rotation around x-axis
    //First node
    answer.at(4, 1) = 0.;
    answer.at(4, 2) = 0;
    answer.at(4, 3) = 0.;
    answer.at(4, 4) = -1.;
    answer.at(4, 5) = 0.;
    answer.at(4, 6) = 0.;
    //Second node
    answer.at(4, 7) = 0.;
    answer.at(4, 8) = 0.;
    answer.at(4, 9) = 0.;
    answer.at(4, 10) = 1.;
    answer.at(4, 11) = 0.;
    answer.at(4, 12) = 0.;

    //Rotation around y-axis
    //First node
    answer.at(5, 1) = 0.;
    answer.at(5, 2) = 0.;
    answer.at(5, 3) = 0.;
    answer.at(5, 4) = 0.;
    answer.at(5, 5) = -1.;
    answer.at(5, 6) = 0.;
    //Second node
    answer.at(5, 7) = 0.;
    answer.at(5, 8) = 0.;
    answer.at(5, 9) =  0.;
    answer.at(5, 10) = 0.;
    answer.at(5, 11) = 1.;
    answer.at(5, 12) = 0.;

    //Rotation around z-axis
    //First node
    answer.at(6, 1) = 0.;
    answer.at(6, 2) = 0.;
    answer.at(6, 3) = 0.;
    answer.at(6, 4) = 0.;
    answer.at(6, 5) = 0.;
    answer.at(6, 6) = -1.;
    //Second node
    answer.at(6, 7) = 0.;
    answer.at(6, 8) = 0.;
    answer.at(6, 9) =  0.;
    answer.at(6, 10) = 0.;
    answer.at(6, 11) = 0.;
    answer.at(6, 12) = 1.;

    return;
}

void
LatticeFrame3d3g::computeBFmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
{
    //Assemble BFmatrix (used to compute forces )
    answer.resize(12, 6);
    answer.zero();
    double tol=1.e-16;

    this->length = computeLength();

    TimeStep *tStep = this->domain->giveEngngModel()->giveCurrentStep();

    FloatArray u;
    this->computeVectorOf(VM_Incremental, tStep, u);
    //    this->computeVectorOf(VM_Total, tStep, u);
    double l1 = this->length*(1.-this->s)/2;
    double l2 = this->length*(1.+this->s)/2;

    //Nx1
    answer.at(1, 1) = -1.;
    answer.at(1, 2) = 0.;
    answer.at(1, 3) = 0.;
    answer.at(1, 4) = 0.;
    answer.at(1, 5) = 0.;
    answer.at(1, 6) = 0.;

    //Fy1
    answer.at(2, 1) = 0.;
    answer.at(2, 2) = -1.;
    answer.at(2, 3) = 0.;
    answer.at(2, 4) = 0.;
    answer.at(2, 5) = 0;
    answer.at(2, 6) = 0.;

    //Fz1
    answer.at(3, 1) = 0.;
    answer.at(3, 2) = 0.;
    answer.at(3, 3) = -1.;
    answer.at(3, 4) = 0.;
    answer.at(3, 5) = 0.;
    answer.at(3, 6) = 0.;

    //Mx1
    answer.at(4, 1) = 0.;
    answer.at(4, 2) = -sin(u.at(5))*l1;
    answer.at(4, 3) = -sin(u.at(6))*l1;
    answer.at(4, 4) = -1.;
    answer.at(4, 5) = 0.;
    answer.at(4, 6) = 0.;

    //My1
    answer.at(5, 1) = (sin( u.at( 4 ) )*sin( u.at( 6 ) )-cos( u.at( 4 ) )*sin( u.at( 5 ) )*cos( u.at( 6 ) ))*l1;
    answer.at(5, 2) = 0.;
    answer.at(5, 3) = cos(u.at(5))*cos(u.at(6))*l1;
    answer.at(5, 4) = 0.;
    answer.at(5, 5) = -1.;
    answer.at(5, 6) = 0.;

    //Mz1
    answer.at(6, 1) = (cos( u.at( 4 ) )*sin( u.at( 6 ) )+sin( u.at( 4 ) )*sin( u.at( 5 ) )*cos( u.at( 6 ) ) )* l1;
    answer.at(6, 2) = -cos(u.at(5))*cos(u.at(6))*l1;
    answer.at(6, 3) = 0.;
    answer.at(6, 4) = 0.;
    answer.at(6, 5) = 0.;
    answer.at(6, 6) = -1.;


    //Nx2
    answer.at(7, 1) = 1.;
    answer.at(7, 2) = 0.;
    answer.at(7, 3) = 0.;
    answer.at(7, 4) = 0.;
    answer.at(7, 5) = 0.;
    answer.at(7, 6) = 0.;

    //Fy2
    answer.at(8, 1) = 0.;
    answer.at(8, 2) = 1.;
    answer.at(8, 3) =  0.;
    answer.at(8, 4) = 0.;
    answer.at(8, 5) = 0;
    answer.at(8, 6) = 0.;

    //Fz2
    answer.at(9, 1) = 0.;
    answer.at(9, 2) = 0.;
    answer.at(9, 3) =  1.;
    answer.at(9, 4) = 0.;
    answer.at(9, 5) = 0.;
    answer.at(9, 6) = 0.;

    //Mx2
    answer.at(10, 1) = 0.;
    answer.at(10, 2) = sin(u.at(11))*l2;
    answer.at(10, 3) = sin(u.at(12))*l2;
    answer.at(10, 4) = 1.;
    answer.at(10, 5) = 0.;
    answer.at(10, 6) = 0.;

    //My2
    answer.at(11, 1) = (sin( u.at( 10 ) )*sin( u.at( 12 ) )-cos( u.at( 10 ) )*sin( u.at( 11 ) )*cos( u.at( 12 ) ))*l2;
    answer.at(11, 2) = 0.;
    answer.at(11, 3) = cos(u.at(11))*cos(u.at(12))*l2;
    answer.at(11, 4) = 0.;
    answer.at(11, 5) = 1.;
    answer.at(11, 6) = 0.;

    //Mz2
    answer.at(12, 1) = (cos( u.at( 10 ) )*sin( u.at( 12 ) )+sin( u.at( 10 ) )*sin( u.at( 11 ) )*cos( u.at( 12 ) ) )*l2;
    answer.at(12, 2) = -cos(u.at(11))*cos(u.at(12))*l2;
    answer.at(12, 3) =  0.;
    answer.at(12, 4) = 0.;
    answer.at(12, 5) = 0.;
    answer.at(12, 6) = 1.;

    return;
}
void
LatticeFrame3d3g::computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer =  static_cast< LatticeCrossSection * >( this->giveCrossSection() )->give3dFrameStiffnessMatrix(rMode, gp, tStep);
}
int
LatticeFrame3d3g::computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    double ksi, n1, n2;

    ksi = lcoords.at(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.resize(3);
    answer.at(1) = n1 * this->giveNode(1)->giveCoordinate(1) + n2 * this->giveNode(2)->giveCoordinate(1);
    answer.at(2) = n1 * this->giveNode(1)->giveCoordinate(2) + n2 * this->giveNode(2)->giveCoordinate(2);
    answer.at(3) = n1 * this->giveNode(1)->giveCoordinate(3) + n2 * this->giveNode(2)->giveCoordinate(3);

    return 1;
}
void
LatticeFrame3d3g::computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
    TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    FloatMatrix d, bi, bj, bjt, dbj, dij, bf, g, b, bt;
    FloatArray u;
    this->length = computeLength();
    answer.resize(12, 12);
    answer.zero();

    this->computeBDmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), bj);
    this->LatticeFrame3d::computeBmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), b);
    //printf("Bmatrix/n");
    // bj.printYourself();
    this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);
    computeBFmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), bf);
    //  printf("BFmatrix/n");
    // bf.printYourself();
    double l1 = this->length*(1.-this->s)/2;
    double l2 = this->length*(1.+this->s)/2;
    this->computeVectorOf(VM_Incremental, tStep, u);
    dbj.beProductOf(d, bj);
    dbj.times(1. / length);
    bjt.beTranspositionOf(bj);
    bt.beTranspositionOf(b);

   // answer.beProductOf(bt, dbj);
    //answer.beProductOf(bf, dbj);
    //Axial 1
    answer.at(1,1)= d.at(1,1);
    answer.at(1,2)= 0.;
    answer.at(1,3)= 0.;
    answer.at(1,4)= 0.;
    answer.at(1,5)= -d.at(1, 1)*l1*sin(u.at(5))*cos(u.at(6));
    answer.at(1,6)= -d.at(1, 1)*l1*cos(u.at(5))*sin(u.at(6));
    answer.at(1,7)= -d.at(1,1);
    answer.at(1,8)= 0.;
    answer.at(1,9)= 0.;
    answer.at(1,10)= 0.;
    answer.at(1,11)= -d.at(1, 1)*l2*sin(u.at(11))*cos(u.at(12));
    answer.at(1,12)= -d.at(1, 1)*l2*cos(u.at(11))*sin(u.at(12));

    //Shear Y 1
    answer.at(2,1)= 0;
    answer.at(2,2)= d.at(2,2);
    answer.at(2,3)= 0.;
    answer.at(2,4)= -d.at(2,2)*(sin(u.at(4))*sin(u.at(6))-cos(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1;
    answer.at(2,5)= -d.at(2, 2)*(-sin(u.at(4))*cos(u.at(5))*cos(u.at(6)))*l1;
    answer.at(2,6)= -d.at(2, 2)*(-cos(u.at(4))*cos(u.at(6))+sin(u.at(4))*sin(u.at(5))*sin(u.at(6)))*l1;
    answer.at(2,7)= 0.;
    answer.at(2,8)= -d.at(2,2);
    answer.at(2,9)= 0.;
    answer.at(2,10)= -d.at(2,2)*(sin(u.at(10))*sin(u.at(12))-cos(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2;
    answer.at(2,11)= -d.at(2, 2)*(-sin(u.at(10))*cos(u.at(11))*cos(u.at(12)))*l2;
    answer.at(2,12)= -d.at(2, 2)*(-cos(u.at(10))*cos(u.at(12))+sin(u.at(10))*sin(u.at(11))*sin(u.at(12)))*l2;


    //Shear Z 1
    answer.at(3,1)= 0;
    answer.at(3,2)= 0.;
    answer.at(3,3)= d.at(3,3);
    answer.at(3,4)= d.at(3,3)*(cos(u.at(4))*sin(u.at(6))+sin(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1;
    answer.at(3,5)= -d.at(3, 3)*(cos(u.at(4))*cos(u.at(5))*cos(u.at(6)))*l1;
    answer.at(3,6)= d.at(3, 3)*(sin(u.at(4))*cos(u.at(6))+cos(u.at(4))*sin(u.at(5))*sin(u.at(6)))*l1;
    answer.at(3,7)= 0.;
    answer.at(3,8)= 0;
    answer.at(3,9)= -d.at(3,3);
    answer.at(3,10)= d.at(3,3)*(cos(u.at(10))*sin(u.at(12))+sin(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2;
    answer.at(3,11)= -d.at(3, 3)*(cos(u.at(10))*cos(u.at(11))*cos(u.at(12)))*l2;
    answer.at(3,12)= d.at(3, 3)*(sin(u.at(10))*cos(u.at(12))+cos(u.at(10))*sin(u.at(11))*sin(u.at(12)))*l2;

    // Mx 1
    answer.at(4,1)= 0;
    answer.at(4,2)= 0.;
    answer.at(4,3)= 0.;
    answer.at(4,4)= d.at(4,4);
    answer.at(4,5)= 0.;
    answer.at(4,6)= 0.;
    answer.at(4,7)= 0.;
    answer.at(4,8)= 0;
    answer.at(4,9)= 0.;
    answer.at(4,10)= -d.at(4,4);
    answer.at(4,11)= 0.;
    answer.at(4,12)= 0.;

    // My 1
    answer.at(5,1)= -d.at(1,1)*(sin(u.at(4))*sin(u.at(6))-cos(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1;
    answer.at(5,2)= 0.;
    answer.at(5,3)= -d.at(3,3)*cos(u.at(5))*cos(u.at(6))*l1;
    answer.at(5,4)= d.at(1,1)*((u.at(7)-u.at(1))+(1-cos(u.at(11))*cos(u.at(12)))*l2+(1-cos(u.at(5))*cos(u.at(6)))*l1)*(cos(u.at(4))*sin(u.at(6))+sin(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1+d.at(3,3)*(-cos(u.at(4))*sin(u.at(6))*l1*cos(u.at(5))*cos(u.at(6))*l1-sin(u.at(4))*sin(u.at(5))*cos(u.at(6))*l1*cos(u.at(5))*cos(u.at(6))*l1);
    answer.at(5,5)= d.at(1,1)*((sin(u.at(4))*sin(u.at(6))-cos(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1*sin(u.at(5))*cos(u.at(6))*l1+((u.at(7)-u.at(1))+(1-cos(u.at(11))*cos(u.at(12)))*l2+(1-cos(u.at(5))*cos(u.at(6)))*l1)*(-cos(u.at(4))*cos(u.at(5))*cos(u.at(6))*l1))+d.at(3,3)*(-(u.at(9)-u.at(3))*sin(u.at(5))*cos(u.at(6))*l1+sin(u.at(10))*sin(u.at(12))*l2*sin(u.at(5))*cos(u.at(6))*l1-cos(u.at(10))*sin(u.at(11))*cos(u.at(12))*l2*sin(u.at(5))*cos(u.at(6))*l1+sin(u.at(4))*sin(u.at(6))*l1*sin(u.at(5))*cos(u.at(6))*l1+cos(u.at(4))*cos(2*u.at(5))*cos(u.at(6))*l1*cos(u.at(6))*l1)+d.at(5,5);
    answer.at(5,6)= d.at(1,1)*((sin(u.at(4))*sin(u.at(6))-cos(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1*cos(u.at(5))*sin(u.at(6))*l1+((u.at(7)-u.at(1))+(1-cos(u.at(11))*cos(u.at(12)))*l2+(1-cos(u.at(5))*cos(u.at(6)))*l1)*(sin(u.at(4))*cos(u.at(6))+cos(u.at(4))*sin(u.at(5))*sin(u.at(6)))*l1)+d.at(3,3)*(-(u.at(9)-u.at(3))*cos(u.at(5))*sin(u.at(6))*l1+sin(u.at(10))*sin(u.at(12))*l2*cos(u.at(5))*sin(u.at(6))*l1-cos(u.at(10))*sin(u.at(11))*cos(u.at(12))*l2*cos(u.at(5))*sin(u.at(6))*l1-sin(u.at(4))*cos(2*u.at(6))*l1*cos(u.at(5))*l1-cos(u.at(4))*sin(u.at(5))*2*sin(u.at(6))*cos(u.at(6))*l1*cos(u.at(5))*l1);
    answer.at(5,7)= d.at(1,1)*(sin(u.at(4))*sin(u.at(6))-cos(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1;
    answer.at(5,8)= 0;
    answer.at(5,9)= d.at(3,3)*cos(u.at(5))*cos(u.at(6))*l1;
    answer.at(5,10)= d.at(3,3)*(-cos(u.at(10))*sin(u.at(12))*l2*cos(u.at(5))*cos(u.at(6))*l1-sin(u.at(10))*sin(u.at(11))*cos(u.at(12))*l2*cos(u.at(5))*cos(u.at(6))*l1);
    answer.at(5,11)= d.at(1,1)*((sin(u.at(4))*sin(u.at(6))-cos(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1*sin(u.at(11))*cos(u.at(12))*l2)+d.at(3,3)*(cos(u.at(10))*cos(u.at(11))*cos(u.at(12))*l2*cos(u.at(5))*cos(u.at(6))*l1)-d.at(5,5);
    answer.at(5,12)= d.at(1,1)*((sin(u.at(4))*sin(u.at(6))-cos(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1*cos(u.at(11))*sin(u.at(12))*l2)+d.at(3,3)*(-sin(u.at(10))*cos(u.at(12))*l2*cos(u.at(5))*cos(u.at(6))*l1-cos(u.at(10))*sin(u.at(11))*sin(u.at(12))*l2*cos(u.at(5))*cos(u.at(6))*l1);

    // Mz 1
    answer.at(6,1)= -d.at(1, 1)*(cos(u.at(4))*sin(u.at(6))+sin(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1;
    answer.at(6,2)= d.at(2, 2)*cos(u.at(5))*cos(u.at(6))*l1;
    answer.at(6,3)= 0;
    answer.at(6,4)= d.at(1, 1)*((u.at(7)-u.at(1))+(1-cos(u.at(11))*cos(u.at(12)))*l2+(1-cos(u.at(5))*cos(u.at(6)))*l1)*(-sin(u.at(4))*sin(u.at(6))+cos(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1-d.at(2, 2)*(sin(u.at(4))*sin(u.at(6))*l1*cos(u.at(5))*cos(u.at(6))*l1-cos(u.at(4))*sin(u.at(5))*cos(u.at(6))*l1*cos(u.at(5))*cos(u.at(6))*l1);
    answer.at(6,5)= d.at(1,1)*((cos(u.at(4))*sin(u.at(6))+sin(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1*sin(u.at(5))*cos(u.at(6))*l1+((u.at(7)-u.at(1))+(1-cos(u.at(11))*cos(u.at(12)))*l2+(1-cos(u.at(5))*cos(u.at(6)))*l1)*sin(u.at(4))*cos(u.at(5))*cos(u.at(6))*l1)-d.at(2,2)*(-(u.at(8)-u.at(2))*l1*sin(u.at(5))*cos(u.at(6))+cos(u.at(10))*sin(u.at(12))*l2*sin(u.at(5))*cos(u.at(6))*l1+sin(u.at(10))*sin(u.at(11))*cos(u.at(12))*l2*sin(u.at(5))*cos(u.at(6))*l1+cos(u.at(4))*sin(u.at(6))*l1*sin(u.at(5))*cos(u.at(6))*l1-sin(u.at(4))*cos(u.at(5)*2)*cos(u.at(6))*l1*cos(u.at(6))*l1);
    answer.at(6,6)= d.at(1, 1)*((cos(u.at(4))*sin(u.at(6))+sin(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1*cos(u.at(5))*sin(u.at(6))*l1+((u.at(7)-u.at(1))+(1-cos(u.at(11))*cos(u.at(12)))*l2+(1-cos(u.at(5))*cos(u.at(6)))*l1)*(cos(u.at(4))*cos(u.at(6))-sin(u.at(4))*sin(u.at(5))*sin(u.at(6)))*l1)- d.at(2, 2)*((-u.at(8)+u.at(2))*cos(u.at(5))*sin(u.at(6))*l1+l1*l2*cos(u.at(10))*sin(u.at(12))*cos(u.at(5))*sin(u.at(6))+l1*l2*sin(u.at(10))*sin(u.at(11))*cos(u.at(12))*cos(u.at(5))*sin(u.at(6))-cos(u.at(4))*cos(u.at(5))*cos(2*u.at(6))*l1*l1+sin(u.at(4))*sin(u.at(5))*cos(u.at(5))*cos(u.at(6))*sin(u.at(6))*2*l1*l1)+d.at(6, 6);
    answer.at(6,7)= d.at(1, 1)*(cos(u.at(4))*sin(u.at(6))+sin(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1;;
    answer.at(6,8)= -d.at(2, 2)*cos(u.at(5))*cos(u.at(6))*l1;
    answer.at(6,9)= 0.;
    answer.at(6,10)= -d.at(2,2)*(sin(u.at(10))*sin(u.at(12))*l2*cos(u.at(5))*cos(u.at(6))*l1-cos(u.at(10))*sin(u.at(11))*cos(u.at(12))*l2*cos(u.at(5))*cos(u.at(6))*l1);
    answer.at(6,11)= d.at(1, 1)*(cos(u.at(4))*sin(u.at(6))+sin(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1*sin(u.at(11))*cos(u.at(12))*l2-d.at(2,2)*(-sin(u.at(10))*cos(u.at(11))*cos(u.at(12))*l2*cos(u.at(5))*cos(u.at(6))*l1);
    answer.at(6,12)= d.at(1, 1)*((cos(u.at(4))*sin(u.at(6))+sin(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1*cos(u.at(11))*sin(u.at(12))*l2)-d.at(2, 2)*(-cos(u.at(10))*cos(u.at(12))*l2*cos(u.at(5))*cos(u.at(6))*l1+sin(u.at(10))*sin(u.at(11))*sin(u.at(12))*l2*cos(u.at(5))*cos(u.at(6))*l1)-d.at(6, 6);

    //Axial 2
    answer.at(7,1)= -d.at(1,1);
    answer.at(7,2)= 0.;
    answer.at(7,3)= 0.;
    answer.at(7,4)= 0.;
    answer.at(7,5)= d.at(1, 1)*l1*sin(u.at(5))*cos(u.at(6));
    answer.at(7,6)= d.at(1, 1)*l1*cos(u.at(5))*sin(u.at(6));
    answer.at(7,7)= d.at(1,1);
    answer.at(7,8)= 0.;
    answer.at(7,9)= 0.;
    answer.at(7,10)= 0.;
    answer.at(7,11)= d.at(1, 1)*l2*sin(u.at(11))*cos(u.at(12));
    answer.at(7,12)= d.at(1, 1)*l2*cos(u.at(11))*sin(u.at(12));

    //Shear Y 2
    answer.at(8,1)= 0;
    answer.at(8,2)= -d.at(2,2);
    answer.at(8,3)= 0.;
    answer.at(8,4)= d.at(2,2)*(sin(u.at(4))*sin(u.at(6))-cos(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1;
    answer.at(8,5)= d.at(2, 2)*(-sin(u.at(4))*cos(u.at(5))*cos(u.at(6)))*l1;
    answer.at(8,6)= d.at(2, 2)*(-cos(u.at(4))*cos(u.at(6))+sin(u.at(4))*sin(u.at(5))*sin(u.at(6)))*l1;
    answer.at(8,7)= 0.;
    answer.at(8,8)= d.at(2,2);
    answer.at(8,9)= 0.;
    answer.at(8,10)= d.at(2,2)*(sin(u.at(10))*sin(u.at(12))-cos(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2;
    answer.at(8,11)= d.at(2, 2)*(-sin(u.at(10))*cos(u.at(11))*cos(u.at(12)))*l2;
    answer.at(8,12)= d.at(2, 2)*(-cos(u.at(10))*cos(u.at(12))+sin(u.at(10))*sin(u.at(11))*sin(u.at(12)))*l2;

    //Shear Z 2
    answer.at(9,1)= 0;
    answer.at(9,2)= 0.;
    answer.at(9,3)= -d.at(3,3);
    answer.at(9,4)= d.at(3,3)*(cos(u.at(4))*sin(u.at(6))+sin(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1;
    answer.at(9,5)= d.at(3, 3)*(cos(u.at(4))*cos(u.at(5))*cos(u.at(6)))*l1;
    answer.at(9,6)= -d.at(3, 3)*(sin(u.at(4))*cos(u.at(6))+cos(u.at(4))*sin(u.at(5))*sin(u.at(6)))*l1;
    answer.at(9,7)= 0.;
    answer.at(9,8)= 0;
    answer.at(9,9)= d.at(3,3);
    answer.at(9,10)= d.at(3,3)*(cos(u.at(10))*sin(u.at(12))+sin(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2;
    answer.at(9,11)= d.at(3, 3)*(cos(u.at(10))*cos(u.at(11))*cos(u.at(12)))*l2;
    answer.at(9,12)= -d.at(3, 3)*(sin(u.at(10))*cos(u.at(12))+cos(u.at(10))*sin(u.at(11))*sin(u.at(12)))*l2;

    // Mx 2
    answer.at(10,1)= 0;
    answer.at(10,2)= 0.;
    answer.at(10,3)= 0.;
    answer.at(10,4)= -d.at(4,4);
    answer.at(10,5)= 0.;
    answer.at(10,6)= 0.;
    answer.at(10,7)= 0.;
    answer.at(10,8)= 0;
    answer.at(10,9)= 0.;
    answer.at(10,10)= d.at(4,4);
    answer.at(10,11)= 0.;
    answer.at(10,12)= 0.;

    // My 2
    answer.at(11,1)= -d.at(1,1)*(sin(u.at(10))*sin(u.at(12))-cos(u.at(12))*sin(u.at(11))*cos(u.at(12)))*l2;
    answer.at(11,2)= 0.;
    answer.at(11,3)= -d.at(3,3)*cos(u.at(11))*cos(u.at(12))*l2;
    answer.at(11,4)= d.at(3,3)*(-cos(u.at(4))*sin(u.at(6))*l1*cos(u.at(11))*cos(u.at(12))*l2-sin(u.at(4))*sin(u.at(5))*cos(u.at(6))*l1*cos(u.at(11))*cos(u.at(12))*l2);;
    answer.at(11,5)= d.at(1,1)*((sin(u.at(10))*sin(u.at(12))-cos(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2*sin(u.at(5))*cos(u.at(6))*l1)+d.at(3,3)*(cos(u.at(4))*cos(u.at(5))*cos(u.at(6))*l1*cos(u.at(11))*cos(u.at(12))*l2)-d.at(5,5);
    answer.at(11,6)= d.at(1,1)*((sin(u.at(10))*sin(u.at(12))-cos(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2*cos(u.at(5))*sin(u.at(6))*l1)+d.at(3,3)*(-sin(u.at(4))*cos(u.at(6))*l1*cos(u.at(11))*cos(u.at(12))*l2-cos(u.at(4))*sin(u.at(5))*sin(u.at(6))*l1*cos(u.at(11))*cos(u.at(12))*l2);
    answer.at(11,7)= d.at(1,1)*(sin(u.at(10))*sin(u.at(12))-cos(u.at(12))*sin(u.at(11))*cos(u.at(12)))*l2;
    answer.at(11,8)= 0;
    answer.at(11,9)= d.at(3,3)*cos(u.at(11))*cos(u.at(12))*l2;
    answer.at(11,10)= d.at(1,1)*((u.at(7)-u.at(1))+(1-cos(u.at(11))*cos(u.at(12)))*l2+(1-cos(u.at(5))*cos(u.at(6)))*l1)*(cos(u.at(10))*sin(u.at(12))+sin(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2+d.at(3,3)*(-cos(u.at(10))*sin(u.at(12))*l2*cos(u.at(11))*cos(u.at(12))*l2-sin(u.at(10))*sin(u.at(11))*cos(u.at(12))*l2*cos(u.at(11))*cos(u.at(12))*l2);;
    answer.at(11,11)= d.at(1,1)*((sin(u.at(10))*sin(u.at(12))-cos(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2*sin(u.at(11))*cos(u.at(12))*l2+((u.at(7)-u.at(1))+(1-cos(u.at(11))*cos(u.at(12)))*l2+(1-cos(u.at(5))*cos(u.at(6)))*l1)*(-cos(u.at(10))*cos(u.at(11))*cos(u.at(12))*l2))+d.at(3,3)*(-(u.at(9)-u.at(3))*sin(u.at(11))*cos(u.at(12))*l2+sin(u.at(10))*sin(u.at(12))*l2*sin(u.at(11))*cos(u.at(12))*l2+cos(u.at(10))*cos(2*u.at(11))*cos(u.at(12))*l2*cos(u.at(12))*l2 +sin(u.at(4))*sin(u.at(6))*l1*sin(u.at(11))*cos(u.at(12))*l2-cos(u.at(4))*cos(2*u.at(5))*cos(u.at(6))*l1*sin(u.at(11))*cos(u.at(12))*l2)+d.at(5,5);
    answer.at(11,12)= d.at(1,1)*((sin(u.at(10))*sin(u.at(12))-cos(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2*cos(u.at(11))*sin(u.at(12))*l2+((u.at(7)-u.at(1))+(1-cos(u.at(11))*cos(u.at(12)))*l2+(1-cos(u.at(5))*cos(u.at(6)))*l1)*(sin(u.at(10))*cos(u.at(12))+cos(u.at(10))*sin(u.at(11))*sin(u.at(12)))*l2)+d.at(3,3)*(-(u.at(9)-u.at(3))*cos(u.at(11))*sin(u.at(12))*l2-sin(u.at(10))*cos(2*u.at(12))*cos(u.at(11))*l2*l2-cos(u.at(10))*sin(u.at(11))*2*cos(u.at(12))*sin(u.at(12))*l2*cos(u.at(11))*l2+sin(u.at(4))*sin(u.at(6))*l1*cos(u.at(11))*sin(u.at(12))*l2-cos(u.at(4))*sin(u.at(5))*2*cos(u.at(6))*l1*cos(u.at(11))*sin(u.at(12))*l2);
    // Mz 2
    answer.at(12,1)= -d.at(1, 1)*(cos(u.at(10))*sin(u.at(12))+sin(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2;
    answer.at(12,2)= d.at(2, 2)*cos(u.at(11))*cos(u.at(12))*l2;
    answer.at(12,3)= 0.;
    answer.at(12,4)= -d.at(2,2)*(sin(u.at(4))*sin(u.at(6))*l1*cos(u.at(11))*cos(u.at(12))*l2-cos(u.at(4))*sin(u.at(5))*cos(u.at(6))*l1*cos(u.at(11))*cos(u.at(12))*l2);
    answer.at(12,5)= d.at(1, 1)*(cos(u.at(10))*sin(u.at(12))+sin(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2*sin(u.at(5))*cos(u.at(6))*l1-d.at(2,2)*(-sin(u.at(4))*cos(u.at(5))*cos(u.at(6))*l1*cos(u.at(11))*cos(u.at(12))*l2);
    answer.at(12,6)= d.at(1,1)*((cos(u.at(10))*sin(u.at(12))+sin(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2*cos(u.at(5))*sin(u.at(6))*l1)-d.at(2,2)*(-cos(u.at(4))*cos(u.at(6))*l1*cos(u.at(11))*cos(u.at(12))*l2+sin(u.at(4))*sin(u.at(5))*sin(u.at(6))*l1*cos(u.at(11))*cos(u.at(12))*l2)-d.at(6,6);
    answer.at(12,7)= d.at(1, 1)*(cos(u.at(10))*sin(u.at(12))+sin(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2;
    answer.at(12,8)= -d.at(2, 2)*cos(u.at(11))*cos(u.at(12))*l2;
    answer.at(12,9)= 0.;
    answer.at(12,10)= d.at(1, 1)*((u.at(7)-u.at(1))+(1-cos(u.at(11))*cos(u.at(12)))*l2+(1-cos(u.at(5))*cos(u.at(6)))*l1)*(-sin(u.at(10))*sin(u.at(12))+cos(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2-d.at(2, 2)*(sin(u.at(10))*sin(u.at(12))*l2*cos(u.at(11))*cos(u.at(12))*l2-cos(u.at(10))*sin(u.at(11))*cos(u.at(12))*l2*cos(u.at(11))*cos(u.at(12))*l2);
    answer.at(12,11)= d.at(1,1)*((cos(u.at(10))*sin(u.at(12))+sin(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2*sin(u.at(11))*cos(u.at(12))*l2+((u.at(7)-u.at(1))+(1-cos(u.at(11))*cos(u.at(12)))*l2+(1-cos(u.at(5))*cos(u.at(6)))*l1)*sin(u.at(10))*cos(u.at(11))*cos(u.at(12))*l2)-d.at(2,2)*(-(u.at(8)-u.at(2))*l2*sin(u.at(11))*cos(u.at(12))+cos(u.at(10))*sin(u.at(12))*l2*sin(u.at(11))*cos(u.at(12))*l2-sin(u.at(10))*cos(2*u.at(11))*cos(u.at(12))*l2*cos(u.at(12))*l2+cos(u.at(4))*sin(u.at(6))*l1*sin(u.at(11))*cos(u.at(12))*l2+sin(u.at(4))*sin(u.at(5))*cos(u.at(6))*l1*sin(u.at(11))*cos(u.at(12))*l2);
    answer.at(12,12)= d.at(1,1)*((cos(u.at(10))*sin(u.at(12))+sin(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2*(cos(u.at(11))*sin(u.at(12))*l2)+((u.at(7)-u.at(1))+(1-cos(u.at(11))*cos(u.at(12)))*l2+(1-cos(u.at(5))*cos(u.at(6)))*l1)*(cos(u.at(10))*cos(u.at(12))-sin(u.at(10))*sin(u.at(11))*sin(u.at(12)))*l2)-d.at(2,2)*(-(u.at(8)-u.at(2))*cos(u.at(11))*sin(u.at(12))*l2-cos(u.at(10))*cos(2*u.at(12))*l2*cos(u.at(11))*l2+sin(u.at(10))*sin(u.at(11))*2*cos(u.at(12))*sin(u.at(12))*l2*cos(u.at(11))*l2+cos(u.at(4))*sin(u.at(6))*l1*cos(u.at(11))*sin(u.at(12))*l2+sin(u.at(4))*sin(u.at(5))*cos(u.at(6))*l1*cos(u.at(11))*sin(u.at(12))*l2)+d.at(6,6);
   //
    answer.times(1. / length);
    printf("answerSM/n");
    answer.printYourself();
    return;
}
double LatticeFrame3d3g::giveArea() {
    FloatArray lc(1);
    return this->giveCrossSection()->give(CS_Area, lc, this);
}

double LatticeFrame3d3g::giveIy() {
    FloatArray lc(1);
    return this->giveCrossSection()->give(CS_InertiaMomentY, lc, this);
}

double LatticeFrame3d3g::giveIz() {
    FloatArray lc(1);
    return this->giveCrossSection()->give(CS_InertiaMomentZ, lc, this);
}

double LatticeFrame3d3g::giveIk() {
    FloatArray lc(1);
    return this->giveCrossSection()->give(CS_TorsionConstantX, lc, this);
}

double LatticeFrame3d3g::giveShearAreaY() {
    FloatArray lc(1);
    return this->giveCrossSection()->give(CS_ShearAreaY, lc, this);
}

double LatticeFrame3d3g::giveShearAreaZ() {
    FloatArray lc(1);
    return this->giveCrossSection()->give(CS_ShearAreaZ, lc, this);
}

void
LatticeFrame3d3g :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step tStep. The nature of these strains depends
// on the element's type.
{
    FloatMatrix b;
    FloatArray u;
    double l1 = this->length*(1.-this->s)/2;
    double l2 = this->length*(1.+this->s)/2;
    this->computeVectorOf(VM_Incremental, tStep, u);
    LatticeMaterialStatus *lmatStat = dynamic_cast< LatticeMaterialStatus * >( integrationRulesArray [ 0 ]->getIntegrationPoint(0)->giveMaterialStatus() );
    auto strain =lmatStat ->giveLatticeStrain();
    this->LatticeFrame3d3g::computeBDmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), b);
    //answer.beProductOf(b, u);
    answer.resize(6);
    answer.at(1)= (1-cos(u.at(11))*cos(u.at(12)))*l2+(1-cos(u.at(5))*cos(u.at(6)))*l1+u.at(7)-u.at(1);
    answer.at(2)= -(cos(u.at(10))*sin(u.at(12))+sin(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2-(cos(u.at(4))*sin(u.at(6))+sin(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1+u.at(8)-u.at(2);
    answer.at(3)= -(sin(u.at(10))*sin(u.at(12))-cos(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2-(sin(u.at(4))*sin(u.at(6))-cos(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1+u.at(9)-u.at(3);
    answer.at(4)= u.at(10)-u.at(4);
    answer.at(5)= u.at(11)-u.at(5);
    answer.at(6)= u.at(12)-u.at(6);
    answer.times(1./this->length);
    printf("STRAIN/n");
    answer.printYourself();
    answer += strain;
}
//
void
LatticeFrame3d3g::giveInternalForcesVector(FloatArray &answer,
    TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix b, bt, bf;
    FloatArray u, stress, strain;
    this->computeVectorOf(VM_Incremental, tStep, u);
    this->length   = computeLength();
    GaussPoint *gp = this->integrationRulesArray[0]->getIntegrationPoint( 0 );
    this->LatticeFrame3d3g::computeBFmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), bf);
    //    bt.beTranspositionOf( b );
    // Total stress
    this->LatticeFrame3d3g::computeStrainVector( strain, gp, tStep );
    this->computeStressVector( stress, strain, integrationRulesArray[0]->getIntegrationPoint( 0 ), tStep );
    //   totalInternalForces.beProductOf( bf, stress );
    // Old stresses
    LatticeMaterialStatus *lmatStat = dynamic_cast<LatticeMaterialStatus *>( integrationRulesArray[0]->getIntegrationPoint( 0 )->giveMaterialStatus() );
    auto oldStress= lmatStat->giveLatticeStress();

    auto oldInternalForces = lmatStat->giveInternalForces();
    double l1 = this->length*(1.-this->s)/2;
    double l2 = this->length*(1.+this->s)/2;
    FloatArray incrementalStress;
    incrementalStress.beDifferenceOf( stress, oldStress );
    //printf("incrementalStress/n");
    //incrementalStress.printYourself();
    answer.resize(12);
    answer.at(1)= -incrementalStress.at(1),
    answer.at(2)= -incrementalStress.at(2);
    answer.at(3)= -incrementalStress.at(3);
    answer.at(4)= -incrementalStress.at(4);
    answer.at(5)= incrementalStress.at(1)*(sin(u.at(4))*sin(u.at(6))-cos(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1+incrementalStress.at(3)*(cos(u.at(5))*cos(u.at(6)))*l1-incrementalStress.at(5);
    answer.at(6)= incrementalStress.at(1)*(cos(u.at(4))*sin(u.at(6))+sin(u.at(4))*sin(u.at(5))*cos(u.at(6)))*l1-incrementalStress.at(2)*(cos(u.at(5))*cos(u.at(6)))*l1-incrementalStress.at(6);
    answer.at(7)= incrementalStress.at(1);
    answer.at(8)= incrementalStress.at(2);
    answer.at(9)= incrementalStress.at(3);
    answer.at(10)= incrementalStress.at(4);
    answer.at(11)= incrementalStress.at(1)*(sin(u.at(10))*sin(u.at(12))-cos(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2+incrementalStress.at(3)*(cos(u.at(11))*cos(u.at(12)))*l2+incrementalStress.at(5);
    answer.at(12)= incrementalStress.at(1)*(cos(u.at(10))*sin(u.at(12))+sin(u.at(10))*sin(u.at(11))*cos(u.at(12)))*l2-incrementalStress.at(2)*(cos(u.at(11))*cos(u.at(12)))*l2+incrementalStress.at(6);        //answer.beProductOf( bf, stress );
    answer += oldInternalForces;
    printf("answerFORCE/n");
    answer.printYourself();
    lmatStat->letTempInternalForcesBe(answer);
}
bool
LatticeFrame3d3g::computeGtoLRotationMatrix(FloatMatrix &answer)
{
    FloatMatrix lcs;
    answer.resize(12, 12);
    answer.zero();

    this->LatticeFrame3d3g::giveLocalCoordinateSystem(lcs);
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(i + 3, j + 3) = lcs.at(i, j);
            answer.at(i + 6, j + 6) = lcs.at(i, j);
            answer.at(i + 9, j + 9) = lcs.at(i, j);
        }
    }

    return 1;
}
int
LatticeFrame3d3g::giveLocalCoordinateSystem(FloatMatrix &answer)
{
    FloatArray lx, ly, lz, help(3);
    FloatArray coordA, coordB;
    FloatArray uA(6),uAIncr(6), uB(6),uBIncr(6);
    IntArray dofid = {1,2,3,4,5,6};

    TimeStep *tStep = this->domain->giveEngngModel()->giveCurrentStep();

    Node *nodeA, *nodeB;
    nodeA = this->giveNode(1);
    nodeB = this->giveNode(2);


    coordA = nodeA->giveCoordinates();
    nodeA->giveUnknownVector(uA,dofid,VM_Total,tStep,false);
    nodeA->giveUnknownVector(uAIncr,dofid,VM_Incremental,tStep,false);
    //    nodeA->giveUnknownVectorOfType(uA, DisplacementVector, VM_Total, tStep);
    // nodeA->giveUnknownVectorOfType(uAIncr, DisplacementVector, VM_Incremental, tStep);
    for(int i=1;i<=3;i++){
      coordA.at(i) += uA.at(i)-uAIncr.at(i);
    }

    coordB = nodeB->giveCoordinates();
    nodeB->giveUnknownVector(uB,dofid,VM_Total,tStep,false);
    nodeB->giveUnknownVector(uBIncr,dofid,VM_Incremental,tStep,false);

    //      nodeB->giveUnknownVector(uB,dofid,VM_Total,tStep,false);
    for(int i=1;i<=3;i++){
      coordB.at(i) += uB.at(i)-uBIncr.at(i);
    }

    /* uA.at(1) = nodeA->giveUpdatedCoordinate(1,tStep,1.); */
    /* uA.at(2) = nodeA->giveUpdatedCoordinate(2,tStep,1.); */
    /* uA.at(3) = nodeA->giveUpdatedCoordinate(3,tStep,1.); */

    /* uB.at(1) = nodeB->giveUpdatedCoordinate(1,tStep,1.); */
    /* uB.at(2) = nodeB->giveUpdatedCoordinate(2,tStep,1.); */
    /* uB.at(3) = nodeB->giveUpdatedCoordinate(3,tStep,1.); */


    lx.beDifferenceOf(coordB, coordA );
    lx.normalize();

    if ( this->referenceNode ) {
        Node *refNode = this->giveDomain()->giveNode(this->referenceNode);
        help.beDifferenceOf(refNode->giveCoordinates(), nodeA->giveCoordinates() );

        lz.beVectorProductOf(lx, help);
        lz.normalize();
    } else if ( this->zaxis.giveSize() > 0 ) {
        lz = this->zaxis;
        lz.add(lz.dotProduct(lx), lx);
        lz.normalize();
    } else {
        FloatMatrix rot(3, 3);
        double theta = referenceAngle * M_PI / 180.0;

        rot.at(1, 1) = cos(theta) + pow(lx.at(1), 2) * ( 1 - cos(theta) );
        rot.at(1, 2) = lx.at(1) * lx.at(2) * ( 1 - cos(theta) ) - lx.at(3) * sin(theta);
        rot.at(1, 3) = lx.at(1) * lx.at(3) * ( 1 - cos(theta) ) + lx.at(2) * sin(theta);

        rot.at(2, 1) = lx.at(2) * lx.at(1) * ( 1 - cos(theta) ) + lx.at(3) * sin(theta);
        rot.at(2, 2) = cos(theta) + pow(lx.at(2), 2) * ( 1 - cos(theta) );
        rot.at(2, 3) = lx.at(2) * lx.at(3) * ( 1 - cos(theta) ) - lx.at(1) * sin(theta);

        rot.at(3, 1) = lx.at(3) * lx.at(1) * ( 1 - cos(theta) ) - lx.at(2) * sin(theta);
        rot.at(3, 2) = lx.at(3) * lx.at(2) * ( 1 - cos(theta) ) + lx.at(1) * sin(theta);
        rot.at(3, 3) = cos(theta) + pow(lx.at(3), 2) * ( 1 - cos(theta) );

        help.at(3) = 1.0;         // up-vector
        // here is ly is used as a temp var
        if ( fabs(lx.dotProduct(help) ) > 0.999 ) {  // Check if it is vertical
            ly = {
                0., 1., 0.
            };
        } else {
            ly.beVectorProductOf(lx, help);
        }
        lz.beProductOf(rot, ly);
        lz.normalize();
    }

    ly.beVectorProductOf(lz, lx);
    ly.normalize();

    answer.resize(3, 3);
    answer.zero();
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = lx.at(i);
        answer.at(2, i) = ly.at(i);
        answer.at(3, i) = lz.at(i);
    }

    return 1;
}

 

 void
LatticeFrame3d3g::giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, D_w, R_u, R_v, R_w
    };
}

void
LatticeFrame3d3g::initializeFrom(InputRecord &ir)
{
    LatticeFrame3d::initializeFrom(ir);



}
double
LatticeFrame3d3g::computeLength()
{
    double dx, dy, dz;
    Node *nodeA, *nodeB;

    if ( length == 0. ) {
        nodeA   = this->giveNode(1);
        nodeB   = this->giveNode(2);
        dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
        dz      = nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3);
        length  = sqrt(dx * dx + dy * dy + dz * dz);
    }

    return length;
}

double
LatticeFrame3d3g::computeCurrentLength()
{

  TimeStep *tStep = this->domain->giveEngngModel()->givePreviousStep();
  if(tStep->isTheFirstStep()){
    tStep = this->domain->giveEngngModel()->giveCurrentStep();
  }


  double dx, dy, dz;
  Node *nodeA, *nodeB;
  double currentLength;
  FloatArray uA(6), uB(6), coordA(3), coordB(3);
  IntArray dofid = {1,2,3,4,5,6};

  nodeA   = this->giveNode(1);
  coordA = nodeA->giveCoordinates();
  nodeA->giveUnknownVector(uA,dofid,VM_Total,tStep,false);
  for(int i=1;i<=3;i++){
    coordA.at(i) += uA.at(i);
  }

  nodeB   = this->giveNode(2);
  coordB = nodeB->giveCoordinates();
  nodeB->giveUnknownVector(uB,dofid,VM_Total,tStep,false);
  for(int i=1;i<=3;i++){
    coordB.at(i) += uB.at(i);
  }

  dx      = coordB.at(1) - coordA.at(1);
  dy      = coordB.at(2) - coordA.at(2);
  dz      = coordB.at(3) - coordA.at(3);
  currentLength  = sqrt(dx * dx + dy * dy + dz * dz);

  return currentLength;
}

void
LatticeFrame3d3g::computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    double density = static_cast< LatticeCrossSection * >( this->giveCrossSection() )->give('d', gp);
    double halfMass = density * computeVolumeAround(gp) / 2.;
    answer.resize(12, 12);
    answer.zero();
    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = halfMass;
    answer.at(7, 7) = answer.at(8, 8) = answer.at(9, 9) = halfMass;
}


} // end namespace oofem
