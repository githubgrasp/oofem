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
#include "latticeframe3dg.h"
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
REGISTER_Element(LatticeFrame3dg);

LatticeFrame3dg::LatticeFrame3dg(int n, Domain *aDomain) : LatticeFrame3d(n, aDomain)
{
    numberOfDofMans = 2;
}

LatticeFrame3dg::~LatticeFrame3dg()
{}

 
void
LatticeFrame3dg::computeBDmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
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
        answer.at(1, 5) =  l1*(1. - cos(u.at(5)))/u.at(5);
    }
    if ( fabs(u.at(6)) <= tol ) {
        answer.at(1, 6) =  0.;
    } else {
        answer.at(1, 6) =  l1*(1. - cos(u.at(6)))/u.at(6);
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
        answer.at(1, 11) =  l2*(1.-cos(u.at(11)))/u.at(11);
    }
    if ( fabs(u.at(12)) <= tol ) {
        answer.at(1, 12) =  0.;
    } else {
        answer.at(1, 12) =  l2*(1.-cos(u.at(12)))/u.at(12);
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
        answer.at( 2, 6 ) = -sin( u.at( 6 ) ) * l1/u.at(6);
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
        answer.at( 2, 12 ) = -sin( u.at( 12 ) ) * l2/u.at(12);
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
        answer.at( 3, 5 ) = sin( u.at( 5 ) ) * l1/u.at(5);
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
        answer.at( 3, 11 ) = sin( u.at( 11 ) ) *l2/u.at(11);
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
LatticeFrame3dg::computeBFmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
{
    //Assemble BFmatrix (used to compute forces )
    answer.resize(12, 6);
    answer.zero();

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
    answer.at(5, 1) = sin(u.at(5))*l1;
    answer.at(5, 2) = 0.;
    answer.at(5, 3) = cos(u.at(5))*l1;
    answer.at(5, 4) = 0.;
    answer.at(5, 5) = -1.;
    answer.at(5, 6) = 0.;

    //Mz1
    answer.at(6, 1) = sin(u.at(6))*l1;
    answer.at(6, 2) = -cos(u.at(6))*l1;
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
    answer.at(11, 1) = sin(u.at(11))*l2;
    answer.at(11, 2) = 0.;
    answer.at(11, 3) = cos(u.at(11))*l2;
    answer.at(11, 4) = 0.;
    answer.at(11, 5) = 1.;
    answer.at(11, 6) = 0.;

    //Mz2
    answer.at(12, 1) = sin(u.at(12))*l2;
    answer.at(12, 2) = -cos(u.at(12))*l2;
    answer.at(12, 3) =  0.;
    answer.at(12, 4) = 0.;
    answer.at(12, 5) = 0.;
    answer.at(12, 6) = 1.;

    return;
}
void
LatticeFrame3dg::computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer =  static_cast< LatticeCrossSection * >( this->giveCrossSection() )->give3dFrameStiffnessMatrix(rMode, gp, tStep);
}
void
LatticeFrame3dg::computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
    TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    FloatMatrix d, bi, bj, bjt, dbj, dij, bf;

    this->length = computeLength();

    answer.resize(12, 12);
    answer.zero();
    this->computeBDmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), bj);
    this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);
    computeBFmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), bf);

    dbj.beProductOf(d, bj);
    dbj.times(1. / length);
    bjt.beTranspositionOf(bj);
    answer.beProductOf(bf, dbj);
    //printf("answer/n");
    //answer.printYourself();
    return;
}

void
LatticeFrame3dg :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step tStep. The nature of these strains depends
// on the element's type.
{
    FloatMatrix b;
    FloatArray u;
    FloatArray incrementalStrain;

    LatticeMaterialStatus *lmatStat = dynamic_cast< LatticeMaterialStatus * >( integrationRulesArray [ 0 ]->getIntegrationPoint(0)->giveMaterialStatus() );
    auto strain =lmatStat ->giveLatticeStrain();
    this->LatticeFrame3dg::computeBDmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), b);
    this->computeVectorOf(VM_Incremental, tStep, u);
    //this->computeVectorOf(VM_Total, tStep, u);
    answer.beProductOf(b, u);
    answer.times(1./this->length);
    answer += strain;
}

void
LatticeFrame3dg::giveInternalForcesVector(FloatArray &answer,
    TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix b, bt, bf;
    FloatArray u, stress, strain;

    this->length   = computeLength();
    GaussPoint *gp = this->integrationRulesArray[0]->getIntegrationPoint( 0 );
    this->LatticeFrame3dg::computeBFmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), bf);
    //    bt.beTranspositionOf( b );
    // Total stress
    this->LatticeFrame3dg::computeStrainVector( strain, gp, tStep );
    this->computeStressVector( stress, strain, integrationRulesArray[0]->getIntegrationPoint( 0 ), tStep );
    //   totalInternalForces.beProductOf( bf, stress );
    // Old stresses
    LatticeMaterialStatus *lmatStat = dynamic_cast<LatticeMaterialStatus *>( integrationRulesArray[0]->getIntegrationPoint( 0 )->giveMaterialStatus() );
    auto oldStress= lmatStat->giveLatticeStress();

    auto oldInternalForces = lmatStat->giveInternalForces();

    FloatArray incrementalStress;
    incrementalStress.beDifferenceOf( stress, oldStress );

    answer.beProductOf( bf, incrementalStress );
    //answer.beProductOf( bf, stress );
    answer += oldInternalForces;

    lmatStat->letTempInternalForcesBe(answer);
}

bool
LatticeFrame3dg::computeGtoLRotationMatrix(FloatMatrix &answer)
{
    FloatMatrix lcs;
    answer.resize(12, 12);
    answer.zero();

    this->LatticeFrame3dg::giveLocalCoordinateSystem(lcs);
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
LatticeFrame3dg::giveLocalCoordinateSystem(FloatMatrix &answer)
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
LatticeFrame3dg::giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, D_w, R_u, R_v, R_w
    };
}

void
LatticeFrame3dg::initializeFrom(InputRecord &ir)
{
    LatticeStructuralElement::initializeFrom(ir);

    referenceNode = 0;
    referenceAngle = 0;
    this->zaxis.clear();
    if ( ir.hasField(_IFT_LatticeFrame3d_zaxis) ) {
        IR_GIVE_FIELD(ir, this->zaxis, _IFT_LatticeFrame3dg_zaxis);
    } else if ( ir.hasField(_IFT_LatticeFrame3d_refnode) ) {
        IR_GIVE_FIELD(ir, referenceNode, _IFT_LatticeFrame3dg_refnode);
        if ( referenceNode == 0 ) {
            OOFEM_WARNING("wrong reference node specified. Using default orientation.");
        }
    } else if ( ir.hasField(_IFT_LatticeFrame3d_refangle) ) {
        IR_GIVE_FIELD(ir, referenceAngle, _IFT_LatticeFrame3dg_refangle);
    } else {
        throw ValueInputException(ir, _IFT_LatticeFrame3d_zaxis, "axis, reference node, or angle not set");
    }

    this->s = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, s, _IFT_LatticeFrame3d_s);

}


double
LatticeFrame3dg::computeLength()
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
LatticeFrame3dg::computeCurrentLength()
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
LatticeFrame3dg::computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
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
