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
  LatticeFrame3dg::computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer,  TimeStep *tStep)
// Returns the strain matrix of the receiver.
{
    //Assemble Bmatrix (used to compute strains and rotations)
    answer.resize(6, 12);
    answer.zero();

    this->length = computeLength();
    FloatArray u;
    this->computeVectorOf(VM_Total, tStep, u);
     //printf("u/n");
     //u.printYourself();
    //Normal displacement jump in x-direction
    //First node
    answer.at(1, 1) = -1.;
    answer.at(1, 2) = 0.;
    answer.at(1, 3) = 0.;
    answer.at(1, 4) = 0.;
    answer.at(1, 5) =  sin(u.at(5))*this->length*(1.-this->s)/2.;
    answer.at(1, 6) =  sin(u.at(6))*this->length*(1.-this->s)/2.;
    //Second node
    answer.at(1, 7) = 1.;
    answer.at(1, 8) = 0.;
    answer.at(1, 9) = 0.;
    answer.at(1, 10) = 0.;
    answer.at(1, 11) =  sin(u.at(10))*this->length*(1.+this->s)/2.;
    answer.at(1, 12) =  sin(u.at(12))*this->length*(1.+this->s)/2.;

    //Shear displacement jump in y-plane
    //first node
    answer.at(2, 1) = 0.;
    answer.at(2, 2) = -1.;
    answer.at(2, 3) =  0.;
    answer.at(2, 4) = 0.;
    answer.at(2, 5) = 0;
    answer.at(2, 6) =  -cos(u.at(6))*this->length*(1.-this->s)/2.;
    //Second node
    answer.at(2, 7) = 0.;
    answer.at(2, 8) = 1.;
    answer.at(2, 9) =  0.;
    answer.at(2, 10) = 0.;
    answer.at(2, 11) = 0;
    answer.at(2, 12) = -cos(u.at(12))*this->length*(1.+this->s)/2.;

    //Shear displacement jump in z-plane
    //first node
    answer.at(3, 1) = 0.;
    answer.at(3, 2) = 0.;
    answer.at(3, 3) = -1.;
    answer.at(3, 4) = 0.;
    answer.at(3, 5) = cos(u.at(5))*this->length*(1.-this->s)/2.;
    answer.at(3, 6) = 0.;
    //Second node
    answer.at(3, 7) = 0.;
    answer.at(3, 8) = 0.;
    answer.at(3, 9) =  1.;
    answer.at(3, 10) = 0.;
    answer.at(3, 11) = cos(u.at(10))*this->length*(1.+this->s)/2.;
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
LatticeFrame3dg::computeBFmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer,  TimeStep *tStep)
{
    //Assemble BFmatrix (used to compute forces )
    answer.resize(12, 6);
    answer.zero();

    this->length = computeLength();
    FloatArray u;
    this->computeVectorOf(VM_Total, tStep, u);

    //First Nx1
    answer.at(1, 1) = -1.;
    answer.at(1, 2) = 0.;
    answer.at(1, 3) = 0.;
    answer.at(1, 4) = 0.;
    answer.at(1, 5) = 0.;
    answer.at(1, 6) = 0.;

    //Second Nx2
    answer.at(7, 1) = 1.;
    answer.at(7, 2) = 0.;
    answer.at(7, 3) = 0.;
    answer.at(7, 4) = 0.;
    answer.at(7, 5) = 0.;
    answer.at(7, 6) = 0.;

    //Shear Y
    //first node Fy1
    answer.at(2, 1) = 0.;
    answer.at(2, 2) = -1.;
    answer.at(2, 3) = 0.;
    answer.at(2, 4) = 0.;
    answer.at(2, 5) = 0;
    answer.at(2, 6) = 0.;
    //Second node Fy2
    answer.at(8, 1) = 0.;
    answer.at(8, 2) = 1.;
    answer.at(8, 3) =  0.;
    answer.at(8, 4) = 0.;
    answer.at(8, 5) = 0;
    answer.at(8, 6) = 0.;

    //Shear Z
    //first node Fz1
    answer.at(3, 1) = 0.;
    answer.at(3, 2) = 0.;
    answer.at(3, 3) = -1.;
    answer.at(3, 4) = 0.;
    answer.at(3, 5) = 0.;
    answer.at(3, 6) = 0.;
    //Second node Fz2
    answer.at(9, 1) = 0.;
    answer.at(9, 2) = 0.;
    answer.at(9, 3) =  1.;
    answer.at(9, 4) = 0.;
    answer.at(9, 5) = 0.;
    answer.at(9, 6) = 0.;

    //Torsion  x
    //First node Mx1
    answer.at(4, 1) = 0.;
    answer.at(4, 2) = -sin(u.at(5))*this->length*(1.+this->s)/2.;
    answer.at(4, 3) = -sin(u.at(6))*this->length*(1.+this->s)/2.;
    answer.at(4, 4) = -1.;
    answer.at(4, 5) = 0.;
    answer.at(4, 6) = 0.;
    //Second node Mx2
    answer.at(10, 1) = 0.;
    answer.at(10, 2) = sin(u.at(10))*this->length*(1.+this->s)/2.;
    answer.at(10, 3) = sin(u.at(12))*this->length*(1.+this->s)/2.;
    answer.at(10, 4) = 1.;
    answer.at(10, 5) = 0.;
    answer.at(10, 6) = 0.;

    //Moment around y-axis
    //First node my1
    answer.at(5, 1) = sin(u.at(5))*this->length*(1.+this->s)/2.;
    answer.at(5, 2) = 0.;
    answer.at(5, 3) = cos(u.at(6))*this->length*(1.+this->s)/2.;
    answer.at(5, 4) = 0.;
    answer.at(5, 5) = -1.;
    answer.at(5, 6) = 0.;
    //Second node My2
    answer.at(11, 1) = sin(u.at(10))*this->length*(1.+this->s)/2.;
    answer.at(11, 2) = 0.;
    answer.at(11, 3) = cos(u.at(12))*this->length*(1.+this->s)/2.;
    answer.at(11, 4) = 0.;
    answer.at(11, 5) = 1.;
    answer.at(11, 6) = 0.;

    //Moment around z-axis
    //First node Mz1
    answer.at(6, 1) = sin(u.at(5))*this->length*(1.+this->s)/2.;
    answer.at(6, 2) = -cos(u.at(6))*this->length*(1.+this->s)/2.;
    answer.at(6, 3) = 0.;
    answer.at(6, 4) = 0.;
    answer.at(6, 5) = 0.;
    answer.at(6, 6) = -1.;
    //Second node Mz2
    answer.at(12, 1) = sin(u.at(10))*this->length*(1.+this->s)/2.;
    answer.at(12, 2) = -cos(u.at(12))*this->length*(1.+this->s)/2.;
    answer.at(12, 3) =  0.;
    answer.at(12, 4) = 0.;
    answer.at(12, 5) = 0.;
    answer.at(12, 6) = 1.;

    return;
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
    this->LatticeFrame3dg::computeBFmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), bf, tStep);
    this->LatticeFrame3dg::computeBmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), bj, tStep);
    this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

    dbj.beProductOf(d, bj);
    dbj.times(1. / length);
    bjt.beTranspositionOf(bj);
    answer.beProductOf(bf, dbj);
   // printf("answer/n");
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

    if ( !this->LatticeFrame3dg::isActivated(tStep) ) {
        answer.resize(StructuralMaterial :: giveSizeOfVoigtSymVector(gp->giveMaterialMode() ) );
        answer.zero();
        return;
    }

    this->LatticeFrame3dg::computeBmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), b, tStep);

    this->computeVectorOf(VM_Total, tStep, u);
    answer.beProductOf(b, u);

    answer.times(1./this->length);
   // printf("Strain/n");
    //answer.printYourself();
    }

void
LatticeFrame3dg::giveInternalForcesVector(FloatArray &answer,
    TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix b, bt, bf;
    FloatArray u, stress, strain;

    this->length = computeLength();

    this->computeVectorOf(VM_Total, tStep, u);

    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.clear();

    this->LatticeFrame3dg::computeBFmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), bf, tStep);

    if ( useUpdatedGpRecord == 1 ) {
        LatticeMaterialStatus *lmatStat = dynamic_cast< LatticeMaterialStatus * >( integrationRulesArray [ 0 ]->getIntegrationPoint(0)->giveMaterialStatus() );
        stress = lmatStat->giveLatticeStress();
    } else {
        if ( !this->isActivated(tStep) ) {
            strain.zero();
        }
        this->LatticeFrame3dg::computeStrainVector( strain, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep );
        this->computeStressVector(stress, strain, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);
    }
   // printf("Stress/n");
    //stress.printYourself();

    answer.beProductOf(bf, stress);
    //printf("Force/n");
    //answer.printYourself();
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}
bool
LatticeFrame3dg::computeGtoLRotationMatrix(FloatMatrix &answer,  TimeStep *tStep)
{
    FloatMatrix lcs;
    answer.resize(12, 12);
    answer.zero();

    this->giveLocalCoordinateSystem(lcs, tStep);
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
LatticeFrame3dg::giveLocalCoordinateSystem(FloatMatrix &answer,  TimeStep *tStep)
{
    FloatArray lx, ly, lz, help(3);
    FloatArray coord;
    FloatArray uA, uB, u;

    Node *nodeA, *nodeB;
    this->computeVectorOf(VM_Total, tStep, u);

    nodeA = this->giveNode(1);
    coord = nodeA->giveCoordinates();
    uA =  { u.at(1), u.at(2), u.at(3) } ;
    uA += coord;

    nodeB = this->giveNode(2);
    nodeB->giveCoordinates();
    coord = nodeB->giveCoordinates();
    uB =  { u.at(7), u.at(8), u.at(9) } ;
    uB += coord;

    lx.beDifferenceOf(uB, uA );
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

} // end namespace oofem
