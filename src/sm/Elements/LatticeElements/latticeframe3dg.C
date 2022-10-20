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
    answer.at(1, 11) =  sin(u.at(5))*this->length*(1.+this->s)/2.;
    answer.at(1, 12) =  sin(u.at(6))*this->length*(1.+this->s)/2.;

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
    answer.at(2, 12) = -cos(u.at(6))*this->length*(1.+this->s)/2.;

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
    answer.at(3, 11) = cos(u.at(5))*this->length*(1.+this->s)/2.;
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
LatticeFrame3dg::computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                       TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    FloatMatrix d, bi, bj, bjt, dbj, dij;

    this->length = computeLength();

    answer.resize(12, 12);
    answer.zero();
    this->computeBmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), bj, tStep);
    this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

    dbj.beProductOf(d, bj);
    dbj.times(1. / length);
    bjt.beTranspositionOf(bj);
    answer.beProductOf(bjt, dbj);

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

    LatticeMaterialStatus *lmatStat = dynamic_cast< LatticeMaterialStatus * >( integrationRulesArray [ 0 ]->getIntegrationPoint(0)->giveMaterialStatus() );
    auto Strain =lmatStat ->giveLatticeStrain();

    this->computeBmatrixAt( integrationRulesArray[0]->getIntegrationPoint( 0 ), b, tStep );
    this->computeVectorOf(VM_Incremental, tStep, u);

    answer.beProductOf(b, u);
    answer.times(1./this->length);

    answer +=Strain;
}

void
LatticeFrame3dg::giveInternalForcesVector(FloatArray &answer,
                                         TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix b, bt;
    FloatArray u, stress, strain;
    FloatArray incrementalStress;
    FloatArray incrementalInternalForces;
    FloatArray oldInternalForces;

    this->length = computeLength();
    GaussPoint *gp = this->integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    //this->computeVectorOf( VM_Total, tStep, u );
    //   answer.clear();

    this->computeBmatrixAt( integrationRulesArray[0]->getIntegrationPoint( 0 ), b, tStep );
    bt.beTranspositionOf( b );

    //Total stress
    this->LatticeFrame3dg::computeStrainVector(strain, gp, tStep);
    this->computeStressVector( stress, strain, integrationRulesArray[0]->getIntegrationPoint( 0 ), tStep );

    //Old stresses
    LatticeMaterialStatus *lmatStat = dynamic_cast<LatticeMaterialStatus *>( integrationRulesArray[0]->getIntegrationPoint( 0 )->giveMaterialStatus() );
    auto oldStress                          = lmatStat->giveLatticeStress();
    oldInternalForces.beProductOf(bt, oldStress);

    incrementalStress.beDifferenceOf(oldStress, stress);
    answer.beProductOf(bt, incrementalStress);

     answer +=oldInternalForces;



        return;

}
} // end namespace oofem
