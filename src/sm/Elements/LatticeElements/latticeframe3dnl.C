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
 *               Copyright (C) 1993 - 2023   Borek Patzak
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
#include "latticeframe3dnl.h"
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
    REGISTER_Element(LatticeFrame3dNL);

    LatticeFrame3dNL::LatticeFrame3dNL(int n, Domain *aDomain) : LatticeFrame3d(n, aDomain)
    {
        numberOfDofMans = 2;
    }

    LatticeFrame3dNL::~LatticeFrame3dNL()
    {}
    void
    LatticeFrame3dNL::computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
    //Returns the strain matrix of the receiver.
    {
        FloatArray u;
        TimeStep *tStep = this->domain->giveEngngModel()->giveCurrentStep();
        this->computeVectorOf(VM_Total, tStep, u);

        this->length = computeLength();
        double l1 = this->length * ( 1. + this->s ) / 2;
        double l2 = this->length * ( 1. - this->s ) / 2;
        answer.resize(6, 12);
        answer.zero();
        double cx1 = ( cos(u.at(5) ) * cos(u.at(6) ) ) * l1;
        double cy1 = ( cos(u.at(4) ) * sin(u.at(6) ) + sin(u.at(4) ) * sin(u.at(5) ) * cos(u.at(6) ) ) * l1;
        double cz1 = ( sin(u.at(4) ) * sin(u.at(6) ) - cos(u.at(4) ) * sin(u.at(5) ) * cos(u.at(6) ) ) * l1;


        double cx2 = ( cos(u.at(11) ) * cos(u.at(12) ) ) * l2;
        double cy2 = ( cos(u.at(10) ) * sin(u.at(12) ) + sin(u.at(10) ) * sin(u.at(11) ) * cos(u.at(12) ) ) * l2;
        double cz2 = ( sin(u.at(10) ) * sin(u.at(12) ) - cos(u.at(10) ) * sin(u.at(11) ) * cos(u.at(12) ) ) * l2;


        //Normal displacement jump in x-direction
        //First node
        answer.at(1, 1) = -1.;
        answer.at(1, 2) = 0.;
        answer.at(1, 3) = 0.;
        answer.at(1, 4) = 0.;
        answer.at(1, 5) = -cz1;
        answer.at(1, 6) = cy1;
        //Second node
        answer.at(1, 7) = 1.;
        answer.at(1, 8) = 0.;
        answer.at(1, 9) = 0.;
        answer.at(1, 10) = 0.;
        answer.at(1, 11) = -cz2;
        answer.at(1, 12) = cy2;

        //Shear displacement jump in y-plane
        //first node
        answer.at(2, 1) = 0.;
        answer.at(2, 2) = -1.;
        answer.at(2, 3) =  0.;
        answer.at(2, 4) = cz1;
        answer.at(2, 5) = 0;
        answer.at(2, 6) = -cx1;
        //Second node
        answer.at(2, 7) = 0.;
        answer.at(2, 8) = 1.;
        answer.at(2, 9) =  0.;
        answer.at(2, 10) = cz2;
        answer.at(2, 11) = 0;
        answer.at(2, 12) = -cx2;

        //Shear displacement jump in z-plane
        //first node
        answer.at(3, 1) = 0.;
        answer.at(3, 2) = 0.;
        answer.at(3, 3) = -1.;
        answer.at(3, 4) = -cy1;
        answer.at(3, 5) = cx1;
        answer.at(3, 6) = 0.;
        //Second node
        answer.at(3, 7) = 0.;
        answer.at(3, 8) = 0.;
        answer.at(3, 9) =  1.;
        answer.at(3, 10) = -cy2;
        answer.at(3, 11) = cx2;
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
    LatticeFrame3dNL::computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                             TimeStep *tStep)
    {
        FloatMatrix d, bt, db, b;
        FloatArray u;

        this->computeVectorOf(VM_Total, tStep, u);
        this->length = computeLength();

        answer.resize(12, 12);
        answer.zero();
        this->computeBmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), b);
        this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

        db.beProductOf(d, b);
        db.times(1. / length);
        bt.beTranspositionOf(b);
        answer.beProductOf(bt, db);

        return;
    }

    void
    LatticeFrame3dNL::computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
    {
        FloatArray u;
        this->computeVectorOf(VM_Total, tStep, u);

        this->length   = computeLength();
        double l1 = this->length * ( 1. + this->s ) / 2;
        double l2 = this->length * ( 1. - this->s ) / 2;
        LatticeMaterialStatus *lmatStat = dynamic_cast < LatticeMaterialStatus * > ( integrationRulesArray [ 0 ]->getIntegrationPoint(0)->giveMaterialStatus() );
        auto strain = lmatStat->giveLatticeStrain();

        double cx1 = ( cos(u.at(5) ) * cos(u.at(6) ) ) * l1;
        double cy1 = ( cos(u.at(4) ) * sin(u.at(6) ) + sin(u.at(4) ) * sin(u.at(5) ) * cos(u.at(6) ) ) * l1;
        double cz1 = ( sin(u.at(4) ) * sin(u.at(6) ) - cos(u.at(4) ) * sin(u.at(5) ) * cos(u.at(6) ) ) * l1;


        double cx2 = ( cos(u.at(11) ) * cos(u.at(12) ) ) * l2;
        double cy2 = ( cos(u.at(10) ) * sin(u.at(12) ) + sin(u.at(10) ) * sin(u.at(11) ) * cos(u.at(12) ) ) * l2;
        double cz2 = ( sin(u.at(10) ) * sin(u.at(12) ) - cos(u.at(10) ) * sin(u.at(11) ) * cos(u.at(12) ) ) * l2;
        //
        answer.resize(6);
        answer.at(1) = u.at(7) - u.at(1) - cx2 - cx1 + l1 + l2;
        answer.at(2) = u.at(8) - u.at(2) - cy2 - cy1;
        answer.at(3) = u.at(9) - u.at(3) - cz2 - cz1;
        answer.at(4) = u.at(10) - u.at(4);
        answer.at(5) = u.at(11) - u.at(5);
        answer.at(6) = u.at(12) - u.at(6);
        answer.times(1. / this->length);
    }
    //
    void
    LatticeFrame3dNL::giveInternalForcesVector(FloatArray &answer,
                                               TimeStep *tStep, int useUpdatedGpRecord)
    {
        FloatMatrix b, bt, bf, d;
        FloatArray u, stress, strain;
        this->computeVectorOf(VM_Total, tStep, u);
        this->length   = computeLength();
        GaussPoint *gp = this->integrationRulesArray [ 0 ]->getIntegrationPoint(0);

        this->computeStrainVector(strain, gp, tStep);
        this->computeStressVector(stress, strain, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

        LatticeMaterialStatus *lmatStat = dynamic_cast < LatticeMaterialStatus * > ( integrationRulesArray [ 0 ]->getIntegrationPoint(0)->giveMaterialStatus() );
        auto oldStress = lmatStat->giveLatticeStress();

        auto oldInternalForces = lmatStat->giveInternalForces();
        double l1 = this->length * ( 1. + this->s ) / 2;
        double l2 = this->length * ( 1. - this->s ) / 2;

        double cx1 = ( cos(u.at(5) ) * cos(u.at(6) ) ) * l1;
        double cy1 = ( cos(u.at(4) ) * sin(u.at(6) ) + sin(u.at(4) ) * sin(u.at(5) ) * cos(u.at(6) ) ) * l1;
        double cz1 = ( sin(u.at(4) ) * sin(u.at(6) ) - cos(u.at(4) ) * sin(u.at(5) ) * cos(u.at(6) ) ) * l1;


        double cx2 = ( cos(u.at(11) ) * cos(u.at(12) ) ) * l2;
        double cy2 = ( cos(u.at(10) ) * sin(u.at(12) ) + sin(u.at(10) ) * sin(u.at(11) ) * cos(u.at(12) ) ) * l2;
        double cz2 = ( sin(u.at(10) ) * sin(u.at(12) ) - cos(u.at(10) ) * sin(u.at(11) ) * cos(u.at(12) ) ) * l2;

        answer.resize(12);
        answer.at(1) = -stress.at(1);
        answer.at(2) = -stress.at(2);
        answer.at(3) = -stress.at(3);
        answer.at(4) =  stress.at(2) * cz1 - stress.at(3) * cy1 - stress.at(4);
        answer.at(5) = -stress.at(1) * cz1 + stress.at(3) * cx1 - stress.at(5);
        answer.at(6) = stress.at(1) * cy1 - stress.at(2) * cx1 - stress.at(6);
        answer.at(7) = stress.at(1);
        answer.at(8) = stress.at(2);
        answer.at(9) = stress.at(3);
        answer.at(10) =  stress.at(2) * cz2 - stress.at(3) * cy2 + stress.at(4);
        answer.at(11) = -stress.at(1) * cz2 + stress.at(3) * cx2 + stress.at(5);
        answer.at(12) = stress.at(1) * cy2 - stress.at(2) * cx2 + stress.at(6);

        lmatStat->letTempInternalForcesBe(answer);
    }
} // end namespace oofem
