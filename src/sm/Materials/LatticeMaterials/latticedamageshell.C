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

#include "latticedamageshell.h"
#include "latticelinearelastic.h"
#include "gausspoint.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "mathfem.h"
#include "Elements/LatticeElements/latticestructuralelement.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(LatticeDamageShell);

// ---------------------------------------------------------------------------
// LatticeDamageShellStatus
// ---------------------------------------------------------------------------

LatticeDamageShellStatus :: LatticeDamageShellStatus(GaussPoint *g) :
    LatticeDamageStatus(g)
{ }


void
LatticeDamageShellStatus :: initTempStatus()
{
    LatticeDamageStatus :: initTempStatus();
    tempKappaShear  = kappaShear;
    tempDamageShear = damageShear;
}


void
LatticeDamageShellStatus :: updateYourself(TimeStep *tStep)
{
    LatticeDamageStatus :: updateYourself(tStep);
    kappaShear  = tempKappaShear;
    damageShear = tempDamageShear;
}


void
LatticeDamageShellStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    LatticeDamageStatus :: printOutputAt(file, tStep);
    fprintf(file, "kappaShear %.8e damageShear %.8e\n", kappaShear, damageShear);
}


void
LatticeDamageShellStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    LatticeDamageStatus :: saveContext(stream, mode);
    if ( !stream.write(kappaShear) ) { THROW_CIOERR(CIO_IOERR); }
    if ( !stream.write(damageShear) ) { THROW_CIOERR(CIO_IOERR); }
}


void
LatticeDamageShellStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    LatticeDamageStatus :: restoreContext(stream, mode);
    if ( !stream.read(kappaShear) ) { THROW_CIOERR(CIO_IOERR); }
    if ( !stream.read(damageShear) ) { THROW_CIOERR(CIO_IOERR); }
}


// ---------------------------------------------------------------------------
// LatticeDamageShell
// ---------------------------------------------------------------------------

MaterialStatus *
LatticeDamageShell :: CreateStatus(GaussPoint *gp) const
{
    return new LatticeDamageShellStatus(gp);
}


void
LatticeDamageShell :: initializeFrom(InputRecord &ir)
{
    LatticeDamage :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, wfShear, _IFT_LatticeDamageShell_wfShear);

    // Default peak shear strain equals the normal peak strain
    e0MeanShear = e0Mean;
    IR_GIVE_OPTIONAL_FIELD(ir, e0MeanShear, _IFT_LatticeDamageShell_e0Shear);
}


void
LatticeDamageShell :: performDamageEvaluation(GaussPoint *gp, FloatArrayF< 6 > &reducedStrain) const
{
    // omega1 via base class
    LatticeDamage :: performDamageEvaluation(gp, reducedStrain);

    // omega2: out-of-plane shear (comp 3), independent fracture energy
    auto *status = static_cast< LatticeDamageShellStatus * >( this->giveStatus(gp) );

    double le      = static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveLength();
    double e0s     = this->give(e0_ID, gp) * e0MeanShear;
    double eNormals = this->give(eNormal_ID, gp) * eNormalMean;

    double shearEquivStrain = fabs(reducedStrain.at(3));
    double fShear = shearEquivStrain - status->giveKappaShear();

    double tempKappaShear, omegaShear;
    if ( fShear <= 0. ) {
        tempKappaShear = status->giveKappaShear();
        omegaShear     = status->giveDamageShear();
    } else {
        tempKappaShear = shearEquivStrain;
        omegaShear     = computeDamageParamExplicit(tempKappaShear, e0s, wfShear, eNormals, le);
    }

    status->setTempKappaShear(tempKappaShear);
    status->setTempDamageShear(omegaShear);
}


FloatArrayF< 6 >
LatticeDamageShell :: giveLatticeStress3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep)
{
    // Base class applies omega1 to all components; correct comp 3 here to use omega2.
    auto stress = LatticeDamage :: giveLatticeStress3d(strain, gp, tStep);

    auto *status = static_cast< LatticeDamageShellStatus * >( this->giveStatus(gp) );
    double omega2 = min(status->giveTempDamageShear(), 0.99999);
    auto elastic = LatticeLinearElastic :: give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);
    const auto &reducedStrain = status->giveTempReducedLatticeStrain();
    stress.at(3) = elastic.at(3, 3) * reducedStrain.at(3) * ( 1. - omega2 );

    return stress;
}


FloatMatrixF< 6, 6 >
LatticeDamageShell :: give3dLatticeStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto elastic = LatticeLinearElastic :: give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    if ( mode == ElasticStiffness ) {
        return elastic;
    } else if ( mode == SecantStiffness || mode == TangentStiffness ) {
        auto *status = static_cast< LatticeDamageShellStatus * >( this->giveStatus(gp) );
        double omega1 = min(status->giveTempDamage(),      0.99999);
        double omega2 = min(status->giveTempDamageShear(), 0.99999);

        auto result = elastic;
        for ( int i : { 1, 2, 4, 5, 6 } ) {  // comps 1,2,4-6 -> omega1
            result.at(i, i) *= ( 1. - omega1 );
        }
        result.at(3, 3) *= ( 1. - omega2 );  // comp 3 -> omega2
        return result;
    } else {
        OOFEM_ERROR("Unsupported stiffness mode\n");
        return elastic;
    }
}

} // end namespace oofem
