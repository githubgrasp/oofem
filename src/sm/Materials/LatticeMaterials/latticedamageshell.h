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

#ifndef latticedamageshell_h
#define latticedamageshell_h

#include "latticedamage.h"

///@name Input fields for LatticeDamageShell
//@{
#define _IFT_LatticeDamageShell_Name "latticedamageshell"
#define _IFT_LatticeDamageShell_wfShear "wfshear"
#define _IFT_LatticeDamageShell_e0Shear "e0shear"
//@}

namespace oofem {

/**
 * Status for LatticeDamageShell; extends LatticeDamageStatus with a
 * second damage variable (omega2) driven by the out-of-plane shear components.
 */
class LatticeDamageShellStatus : public LatticeDamageStatus
{
protected:
    double kappaShear = 0.;
    double tempKappaShear = 0.;
    double damageShear = 0.;
    double tempDamageShear = 0.;

public:
    LatticeDamageShellStatus(GaussPoint *g);

    double giveKappaShear() const { return kappaShear; }
    double giveTempKappaShear() const { return tempKappaShear; }
    void   setTempKappaShear(double v) { tempKappaShear = v; }

    double giveDamageShear() const { return damageShear; }
    double giveTempDamageShear() const { return tempDamageShear; }
    void   setTempDamageShear(double v) { tempDamageShear = v; }

    const char *giveClassName() const override { return "LatticeDamageShellStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;
    void printOutputAt(FILE *file, TimeStep *tStep) const override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * Lattice damage for shell elements.
 * Two independent variables: omega1 for comps {1,2,4-6}, omega2 for comp 3 (out-of-plane shear).
 */
class LatticeDamageShell : public LatticeDamage
{
protected:
    double wfShear = 0.;
    double e0MeanShear = 0.;

public:
    LatticeDamageShell(int n, Domain *d) : LatticeDamage(n, d) { }

    const char *giveInputRecordName() const override { return _IFT_LatticeDamageShell_Name; }
    const char *giveClassName() const override { return "LatticeDamageShell"; }

    void initializeFrom(InputRecord &ir) override;

    void performDamageEvaluation(GaussPoint *gp, FloatArrayF< 6 > &reducedStrain) const override;

    FloatArrayF< 6 >giveLatticeStress3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) override;

    FloatMatrixF< 6, 6 >give3dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
};

} // end namespace oofem
#endif
