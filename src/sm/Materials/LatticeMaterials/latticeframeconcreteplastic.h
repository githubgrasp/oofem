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

#ifndef latticeframeconcreteplastic_h
#define latticeframeconcreteplastic_h

#include "latticestructuralmaterial.h"
#include "cltypes.h"
#include "randommaterialext.h"
#include "strainvector.h"
#include "stressvector.h"
#include "latticematstatus.h"

///@name Input fields for LatticeFrameConcrePlastic
//@{
#define _IFT_LatticeFrameConcretePlastic_Name "LatticeFrameConcretePlastic"
#define _IFT_LatticeFrameConcretePlastic_talpha "talpha"
#define _IFT_LatticeFrameConcretePlastic_e "e"
#define _IFT_LatticeFrameConcretePlastic_n "n"
#define _IFT_LatticeFrameConcretePlastic_nx0 "nx0"
#define _IFT_LatticeFrameConcretePlastic_mx0 "mx0"
#define _IFT_LatticeFrameConcretePlastic_my0 "my0"
#define _IFT_LatticeFrameConcretePlastic_mz0 "mz0"
#define _IFT_LatticeFrameConcretePlastic_tol "tol"
#define _IFT_LatticeFrameConcretePlastic_iter "iter"
#define _IFT_LatticeFrameConcretePlastic_sub "sub"
#define _IFT_LatticeFrameConcretePlastic_plastic "plastic"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to LatticeFrameConcretePlastic.
 * @authors: Gumaa Abdelrhim, Peter Grassl
 */

class LatticeFrameConcretePlasticStatus : public LatticeMaterialStatus
{
public:

    enum state_flag_values {
        LatticeFrameConcretePlastic_Elastic,
        LatticeFrameConcretePlastic_Unloading,
        LatticeFrameConcretePlastic_Plastic,
    };


    enum LatticeFrameConcretePlastic_ReturnResult {
        RR_NotConverged,
        RR_Converged
    };


protected:

    int tempReturnResult = LatticeFrameConcretePlasticStatus::RR_NotConverged;



public:

    /// Constructor
    LatticeFrameConcretePlasticStatus(int n, Domain *d, GaussPoint *g);


    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "LatticeFrameConcretePlasticStatus"; }

    void letTempReturnResultBe(const int result) { tempReturnResult = result; }

    int giveTempReturnResult() const { return tempReturnResult; }
};


/**
 * This class implements a local random linear elastic model for lattice elements.
 */
class LatticeFrameConcretePlastic : public LatticeStructuralMaterial

{
protected:

    ///Normal modulus
    double e;

    ///Ratio of shear and normal modulus
    double nu;

    ///maximum axial force in x-axis x-axis nx0
    double nx0;

    ///maximum  bending moment about x-axis mx0
    double mx0;

    ///maximum  bending moment about x-axis my0
    double my0;

    ///maximum  bending moment about x-axis mz0
    double mz0;

    /// yield tolerance
    double yieldTol;

    /// maximum number of iterations for stress return
    double newtonIter;

    ///number Of SubIncrements
    double numberOfSubIncrements;

    ///plastic flag
    double plasticFlag;

    enum LatticeFrameConcretePlastic_ReturnResult { RR_NotConverged, RR_Converged };
    //   mutable LatticeFrameConcretePlastic_ReturnResult returnResult = RR_NotConverged; /// FIXME: This must be removed. Not thread safe. Shouldn't be stored at all.

    double initialYieldStress = 0.;

    //

public:
    LatticeFrameConcretePlastic(int n, Domain *d) : LatticeStructuralMaterial(n, d) { };

    FloatArrayF< 4 >computeFVector(const FloatArrayF< 4 > &sigma, GaussPoint *gp, TimeStep *tStep) const;

    FloatMatrixF< 4, 4 >computeDMMatrix(const FloatArrayF< 4 > &sigma, GaussPoint *gp, TimeStep *tStep) const;

    FloatArrayF< 6 >giveThermalDilatationVector(GaussPoint *gp,  TimeStep *tStep) const override;

    FloatArrayF< 6 >giveReducedLatticeStrain(GaussPoint *gp, TimeStep *tStep) const;

    virtual FloatArrayF< 6 >giveReducedStrain(GaussPoint *gp, TimeStep *tStep) const;

    FloatArrayF< 6 >performPlasticityReturn(GaussPoint *gp, const FloatArrayF< 6 > &reducedStrain, TimeStep *tStep) const;

    void performRegularReturn(FloatArrayF< 4 > &stress, double yieldValue, GaussPoint *gp, TimeStep *tStep) const;

    double computeYieldValue(const FloatArrayF< 4 > &sigma, GaussPoint *gp, TimeStep *tStep) const;

    FloatMatrixF< 5, 5 >computeJacobian(const FloatArrayF< 4 > &sigma, const double deltaLambda, GaussPoint *gp, TimeStep *tStep) const;

    FloatArrayF< 6 >giveFrameForces3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) override;

    const char *giveInputRecordName() const override { return _IFT_LatticeFrameConcretePlastic_Name; }

    const char *giveClassName() const override { return "LatticeFrameConcretePlastic"; }

    void initializeFrom(InputRecord &ir) override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    Interface *giveInterface(InterfaceType) override;

    FloatMatrixF< 6, 6 >give3dFrameStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    bool hasMaterialModeCapability(MaterialMode mode) const override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    MaterialStatus *giveStatus(GaussPoint *gp) const override;

protected:
};
} // end namespace oofem

#endif
