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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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

#ifndef latticeplasticitydamage_h
#define latticeplasticitydamage_h

#include "latticelinearelastic.h"
#include "latticematstatus.h"

///@name Input fields for LatticePlasticityDamage
//@{
#define _IFT_LatticePlasticityDamage_Name "latticeplastdam"
#define _IFT_LatticePlasticityDamage_tol "tol"
#define _IFT_LatticePlasticityDamage_iter "iter"
#define _IFT_LatticePlasticityDamage_sub "sub"
#define _IFT_LatticePlasticityDamage_ft "ft"
#define _IFT_LatticePlasticityDamage_fc "fc"
#define _IFT_LatticePlasticityDamage_angle1 "angle1"
#define _IFT_LatticePlasticityDamage_angle2 "angle2"
#define _IFT_LatticePlasticityDamage_flow "flow"
#define _IFT_LatticePlasticityDamage_stype "stype"
#define _IFT_LatticePlasticityDamage_wf "wf"
#define _IFT_LatticePlasticityDamage_ft1 "ft1"
#define _IFT_LatticePlasticityDamage_wf1 "wf1"
#define _IFT_LatticePlasticityDamage_ahard "ahard"
#define _IFT_LatticePlasticityDamage_damage "damage"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to LatticePlasticityDamage.
 * @author: Peter Grassl
 */
class LatticePlasticityDamageStatus : public LatticeMaterialStatus
{
protected:

    double kappaP = 0.;

    double tempKappaP = 0.;

    double kappaDOne = 0., kappaDTwo = 0.;//, kappaDThree;

    double tempKappaDOne = 0., tempKappaDTwo = 0.;//, tempKappaDThree;

    double damage = 0.;

    double tempDamage = 0.;

    double dissipation = 0.;
    
    double tempDissipation = 0.;

    double tensionDissipation = 0.;
    
    double tempTensionDissipation = 0.;

    double compressionDissipation = 0.;
    
    double tempCompressionDissipation = 0.;

    double shearDissipation = 0.;
    
    double tempShearDissipation = 0.;

    double deltaDissipation = 0.;
    
    double tempDeltaDissipation = 0.;

    double tensionDeltaDissipation = 0.;
    
    double tempTensionDeltaDissipation = 0.;

    double compressionDeltaDissipation = 0.;
    
    double tempCompressionDeltaDissipation = 0.;

    double shearDeltaDissipation = 0.;
    
    double tempShearDeltaDissipation = 0.;

    int deletionFlag = 0;
    int tempDeletionFlag = 1.;
    
    //double e0 = 0.;

    int compressionFlag = 0;


public:

    /// Constructor
    LatticePlasticityDamageStatus(int n, Domain *d, GaussPoint *g);

    double giveKappaP() const { return kappaP; }

    double giveTempKappaP() const { return tempKappaP; }

    double giveKappaDOne() const { return kappaDOne; }
    double giveKappaDTwo() const { return kappaDTwo; }

    double giveTempKappaDOne() const { return tempKappaDOne; }
    double giveTempKappaDTwo() const { return tempKappaDTwo; }

    void   setTempKappaP(double newKappa) { tempKappaP = newKappa; }

    void   setTempKappaDOne(double newKappa) { tempKappaDOne = newKappa; }

    void   setTempKappaDTwo(double newKappa) { tempKappaDTwo = newKappa; }

    double giveDamage() const { return damage; }

    double giveTempDamage() const { return tempDamage; }

    void   setTempDamage(double newDamage) { tempDamage = newDamage; }

    double giveDissipation() const override {return dissipation;}
    double giveTensionDissipation() const {return tensionDissipation;}
    double giveCompressionDissipation() const {return compressionDissipation;}
    double giveShearDissipation() const {return shearDissipation;}

    double giveDeltaDissipation() const override {return deltaDissipation;}
    double giveTensionDeltaDissipation() const {return tensionDeltaDissipation;}
    double giveCompressionDeltaDissipation() const {return compressionDeltaDissipation;}
    double giveShearDeltaDissipation() const {return shearDeltaDissipation;}
    
    void setTempDissipation(double newDiss) {this->tempDissipation = newDiss;}
    void setTempTensionDissipation(double newDiss) {this->tempTensionDissipation = newDiss;}
    void setTempCompressionDissipation(double newDiss) {this->tempCompressionDissipation = newDiss;}
    void setTempShearDissipation(double newDiss) {this->tempShearDissipation = newDiss;}

    void setTempDeltaDissipation(double newDiss) {this->tempDeltaDissipation = newDiss;}
    void setTempTensionDeltaDissipation(double newDiss) {this->tempTensionDeltaDissipation = newDiss;}
    void setTempCompressionDeltaDissipation(double newDiss) {this->tempCompressionDeltaDissipation = newDiss;}
    void setTempShearDeltaDissipation(double newDiss) {this->tempShearDeltaDissipation = newDiss;}

    void setTempDeletionFlag(int newFlag) {this->tempDeletionFlag = newFlag;}
    int giveDeletionFlag() const{ return this->deletionFlag;}
    int giveTempDeletionFlag() const{ return this->tempDeletionFlag;}
    
    int giveCompressionFlag() const { return compressionFlag; }

    void setCompressionFlag(int flag) { compressionFlag = flag; }

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "LatticePlasticityDamageStatus"; }

    void initTempStatus() override;

    void updateYourself(TimeStep *) override;

    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * This class implements a local random plasticity damage model for concrete for lattice elements.
 */
class LatticePlasticityDamage : public LatticeLinearElastic
{
protected:

    enum LatticePlasticityDamage_ReturnResult {
        RR_Unknown,
        RR_NotConverged,
        RR_Converged
    };

    double initialYieldStress = 0.;

    /// tensile strength
    double ft = 0.;
    /// compressive strength
    double fc = 0.;
    /// frictional angle of the yield surface
    double frictionAngleOne = 0.;
    /// frictional angle of the yield surface
    double frictionAngleTwo = 0.;
    /// frictional angle of the plastic potential
    double flowAngleOne = 0.;

    /// frictional angle of the plastic potential
    double flowAngleTwo = 0.;

    /// determines the softening -> corresponds to crack opening (not strain) when tension stress vanishes
    double wf = 0.;

    /// softening type determines the type of softening. 0 is exponential and 1 is bilinear.
    int softeningType = 0;

    /// ratio of tensile stress value for bilinear stress-crack opening curve
    double ftOneRatio = 0.;

    /// crack opening value for bilinear stress-crack opening curve
    double wfOne = 0.;

    ///hardening parameter
    double aHard = 0.;

    /// yield tolerance
    double yieldTol = 0.;

    /// maximum number of iterations for stress return
    int newtonIter = 0;
    int numberOfSubIncrements = 0;

    ///damageFlag
    int damageFlag = 0;

    virtual double giveTensileStrength(GaussPoint *gp, TimeStep *tStep) const { return this->give(ft_strength, gp) * this->ft; }

    virtual double giveCompressiveStrength(GaussPoint *gp, TimeStep *tStep) const { return this->give(fc_strength, gp) * this->fc; }

public:

    /// Constructor
    LatticePlasticityDamage(int n, Domain *d);

    const char *giveInputRecordName() const override { return _IFT_LatticePlasticityDamage_Name; }
    const char *giveClassName() const override { return "LatticePlasticityDamage"; }


    void initializeFrom(InputRecord &ir) override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    double give(int aProperty, GaussPoint *gp) const override;

    FloatMatrixF< 6, 6 >give3dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    bool hasMaterialModeCapability(MaterialMode mode) const override;


    double computeDeltaDissipation(const double omega,
				   const FloatArray &reducedStrain,
				   GaussPoint *gp,
				   TimeStep *atTime) const;
    

    FloatArrayF< 3 >computeFVector(const FloatArrayF< 3 > &sigma,
                                   const double deltaLambda,
                                   const int surface,
                                   GaussPoint *gp,
                                   TimeStep *tStep) const;

    FloatArrayF< 3 >computeMVector(const FloatArrayF< 3 > &sigma,
                                   const double deltaLambda,
                                   const int surface,
                                   GaussPoint *gp,
                                   TimeStep *tStep) const;

    FloatMatrixF< 3, 3 >computeDMMatrix(const FloatArrayF< 3 > &sigma,
                                        const double deltaLambda,
                                        const int surface,
                                        GaussPoint *gp, TimeStep *tStep) const;

    int checkTransition(const FloatArrayF< 3 > &sigma, const double tempKappa, GaussPoint *gp, TimeStep *tStep) const;


    FloatMatrixF< 4, 4 >computeJacobian(const FloatArrayF< 3 > &sigma,
                                        const double tempKappa,
                                        const double deltaLambda,
                                        const int surface,
                                        GaussPoint *gp,
                                        TimeStep *tStep) const;

    virtual double computeDamageParam(double kappaOne, double kappaTwo, GaussPoint *gp, TimeStep *tStep) const;

    FloatArrayF< 6 >giveLatticeStress3d(const FloatArrayF< 6 > &jump, GaussPoint *gp, TimeStep *tStep) override;

    FloatArrayF< 6 >performPlasticityReturn(GaussPoint *gp,
                                            const FloatArrayF< 6 > &reducedStrain,
                                            TimeStep *tStep) const;

    void performDamageEvaluation(GaussPoint *gp,
                                 FloatArrayF< 6 > &reducedStrain,
                                 TimeStep *tStep) const;

    double performRegularReturn(FloatArrayF< 3 > &stress,
                                LatticePlasticityDamage_ReturnResult &returnResult,
                                double yieldValue,
                                const int surface,
                                GaussPoint *gp,
                                TimeStep *tStep) const;

    double computeYieldValue(const FloatArrayF< 3 > &sigma,
                             const double tempKappa,
                             const int surface,
                             GaussPoint *gp,
                             TimeStep *tStep) const;

    double computeHardening(const double kappa,
                            GaussPoint *gp) const;


    double computeDHardeningDKappa(const double kappa,
                                   GaussPoint *gp) const;
    double computeDDHardeningDDKappa(const double kappa,
                                     GaussPoint *gp) const;

    double computeDuctilityMeasure(FloatArray &stress, double ductilityParameter) const;

    double computeYieldStress(double kappaP, GaussPoint *gp);
    const

    double computeEquivalentStress(const FloatArray &tempSigma) const;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    virtual FloatArrayF< 6 >giveReducedStrain(GaussPoint *gp, TimeStep *tStep) const;


protected:

    int giveIPValue(FloatArray &answer,
                    GaussPoint *gp,
                    InternalStateType type,
                    TimeStep *atTime) override;

    int giveIPValueSize(InternalStateType type,
                        GaussPoint *gp);

    int giveIntVarCompFullIndx(IntArray &answer,
                               InternalStateType type,
                               MaterialMode mmode);

    InternalStateValueType giveIPValueType(InternalStateType type);
};
} // end namespace oofem

#endif
