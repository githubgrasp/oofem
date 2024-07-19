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
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef CDPM2F_h
#define CDPM2F_h

#include "sm/Materials/ConcreteMaterials/concretedpm2.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "mathfem.h"

#define _IFT_CDPM2F_Name "cdpm2f"
#define _IFT_CDPM2F_Lf "lf"
#define _IFT_CDPM2F_Vf0 "vf0"
#define _IFT_CDPM2F_Df "df"
#define _IFT_CDPM2F_Tau0 "tau0"
#define _IFT_CDPM2F_Beta "beta"
#define _IFT_CDPM2F_f "f"
#define _IFT_CDPM2F_Ef "ef"
#define _IFT_CDPM2F_Sm "sm"
#define _IFT_CDPM2F_Alpha "alpha"
#define _IFT_CDPM2F_xi "xi"
//@}

namespace oofem {
/**
 * This class implements the material status associated to CDPM2F.
 * @author Chao Zhou, Peter Grassl
 */
class CDPM2FStatus : public ConcreteDPM2Status
{
public:
    /// Constructor
    CDPM2FStatus(GaussPoint *gp);
};


//   ********************************
//   *** CLASS CDPM2F   ***
//   ********************************

/**
 * This class implements an extension of CDPM2 for fibre reinforced materials which exhibit strain hardening.
 * The fibre equation used are described in "CDPM2F: A damage-plasticity approach to modelling the failure of engineered cementitious composites" by Chao Zhou, Antoine Merlot and Peter Grassl
 * The fibre tensile response is based on the paper "Crack bridging in fibre reinforced cementitious composites with slip-hardening interfaces." by Lin and Li. 1997
 * @author Chao Zhou, Peter Grassl
 */
class CDPM2F : public ConcreteDPM2
{
public:

protected:
    /// length of fibre
    double lf = 0.;
    /// base volume fraction of fibres
    double vf0 = 0.;
    double vf = 0.;
    double vfm = 0.;
    double z = 0.;
    /// diameter of fibre
    double df = 0.;
    /// Young's modulus of fibre
    double ef = 0.;
    /// initial bond strength of interface between fibre and matrix
    double tau0 = 0.;
    /// hardening parameter of interface between fibre and matrix
    double beta = 0.;
    ///snubbing factor
    double f = 0.;
    double eta = 0.;
    double g = 0.;
    double s0 = 0.;
    double omega = 0.;
    double k = 0.;
    double lambda = 0.;
    double deltaStar = 0.;
    double c = 0.;
    double deltaCu = 0.;
    ///saturated crack spacing
    double sm = 0.;
    ///parameter related to dispersion (smaller or equal to 1)
    double alpha = 0.;
    double alphaMin = 0.;
    double ap = 0.;
    ///parameter affecting the hardening response.
    double xi = 0.;
    double deltaUl = 0.;
    double stressCu = 0.;


public:
    /// Constructor
    CDPM2F(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;

    const char *giveClassName() const override { return "CDPM2F"; }
    const char *giveInputRecordName() const override { return _IFT_CDPM2F_Name; }

    double  computeFibreStress(double delta) const;

    double computeCrackOpening(double epscr, const double le) const;

    double computeMatrixStress(double delta) const;

    double computeStressResidual(double equivStrain, double omega, double kappaOne, double kappaTwo, double le) const;

    /// Compute damage parameter in tension.
    virtual double computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double omegaOld, double rateFactor) const override;
};
} //end namespace oofem
#endif
