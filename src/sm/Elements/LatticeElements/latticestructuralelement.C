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

#include "sm/Elements/LatticeElements/latticestructuralelement.h"
#include "gausspoint.h"
#include "crosssection.h"

namespace oofem {
LatticeStructuralElement :: LatticeStructuralElement(int n, Domain *aDomain) : StructuralElement(n, aDomain)
{ }

void
LatticeStructuralElement :: initializeFrom(InputRecord &ir)
{
    StructuralElement :: initializeFrom(ir);
}


void LatticeStructuralElement :: giveSectionScaleFactors3d(FloatArray &q, GaussPoint *gp)
{


    q.resize(6);
    double A  = this->giveArea(gp);
    double Aq1 = this->giveShearArea1(gp);
    double Aq2 = this->giveShearArea2(gp);
    double I1 = this->giveI1(gp);
    double I2 = this->giveI2(gp);
    double J = this->giveJ(gp);

    q.at(1) = A;
    q.at(2) = Aq1;
    q.at(3) = Aq2;
    q.at(4) = J;
    q.at(5) = I1;
    q.at(6) = I2;
}

void LatticeStructuralElement :: giveSectionScaleFactors2d(FloatArray &q, GaussPoint *gp)
{
    q.resize(3);
    double A  = this->giveArea(gp);
    double I2 = this->giveI2(gp);


    q.at(1) = A;
    q.at(2) = A;
    q.at(3) = I2;
} 

void LatticeStructuralElement :: convertStressToResultants3d(FloatArray &S, const FloatArray &sigma, GaussPoint *gp)
{

    S = sigma;

    auto *mat = static_cast<LatticeStructuralMaterial*>(this->giveCrossSection()->giveMaterial(gp));
    if ( mat->giveLatticeResponseType() == LatticeStructuralMaterial::LRT_StressBased ) {
        FloatArray q;
        this->giveSectionScaleFactors3d(q, gp);
        for (int i = 1; i <= 6; ++i) {
            S.at(i) *= q.at(i);
        }
    }
}

void LatticeStructuralElement :: convertStressToResultants2d(FloatArray &S, const FloatArray &sigma, GaussPoint *gp)
{
    S = sigma;

    auto *mat = static_cast<LatticeStructuralMaterial*>(this->giveCrossSection()->giveMaterial(gp));
    if ( mat->giveLatticeResponseType() == LatticeStructuralMaterial::LRT_StressBased ) {
        FloatArray q;
        this->giveSectionScaleFactors2d(q, gp);

        for (int i = 1; i <= 3; ++i) {
            S.at(i) *= q.at(i);
        }
    }
}




 
void LatticeStructuralElement :: convertTangentToResultantTangent3d(FloatMatrix &DS, const FloatMatrix &Dsig, GaussPoint *gp)
{

    DS = Dsig;

    auto *mat = static_cast<LatticeStructuralMaterial*>(this->giveCrossSection()->giveMaterial(gp));
    if ( mat->giveLatticeResponseType() == LatticeStructuralMaterial::LRT_StressBased ) {
        FloatArray q;
        this->giveSectionScaleFactors3d(q, gp);

        // DS = Q * Dsig. Scale rows
        for (int i = 1; i <= 6; ++i) {
            for (int j = 1; j <= 6; ++j) {
                DS.at(i,j) *= q.at(i);
            }
        }
    }
}


void LatticeStructuralElement :: convertTangentToResultantTangent2d(FloatMatrix &DS, const FloatMatrix &Dsig, GaussPoint *gp)
{
    DS = Dsig;

    auto *mat = static_cast<LatticeStructuralMaterial*>(this->giveCrossSection()->giveMaterial(gp));
    if ( mat->giveLatticeResponseType() == LatticeStructuralMaterial::LRT_StressBased ) {
        FloatArray q;
        this->giveSectionScaleFactors2d(q, gp);
        // DS = Q * Dsig. Scale rows
        for (int i = 1; i <= 3; ++i) {
            for (int j = 1; j <= 3; ++j) {
                DS.at(i,j) *= q.at(i);
            }
        }
    }
    }
 
 
void
LatticeStructuralElement :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step tStep. The nature of these strains depends
// on the element's type.
{
    FloatMatrix b;
    FloatArray u;

    if ( !this->isActivated(tStep) ) {
        this->computeBmatrixAt(gp, b);
        answer.resize(b.giveNumberOfRows());
        answer.zero();
        return;
    }

    this->computeBmatrixAt(gp, b);
    this->computeVectorOf(VM_Total, tStep, u);

    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }
    answer.beProductOf(b, u);
    answer.times(1./giveLength());
}


 
void
LatticeStructuralElement :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralElement :: printOutputAt(file, tStep);
}



} // end namespace oofem
