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
    double I1 = this->giveI1(gp);
    double I2 = this->giveI2(gp);
    double Ip = this->giveIp(gp);

    q.at(1) = A;
    q.at(2) = A;
    q.at(3) = A;
    q.at(4) = Ip;
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
    FloatArray q;
    this->giveSectionScaleFactors3d(q, gp);

    S = sigma;
    for (int i = 1; i <= 6; ++i) {
        S.at(i) *= q.at(i);
    }
}

void LatticeStructuralElement :: convertStressToResultants2d(FloatArray &S, const FloatArray &sigma, GaussPoint *gp)
{
    FloatArray q;
    this->giveSectionScaleFactors2d(q, gp);

    S = sigma;
    for (int i = 1; i <= 3; ++i) {
        S.at(i) *= q.at(i);
    }
}




 
void LatticeStructuralElement :: convertTangentToResultantTangent3d(FloatMatrix &DS, const FloatMatrix &Dsig, GaussPoint *gp)
{
    FloatArray q;
    this->giveSectionScaleFactors3d(q, gp);

    DS = Dsig;
    // DS = Q * Dsig. Scale rows
    for (int i = 1; i <= 6; ++i) {
        for (int j = 1; j <= 6; ++j) {
            DS.at(i,j) *= q.at(i);
        }
    }
}


void LatticeStructuralElement :: convertTangentToResultantTangent2d(FloatMatrix &DS, const FloatMatrix &Dsig, GaussPoint *gp)
{
    FloatArray q;
    this->giveSectionScaleFactors2d(q, gp);

    DS = Dsig;
    // DS = Q * Dsig. Scale rows
    for (int i = 1; i <= 3; ++i) {
        for (int j = 1; j <= 3; ++j) {
            DS.at(i,j) *= q.at(i);
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

    /// FIXME: This output should just be moved to the elements themselves. But, they don't exist yet? / Mikael
    FloatArray forces;
    if ( this->giveClassName() == std::string("LatticeBeam3d") ) {
        this->giveInternalForcesVector(forces, tStep, 0);
        fprintf(file, "LatticeBeam forces = %e %e %e %e %e %e.\n", forces.at(7), forces.at(8), forces.at(9), forces.at(10), forces.at(11), forces.at(12) );
    } else if ( this->giveClassName() == std::string("LatticeBeam3dBoundary") ) {
        this->giveInternalForcesVector(forces, tStep, 0);
        fprintf(file, "LatticeBeam3dBoundary forces = %e %e %e %e %e %e.\n", forces.at(7), forces.at(8), forces.at(9), forces.at(10), forces.at(11), forces.at(12) );
    }
}



} // end namespace oofem
