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
    LatticeFrame3dNL::computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                             TimeStep *tStep)
    {
        FloatMatrix d, bi, bj, bjt, dbj, dij, bf, g, b, bt;
        FloatArray u;
        this->length = computeLength();
        answer.resize(12, 12);
        answer.zero();

        this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

        double l1 = this->length * ( 1. - this->s ) / 2;
        double l2 = this->length * ( 1. + this->s ) / 2;

        this->computeVectorOf(VM_Incremental, tStep, u);

        //Axial 1
        answer.at(1, 1) = d.at(1, 1);
        answer.at(1, 2) = 0.;
        answer.at(1, 3) = 0.;
        answer.at(1, 4) = 0.;
        answer.at(1, 5) = -d.at(1, 1) * l1 * sin( u.at(5) ) * cos( u.at(6) );
        answer.at(1, 6) = -d.at(1, 1) * l1 * cos( u.at(5) ) * sin( u.at(6) );
        answer.at(1, 7) = -d.at(1, 1);
        answer.at(1, 8) = 0.;
        answer.at(1, 9) = 0.;
        answer.at(1, 10) = 0.;
        answer.at(1, 11) = -d.at(1, 1) * l2 * sin( u.at(11) ) * cos( u.at(12) );
        answer.at(1, 12) = -d.at(1, 1) * l2 * cos( u.at(11) ) * sin( u.at(12) );

        //Shear Y 1
        answer.at(2, 1) = 0;
        answer.at(2, 2) = d.at(2, 2);
        answer.at(2, 3) = 0.;
        answer.at(2, 4) = -d.at(2, 2) * ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1;
        answer.at(2, 5) = -d.at(2, 2) * ( -sin( u.at(4) ) * cos( u.at(5) ) * cos( u.at(6) ) ) * l1;
        answer.at(2, 6) = -d.at(2, 2) * ( -cos( u.at(4) ) * cos( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * sin( u.at(6) ) ) * l1;
        answer.at(2, 7) = 0.;
        answer.at(2, 8) = -d.at(2, 2);
        answer.at(2, 9) = 0.;
        answer.at(2, 10) = -d.at(2, 2) * ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2;
        answer.at(2, 11) = -d.at(2, 2) * ( -sin( u.at(10) ) * cos( u.at(11) ) * cos( u.at(12) ) ) * l2;
        answer.at(2, 12) = -d.at(2, 2) * ( -cos( u.at(10) ) * cos( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * sin( u.at(12) ) ) * l2;

        //Shear Z 1
        answer.at(3, 1) = 0;
        answer.at(3, 2) = 0.;
        answer.at(3, 3) = d.at(3, 3);
        answer.at(3, 4) = -d.at(3, 3) * ( -cos( u.at(4) ) * sin( u.at(6) ) - sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1;
        answer.at(3, 5) = -d.at(3, 3) * ( cos( u.at(4) ) * cos( u.at(5) ) * cos( u.at(6) ) ) * l1;
        answer.at(3, 6) = -d.at(3, 3) * ( -sin( u.at(4) ) * cos( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * sin( u.at(6) ) ) * l1;
        answer.at(3, 7) = 0.;
        answer.at(3, 8) = 0;
        answer.at(3, 9) = -d.at(3, 3);
        answer.at(3, 10) = -d.at(3, 3) * ( -cos( u.at(10) ) * sin( u.at(12) ) - sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2;
        answer.at(3, 11) = -d.at(3, 3) * ( cos( u.at(10) ) * cos( u.at(11) ) * cos( u.at(12) ) ) * l2;
        answer.at(3, 12) = -d.at(3, 3) * ( -sin( u.at(10) ) * cos( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * sin( u.at(12) ) ) * l2;

        // Mx 1
        answer.at(4, 1) = 0;
        answer.at(4, 2) = d.at(2, 2) * ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1;
        answer.at(4, 3) = d.at(3, 3) * ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1;
        answer.at(4, 4) = -d.at(2, 2) * ( ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( -( -sin( u.at(4) ) * sin( u.at(6) ) + cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) + ( u.at(8) - u.at(2) - ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) - d.at(3, 3) * ( ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( -( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) + ( u.at(9) - u.at(3) - ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( -sin( u.at(4) ) * sin( u.at(6) ) + cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) + d.at(4, 4);
        answer.at(4, 5) = -d.at(2, 2) * ( ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( -( sin( u.at(4) ) * cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) + ( u.at(8) - u.at(2) - ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( -cos( u.at(4) ) * cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) - d.at(3, 3) * ( ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( -( -cos( u.at(4) ) * cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) + ( u.at(9) - u.at(3) - ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( sin( u.at(4) ) * cos( u.at(5) ) * cos( u.at(6) ) ) * l1 );
        answer.at(4, 6) = -d.at(2, 2) * ( ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( -( cos( u.at(4) ) * cos( u.at(6) ) - sin( u.at(4) ) * sin( u.at(5) ) * sin( u.at(6) ) ) * l1 ) + ( u.at(8) - u.at(2) - ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( sin( u.at(4) ) * cos( u.at(6) ) + cos( u.at(4) ) * sin( u.at(5) ) * sin( u.at(6) ) ) * l1 ) - d.at(3, 3) * ( ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( -( sin( u.at(4) ) * sin( u.at(6) ) + cos( u.at(4) ) * sin( u.at(5) ) * sin( u.at(6) ) ) * l1 ) + ( u.at(9) - u.at(3) - ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( cos( u.at(4) ) * cos( u.at(6) ) - sin( u.at(4) ) * sin( u.at(5) ) * sin( u.at(6) ) ) * l1 );
        answer.at(4, 7) = 0.;
        answer.at(4, 8) = -d.at(2, 2) * ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1;
        answer.at(4, 9) = -d.at(3, 3) * ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1;
        answer.at(4, 10) = -d.at(2, 2) * ( ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( -( -sin( u.at(10) ) * sin( u.at(12) ) + cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 ) ) - d.at(3, 3) * ( ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( -( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 ) ) - d.at(4, 4);
        answer.at(4, 11) = -d.at(2, 2) * ( ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( -( sin( u.at(10) ) * cos( u.at(11) ) * cos( u.at(12) ) ) * l2 ) ) - d.at(3, 3) * ( ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( -( -cos( u.at(10) ) * cos( u.at(11) ) * cos( u.at(12) ) ) * l2 ) );
        answer.at(4, 12) = -d.at(2, 2) * ( ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( -( cos( u.at(10) ) * cos( u.at(12) ) - sin( u.at(10) ) * sin( u.at(11) ) * sin( u.at(12) ) ) * l2 ) ) - d.at(3, 3) * ( ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( -( sin( u.at(10) ) * cos( u.at(12) ) + cos( u.at(10) ) * sin( u.at(11) ) * sin( u.at(12) ) ) * l2 ) );

        // My 1
        answer.at(5, 1) = -d.at(1, 1) * ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1;
        answer.at(5, 2) = 0.;
        answer.at(5, 3) = -d.at(3, 3) * cos( u.at(5) ) * cos( u.at(6) ) * l1;
        answer.at(5, 4) = d.at(1, 1) * ( ( ( u.at(7) - u.at(1) ) + ( 1 - cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( 1 - cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) + d.at(3, 3) * ( ( cos( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( -cos( u.at(4) ) * sin( u.at(6) ) - sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 );
        answer.at(5, 5) = d.at(1, 1) * ( ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( sin( u.at(5) ) * cos( u.at(6) ) ) * l1 + ( ( u.at(7) - u.at(1) ) + ( 1 - cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( 1 - cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( -cos( u.at(4) ) * cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) + d.at(3, 3) * ( ( cos( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( cos( u.at(4) ) * cos( u.at(5) ) * cos( u.at(6) ) ) * l1 + ( ( u.at(9) - u.at(3) ) - ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( -sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) + d.at(5, 5);
        answer.at(5, 6) = d.at(1, 1) * ( ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( cos( u.at(5) ) * sin( u.at(6) ) ) * l1 + ( ( u.at(7) - u.at(1) ) + ( 1 - cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( 1 - cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( sin( u.at(4) ) * cos( u.at(6) ) + cos( u.at(4) ) * sin( u.at(5) ) * sin( u.at(6) ) ) * l1 ) + d.at(3, 3) * ( ( cos( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( -sin( u.at(4) ) * cos( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * sin( u.at(6) ) ) * l1 + ( ( u.at(9) - u.at(3) ) - ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( -cos( u.at(5) ) * sin( u.at(6) ) ) * l1 );
        answer.at(5, 7) = d.at(1, 1) * ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1;
        answer.at(5, 8) = 0;
        answer.at(5, 9) = d.at(3, 3) * cos( u.at(5) ) * cos( u.at(6) ) * l1;
        answer.at(5, 10) =  d.at(3, 3) * ( ( cos( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( -cos( u.at(10) ) * sin( u.at(12) ) - sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 );
        answer.at(5, 11) = d.at(1, 1) * ( ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( sin( u.at(11) ) * cos( u.at(12) ) ) * l2 ) + d.at(3, 3) * ( ( cos( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( cos( u.at(10) ) * cos( u.at(11) ) * cos( u.at(12) ) ) * l2 ) - d.at(5, 5);
        answer.at(5, 12) = d.at(1, 1) * ( ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( cos( u.at(11) ) * sin( u.at(12) ) ) * l2 ) + d.at(3, 3) * ( ( cos( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( -sin( u.at(10) ) * cos( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * sin( u.at(12) ) ) * l2 );

        // Mz 1
        answer.at(6, 1) = -d.at(1, 1) * ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1;
        answer.at(6, 2) = d.at(2, 2) * cos( u.at(5) ) * cos( u.at(6) ) * l1;
        answer.at(6, 3) = 0;
        answer.at(6, 4) = d.at(1, 1) * ( ( ( u.at(7) - u.at(1) ) + ( 1 - cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( 1 - cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( -sin( u.at(4) ) * sin( u.at(6) ) + cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1  )  - d.at(2, 2) * ( ( cos( u.at(5) ) * cos( u.at(6) ) * l1 ) * ( ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) );
        answer.at(6, 5) = d.at(1, 1) * ( ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( sin( u.at(5) ) * cos( u.at(6) ) ) * l1 + ( ( u.at(7) - u.at(1) ) + ( 1 - cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( 1 - cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( sin( u.at(4) ) * cos( u.at(5) ) * cos( u.at(6) ) ) * l1 )  - d.at(2, 2) * ( ( cos( u.at(5) ) * cos( u.at(6) ) * l1 ) * ( -( sin( u.at(4) ) * cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) + ( ( u.at(8) - u.at(2) ) - ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( -sin( u.at(5) ) * cos( u.at(6) ) * l1 ) );
        answer.at(6, 6) = d.at(1, 1) * ( ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( ( cos( u.at(5) ) * sin( u.at(6) ) ) * l1 ) + ( ( u.at(7) - u.at(1) ) + ( 1 - cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( 1 - cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( cos( u.at(4) ) * cos( u.at(6) ) - sin( u.at(4) ) * sin( u.at(5) ) * sin( u.at(6) ) ) * l1 )  - d.at(2, 2) * ( ( cos( u.at(5) ) * cos( u.at(6) ) * l1 ) * ( -( cos( u.at(4) ) * cos( u.at(6) ) - sin( u.at(4) ) * sin( u.at(5) ) * sin( u.at(6) ) ) * l1 ) + ( ( u.at(8) - u.at(2) ) - ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( -cos( u.at(5) ) * sin( u.at(6) ) * l1 ) ) + d.at(6, 6);
        answer.at(6, 7) = d.at(1, 1) * ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1;
        answer.at(6, 8) = -d.at(2, 2) * cos( u.at(5) ) * cos( u.at(6) ) * l1;
        answer.at(6, 9) = 0.;
        answer.at(6, 10) = -d.at(2, 2) * ( ( cos( u.at(5) ) * cos( u.at(6) ) * l1 ) * ( -( -sin( u.at(10) ) * sin( u.at(12) ) + cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 ) );
        answer.at(6, 11) = d.at(1, 1) * ( ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( ( sin( u.at(11) ) * cos( u.at(12) ) ) * l2  ) )  - d.at(2, 2) * ( ( cos( u.at(5) ) * cos( u.at(6) ) * l1 ) * ( -( sin( u.at(10) ) * cos( u.at(11) ) * cos( u.at(12) ) ) * l2 ) );
        answer.at(6, 12) = d.at(1, 1) * ( ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 * ( ( cos( u.at(11) ) * sin( u.at(12) ) ) * l2 ) )  - d.at(2, 2) * ( ( cos( u.at(5) ) * cos( u.at(6) ) * l1 ) * ( -( cos( u.at(10) ) * cos( u.at(12) ) - sin( u.at(10) ) * sin( u.at(11) ) * sin( u.at(12) ) ) * l2 ) ) - d.at(6, 6);

        //Axial 2
        answer.at(7, 1) = -d.at(1, 1);
        answer.at(7, 2) = 0.;
        answer.at(7, 3) = 0.;
        answer.at(7, 4) = 0.;
        answer.at(7, 5) = d.at(1, 1) * l1 * sin( u.at(5) ) * cos( u.at(6) );
        answer.at(7, 6) = d.at(1, 1) * l1 * cos( u.at(5) ) * sin( u.at(6) );
        answer.at(7, 7) = d.at(1, 1);
        answer.at(7, 8) = 0.;
        answer.at(7, 9) = 0.;
        answer.at(7, 10) = 0.;
        answer.at(7, 11) = d.at(1, 1) * l2 * sin( u.at(11) ) * cos( u.at(12) );
        answer.at(7, 12) = d.at(1, 1) * l2 * cos( u.at(11) ) * sin( u.at(12) );

        //Shear Y 2
        answer.at(8, 1) = 0;
        answer.at(8, 2) = -d.at(2, 2);
        answer.at(8, 3) = 0.;
        answer.at(8, 4) = d.at(2, 2) * ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1;
        answer.at(8, 5) = d.at(2, 2) * ( -sin( u.at(4) ) * cos( u.at(5) ) * cos( u.at(6) ) ) * l1;
        answer.at(8, 6) = d.at(2, 2) * ( -cos( u.at(4) ) * cos( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * sin( u.at(6) ) ) * l1;
        answer.at(8, 7) = 0.;
        answer.at(8, 8) = d.at(2, 2);
        answer.at(8, 9) = 0.;
        answer.at(8, 10) = d.at(2, 2) * ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2;
        answer.at(8, 11) = d.at(2, 2) * ( -sin( u.at(10) ) * cos( u.at(11) ) * cos( u.at(12) ) ) * l2;
        answer.at(8, 12) = d.at(2, 2) * ( -cos( u.at(10) ) * cos( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * sin( u.at(12) ) ) * l2;

        //Shear Z 2
        answer.at(9, 1) = 0;
        answer.at(9, 2) = 0.;
        answer.at(9, 3) = -d.at(3, 3);
        answer.at(9, 4) = d.at(3, 3) * ( -cos( u.at(4) ) * sin( u.at(6) ) - sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1;
        answer.at(9, 5) = d.at(3, 3) * ( cos( u.at(4) ) * cos( u.at(5) ) * cos( u.at(6) ) ) * l1;
        answer.at(9, 6) = d.at(3, 3) * ( -sin( u.at(4) ) * cos( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * sin( u.at(6) ) ) * l1;
        answer.at(9, 7) = 0.;
        answer.at(9, 8) = 0;
        answer.at(9, 9) = d.at(3, 3);
        answer.at(9, 10) = d.at(3, 3) * ( -cos( u.at(10) ) * sin( u.at(12) ) - sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2;
        answer.at(9, 11) = d.at(3, 3) * ( cos( u.at(10) ) * cos( u.at(11) ) * cos( u.at(12) ) ) * l2;
        answer.at(9, 12) = d.at(3, 3) * ( -sin( u.at(10) ) * cos( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * sin( u.at(12) ) ) * l2;

        // Mx 2
        answer.at(10, 1) = 0;
        answer.at(10, 2) = -d.at(2, 2) * ( -( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 );
        answer.at(10, 3) = -d.at(3, 3) * ( -( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 );
        answer.at(10, 4) = -d.at(2, 2) * ( ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( -( -sin( u.at(4) ) * sin( u.at(6) ) + cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) ) - d.at(3, 3) * ( ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( -( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) ) - d.at(4, 4);
        answer.at(10, 5) = -d.at(2, 2) * ( ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( -sin( u.at(4) ) * cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) - d.at(3, 3) * ( ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( cos( u.at(4) ) * cos( u.at(5) ) * cos( u.at(6) ) ) * l1 );
        answer.at(10, 6) = -d.at(2, 2) * ( ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( -( cos( u.at(4) ) * cos( u.at(6) ) - sin( u.at(4) ) * sin( u.at(5) ) * sin( u.at(6) ) ) * l1 ) ) - d.at(3, 3) * ( ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( -( sin( u.at(4) ) * cos( u.at(6) ) + cos( u.at(4) ) * sin( u.at(5) ) * sin( u.at(6) ) ) * l1 ) );
        answer.at(10, 7) = 0.;
        answer.at(10, 8) = -d.at(2, 2) * ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2;
        answer.at(10, 9) = -d.at(3, 3) * ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2;
        answer.at(10, 10) = -d.at(2, 2) * ( ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( -( -sin( u.at(10) ) * sin( u.at(12) ) + cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 ) + ( u.at(8) - u.at(2) - ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 ) - d.at(3, 3) * ( ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( -( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 ) + ( u.at(9) - u.at(3) - ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( -sin( u.at(10) ) * sin( u.at(12) ) + cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 ) + d.at(4, 4);
        answer.at(10, 11) = -d.at(2, 2) * ( ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( -sin( u.at(10) ) * cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( u.at(8) - u.at(2) - ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( -cos( u.at(10) ) * cos( u.at(11) ) * cos( u.at(12) ) ) * l2 ) - d.at(3, 3) * ( ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( cos( u.at(10) ) * cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( u.at(9) - u.at(3) - ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( sin( u.at(10) ) * cos( u.at(11) ) * cos( u.at(12) ) ) * l2 );
        answer.at(10, 12) = -d.at(2, 2) * ( ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( -cos( u.at(10) ) * cos( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * sin( u.at(12) ) ) * l2 + ( u.at(8) - u.at(2) - ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( sin( u.at(10) ) * cos( u.at(12) ) + cos( u.at(10) ) * sin( u.at(11) ) * sin( u.at(12) ) ) * l2 ) - d.at(3, 3) * ( ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( -sin( u.at(10) ) * cos( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * sin( u.at(12) ) ) * l2 + ( u.at(9) - u.at(3) - ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( cos( u.at(10) ) * cos( u.at(12) ) - sin( u.at(10) ) * sin( u.at(11) ) * sin( u.at(12) ) ) * l2 );

        // My 2
        answer.at(11, 1) = -d.at(1, 1) * ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(12) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2;
        answer.at(11, 2) = 0.;
        answer.at(11, 3) = -d.at(3, 3) * cos( u.at(11) ) * cos( u.at(12) ) * l2;
        answer.at(11, 4) =  d.at(3, 3) * ( ( cos( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( -cos( u.at(4) ) * sin( u.at(6) ) - sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 );
        answer.at(11, 5) = d.at(1, 1) * ( ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) + d.at(3, 3) * ( ( cos( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( cos( u.at(4) ) * cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) - d.at(5, 5);
        answer.at(11, 6) = d.at(1, 1) * ( ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( cos( u.at(5) ) * sin( u.at(6) ) ) * l1 ) + d.at(3, 3) * ( ( cos( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( -sin( u.at(4) ) * cos( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * sin( u.at(6) ) ) * l1 );
        answer.at(11, 7) = d.at(1, 1) * ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(12) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2;
        answer.at(11, 8) = 0;
        answer.at(11, 9) = d.at(3, 3) * cos( u.at(11) ) * cos( u.at(12) ) * l2;
        answer.at(11, 10) = d.at(1, 1) * ( ( ( u.at(7) - u.at(1) ) + ( 1 - cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( 1 - cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 ) + d.at(3, 3) * ( ( cos( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( -cos( u.at(10) ) * sin( u.at(12) ) - sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 );
        answer.at(11, 11) = d.at(1, 1) * ( ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( sin( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( ( u.at(7) - u.at(1) ) + ( 1 - cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( 1 - cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( -cos( u.at(10) ) * cos( u.at(11) ) * cos( u.at(12) ) ) * l2 ) + d.at(3, 3) * ( ( cos( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( cos( u.at(10) ) * cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( ( u.at(9) - u.at(3) ) - ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( -sin( u.at(11) ) * cos( u.at(12) ) ) * l2 ) + d.at(5, 5);
        answer.at(11, 12) = d.at(1, 1) * ( ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( cos( u.at(11) ) * sin( u.at(12) ) ) * l2 + ( ( u.at(7) - u.at(1) ) + ( 1 - cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( 1 - cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( sin( u.at(11) ) * cos( u.at(12) ) + cos( u.at(10) ) * sin( u.at(11) ) * sin( u.at(12) ) ) * l2 ) + d.at(3, 3) * ( ( cos( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( -sin( u.at(10) ) * cos( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * sin( u.at(12) ) ) * l2 + ( ( u.at(9) - u.at(3) ) - ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( -cos( u.at(11) ) * sin( u.at(12) ) ) * l2 );

        // Mz 2
        answer.at(12, 1) = -d.at(1, 1) * ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2;
        answer.at(12, 2) = d.at(2, 2) * cos( u.at(11) ) * cos( u.at(12) ) * l2;
        answer.at(12, 3) = 0.;
        answer.at(12, 4) = -d.at(2, 2) * ( ( cos( u.at(11) ) * cos( u.at(12) ) * l2 ) * ( -( -sin( u.at(4) ) * sin( u.at(6) ) + cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) );
        answer.at(12, 5) = d.at(1, 1) * ( ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( ( sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) )  - d.at(2, 2) * ( ( cos( u.at(11) ) * cos( u.at(12) ) * l2 ) * ( -( sin( u.at(4) ) * cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) );
        answer.at(12, 6) = d.at(1, 1) * ( ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( ( cos( u.at(5) ) * sin( u.at(6) ) ) * l1 ) )  - d.at(2, 2) * ( ( cos( u.at(11) ) * cos( u.at(12) ) * l2 ) * ( -( cos( u.at(4) ) * cos( u.at(6) ) - sin( u.at(4) ) * sin( u.at(5) ) * sin( u.at(6) ) ) * l1 ) ) - d.at(6, 6);
        answer.at(12, 7) = d.at(1, 1) * ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2;
        answer.at(12, 8) = -d.at(2, 2) * cos( u.at(11) ) * cos( u.at(12) ) * l2;
        answer.at(12, 9) = 0.;
        answer.at(12, 10) = d.at(1, 1) * ( ( ( u.at(7) - u.at(1) ) + ( 1 - cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( 1 - cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( -sin( u.at(10) ) * sin( u.at(12) ) + cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 )  - d.at(2, 2) * ( ( cos( u.at(11) ) * cos( u.at(12) ) * l2 ) * ( -( -sin( u.at(10) ) * sin( u.at(12) ) + cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 ) );
        answer.at(12, 11) = d.at(1, 1) * ( ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( ( sin( u.at(11) ) * cos( u.at(12) ) ) * l2 ) + ( ( u.at(7) - u.at(1) ) + ( 1 - cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( 1 - cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( sin( u.at(10) ) * cos( u.at(11) ) * cos( u.at(12) ) ) * l2 )  - d.at(2, 2) * ( ( cos( u.at(11) ) * cos( u.at(12) ) * l2 ) * ( -(  sin( u.at(10) ) * cos( u.at(11) ) * cos( u.at(12) ) ) * l2 ) + ( ( u.at(8) - u.at(2) ) - ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( -sin( u.at(11) ) * cos( u.at(12) ) * l2 ) );
        answer.at(12, 12) = d.at(1, 1) * ( ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 * ( ( sin( u.at(11) ) * cos( u.at(12) ) ) * l2 ) + ( ( u.at(7) - u.at(1) ) + ( 1 - cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( 1 - cos( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( cos( u.at(10) ) * cos( u.at(12) ) - sin( u.at(10) ) * sin( u.at(11) ) * sin( u.at(12) ) ) * l2 )  - d.at(2, 2) * ( ( cos( u.at(11) ) * cos( u.at(12) ) * l2 ) * ( -( cos( u.at(10) ) * cos( u.at(12) ) - sin( u.at(10) ) * sin( u.at(11) ) * sin( u.at(12) ) ) * l2 ) + ( ( u.at(8) - u.at(2) ) - ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 ) * ( -cos( u.at(11) ) * sin( u.at(12) ) * l2 ) ) + d.at(6, 6);

        answer.times(1. / length);

        return;
    }

    void
    LatticeFrame3dNL::computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
    // Computes the vector containing the strains at the Gauss point gp of
    // the receiver, at time step tStep. The nature of these strains depends
    // on the element's type.
    {
        FloatMatrix b;
        FloatArray u;
        double l1 = this->length * ( 1. - this->s ) / 2;
        double l2 = this->length * ( 1. + this->s ) / 2;
        this->computeVectorOf(VM_Incremental, tStep, u);
        LatticeMaterialStatus *lmatStat = dynamic_cast < LatticeMaterialStatus * > ( integrationRulesArray [ 0 ]->getIntegrationPoint(0)->giveMaterialStatus() );
        auto strain = lmatStat->giveLatticeStrain();

        answer.resize(6);
        answer.at(1) = ( 1 - cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + ( 1 - cos( u.at(5) ) * cos( u.at(6) ) ) * l1 + u.at(7) - u.at(1);
        answer.at(2) = -( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 + u.at(8) - u.at(2);
        answer.at(3) = -( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 + u.at(9) - u.at(3);
        answer.at(4) = u.at(10) - u.at(4);
        answer.at(5) = u.at(11) - u.at(5);
        answer.at(6) = u.at(12) - u.at(6);
        answer.times(1. / this->length);
        // printf("Strain/n");
        // answer.printYourself();
        answer += strain;
    }
    //
    void
    LatticeFrame3dNL::giveInternalForcesVector(FloatArray &answer,
                                               TimeStep *tStep, int useUpdatedGpRecord)
    {
        FloatMatrix b, bt, bf;
        FloatArray u, stress, strain;
        this->computeVectorOf(VM_Incremental, tStep, u);
        this->length   = computeLength();
        GaussPoint *gp = this->integrationRulesArray [ 0 ]->getIntegrationPoint(0);

        // Total stress
        this->LatticeFrame3dNL::computeStrainVector(strain, gp, tStep);
        this->computeStressVector(stress, strain, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

        // Old stresses
        LatticeMaterialStatus *lmatStat = dynamic_cast < LatticeMaterialStatus * > ( integrationRulesArray [ 0 ]->getIntegrationPoint(0)->giveMaterialStatus() );
        auto oldStress = lmatStat->giveLatticeStress();

        auto oldInternalForces = lmatStat->giveInternalForces();
        double l1 = this->length * ( 1. - this->s ) / 2;
        double l2 = this->length * ( 1. + this->s ) / 2;
        FloatArray incrementalStress;
        incrementalStress.beDifferenceOf(stress, oldStress);

        answer.resize(12);
        answer.at(1) = -incrementalStress.at(1),
        answer.at(2) = -incrementalStress.at(2);
        answer.at(3) = -incrementalStress.at(3);
        answer.at(4) = +incrementalStress.at(2) * ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 - incrementalStress.at(3) * ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 - incrementalStress.at(4);
        answer.at(5) = -incrementalStress.at(1) * ( sin( u.at(4) ) * sin( u.at(6) ) - cos( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 + incrementalStress.at(3) * ( cos( u.at(5) ) * cos( u.at(6) ) ) * l1 - incrementalStress.at(5);
        answer.at(6) = incrementalStress.at(1) * ( cos( u.at(4) ) * sin( u.at(6) ) + sin( u.at(4) ) * sin( u.at(5) ) * cos( u.at(6) ) ) * l1 - incrementalStress.at(2) * ( cos( u.at(5) ) * cos( u.at(6) ) ) * l1 - incrementalStress.at(6);
        answer.at(7) = incrementalStress.at(1);
        answer.at(8) = incrementalStress.at(2);
        answer.at(9) = incrementalStress.at(3);
        answer.at(10) = incrementalStress.at(2) * ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - incrementalStress.at(3) * ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 + incrementalStress.at(4);
        answer.at(11) = -incrementalStress.at(1) * ( sin( u.at(10) ) * sin( u.at(12) ) - cos( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 + incrementalStress.at(3) * ( cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + incrementalStress.at(5);
        answer.at(12) = incrementalStress.at(1) * ( cos( u.at(10) ) * sin( u.at(12) ) + sin( u.at(10) ) * sin( u.at(11) ) * cos( u.at(12) ) ) * l2 - incrementalStress.at(2) * ( cos( u.at(11) ) * cos( u.at(12) ) ) * l2 + incrementalStress.at(6);
        answer += oldInternalForces;

        lmatStat->letTempInternalForcesBe(answer);
    }

    int
    LatticeFrame3dNL::giveLocalCoordinateSystem(FloatMatrix &answer)
    {
        FloatArray lx, ly, lz, help(3);
        FloatArray coordA, coordB;
        FloatArray uA(6), uAIncr(6), uB(6), uBIncr(6);
        IntArray dofid = {
            1, 2, 3, 4, 5, 6
        };

        TimeStep *tStep = this->domain->giveEngngModel()->giveCurrentStep();

        Node *nodeA, *nodeB;
        nodeA = this->giveNode(1);
        nodeB = this->giveNode(2);

        //Local coordinate system is determined from the displacement of last step.

        coordA = nodeA->giveCoordinates();
        nodeA->giveUnknownVector(uA, dofid, VM_Total, tStep, false);
        nodeA->giveUnknownVector(uAIncr, dofid, VM_Incremental, tStep, false);
        for (int i = 1; i <= 3; i++) {
            coordA.at(i) += uA.at(i) - uAIncr.at(i);
        }

        coordB = nodeB->giveCoordinates();
        nodeB->giveUnknownVector(uB, dofid, VM_Total, tStep, false);
        nodeB->giveUnknownVector(uBIncr, dofid, VM_Incremental, tStep, false);

        for (int i = 1; i <= 3; i++) {
            coordB.at(i) += uB.at(i) - uBIncr.at(i);
        }

        lx.beDifferenceOf(coordB, coordA);
        lx.normalize();


        if ( this->referenceNode ) {
            Node *refNode = this->giveDomain()->giveNode(this->referenceNode);
            help.beDifferenceOf( refNode->giveCoordinates(), nodeA->giveCoordinates() );

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

            help.at(3) = 1.0;     // up-vector

            // here is ly is used as a temp var
            if ( fabs( lx.dotProduct(help) ) > 0.999 ) {  // Check if it is vertical
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
