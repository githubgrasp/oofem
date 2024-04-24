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
        FloatArray u, uIncr,  un1d;
        this->length = computeLength();
        answer.resize(12, 12);
        answer.zero();

        this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

        double l1 = this->length * ( 1. - this->s ) / 2;
        double l2 = this->length * ( 1. + this->s ) / 2;

        this->computeVectorOf(VM_Incremental, tStep, u);
        this->computeVectorOf(VM_Incremental, tStep, uIncr);
        this->computeVectorOf(VM_Total, tStep, un1d);
        auto un1     = un1d - uIncr;

        double xin1=(cos(un1.at(5))*cos(un1.at(6)))*l1;
        double yin1=(cos(un1.at(4))*sin(un1.at(6))+sin(un1.at(4))*sin(un1.at(5))*cos(un1.at(6)))*l1;
        double zin1=(sin(un1.at(4))*sin(un1.at(6))-cos(un1.at(4))*sin(un1.at(5))*cos(un1.at(6)))*l1;


        double xjn1=(cos(un1.at(11))*cos(un1.at(12)))*l2;
        double yjn1=(cos(un1.at(10))*sin(un1.at(12))+sin(un1.at(10))*sin(un1.at(11))*cos(un1.at(12)))*l2;
        double zjn1=(sin(un1.at(10))*sin(un1.at(12))-cos(un1.at(10))*sin(un1.at(11))*cos(un1.at(12)))*l2;



        //Axial 1
        answer.at(1, 1) = d.at(1, 1);
        answer.at(1, 2) = 0.;
        answer.at(1, 3) = 0.;
        answer.at(1, 4) = 0.;
        answer.at(1, 5) = d.at(1,1)*(-zin1);
        answer.at(1, 6) = -d.at(1,1)*(yin1);
        answer.at(1, 7) = -d.at(1, 1);
        answer.at(1, 8) = 0.;
        answer.at(1, 9) = 0.;
        answer.at(1, 10) = 0.;
        answer.at(1, 11) = d.at(1,1)*(-zjn1);
        answer.at(1, 12) = -d.at(1,1)*(yjn1);

        //Shear Y 1
        answer.at(2, 1) = 0;
        answer.at(2, 2) = d.at(2, 2);
        answer.at(2, 3) = 0.;
        answer.at(2, 4) = -d.at(2,2)*(zin1);
        answer.at(2, 5) = 0;
        answer.at(2, 6) = -d.at(2,2)*(-xin1);
        answer.at(2, 7) = 0.;
        answer.at(2, 8) = -d.at(2,2);
        answer.at(2, 9) = 0.;
        answer.at(2, 10) = -d.at(2,2)*(zjn1);
        answer.at(2, 11) = 0;
        answer.at(2, 12) = -d.at(2,2)*(-xjn1);

        //Shear Z 1
        answer.at(3, 1) = 0;
        answer.at(3, 2) = 0.;
        answer.at(3, 3) = d.at(3, 3);
        answer.at(3, 4) = -d.at(3,3)*(-zin1);
        answer.at(3, 5) = -d.at(3,3)*(xin1);
        answer.at(3, 6) = 0;
        answer.at(3, 7) = 0.;
        answer.at(3, 8) = 0;
        answer.at(3, 9) = -d.at(3, 3);
        answer.at(3, 10) = -d.at(3,3)*(-zjn1);
        answer.at(3, 11) = -d.at(3,3)*(xjn1);
        answer.at(3, 12) = 0;

        // Mx 1
        answer.at(4, 1) = 0;
        answer.at(4, 2) = d.at(2,2)*(-zin1);
        answer.at(4, 3) = -d.at(3,3)*(-yin1);
        answer.at(4, 4) = -d.at(3,3)*(-yin1*yin1)+d.at(2,2)*(zin1*zin1)+d.at(4,4);
        answer.at(4, 5) = -d.at(3,3)*(yin1*xin1);
        answer.at(4, 6) = d.at(2,2)*(-zin1*xin1);
        answer.at(4, 7) = 0.;
        answer.at(4, 8) = d.at(2,2)*(zin1);
        answer.at(4, 9) = -d.at(3,3)*(yjn1);
        answer.at(4, 10) = -d.at(3,3)*(-yin1*yjn1)+d.at(2,2)*(+zin1*zjn1)-d.at(4,4);
        answer.at(4, 11) = -d.at(3,3)*(yin1*xjn1);
        answer.at(4, 12) = d.at(2,2)*(-zin1*xjn1);

        // My 1
        answer.at(5, 1) = -d.at(1,1)*(-zin1);
        answer.at(5, 2) = 0;
        answer.at(5, 3) = d.at(3,3)*(-xin1);
        answer.at(5, 4) = d.at(3,3)*(-xin1*yin1);
        answer.at(5, 5) = -d.at(1,1)*(-zin1*zin1)+d.at(3,3)*(xin1*xin1)+d.at(5,5);
        answer.at(5, 6) = -d.at(1,1)*(+zin1*yin1);
        answer.at(5, 7) = -d.at(1,1)*(zin1);
        answer.at(5, 8) = 0;
        answer.at(5, 9) = d.at(3,3)*(xin1);
        answer.at(5, 10) = d.at(3,3)*(-xin1*yjn1);
        answer.at(5, 11) = -d.at(1,1)*(-zin1*zjn1)+d.at(3,3)*(xin1*xjn1)-d.at(5,5);
        answer.at(5, 12) = -d.at(1,1)*(+zin1*yjn1);


        // Mz 1
        answer.at(6, 1) = d.at(1,1)*(-yin1);
        answer.at(6, 2) = -d.at(2,2)*(-xin1);
        answer.at(6, 3) = 0;
        answer.at(6, 4) = -d.at(2,2)*(xin1*zin1);
        answer.at(6, 5) = d.at(1,1)*(-yin1*zin1);
        answer.at(6, 6) = d.at(1,1)*(yin1*yin1)-d.at(2,2)*(-xin1*xin1)+d.at(6,6);
        answer.at(6, 7) = d.at(1,1)*(yin1);
        answer.at(6, 8) = -d.at(2,2)*(xin1);
        answer.at(6, 9) = 0;
        answer.at(6, 10) = -d.at(2,2)*(xin1*zjn1);
        answer.at(6, 11) = d.at(1,1)*(-yin1*zjn1);
        answer.at(6, 12) = d.at(1,1)*(yin1*yjn1)-d.at(2,2)*(-xin1*xjn1)-d.at(6,6);



        //Axial 2
        answer.at(7, 1) = -d.at(1, 1);
        answer.at(7, 2) = 0.;
        answer.at(7, 3) = 0.;
        answer.at(7, 4) = 0.;
        answer.at(7, 5) = -d.at(1,1)*(-zin1);
        answer.at(7, 6) = d.at(1,1)*(yin1);
        answer.at(7, 7) = d.at(1, 1);
        answer.at(7, 8) = 0.;
        answer.at(7, 9) = 0.;
        answer.at(7, 10) = 0.;
        answer.at(7, 11) = d.at(1,1)*(-zjn1);
        answer.at(7, 12) = d.at(1,1)*(yjn1);

        //Shear Y 2
        answer.at(8, 1) = 0;
        answer.at(8, 2) = -d.at(2, 2);
        answer.at(8, 3) = 0.;
        answer.at(8, 4) = d.at(2,2)*(zin1);
        answer.at(8, 5) = 0;
        answer.at(8, 6) = d.at(2,2)*(-xin1);
        answer.at(8, 7) = 0.;
        answer.at(8, 8) = d.at(2,2);
        answer.at(8, 9) = 0.;
        answer.at(8, 10) = d.at(2,2)*(zjn1);
        answer.at(8, 11) = 0;
        answer.at(8, 12) = d.at(2,2)*(-xjn1);

        //Shear Z 2
        answer.at(9, 1) = 0;
        answer.at(9, 2) = 0.;
        answer.at(9, 3) = -d.at(3, 3);
        answer.at(9, 4) = d.at(3,3)*(-zin1);
        answer.at(9, 5) = d.at(3,3)*(xin1);
        answer.at(9, 6) = 0;
        answer.at(9, 7) = 0.;
        answer.at(9, 8) = 0;
        answer.at(9, 9) = d.at(3, 3);
        answer.at(9, 10) = d.at(3,3)*(-zjn1);
        answer.at(9, 11) = d.at(3,3)*(xjn1);
        answer.at(9, 12) = 0;

        // Mx 2
        answer.at(10, 1) = 0;
        answer.at(10, 2) = -d.at(2,2)*(zjn1);
        answer.at(10, 3) = d.at(3,3)*(yin1);
        answer.at(10, 4) = d.at(3,3)*(yjn1*yin1)-d.at(2,2)*(-zjn1*zin1)-d.at(4,4);
        answer.at(10, 5) = d.at(3,3)*(-yjn1*xin1);
        answer.at(10, 6) = -d.at(2,2)*(zjn1*xin1);
        answer.at(10, 7) = 0.;
        answer.at(10, 8) = -d.at(2,2)*(-zjn1);
        answer.at(10, 9) = d.at(3,3)*(-yjn1);
        answer.at(10, 10) = d.at(3,3)*(yjn1*yjn1)-d.at(2,2)*(-zjn1*zjn1)+d.at(4,4);
        answer.at(10, 11) = d.at(3,3)*(-yjn1*xjn1);
        answer.at(10, 12) = -d.at(2,2)*(+zjn1*xjn1);
        // My 2
        answer.at(11, 1) = d.at(1,1)*(zjn1);
        answer.at(11, 2) = 0;
        answer.at(11, 3) = -d.at(3,3)*(xjn1);
        answer.at(11, 4) = -d.at(3,3)*(xjn1*yin1);
        answer.at(11, 5) = d.at(1,1)*(zjn1*zin1)-d.at(3,3)*(-xjn1*xin1)-d.at(5,5);
        answer.at(11, 6) = d.at(1,1)*(-zjn1*yin1);
        answer.at(11, 7) = d.at(1,1)*(-zjn1);
        answer.at(11, 8) = 0;
        answer.at(11, 9) = -d.at(3,3)*(-xjn1);
        answer.at(11, 10) = -d.at(3,3)*(xjn1*yjn1);
        answer.at(11, 11) = d.at(1,1)*(zjn1*zjn1)-d.at(3,3)*(-xjn1*xjn1)+d.at(5,5);
        answer.at(11, 12) = d.at(1,1)*(-zjn1*yjn1);
        // Mz 2
        answer.at(12, 1) = -yjn1*d.at(1,1);
        answer.at(12, 2) = xjn1*d.at(2,2);
        answer.at(12, 3) = 0;
        answer.at(12, 4) = -xjn1*d.at(2,2)*(zin1);
        answer.at(12, 5) = yjn1*d.at(1,1)*(-zin1);
        answer.at(12, 6) = yjn1*d.at(1,1)*(yin1)-xjn1*d.at(2,2)*(-xin1)-d.at(6,6);
        answer.at(12, 7) = yjn1*d.at(1,1);
        answer.at(12, 8) = -xjn1*d.at(2,2);
        answer.at(12, 9) = 0;
        answer.at(12, 10) = -xjn1*d.at(2,2)*(zjn1);
        answer.at(12, 11) = yjn1*d.at(1,1)*(-zjn1);
        answer.at(12, 12) =  yjn1*d.at(1,1)*(yjn1)-xjn1*d.at(2,2)*(-xjn1)+d.at(6,6);
        answer.times(1. / this->length);
        return;
    }


    void
    LatticeFrame3dNL::computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
    // Computes the vector containing the strains at the Gauss point gp of
    // the receiver, at time step tStep. The nature of these strains depends
    // on the element's type.
    {


        FloatArray  uIncr, un1, un;
        FloatArray pi0ni0, pj0nj0, rinpi0ni0, rjnpj0nj0;
        FloatArray pinnin, pjnnjn, rin1pinnin, rjn1pjnnjn;
        FloatArray deltaPi, deltaPj, deltaP, answerold;
        this->computeVectorOf(VM_Incremental, tStep, uIncr);
        this->computeVectorOf(VM_Total, tStep, un1);
         un     = un1 - uIncr;

        this->length   = computeLength();
        double l1 = this->length * ( 1. - this->s ) / 2;
        double l2 = this->length * ( 1. + this->s ) / 2;
        LatticeMaterialStatus *lmatStat = dynamic_cast < LatticeMaterialStatus * > ( integrationRulesArray [ 0 ]->getIntegrationPoint(0)->giveMaterialStatus() );
        auto strain = lmatStat->giveLatticeStrain();

        //Coordinates of point ni^0
        FloatArrayF< 3 >ni0;
        ni0.at(1) = 0;
        ni0.at(2) = 0;
        ni0.at(3) = 0;
        //Coordinates of point Pi^0
        FloatArrayF< 3 >pi0;
        pi0.at(1) = l1;
        pi0.at(2) = 0;
        pi0.at(3) = 0;
        //Coordinates of point pj^0
        FloatArrayF< 3 >pj0;
        pj0.at(1) = l1;
        pj0.at(2) = 0;
        pj0.at(3) = 0;
        //Coordinates of point nj^0
        FloatArrayF< 3 >nj0;
        nj0.at(1) = l1+l2;
        nj0.at(2) = 0;
        nj0.at(3) = 0;
        //Coordinates of point ni^n
        FloatArrayF< 3 >nin;
        nin.at(1) =  ni0.at(1) + un.at(1);
        nin.at(2) =  ni0.at(2) + un.at(2);
        nin.at(3) =  ni0.at(3) + un.at(3);
        //Coordinates of point nj^n
        FloatArrayF< 3 >njn;
        njn.at(1) =  nj0.at(1) + un.at(7);
        njn.at(2) =  nj0.at(2) + un.at(8);
        njn.at(3) =  nj0.at(3) + un.at(9);
        //Coordinates of point ni^n+1
        FloatArrayF< 3 >nin1;
        nin1.at(1) =  nin.at(1) + uIncr.at(1);
        nin1.at(2) =  nin.at(2) + uIncr.at(2);
        nin1.at(3) =  nin.at(3) + uIncr.at(3);
        //Coordinates of point nj^n+1
        FloatArrayF< 3 >njn1;
        njn1.at(1) =  njn.at(1) + uIncr.at(7);
        njn1.at(2) =  njn.at(2) + uIncr.at(8);
        njn1.at(3) =  njn.at(3) + uIncr.at(9);
        //Rotation matrix (Rj^n) of the first element (Ni0-Pi0) using total displacement at step n
        FloatMatrixF< 3, 3 >rin;
        rin.at(1, 1) = cos( un.at(5) )*cos(un.at(6));
        rin.at(1, 2) = -cos( un.at(5) )*sin(un.at(6));
        rin.at(1, 3) = sin(un.at(5));
        //
        rin.at(2, 1) = cos(un.at(4))*sin(un.at(6))+sin(un.at(4))*sin(un.at(5))*cos(un.at(6));
        rin.at(2, 2) = -sin(un.at(4))*sin(un.at(5))*sin(un.at(6))+cos(un.at(4))*cos(un.at(6));
        rin.at(2, 3) = -sin(un.at(4))*cos(un.at(5));
        //
        rin.at(3, 1) = sin(un.at(4))*sin(un.at(6))-cos(un.at(4))*sin(un.at(5))*cos(un.at(6));
        rin.at(3, 2) = cos(un.at(4))*sin(un.at(5))*sin(un.at(6))+sin(un.at(4))*cos(un.at(6));
        rin.at(3, 3) = cos(un.at(4))*cos(un.at(5));
        //Rotation matrix (Rj^n) of the second element (Nj0-Pj0) using total displacement at step n
        FloatMatrixF< 3, 3 >rjn;
        rjn.at(1, 1) = cos( un.at(11) )*cos(un.at(12));
        rjn.at(1, 2) = -cos( un.at(11) )*sin(un.at(12));
        rjn.at(1, 3) = sin(un.at(11));
        //
        rjn.at(2, 1) = cos(un.at(10))*sin(un.at(12))+sin(un.at(10))*sin(un.at(11))*cos(un.at(12));
        rjn.at(2, 2) = -sin(un.at(10))*sin(un.at(11))*sin(un.at(12))+cos(un.at(10))*cos(un.at(12));
        rjn.at(2, 3) = -sin(un.at(10))*cos(un.at(11));
        //
        rjn.at(3, 1) = sin(un.at(10))*sin(un.at(12))-cos(un.at(10))*sin(un.at(11))*cos(un.at(12));
        rjn.at(3, 2) = cos(un.at(10))*sin(un.at(11))*sin(un.at(12))+sin(un.at(10))*cos(un.at(12));
        rjn.at(3, 3) = cos(un.at(10))*cos(un.at(11));
        //
        //Calculate the coordinates of point Pi^n using equation Pi^n=Ni^0+Ui^n+Ri^n(Pi^0-Ni^0).
        //First we calculate (Pi^0-Ni^0) and multiply it with Ri^n
        pi0ni0.beDifferenceOf(pi0, ni0);
        rinpi0ni0.beProductOf(rin, pi0ni0);
        //We add the result with (Ni^0+ui^n) to get Pi^n
        FloatArrayF< 3 >pin;
        pin.at(1) =  nin.at(1) + rinpi0ni0.at(1);
        pin.at(2) =  nin.at(2) + rinpi0ni0.at(2);
        pin.at(3) =  nin.at(3) + rinpi0ni0.at(3);
        //
        //
        //Calculate the coordinates of point Pj^n using equation Pj^n=Nj^0+uj^n+Rj^n(Pj^0-Nj^0).
        //First we calculate (Pj^0-Nj^0) and multiply it with Rj^n
        pj0nj0.beDifferenceOf(pj0, nj0);
        rjnpj0nj0.beProductOf(rjn, pj0nj0);
        //We add the result with (Nj^0+uj^n) to get Pj^n
        FloatArrayF< 3 >pjn;
        pjn.at(1) =  njn.at(1) + rjnpj0nj0.at(1);
        pjn.at(2) =  njn.at(2) + rjnpj0nj0.at(2);
        pjn.at(3) =  njn.at(3) + rjnpj0nj0.at(3);
        //
        //Rotation matrix (Ri^n+1) of the first element (Nin-Pin) using uIncr. displacement at step n+1 (cos(theta)=1, cos(theta)= theta)
        FloatMatrixF< 3, 3 >rin1;
        rin1.at(1, 1) =  1;
        rin1.at(1, 2) = - uIncr.at(6);
        rin1.at(1, 3) = uIncr.at(5);
        //
        rin1.at(2, 1) = uIncr.at(6)+uIncr.at(4)*uIncr.at(5);
        rin1.at(2, 2) = -uIncr.at(4)*uIncr.at(5)*uIncr.at(6)+1;
        rin1.at(2, 3) = -uIncr.at(4);
        //
        rin1.at(3, 1) = uIncr.at(4)*uIncr.at(6)-uIncr.at(5);
        rin1.at(3, 2) = uIncr.at(5)*uIncr.at(6)+uIncr.at(4);
        rin1.at(3, 3) = 1;
        //Rotation matrix (Rj^n+1) of the first element (Njn-Pjn) using uIncr. displacement at step n+1 (cos(theta)=1, cos(theta)= theta)
        FloatMatrixF< 3, 3 >rjn1;
        rjn1.at(1, 1) = 1;
        rjn1.at(1, 2) = -uIncr.at(12);
        rjn1.at(1, 3) = uIncr.at(11);
        //
        rjn1.at(2, 1) = uIncr.at(12)+uIncr.at(10)*uIncr.at(11);
        rjn1.at(2, 2) = -uIncr.at(10)*uIncr.at(11)*uIncr.at(12)+1;
        rjn1.at(2, 3) = -uIncr.at(10);
        //
        rjn1.at(3, 1) = uIncr.at(10)*uIncr.at(12)-uIncr.at(11);
        rjn1.at(3, 2) = uIncr.at(11)*uIncr.at(12)+uIncr.at(10);
        rjn1.at(3, 3) = 1;
        //

        //Calculate the coordinates of point Pi^n+1 using equation Pi^n=Ni^n+uIncr+Ri^n+1(Pi^n-Ni^n).
        //First we calculate (Pi^n-Ni^n) and multiply it with Ri^n+1
        pinnin.beDifferenceOf(pin, nin);
        rin1pinnin.beProductOf(rin1, pinnin);
        //We add the result with (Ni^n+uIncr) to get Pi^n+1
        FloatArrayF< 3 >pin1;
        pin1.at(1) =  nin1.at(1) + rin1pinnin.at(1);
        pin1.at(2) =  nin1.at(2) + rin1pinnin.at(2);
        pin1.at(3) =  nin1.at(3) + rin1pinnin.at(3);
        //
        //Calculate the coordinates of point Pj^n+1 using equation Pj^n=Nj^n+uIncr+Rj^n+1(Pj^n-Nj^n).
        //First we calculate (Pj^n-Nj^n) and multiply it with Rj^n+1
        pjnnjn.beDifferenceOf(pjn, njn);
        rjn1pjnnjn.beProductOf(rjn1, pjnnjn);
        //We add the result with (Nj^n+uIncr) to get Pj^n+1
        FloatArrayF< 3 >pjn1;
        pjn1.at(1) =  njn1.at(1) + rjn1pjnnjn.at(1);
        pjn1.at(2) =  njn1.at(2) + rjn1pjnnjn.at(2);
        pjn1.at(3) =  njn1.at(3) + rjn1pjnnjn.at(3);
        //
        //Calculate delta (Pi) and delta (Pj)
        deltaPi.beDifferenceOf(pin1, pin);
        deltaPj.beDifferenceOf(pjn1, pjn);
        //Calculate delta (P)
        //
       // FloatArray deltaP,
       deltaP.beDifferenceOf(deltaPj, deltaPi);

        //
        //Calculate delta (theta)
        FloatArrayF< 3 >deltaT;
        deltaT.at(1) = uIncr.at(10)-uIncr.at(4);
        deltaT.at(2) = uIncr.at(11)-uIncr.at(5);
        deltaT.at(3) = uIncr.at(12)-uIncr.at(6);
        //Incremental displacement of point Ni which equal to UiN_n+1 - UiN_n
      //  double ln = sqrt(pow(njn.at(1)-nin.at(1), 2) + pow(njn.at(2)-nin.at(2), 2) + pow(njn.at(3)-nin.at(3), 2) );
        //
        //
        answer.resize(6);
        answer.at(1) = deltaP.at(1);
        answer.at(2) = deltaP.at(2);
        answer.at(3) = deltaP.at(3);
        answer.at(4) = deltaT.at(1);
        answer.at(5) = deltaT.at(2);
        answer.at(6) = deltaT.at(3);
        answer.times(1. / this->length);
//         printf("Strain/n");
//         answer.printYourself();
         answer += strain;
        //
        // FloatMatrix b;
        //  FloatArray u;
        //  double l1 = this->length * ( 1. - this->s ) / 2;
        //  double l2 = this->length * ( 1. + this->s ) / 2;
        //  this->computeVectorOf(VM_Incremental, tStep, u);
        //  LatticeMaterialStatus *lmatStat = dynamic_cast < LatticeMaterialStatus * > ( integrationRulesArray [ 0 ]->getIntegrationPoint(0)->giveMaterialStatus() );
        //  auto strain = lmatStat->giveLatticeStrain();
//
        answerold.resize(6);
        answerold.at(1) = ( 1 - cos( uIncr.at(11) ) * cos( uIncr.at(12) ) ) * l2 + ( 1 - cos( uIncr.at(5) ) * cos( uIncr.at(6) ) ) * l1 + uIncr.at(7) - uIncr.at(1);
        answerold.at(2) = -( cos( uIncr.at(10) ) * sin( uIncr.at(12) ) + sin( uIncr.at(10) ) * sin( uIncr.at(11) ) * cos( uIncr.at(12) ) ) * l2 - ( cos( uIncr.at(4) ) * sin( uIncr.at(6) ) + sin( uIncr.at(4) ) * sin( uIncr.at(5) ) * cos( uIncr.at(6) ) ) * l1 + uIncr.at(8) - uIncr.at(2);
        answerold.at(3) = -( sin( uIncr.at(10) ) * sin( uIncr.at(12) ) - cos( uIncr.at(10) ) * sin( uIncr.at(11) ) * cos( uIncr.at(12) ) ) * l2 - ( sin( uIncr.at(4) ) * sin( uIncr.at(6) ) - cos( uIncr.at(4) ) * sin( uIncr.at(5) ) * cos( uIncr.at(6) ) ) * l1 + uIncr.at(9) - uIncr.at(3);
        answerold.at(4) = uIncr.at(10) - uIncr.at(4);
        answerold.at(5) = uIncr.at(11) - uIncr.at(5);
        answerold.at(6) = uIncr.at(12) - uIncr.at(6);
        answerold.times(1. / this->length);
      //   answerold.times(1. / ln);
//         printf("StrainOld/n");
//         answerold.printYourself();
//        answer += strain;
    }
    //
    void
    LatticeFrame3dNL::giveInternalForcesVector(FloatArray &answer,
                                               TimeStep *tStep, int useUpdatedGpRecord)
    {
        FloatMatrix b, bt, bf, d;
        FloatArray u, stress, strain;

        FloatArray  uIncr, un1;
        FloatArray pi0ni0, pj0nj0, rinpi0ni0, rjnpj0nj0;
        FloatArray pinnin, pjnnjn, rin1pinnin, rjn1pjnnjn;
        FloatArray deltaPi, deltaPj, deltaP, pnsfi, pnsfj, answerold;
        this->computeVectorOf(VM_Incremental, tStep, uIncr);
        this->computeVectorOf(VM_Total, tStep, un1);
        auto un     = un1 - uIncr;


     //   this->computeVectorOf(VM_Incremental, tStep, u);
        this->length   = computeLength();
        GaussPoint *gp = this->integrationRulesArray [ 0 ]->getIntegrationPoint(0);

        // Total stress
        this->computeStrainVector(strain, gp, tStep);
        this->computeStressVector(stress, strain, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

        // Old stresses
        LatticeMaterialStatus *lmatStat = dynamic_cast < LatticeMaterialStatus * > ( integrationRulesArray [ 0 ]->getIntegrationPoint(0)->giveMaterialStatus() );
        auto oldStress = lmatStat->giveLatticeStress();

        auto oldInternalForces = lmatStat->giveInternalForces();
        double l1 = this->length * ( 1. - this->s ) / 2;
        double l2 = this->length * ( 1. + this->s ) / 2;
        FloatArray incrementalStress;
        incrementalStress.beDifferenceOf(stress, oldStress);

        //Coordinates of point ni^0
        FloatArrayF< 3 >ni0;
        ni0.at(1) = 0;
        ni0.at(2) = 0;
        ni0.at(3) = 0;
        //Coordinates of point Pi^0
        FloatArrayF< 3 >pi0;
        pi0.at(1) = l1;
        pi0.at(2) = 0;
        pi0.at(3) = 0;
        //Coordinates of point pj^0
        FloatArrayF< 3 >pj0;
        pj0.at(1) = l1;
        pj0.at(2) = 0;
        pj0.at(3) = 0;
        //Coordinates of point nj^0
        FloatArrayF< 3 >nj0;
        nj0.at(1) = l1+l2;
        nj0.at(2) = 0;
        nj0.at(3) = 0;
        //Coordinates of point ni^n
        FloatArrayF< 3 >nin;
        nin.at(1) =  ni0.at(1) + un.at(1);
        nin.at(2) =  ni0.at(2) + un.at(2);
        nin.at(3) =  ni0.at(3) + un.at(3);
        //Coordinates of point nj^n
        FloatArrayF< 3 >njn;
        njn.at(1) =  nj0.at(1) + un.at(7);
        njn.at(2) =  nj0.at(2) + un.at(8);
        njn.at(3) =  nj0.at(3) + un.at(9);
        //Coordinates of point ni^n+1
        FloatArrayF< 3 >nin1;
        nin1.at(1) =  nin.at(1) + uIncr.at(1);
        nin1.at(2) =  nin.at(2) + uIncr.at(2);
        nin1.at(3) =  nin.at(3) + uIncr.at(3);
        //Coordinates of point nj^n+1
        FloatArrayF< 3 >njn1;
        njn1.at(1) =  njn.at(1) + uIncr.at(7);
        njn1.at(2) =  njn.at(2) + uIncr.at(8);
        njn1.at(3) =  njn.at(3) + uIncr.at(9);
        //Rotation matrix (Rj^n) of the first element (Ni0-Pi0) using total displacement at step n
        FloatMatrixF< 3, 3 >rin;
        rin.at(1, 1) = cos( un.at(5) )*cos(un.at(6));
        rin.at(1, 2) = -cos( un.at(5) )*sin(un.at(6));
        rin.at(1, 3) = sin(un.at(5));
        //
        rin.at(2, 1) = cos(un.at(4))*sin(un.at(6))+sin(un.at(4))*sin(un.at(5))*cos(un.at(6));
        rin.at(2, 2) = -sin(un.at(4))*sin(un.at(5))*sin(un.at(6))+cos(un.at(4))*cos(un.at(6));
        rin.at(2, 3) = -sin(un.at(4))*cos(un.at(5));
        //
        rin.at(3, 1) = sin(un.at(4))*sin(un.at(6))-cos(un.at(4))*sin(un.at(5))*cos(un.at(6));
        rin.at(3, 2) = cos(un.at(4))*sin(un.at(5))*sin(un.at(6))+sin(un.at(4))*cos(un.at(6));
        rin.at(3, 3) = cos(un.at(4))*cos(un.at(5));
        //Rotation matrix (Rj^n) of the second element (Nj0-Pj0) using total displacement at step n
        FloatMatrixF< 3, 3 >rjn;
        rjn.at(1, 1) = cos( un.at(11) )*cos(un.at(12));
        rjn.at(1, 2) = -cos( un.at(11) )*sin(un.at(12));
        rjn.at(1, 3) = sin(un.at(11));
        //
        //
        rjn.at(2, 1) = cos(un.at(10))*sin(un.at(12))+sin(un.at(10))*sin(un.at(11))*cos(un.at(12));
        rjn.at(2, 2) = -sin(un.at(10))*sin(un.at(11))*sin(un.at(12))+cos(un.at(10))*cos(un.at(12));
        rjn.at(2, 3) = -sin(un.at(10))*cos(un.at(11));
        //
        rjn.at(3, 1) = sin(un.at(10))*sin(un.at(12))-cos(un.at(10))*sin(un.at(11))*cos(un.at(12));
        rjn.at(3, 2) = cos(un.at(10))*sin(un.at(11))*sin(un.at(12))+sin(un.at(10))*cos(un.at(12));
        rjn.at(3, 3) = cos(un.at(10))*cos(un.at(11));
        //
        //Calculate the coordinates of point Pi^n using equation Pi^n=Ni^0+Ui^n+Ri^n(Pi^0-Ni^0).
        //First we calculate (Pi^0-Ni^0) and multiply it with Ri^n
        pi0ni0.beDifferenceOf(pi0, ni0);
        rinpi0ni0.beProductOf(rin, pi0ni0);
        //We add the result with (Ni^0+ui^n) to get Pi^n
        FloatArrayF< 3 >pin;
        pin.at(1) =  nin.at(1) + rinpi0ni0.at(1);
        pin.at(2) =  nin.at(2) + rinpi0ni0.at(2);
        pin.at(3) =  nin.at(3) + rinpi0ni0.at(3);
        //
        //
        //Calculate the coordinates of point Pj^n using equation Pj^n=Nj^0+uj^n+Rj^n(Pj^0-Nj^0).
        //First we calculate (Pj^0-Nj^0) and multiply it with Rj^n
        pj0nj0.beDifferenceOf(pj0, nj0);
        rjnpj0nj0.beProductOf(rjn, pj0nj0);
        //We add the result with (Nj^0+uj^n) to get Pj^n
        FloatArrayF< 3 >pjn;
        pjn.at(1) =  nj0.at(1) +un.at(7)+ rjnpj0nj0.at(1);
        pjn.at(2) =  nj0.at(2) +un.at(8)+ rjnpj0nj0.at(2);
        pjn.at(3) =  nj0.at(3) +un.at(9)+ rjnpj0nj0.at(3);
        //
        //Rotation matrix (Ri^n+1) of the first element (Nin-Pin) using uIncr. displacement at step n+1 (cos(theta)=1, cos(theta)= theta)
        FloatMatrixF< 3, 3 >rin1;
        rin1.at(1, 1) =  1;
        rin1.at(1, 2) = - uIncr.at(6);
        rin1.at(1, 3) = uIncr.at(5);
        //
        rin1.at(2, 1) = uIncr.at(6)+uIncr.at(4)*uIncr.at(5);
        rin1.at(2, 2) = -uIncr.at(4)*uIncr.at(5)*uIncr.at(6)+1;
        rin1.at(2, 3) = -uIncr.at(4);
        //
        rin1.at(3, 1) = uIncr.at(4)*uIncr.at(6)-uIncr.at(5);
        rin1.at(3, 2) = uIncr.at(5)*uIncr.at(6)+uIncr.at(4);
        rin1.at(3, 3) = 1;
        //Rotation matrix (Rj^n+1) of the first element (Njn-Pjn) using uIncr. displacement at step n+1 (cos(theta)=1, cos(theta)= theta)
        FloatMatrixF< 3, 3 >rjn1;
        rjn1.at(1, 1) = 1;
        rjn1.at(1, 2) = -uIncr.at(12);
        rjn1.at(1, 3) = uIncr.at(11);
        //
        rjn1.at(2, 1) = uIncr.at(12)+uIncr.at(10)*uIncr.at(11);
        rjn1.at(2, 2) = -uIncr.at(10)*uIncr.at(11)*uIncr.at(12)+1;
        rjn1.at(2, 3) = -uIncr.at(10);
        //
        rjn1.at(3, 1) = uIncr.at(10)*uIncr.at(12)-uIncr.at(11);
        rjn1.at(3, 2) = uIncr.at(11)*uIncr.at(12)+uIncr.at(10);
        rjn1.at(3, 3) = 1;
        //
        //Calculate the coordinates of point Pi^n+1 using equation Pi^n=Ni^n+uIncr+Ri^n+1(Pi^n-Ni^n).
        //First we calculate (Pi^n-Ni^n) and multiply it with Ri^n+1
        pinnin.beDifferenceOf(pin, nin);
        rin1pinnin.beProductOf(rin1, pinnin);
        //We add the result with (Ni^n+uIncr) to get Pi^n+1
        FloatArrayF< 3 >pin1;
        pin1.at(1) =  nin1.at(1) + rin1pinnin.at(1);
        pin1.at(2) =  nin1.at(2) + rin1pinnin.at(2);
        pin1.at(3) =  nin1.at(3) + rin1pinnin.at(3);
        //
        //Calculate the coordinates of point Pj^n+1 using equation Pj^n=Nj^n+uIncr+Rj^n+1(Pj^n-Nj^n).
        //First we calculate (Pj^n-Nj^n) and multiply it with Rj^n+1
        pjnnjn.beDifferenceOf(pjn, njn);
        rjn1pjnnjn.beProductOf(rjn1, pjnnjn);
        //We add the result with (Nj^n+uIncr) to get Pj^n+1
        FloatArrayF< 3 >pjn1;
        pjn1.at(1) =  njn1.at(1) + rjn1pjnnjn.at(1);
        pjn1.at(2) =  njn1.at(2) + rjn1pjnnjn.at(2);
        pjn1.at(3) =  njn1.at(3) + rjn1pjnnjn.at(3);
        //
        //Calculate Delta(Mi)={mi1, mi2, mi3}= -(Pi^n+1 - Ni^n+1)xDelta(Sf)-Delta(Sm)
        //Calculate (Pi^n+1-Ni^n+1) and multiply it with Delta(Sf)
        FloatArrayF< 3 >pin1nin1;
        pin1nin1.at(1)=pin1.at(1)-nin1.at(1);
        pin1nin1.at(2)=pin1.at(2)-nin1.at(2);
        pin1nin1.at(3)=pin1.at(3)-nin1.at(3);
        //
        FloatArrayF<3>DeltaSfi;
        DeltaSfi.at(1)=incrementalStress.at(1);
        DeltaSfi.at(2)=incrementalStress.at(2);
        DeltaSfi.at(3)=incrementalStress.at(3);
        //
        pnsfi.beVectorProductOf(  pin1nin1, DeltaSfi);
        //
        //Calculate Delta(Mj)={mj1, mj2, mj3}= -(Pj^n+1 - Nj^n+1)xDelta(Sf)+Delta(Sm)
        //Calculate (Pj^n+1-Nj^n+1) and multiply it with Delta(Sf)
        FloatArrayF< 3 >pjn1njn1;
        pjn1njn1.at(1)=pjn1.at(1)-njn1.at(1);
        pjn1njn1.at(2)=pjn1.at(2)-njn1.at(2);
        pjn1njn1.at(3)=pjn1.at(3)-njn1.at(3);
        //
        FloatArrayF<3>DeltaSfj;
        DeltaSfj.at(1)=incrementalStress.at(1);
        DeltaSfj.at(2)=incrementalStress.at(2);
        DeltaSfj.at(3)=incrementalStress.at(3);
        //
        pnsfj.beVectorProductOf( pjn1njn1, -DeltaSfj);
        //
        //
        answer.resize(12);
        answer.at(1) = -incrementalStress.at(1);
        answer.at(2) = -incrementalStress.at(2);
        answer.at(3) = -incrementalStress.at(3);
        answer.at(4) = -pnsfi.at(1)-incrementalStress.at(4);
        answer.at(5) = -pnsfi.at(2)-incrementalStress.at(5);
        answer.at(6) = -pnsfi.at(3)-incrementalStress.at(6);
        answer.at(7) = incrementalStress.at(1);
        answer.at(8) = incrementalStress.at(2);
        answer.at(9) = incrementalStress.at(3);
        answer.at(10) = -pnsfj.at(1)+incrementalStress.at(4);
        answer.at(11) = -pnsfj.at(2)+incrementalStress.at(5);
        answer.at(12) = -pnsfj.at(3)+incrementalStress.at(6);
        answer += oldInternalForces;
        printf("Force/n");
        answer.printYourself();

       lmatStat->letTempInternalForcesBe(answer);

//        answer.resize(12);
//        answer.at(1) = -incrementalStress.at(1);
//        answer.at(2) = -incrementalStress.at(2);
//        answer.at(3) = -incrementalStress.at(3);
//        answer.at(4) = +incrementalStress.at(2) * ( sin( uIncr.at(4) ) * sin( uIncr.at(6) ) - cos( uIncr.at(4) ) * sin( uIncr.at(5) ) * cos( uIncr.at(6) ) ) * l1 - incrementalStress.at(3) * ( cos( uIncr.at(4) ) * sin( uIncr.at(6) ) + sin( uIncr.at(4) ) * sin( uIncr.at(5) ) * cos( uIncr.at(6) ) ) * l1 - incrementalStress.at(4);
//        answer.at(5) = -incrementalStress.at(1) * ( sin( uIncr.at(4) ) * sin( uIncr.at(6) ) - cos( uIncr.at(4) ) * sin( uIncr.at(5) ) * cos( uIncr.at(6) ) ) * l1 + incrementalStress.at(3) * ( cos( uIncr.at(5) ) * cos( uIncr.at(6) ) ) * l1 - incrementalStress.at(5);
//        answer.at(6) = incrementalStress.at(1) * ( cos( uIncr.at(4) ) * sin( uIncr.at(6) ) + sin( uIncr.at(4) ) * sin( uIncr.at(5) ) * cos( uIncr.at(6) ) ) * l1 - incrementalStress.at(2) * ( cos( uIncr.at(5) ) * cos( uIncr.at(6) ) ) * l1 - incrementalStress.at(6);
//        answer.at(7) = incrementalStress.at(1);
//        answer.at(8) = incrementalStress.at(2);
//        answer.at(9) = incrementalStress.at(3);
//        answer.at(10) = incrementalStress.at(2) * ( sin( uIncr.at(10) ) * sin( uIncr.at(12) ) - cos( uIncr.at(10) ) * sin( uIncr.at(11) ) * cos( uIncr.at(12) ) ) * l2 - incrementalStress.at(3) * ( cos( uIncr.at(10) ) * sin( uIncr.at(12) ) + sin( uIncr.at(10) ) * sin( uIncr.at(11) ) * cos( uIncr.at(12) ) ) * l2 + incrementalStress.at(4);
//        answer.at(11) = -incrementalStress.at(1) * ( sin( uIncr.at(10) ) * sin( uIncr.at(12) ) - cos( uIncr.at(10) ) * sin( uIncr.at(11) ) * cos( uIncr.at(12) ) ) * l2 + incrementalStress.at(3) * ( cos( uIncr.at(11) ) * cos( uIncr.at(12) ) ) * l2 + incrementalStress.at(5);
//        answer.at(12) = incrementalStress.at(1) * ( cos( uIncr.at(10) ) * sin( uIncr.at(12) ) + sin( uIncr.at(10) ) * sin( uIncr.at(11) ) * cos( uIncr.at(12) ) ) * l2 - incrementalStress.at(2) * ( cos( uIncr.at(11) ) * cos( uIncr.at(12) ) ) * l2 + incrementalStress.at(6);
//        answer += oldInternalForces;
//        printf("OldForce/n");
//        answer.printYourself();
//        lmatStat->letTempInternalForcesBe(answer);
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
