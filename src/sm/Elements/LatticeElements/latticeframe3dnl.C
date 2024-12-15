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
        rotationMatrixOne.resize(3,3);
        rotationMatrixOne.zero();
        rotationMatrixOne.at(1,1) = 1.;
        rotationMatrixOne.at(2,2) = 1.;
        rotationMatrixOne.at(3,3) = 1.;

        rotationMatrixTwo.resize(3,3);
        rotationMatrixTwo.zero();
        rotationMatrixTwo.at(1,1) = 1.;
        rotationMatrixTwo.at(2,2) = 1.;
        rotationMatrixTwo.at(3,3) = 1.;


        numberOfDofMans = 2;
    }

    LatticeFrame3dNL::~LatticeFrame3dNL()
    {}
    
    void
    LatticeFrame3dNL::computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
    //Returns the strain matrix of the receiver.
    {		

      answer.resize(6, 12);
      answer.zero();

        FloatArray coordA(3), coordB(3), l1(3), l2(3), help(3);       	
        coordA = giveNode(1)->giveCoordinates();
        coordB = giveNode(2)->giveCoordinates();

	help = coordB-coordA;
	l1 = help;
	l1.times( (1. + this->s ) / 2.);
	l2 = help;
	l2.times( (1. - this->s ) / 2.);

	
	FloatArray c1(3);
	c1.beProductOf(this->tempRotationMatrixOne,l1);
	
	FloatArray c2(3);
	c2.beProductOf(this->tempRotationMatrixTwo,l2);

        //Normal displacement jump in x-direction
        //First node
        answer.at(1, 1) = -1.;
        answer.at(1, 2) = 0.;
        answer.at(1, 3) = 0.;
        answer.at(1, 4) = 0.;
        answer.at(1, 5) = -c1.at(3);
        answer.at(1, 6) = c1.at(2);
        //Second node
        answer.at(1, 7) = 1.;
        answer.at(1, 8) = 0.;
        answer.at(1, 9) = 0.;
        answer.at(1, 10) = 0.;
        answer.at(1, 11) = -c2.at(3);
        answer.at(1, 12) = c2.at(2);

        //Shear displacement jump in y-plane
        //first node
        answer.at(2, 1) = 0.;
        answer.at(2, 2) = -1.;
        answer.at(2, 3) =  0.;
        answer.at(2, 4) = c1.at(3);
        answer.at(2, 5) = 0;
        answer.at(2, 6) = -c1.at(1);
        //Second node
        answer.at(2, 7) = 0.;
        answer.at(2, 8) = 1.;
        answer.at(2, 9) =  0.;
        answer.at(2, 10) = c2.at(3);
        answer.at(2, 11) = 0;
        answer.at(2, 12) = -c2.at(1);

        //Shear displacement jump in z-plane
        //first node
        answer.at(3, 1) = 0.;
        answer.at(3, 2) = 0.;
        answer.at(3, 3) = -1.;
        answer.at(3, 4) = -c1.at(2);
        answer.at(3, 5) = c1.at(1);
        answer.at(3, 6) = 0.;
        //Second node
        answer.at(3, 7) = 0.;
        answer.at(3, 8) = 0.;
        answer.at(3, 9) =  1.;
        answer.at(3, 10) = -c2.at(2);
        answer.at(3, 11) = c2.at(1);
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


void LatticeFrame3dNL::updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
    LatticeFrame3d::updateYourself(tStep);

    rotationMatrixOne = tempRotationMatrixOne;
    rotationMatrixTwo = tempRotationMatrixTwo;
    
}

void
LatticeFrame3dNL::updateRotationMatrices(TimeStep *tStep)
{
    // test if not previously done
    /* if ( tStep->giveSolutionStateCounter() == tempRotationCounter ) { */
    /*        return; */
    /* } */
  //  @TODO: One should only compute this once per step. But this does not work. 

  
          FloatMatrix  localR( 3, 3 ), help( 3, 3 ), globalRTotal( 3, 3 ), globalR( 3, 3 ),  globalRInc( 3, 3 ) , globalRIncX( 3, 3 ) , globalRIncY( 3, 3 ) , globalRIncZ( 3, 3 ) ;

        FloatArray uInc(12);

	this->computeVectorOf(VM_Incremental, tStep, uInc);

	
	FloatMatrix r(12, 12), rT(12, 12);

	Node 1
	FloatArray rotationOne(3);
        rotationOne.at(1) = uInc.at(4);
        rotationOne.at(2) = uInc.at(5);
        rotationOne.at(3) = uInc.at(6);

	computeGlobalRotationMatrix(globalRInc,rotationOne);
	this->tempRotationMatrixOne.beProductOf( globalRInc, this->rotationMatrixOne);

	Node 2
	FloatArray rotationTwo(3);
        rotationTwo.at(1) = uInc.at(10);
        rotationTwo.at(2) = uInc.at(11);
        rotationTwo.at(3) = uInc.at(12);

	computeGlobalRotationMatrix(globalRInc,rotationTwo);       
        this->tempRotationMatrixTwo.beProductOf( globalRInc, this->rotationMatrixTwo );

	tempRotationCounter = tStep->giveSolutionStateCounter();	
}


 
void
LatticeFrame3dNL::initForNewStep()
initializes receiver to new time step or can be used
if current time step must be restarted
{
    LatticeFrame3d::initForNewStep();

    this->tempRotationMatrixOne = this->rotationMatrixOne;

    this->tempRotationMatrixTwo = this->rotationMatrixTwo;
    
}
    
    void
      LatticeFrame3dNL::computeGlobalRotationMatrix(FloatMatrix &answer, FloatArray &rotation)
    {
      answer.resize(3,3);
      answer.zero();
      double thetaX, thetaY,thetaZ;    
      thetaX = rotation.at(1);
      thetaY = rotation.at(2);
      thetaZ = rotation.at(3);
      
      First x, then y and finally z
      answer.at(1,1) = cos(thetaY) * cos(thetaZ);
      answer.at(1,2) = cos(thetaZ) * sin(thetaY) * sin(thetaX) - sin(thetaZ) * cos(thetaX);
      answer.at(1,3) = cos(thetaZ) * sin(thetaY) * cos(thetaX) + sin(thetaZ) * sin(thetaX);
      
      answer.at(2,1) = cos(thetaY) * sin(thetaZ);
      answer.at(2,2) = sin(thetaZ) * sin(thetaY) * sin(thetaX) + cos(thetaZ) * cos(thetaX);
      answer.at(2,3) = sin(thetaZ) * sin(thetaY) * cos(thetaX) - cos(thetaZ) * sin(thetaX);
      
      answer.at(3,1) = - sin(thetaY);
      answer.at(3,2) = cos(thetaY) * sin(thetaX);
      answer.at(3,3) = cos(thetaY) * cos(thetaX);
      
      return;
    }



    void
    LatticeFrame3dNL::computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                              TimeStep *tStep)
    {

      updateRotationMatrices(tStep);

        FloatMatrix d, bt, db, b;
        this->length = computeLength();

        answer.resize(12, 12);
        answer.zero();
        this->computeBmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), b);

        this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

        Rotate constitutive stiffness matrix
        FloatMatrix r(6, 6), rT(6, 6), dR(6, 6), rTDR(6, 6);
        computeGtoLStrainRotationMatrix(r);
        rT.beTranspositionOf(r);

        dR.beProductOf(d, r);
        rTDR.beProductOf(rT, dR);

        db.beProductOf(rTDR, b);
        db.times(1. / length);
        bt.beTranspositionOf(b);
        answer.beProductOf(bt, db);

        return;
    }

    void
    LatticeFrame3dNL::computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
    {
      updateRotationMatrices(tStep);
      
        //Compute normalised displacement jumps in global coordinate system first and then rotate it to the local coordiante system.

        FloatArray u;
        this->computeVectorOf(VM_Total, tStep, u);

        FloatArray coordA(3), coordB(3), l1(3), l2(3), help(3);       	
        coordA = giveNode(1)->giveCoordinates();
        coordB = giveNode(2)->giveCoordinates();

	help = coordB-coordA;	
	l1 = help;
	l1.times(( 1. + this->s ) / 2.);
	l2 = help;
	l2.times(( 1. - this->s ) / 2.);

        this->length = computeLength();

	FloatArray c1(3);
	c1.beProductOf(this->tempRotationMatrixOne,l1);
	
	FloatArray c2(3);
	c2.beProductOf(this->tempRotationMatrixTwo,l2);

        answer.resize(6);
        answer.at(1) = u.at(7) - u.at(1) - c1.at(1) - c2.at(1) + l1.at(1) + l2.at(1);
        answer.at(2) = u.at(8) - u.at(2) - c1.at(2) - c2.at(2) + l1.at(2) + l2.at(2);
        answer.at(3) = u.at(9) - u.at(3) - c1.at(3) - c2.at(3) + l1.at(3) + l2.at(3);
        answer.at(4) = u.at(10) - u.at(4);
        answer.at(5) = u.at(11) - u.at(5);
        answer.at(6) = u.at(12) - u.at(6);

        answer.times(1. / this->length);

        //Rotate strain vector to local coordinate system
        FloatMatrix rotMatrix(6, 6);
        computeGtoLStrainRotationMatrix(rotMatrix);
        answer.rotatedWith(rotMatrix, 'n');

    }


    bool
    LatticeFrame3dNL::computeGtoLStrainRotationMatrix(FloatMatrix &answer)
    {
        FloatMatrix lcs;
        answer.resize(6, 6);
        answer.zero();

        this->giveLocalCoordinateSystem(lcs);
        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 3; j++ ) {
                answer.at(i, j) = lcs.at(i, j);
                answer.at(i + 3, j + 3) = lcs.at(i, j);
            }
        }

        return 1;
    }


    bool
    LatticeFrame3dNL::computeGtoLRotationMatrix(FloatMatrix &answer)
    {
        return false;
    }


    void
    LatticeFrame3dNL::giveInternalForcesVector(FloatArray &answer,
                                                TimeStep *tStep, int useUpdatedGpRecord)
    {
      updateRotationMatrices(tStep);
      
        FloatMatrix b, bt, bf, d;
        FloatArray u, stress, strain;
        this->computeVectorOf(VM_Total, tStep, u);
        this->length   = computeLength();
        GaussPoint *gp = this->integrationRulesArray [ 0 ]->getIntegrationPoint(0);

        FloatArray coordA(3), coordB(3), l1(3), l2(3), help(3);       	
        coordA = giveNode(1)->giveCoordinates();
        coordB = giveNode(2)->giveCoordinates();

	help = coordB-coordA;	
	l1 = help;
	l1.times(( 1. + this->s ) / 2.);
	l2 = help;
	l2.times((1. - this->s ) / 2.);

        this->length = computeLength();

        if ( useUpdatedGpRecord == 1 ) {
	  LatticeMaterialStatus *lmatStat = dynamic_cast < LatticeMaterialStatus * >( this->giveMaterial()->giveStatus(gp));
	  stress = lmatStat->giveLatticeStress();
        } else {
            if ( !this->isActivated(tStep) ) {
                strain.zero();
           }
        this->computeStrainVector(strain, gp, tStep);
        this->computeStressVector(stress, strain, gp, tStep);
    }


    FloatArray c1(3);
    c1.beProductOf(this->tempRotationMatrixOne,l1);

    FloatArray c2(3);
    c2.beProductOf(this->tempRotationMatrixTwo,l2);

    FloatMatrix rotMatrix(6, 6);
    computeGtoLStrainRotationMatrix(rotMatrix);
    stress.rotatedWith(rotMatrix, 't');
    
    answer.resize(12);
    answer.at(1) = -stress.at(1);
    answer.at(2) = -stress.at(2);
    answer.at(3) = -stress.at(3);
    answer.at(4) =  stress.at(2) * c1.at(3) - stress.at(3) * c1.at(2) - stress.at(4);
    answer.at(5) = -stress.at(1) * c1.at(3) + stress.at(3) * c1.at(1) - stress.at(5);
    answer.at(6) = stress.at(1) * c1.at(2) - stress.at(2) * c1.at(1) - stress.at(6);
    answer.at(7) = stress.at(1);
    answer.at(8) = stress.at(2);
    answer.at(9) = stress.at(3);
    answer.at(10) =  stress.at(2) * c2.at(3) - stress.at(3) * c2.at(2) + stress.at(4);
    answer.at(11) = -stress.at(1) * c2.at(3) + stress.at(3) * c2.at(1) + stress.at(5);
    answer.at(12) = stress.at(1) * c2.at(2) - stress.at(2) * c2.at(1) + stress.at(6);

    }
} // end namespace oofem
