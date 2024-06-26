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

#include "domain.h"
#include "../sm/Elements/Interfaces/bondlink1d.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"
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
#include "../sm/Elements/structuralelement.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"
#include "../sm/Materials/structuralmaterial.h"
#include "sm/CrossSections/structuralinterfacecrosssection.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(BondLink1d);

BondLink1d::BondLink1d(int n, Domain *aDomain) : StructuralElement(n, aDomain)
{
    numberOfDofMans = 2;
}


double
BondLink1d::computeVolumeAround(GaussPoint *aGaussPoint)
{
    //Returns artifical volume (bond area times bond length) so that general parts of post processing work (dissipated energy, etc.)
    //Later on this artificial volume is divided by the bond length again so that the correct bond area is used.
    return pow(this->bondLength, 2.) * this->bondDiameter * M_PI;
}


void
BondLink1d::computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the strain matrix of the receiver.
{

    //Assemble Bmatrix based on three rigid arm components
    answer.resize(1, 2);
    answer.zero();

    //Normal displacement jump in x-direction
    //First node
    answer.at(1, 1) = 1.;
    //Second node
    answer.at(1, 2) = -1.;

    return;
}

void
BondLink1d::computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                        TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    FloatMatrix d, b, bt, db;
    FloatArray u, slip;

    this->computeVectorOf(VM_Total, tStep, u);

    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    answer.clear();

    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    this->computeBmatrixAt(gp, b);
    bt.beTranspositionOf(b);

    if ( !this->isActivated(tStep) ) {
        slip.resize( StructuralMaterial::giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        slip.zero();
    }
    slip.beProductOf(b, u);

    answer.resize(2, 2);
    answer.zero();

    this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

    db.beProductOf(d, b);
    answer.beProductOf(bt, db);

    //Introduce integration of bond strength
    double area = this->computeVolumeAround(gp) / this->giveLength();
    answer.times(area);


    return;
}

void BondLink1d::computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    integrationRulesArray.resize(1);
    integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
    integrationRulesArray [ 0 ]->SetUpPointsOnLine(1, _1dMat);
}


int
BondLink1d::computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{

  if(geometryFlag ==0){
    answer.resize(3);



  //coordinates of the two nodes
  Node *nodeA, *nodeB;
  FloatArray coordsA(3), coordsB(3);

    //Order of nodes. Important, because continuum node does not have rotational DOFs.
    //Beam node
    nodeA  = this->giveNode(1);
    //Continuum node
    nodeB  = this->giveNode(2);


    //Calculate components of distance from continuum node to lattice node.
    for ( int i = 0; i < 3; i++ ) {
      answer.at(i + 1) =  (nodeA->giveCoordinate(i + 1)+nodeB->giveCoordinate(i + 1))/2.;
    }
    this->globalCentroid = answer;
    geometryFlag = 1;
    
  }
  else{
    answer = this->globalCentroid;
  }
  
    return 1;
}
 


double
BondLink1d::giveBondLength()
{
    return this->bondLength;
}



double
BondLink1d::giveBondDiameter()
{
    return this->bondDiameter;
}


void
BondLink1d::giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = { D_u };
}

void
BondLink1d::initializeFrom(InputRecord &ir)
{
    // first call parent
    StructuralElement::initializeFrom(ir);

    IR_GIVE_FIELD(ir, this->bondLength, _IFT_BondLink1d_length);

    IR_GIVE_FIELD(ir, this->bondDiameter, _IFT_BondLink1d_diameter);

    IR_GIVE_FIELD(ir, this->directionVector, _IFT_BondLink1d_dirvector);
}



void
BondLink1d::saveContext(DataStream &stream, ContextMode mode)
{
    StructuralElement::saveContext(stream, mode);
}


void
BondLink1d::restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralElement::restoreContext(stream, mode);
}


void
BondLink1d::giveInternalForcesVector(FloatArray &answer,
                                          TimeStep *tStep, int useUpdatedGpRecord)
//
// returns nodal representation of real internal forces - necessary only for
// non-linear analysis.
// if useGpRecord == 1 then data stored in gp->giveStressVector() are used
// instead computing stressVector through this->ComputeStressVector();
// this must be done after you want internal forces after element->updateYourself()
// has been called for the same time step.
//
{
    FloatMatrix b;
    FloatArray u, stress, strain;

    // This function can be quite costly to do inside the loops when one has many slave dofs.
    this->computeVectorOf(VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.clear();

    for ( GaussPoint *gp : * this->giveDefaultIntegrationRulePtr() ) {
        this->computeBmatrixAt(gp, b);

        if ( useUpdatedGpRecord == 1 ) {
            auto status = gp->giveMaterialStatus();
            StructuralMaterialStatus *matStat = dynamic_cast< StructuralMaterialStatus * >( status );
            if ( matStat ) {
                stress = matStat->giveStressVector();
            } else {
                StructuralInterfaceMaterialStatus *ms = static_cast< StructuralInterfaceMaterialStatus * >( status );
                stress = ms->giveTraction();
            }
        } else {
            if ( !this->isActivated(tStep) ) {
                strain.resize( StructuralMaterial::giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
                strain.zero();
            }
            strain.beProductOf(b, u);
            this->computeStressVector(stress, strain, gp, tStep);
        }

        // updates gp stress and strain record  acording to current
        // increment of displacement
        if ( stress.giveSize() == 0 ) {
            break;
        }

        // now every gauss point has real stress vector
        // compute nodal representation of internal forces using f = B^T*Sigma dV
        double dV = this->computeVolumeAround(gp) / this->giveLength();
        ;
        if ( stress.giveSize() == 6 ) {
            // It may happen that e.g. plane strain is computed
            // using the default 3D implementation. If so,
            // the stress needs to be reduced.
            // (Note that no reduction will take place if
            //  the simulation is actually 3D.)
            FloatArray stressTemp;
            StructuralMaterial::giveReducedSymVectorForm( stressTemp, stress, gp->giveMaterialMode() );
            answer.plusProduct(b, stressTemp, dV);
        } else {
            answer.plusProduct(b, stress, dV);
        }
    }

    // if inactive update state, but no contribution to global system
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}


double
BondLink1d::giveLength()
{
    //returns the bond length, not the length between the two nodes, which should be zero
    return this->bondLength;
}

void
BondLink1d::computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer = static_cast< StructuralInterfaceCrossSection * >( this->giveCrossSection() )->give1dStiffnessMatrix_Eng(rMode, gp, tStep);
}

void
BondLink1d::computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
  answer.resize(1);
  answer.at(1) = static_cast< StructuralInterfaceCrossSection * >( this->giveCrossSection() )->giveEngTraction_1d(strain.at(1), gp, tStep);
}



#ifdef __OOFEG

void
BondLink1d::drawYourself(oofegGraphicContext &gc, TimeStep *tStep)
{
    OGC_PlotModeType mode = gc.giveIntVarPlotMode();

    if ( mode == OGC_rawGeometry ) {
        this->drawRawGeometry(gc, tStep);
    } else if ( mode == OGC_deformedGeometry ) {
        this->drawDeformedGeometry(gc, tStep, DisplacementVector);
    } else if ( mode == OGC_eigenVectorGeometry ) {
        this->drawDeformedGeometry(gc, tStep, EigenVector);
    } else if ( mode == OGC_scalarPlot ) {
        this->drawScalar(gc, tStep);
    } else if ( mode == OGC_elemSpecial ) {
        this->drawSpecial(gc, tStep);
    } else {
        OOFEM_ERROR("unsupported mode");
    }
}



void BondLink1d::drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;

    WCRec p[ 2 ];  /* points */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);

    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}

void BondLink1d::drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    double defScale = gc.getDefScale();

    WCRec p[ 2 ];  /* points */

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);

    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, defScale);

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, defScale);

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}

#endif
} // end namespace oofem
