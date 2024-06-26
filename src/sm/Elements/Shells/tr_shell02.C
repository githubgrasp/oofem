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

#include "sm/Elements/Shells/tr_shell02.h"
#include "fei2dtrlin.h"
#include "contextioerr.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "node.h"
#include "load.h"

#ifdef __OOFEG
 #include "node.h"
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif


namespace oofem {
REGISTER_Element(TR_SHELL02);

IntArray TR_SHELL02 :: loc_plate = {3, 4, 5, 9, 10, 11, 15, 16, 17};
IntArray TR_SHELL02 :: loc_membrane = {1, 2, 6, 7, 8, 12, 13, 14, 18};

TR_SHELL02 :: TR_SHELL02(int n, Domain *aDomain) : StructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this), ZZErrorEstimatorInterface(this), SpatialLocalizerInterface(this),
    plate(std::make_unique<DKTPlate3d>(n, aDomain)),
    membrane(std::make_unique<TrPlanestressRotAllman3d>(n, aDomain))
{
    numberOfDofMans = 3;
}


void
TR_SHELL02 :: initializeFrom(InputRecord &ir)
{
    // proc tady neni return = this...   ??? termitovo
    StructuralElement :: initializeFrom(ir);
    plate->initializeFrom(ir);
    membrane->initializeFrom(ir);
}

void
TR_SHELL02 :: postInitialize()
{
    StructuralElement :: postInitialize();

    if ( plate->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints() != membrane->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints() ) {
        OOFEM_ERROR("incompatible integration rules detected");
    }
}

void
TR_SHELL02 :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    StructuralElement :: updateLocalNumbering(f);
    plate->updateLocalNumbering(f);
    membrane->updateLocalNumbering(f);
}

void TR_SHELL02 :: setCrossSection(int csIndx)
{
    StructuralElement :: setCrossSection(csIndx);
    plate->setCrossSection(csIndx);
    membrane->setCrossSection(csIndx);
}



void
TR_SHELL02 :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode, TimeStep *tStep)
//
// returns characteristics vector of receiver accordind to mtrx
//
{
    FloatArray aux;

    answer.resize(18);
    answer.zero();

    plate->giveCharacteristicVector(aux, mtrx, mode, tStep);
    if ( !aux.isEmpty() ) answer.assemble(aux, loc_plate);

    membrane->giveCharacteristicVector(aux, mtrx, mode, tStep);
    if ( !aux.isEmpty() ) answer.assemble(aux, loc_membrane);
}

void
TR_SHELL02 :: giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver accordind to mtrx
//
{
    FloatMatrix aux;

    answer.resize(18, 18);
    answer.zero();

    plate->giveCharacteristicMatrix(aux, mtrx, tStep);
    if ( aux.isNotEmpty() ) answer.assemble(aux, loc_plate);

    membrane->giveCharacteristicMatrix(aux, mtrx, tStep);
    if ( aux.isNotEmpty() ) answer.assemble(aux, loc_membrane);
}

bool
TR_SHELL02 :: giveRotationMatrix(FloatMatrix &answer)
{
    FloatMatrix aux1, aux2;
    int ncol;

    bool t1 = plate->giveRotationMatrix(aux1);
    bool t2 =  membrane->giveRotationMatrix(aux2);

    if ( t1 != t2 ) {
        OOFEM_ERROR("Transformation demand mismatch");
    }

    if ( t1 ) {
        ncol = aux1.giveNumberOfColumns();
        answer.resize(18, ncol);

        for ( int i = 1; i <= 9; i++ ) { // row index
            for ( int j = 1; j <= ncol; j++ ) {
                answer.at(loc_plate.at(i), j) = aux1.at(i, j);
            }
        }

        for ( int i = 1; i <= 9; i++ ) { // row index
            for ( int j = 1; j <= ncol; j++ ) {
                answer.at(loc_membrane.at(i), j) = aux2.at(i, j);
            }
        }
    }

    return t1;
}

void
TR_SHELL02 :: updateInternalState(TimeStep *tStep)
// Updates the receiver at end of step.
{
    plate->updateInternalState(tStep);
    membrane->updateInternalState(tStep);
}

void
TR_SHELL02 :: updateYourself(TimeStep *tStep)
{
    StructuralElement :: updateYourself(tStep);

    plate->updateYourself(tStep);
    membrane->updateYourself(tStep);
}


Interface *
TR_SHELL02 :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    } else if ( interface == ZZErrorEstimatorInterfaceType ) {
        return static_cast< ZZErrorEstimatorInterface * >(this);
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    }


    return NULL;
}

double
TR_SHELL02 :: computeVolumeAround(GaussPoint *gp)
{
    return plate->computeVolumeAround(gp);
}

void
TR_SHELL02 :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
{
        OOFEM_ERROR("This function is not implemented yet.");
}

int
TR_SHELL02 :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_ShellForceTensor || type == IST_ShellStrainTensor ||
        type == IST_ShellMomentTensor || type == IST_CurvatureTensor ) {
        FloatArray aux;
        GaussPoint *membraneGP = membrane->giveDefaultIntegrationRulePtr()->getIntegrationPoint(gp->giveNumber() - 1);
        GaussPoint *plateGP = plate->giveDefaultIntegrationRulePtr()->getIntegrationPoint(gp->giveNumber() - 1);

        plate->giveIPValue(answer, plateGP, type, tStep);
        membrane->giveIPValue(aux, membraneGP, type, tStep);
        answer.add(aux);
        return 1;
    } else {
        return StructuralElement :: giveIPValue(answer, gp, type, tStep);
    }
}


//
// The element interface required by NodalAveragingRecoveryModel
//
void
TR_SHELL02 :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                         InternalStateType type, TimeStep *tStep)
{
    this->giveIPValue(answer, NULL, type, tStep);
}




void
TR_SHELL02 :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
    FloatArray v, aux;
    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    fprintf( file, "element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );

    for ( GaussPoint *gp: *iRule ) {
        fprintf( file, "  GP %2d.%-2d :", iRule->giveNumber(), gp->giveNumber() );
        GaussPoint *membraneGP = membrane->giveDefaultIntegrationRulePtr()->getIntegrationPoint(gp->giveNumber() - 1);
        // Strain - Curvature
        plate->giveIPValue(v, gp, IST_ShellStrainTensor, tStep);
        membrane->giveIPValue(aux, membraneGP, IST_ShellStrainTensor, tStep);
        v.add(aux);

        fprintf(file, "  strains    ");
        for ( auto &val : v ) fprintf(file, " %.4e", val);

        plate->giveIPValue(v, gp, IST_CurvatureTensor, tStep);
        membrane->giveIPValue(aux, membraneGP, IST_CurvatureTensor, tStep);
        v.add(aux);

        fprintf(file, "\n              curvatures ");
        for ( auto &val : v ) fprintf(file, " %.4e", val);

        // Forces - Moments
        plate->giveIPValue(v, gp, IST_ShellForceTensor, tStep);
        membrane->giveIPValue(aux, membraneGP, IST_ShellForceTensor, tStep);
        v.add(aux);

        fprintf(file, "\n              stresses   ");
        for ( auto &val : v ) fprintf(file, " %.4e", val);

        plate->giveIPValue(v, gp, IST_ShellMomentTensor, tStep);
        membrane->giveIPValue(aux, membraneGP, IST_ShellMomentTensor, tStep);
        v.add(aux);

        fprintf(file, "\n              moments    ");
        for ( auto &val : v ) fprintf(file, " %.4e", val);

        fprintf(file, "\n");
    }
}


void
TR_SHELL02 :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralElement :: saveContext(stream, mode);
    this->plate->saveContext(stream, mode);
    this->membrane->saveContext(stream, mode);
}


void
TR_SHELL02 :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralElement :: restoreContext(stream, mode);
    this->plate->restoreContext(stream, mode);
    this->membrane->restoreContext(stream, mode);
}


IntegrationRule *
TR_SHELL02 :: ZZErrorEstimatorI_giveIntegrationRule()
{
    if ( !this->compositeIR ) {
        this->compositeIR = std::make_unique<GaussIntegrationRule>(1, this, 1, 12);
        this->compositeIR->SetUpPointsOnTriangle(plate->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints(), _3dShell);
    }
    return this->compositeIR.get();
}

void
TR_SHELL02 :: ZZErrorEstimatorI_computeLocalStress(FloatArray &answer, FloatArray &sig)
{
    // sig is global ShellForceMomentTensor
    FloatMatrix globTensor(3, 3);
    const FloatMatrix *GtoLRotationMatrix = plate->computeGtoLRotationMatrix();
    FloatMatrix LtoGRotationMatrix;

    answer.resize(8); // reduced, local form
    LtoGRotationMatrix.beTranspositionOf(* GtoLRotationMatrix);

    // Forces
    globTensor.at(1, 1) = sig.at(1);  //sxForce
    globTensor.at(1, 2) = sig.at(6);  //qxyForce
    globTensor.at(1, 3) = sig.at(5);  //qxzForce

    globTensor.at(2, 1) = sig.at(6);  //qxyForce
    globTensor.at(2, 2) = sig.at(2);  //syForce
    globTensor.at(2, 3) = sig.at(4);  //syzForce

    globTensor.at(3, 1) = sig.at(5);  //qxzForce
    globTensor.at(3, 2) = sig.at(4);  //syzForce
    globTensor.at(3, 3) = sig.at(3);  //szForce

    globTensor.rotatedWith(LtoGRotationMatrix);
    // Forces: now globTensoris transformed into local c.s

    // answer should be in reduced, local  form
    answer.at(1) = globTensor.at(1, 1); //sxForce
    answer.at(2) = globTensor.at(2, 2); //syForce
    answer.at(3) = globTensor.at(1, 2); //qxyForce
    answer.at(7) = globTensor.at(2, 3); //syzForce
    answer.at(8) = globTensor.at(1, 3); //qxzForce


    // Moments:
    globTensor.at(1, 1) = sig.at(7);  //mxForce
    globTensor.at(1, 2) = sig.at(12); //mxyForce
    globTensor.at(1, 3) = sig.at(11); //mxzForce

    globTensor.at(2, 1) = sig.at(12); //mxyForce
    globTensor.at(2, 2) = sig.at(8);  //myForce
    globTensor.at(2, 3) = sig.at(10); //myzForce

    globTensor.at(3, 1) = sig.at(11); //mxzForce
    globTensor.at(3, 2) = sig.at(10); //myzForce
    globTensor.at(3, 3) = sig.at(9);  //mzForce

    globTensor.rotatedWith(LtoGRotationMatrix);
    // now globTensoris transformed into local c.s

    answer.at(4)  = globTensor.at(1, 1); //mxForce
    answer.at(5)  = globTensor.at(2, 2); //myForce
    answer.at(6) = globTensor.at(1, 2); //mxyForce
}


void
TR_SHELL02 :: SpatialLocalizerI_giveBBox(FloatArray &bb0, FloatArray &bb1)
{
    FloatArray lt3, gt3; // global vector in the element thickness direction of lenght thickeness/2
    const FloatMatrix *GtoLRotationMatrix = plate->computeGtoLRotationMatrix();

    // setup vector in the element local cs. perpendicular to element plane of thickness/2 length
    lt3 = {0., 0., 1.}; //this->giveCrossSection()->give(CS_Thickness)/2.0; // HUHU
    // transform it to globa cs
    gt3.beTProductOf(* GtoLRotationMatrix, lt3);

    // use gt3 to construct element bounding box respecting true element volume

    FloatArray _c;

    for ( int i = 1; i <= this->giveNumberOfNodes(); ++i ) {
        const auto &coordinates = this->giveNode(i)->giveCoordinates();

        _c = coordinates;
        _c.add(gt3);
        if ( i == 1 ) {
            bb0 = bb1 = _c;
        } else {
            bb0.beMinOf(bb0, _c);
            bb1.beMaxOf(bb1, _c);
        }

        _c = coordinates;
        _c.subtract(gt3);
        bb0.beMinOf(bb0, _c);
        bb1.beMaxOf(bb1, _c);
    }
}

//
// io routines
//
#ifdef __OOFEG

void
TR_SHELL02 :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec p [ 3 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( this->giveMaterial()->isActivated(tStep) ) {
        EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
        EASValsSetColor( gc.getElementColor() );
        EASValsSetEdgeColor( gc.getElementEdgeColor() );
        EASValsSetEdgeFlag(true);
        EASValsSetFillStyle(FILL_SOLID);
        EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
        p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
        p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
        p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
        p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
        p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
        p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
        p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
        p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
        p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveCoordinate(3);

        go =  CreateTriangle3D(p);
        EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
        EGAttachObject(go, ( EObjectP ) this);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}

void
TR_SHELL02 :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    WCRec p [ 3 ];
    GraphicObj *go;
    double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( this->giveMaterial()->isActivated(tStep) ) {
        EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
        EASValsSetColor( gc.getDeformedElementColor() );
        EASValsSetEdgeColor( gc.getElementEdgeColor() );
        EASValsSetEdgeFlag(true);
        EASValsSetFillStyle(FILL_SOLID);
        EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
        p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
        p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
        p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, defScale);
        p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
        p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
        p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, defScale);
        p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale);
        p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale);
        p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(3, tStep, defScale);

        go =  CreateTriangle3D(p);
        EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}

void
TR_SHELL02 :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, indx, result = 0;
    WCRec p [ 3 ];
    GraphicObj *tr;
    FloatArray v1, v2, v3;
    double s [ 3 ], defScale;
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( !this->giveMaterial()->isActivated(tStep) ) {
        return;
    }

    if ( gc.giveIntVarMode() == ISM_recovered ) {
        result += this->giveInternalStateAtNode(v1, gc.giveIntVarType(), gc.giveIntVarMode(), 1, tStep);
        result += this->giveInternalStateAtNode(v2, gc.giveIntVarType(), gc.giveIntVarMode(), 2, tStep);
        result += this->giveInternalStateAtNode(v3, gc.giveIntVarType(), gc.giveIntVarMode(), 3, tStep);
    } else if ( gc.giveIntVarMode() == ISM_local ) {
        double tot_w = 0.;
        FloatArray a, v;
        for ( GaussPoint *gp: *plate->giveDefaultIntegrationRulePtr() ) {
            this->giveIPValue(a, gp, IST_ShellMomentTensor, tStep);
            v.add(gp->giveWeight(), a);
            tot_w += gp->giveWeight();
        }
        v.times(1. / tot_w);
        v1 = v;
        v2 = v;
        v3 = v;
    }

    indx = gc.giveIntVarIndx();

    s [ 0 ] = v1.at(indx);
    s [ 1 ] = v2.at(indx);
    s [ 2 ] = v3.at(indx);
    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( gc.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 3; i++ ) {
            if ( gc.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = gc.getDefScale();
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(3, tStep, defScale);
            } else {
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
            }
        }
        //     //EASValsSetColor(gc.getYieldPlotColor(ratio));
        gc.updateFringeTableMinMax(s, 3);
        tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
        EGWithMaskChangeAttributes(LAYER_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    }
}

#endif
} // end namespace oofem
