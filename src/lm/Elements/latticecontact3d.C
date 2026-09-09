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
 *               Copyright (C) 1993 - 2026   Borek Patzak
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
#include "Elements/latticecontact3d.h"
#include "Materials/latticematstatus.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "mathfem.h"
#include "Elements/latticestructuralelement.h"
#include "parametermanager.h"
#include "paramkey.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"
#include "sm/Materials/structuralmaterial.h"
#include "lm/CrossSections/latticecrosssection.h"

namespace oofem {
REGISTER_Element(LatticeContact3d);

ParamKey LatticeContact3d::IPK_LatticeContact3d_area("area");
ParamKey LatticeContact3d::IPK_LatticeContact3d_normal("normal");
ParamKey LatticeContact3d::IPK_LatticeContact3d_polycoords("polycoords");

LatticeContact3d :: LatticeContact3d(int n, Domain *aDomain) : LatticeStructuralElement(n, aDomain)
{
    numberOfDofMans = 2;
    geometryFlag = 0;
}

LatticeContact3d :: ~LatticeContact3d()
{}


double
LatticeContact3d :: computeVolumeAround(GaussPoint *aGaussPoint)
{
    // Artificial volume (contact area times unit length) so mass/post-processing work.
    return this->area;
}


double
LatticeContact3d :: giveLength()
{
    // Constitutive input is a displacement jump, not a strain, so no physical length enters.
    return 1.;
}


void
LatticeContact3d :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Maps the local translational DOFs of both nodes to the displacement jump.
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    // Three translational jump components (normal + two shears) from the six DOFs.
    answer.resize(3, 6);
    answer.zero();

    // Normal jump (local x)
    answer.at(1, 1) = -1.;
    answer.at(1, 4) =  1.;

    // Shear jump (local y)
    answer.at(2, 2) = -1.;
    answer.at(2, 5) =  1.;

    // Shear jump (local z)
    answer.at(3, 3) = -1.;
    answer.at(3, 6) =  1.;

    return;
}


void
LatticeContact3d :: giveGpCoordinates(FloatArray &coords, GaussPoint *gp)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }
    coords = this->globalCentroid;
}


void
LatticeContact3d :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    double density = static_cast< LatticeCrossSection * >( this->giveCrossSection() )->give('d', gp);
    double halfMass = density * computeVolumeAround(gp) / 2.;
    answer.resize(6, 6);
    answer.zero();
    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = halfMass;
    answer.at(4, 4) = answer.at(5, 5) = answer.at(6, 6) = halfMass;
}


void
LatticeContact3d :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    FloatMatrix d, b, bt, db;

    answer.clear();
    if ( !this->isActivated(tStep) ) {
        answer.resize(6, 6);
        answer.zero();
        return;
    }

    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    this->computeBmatrixAt(gp, b);
    bt.beTranspositionOf(b);

    this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);

    db.beProductOf(d, b);
    answer.beProductOf(bt, db);

    // Scale by contact area (traction per unit jump, no length integration).
    answer.times(this->area);
}


void
LatticeContact3d :: computeGaussPoints()
{
    integrationRulesArray.resize(1);
    integrationRulesArray [ 0 ].reset(new GaussIntegrationRule(1, this, 1, 3) );
    integrationRulesArray [ 0 ]->SetUpPointsOnLine(1, _3dLattice);
}


bool
LatticeContact3d :: computeGtoLRotationMatrix(FloatMatrix &answer)
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

    return true;
}


int
LatticeContact3d :: giveLocalCoordinateSystem(FloatMatrix &answer)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer = this->localCoordinateSystem;

    return 1;
}


void
LatticeContact3d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, D_w
    };
}


void
LatticeContact3d :: initializeFrom(const std::shared_ptr<InputRecord> &ir, int priority)
{
    ParameterManager &ppm = this->giveDomain()->elementPPM;
    LatticeStructuralElement :: initializeFrom(ir, priority);

    PM_UPDATE_PARAMETER(area, ppm, ir, this->number, IPK_LatticeContact3d_area, priority);
    PM_UPDATE_PARAMETER(normalVector, ppm, ir, this->number, IPK_LatticeContact3d_normal, priority);
    PM_UPDATE_PARAMETER(polygonCoords, ppm, ir, this->number, IPK_LatticeContact3d_polycoords, priority);
    numberOfPolygonVertices = (int) ( polygonCoords.giveSize() / 3 );
}


void
LatticeContact3d :: postInitialize()
{
    ParameterManager &ppm = this->giveDomain()->elementPPM;
    LatticeStructuralElement :: postInitialize();
    numberOfPolygonVertices = (int) ( polygonCoords.giveSize() / 3 );
    // Area and normal come from the polygon (Voronoi facet) when given; otherwise
    // 'area' is required and the normal falls back to input/node line.
    if ( numberOfPolygonVertices < 3 ) {
        PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_LatticeContact3d_area);
    }
}


int
LatticeContact3d :: computeGlobalCoordinates(Coordinates &answer, const FloatArray &lcoords)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer = this->globalCentroid;

    return 1;
}


void
LatticeContact3d :: computeGeometryProperties()
{
    Node *nodeA = this->giveNode(1);
    Node *nodeB = this->giveNode(2);

    FloatArray coordsA(3), coordsB(3);
    for ( int i = 1; i <= 3; i++ ) {
        coordsA.at(i) = nodeA->giveCoordinate(i);
        coordsB.at(i) = nodeB->giveCoordinate(i);
    }

    FloatArray normal(3), centroid(3);
    for ( int i = 1; i <= 3; i++ ) {
        centroid.at(i) = 0.5 * ( coordsA.at(i) + coordsB.at(i) );
    }

    if ( this->numberOfPolygonVertices >= 3 ) {
        // Contact-facet polygon (Voronoi facet, as lattice3d): Newell's method
        // gives the area-weighted plane normal in one pass; |sum| = 2 * area.
        FloatArray nw(3);
        nw.zero();
        centroid.zero();
        const int nv = this->numberOfPolygonVertices;
        for ( int k = 0; k < nv; k++ ) {
            const int kn = ( k + 1 ) % nv;
            const double xi = polygonCoords.at(3 * k + 1),  yi = polygonCoords.at(3 * k + 2),  zi = polygonCoords.at(3 * k + 3);
            const double xj = polygonCoords.at(3 * kn + 1), yj = polygonCoords.at(3 * kn + 2), zj = polygonCoords.at(3 * kn + 3);
            nw.at(1) += ( yi - yj ) * ( zi + zj );
            nw.at(2) += ( zi - zj ) * ( xi + xj );
            nw.at(3) += ( xi - xj ) * ( yi + yj );
            centroid.at(1) += xi;
            centroid.at(2) += yi;
            centroid.at(3) += zi;
        }
        const double twiceArea = nw.computeNorm();
        if ( twiceArea < 1e-20 ) {
            OOFEM_ERROR("LatticeContact3d: degenerate contact polygon (zero area).");
        }
        this->area = 0.5 * twiceArea;
        normal = nw;
        normal.times(1.0 / twiceArea);
        centroid.times(1.0 / nv);
    } else {
        // Explicit: normal user-supplied when given, otherwise the node line.
        if ( this->normalVector.giveSize() == 3 && this->normalVector.computeNorm() > 1e-12 ) {
            normal = this->normalVector;
        } else {
            normal.beDifferenceOf(coordsB, coordsA);
            if ( normal.computeNorm() < 1e-12 ) {
                OOFEM_ERROR("LatticeContact3d: contact normal is undefined; supply 'normal' or 'polycoords' for coincident nodes.");
            }
        }
        normal.normalize();
    }

    // Orient along node1->node2 when the nodes are distinct, so the compression
    // sign convention is independent of polygon winding.
    FloatArray axis(3);
    axis.beDifferenceOf(coordsB, coordsA);
    if ( axis.computeNorm() > 1e-12 && normal.dotProduct(axis) < 0. ) {
        normal.times(-1.0);
    }

    // Two axes spanning the plane orthogonal to the normal.
    FloatArray s(3), t(3);
    if ( normal.at(1) == 0 ) {
        s.at(1) = 0.;
        s.at(2) = normal.at(3);
        s.at(3) = -normal.at(2);
    } else if ( normal.at(2) == 0 ) {
        s.at(1) = normal.at(3);
        s.at(2) = 0.;
        s.at(3) = -normal.at(1);
    } else {
        s.at(1) = normal.at(2);
        s.at(2) = -normal.at(1);
        s.at(3) = 0.;
    }
    s.normalize();

    t.beVectorProductOf(normal, s);
    t.normalize();

    this->localCoordinateSystem.resize(3, 3);
    this->localCoordinateSystem.zero();
    for ( int i = 1; i <= 3; i++ ) {
        this->localCoordinateSystem.at(1, i) = normal.at(i);
        this->localCoordinateSystem.at(2, i) = s.at(i);
        this->localCoordinateSystem.at(3, i) = t.at(i);
    }

    this->globalCentroid = centroid;

    this->geometryFlag = 1;
}


void
LatticeContact3d :: saveContext(DataStream &stream, ContextMode mode)
{
    LatticeStructuralElement :: saveContext(stream, mode);
}


void
LatticeContact3d :: restoreContext(DataStream &stream, ContextMode mode)
{
    LatticeStructuralElement :: restoreContext(stream, mode);
}


void
LatticeContact3d :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix b, bt;
    FloatArray u, stress(3), strain(3);

    this->computeVectorOf(VM_Total, tStep, u);
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    answer.clear();

    for ( GaussPoint *gp: * this->giveDefaultIntegrationRulePtr() ) {
        LatticeMaterialStatus *matStat = static_cast< LatticeMaterialStatus * >( gp->giveMaterialStatus() );
        this->computeBmatrixAt(gp, b);
        bt.beTranspositionOf(b);

        if ( useUpdatedGpRecord == 1 ) {
            // Status keeps the full six-component stress; use the translational part.
            FloatArray fullStress = matStat->giveLatticeStress();
            stress.resize(3);
            stress.at(1) = fullStress.at(1);
            stress.at(2) = fullStress.at(2);
            stress.at(3) = fullStress.at(3);
        } else {
            if ( !this->isActivated(tStep) ) {
                strain.resize(3);
                strain.zero();
            } else {
                strain.beProductOf(b, u);
            }
            this->computeStressVector(stress, strain, gp, tStep);
        }

        if ( stress.giveSize() == 0 ) {
            break;
        }

        // f = B^T * traction * area
        answer.beProductOf(bt, stress);
        answer.times(this->area);
    }

    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}


void
LatticeContact3d :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer = static_cast< LatticeCrossSection * >( this->giveCrossSection() )->giveContactStiffnessMatrix(rMode, gp, tStep);
}


void
LatticeContact3d :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    answer = static_cast< LatticeCrossSection * >( this->giveCrossSection() )->giveLatticeContactStress(strain, gp, tStep);
}
} // end namespace oofem
