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
#include "lattice3d.h"
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
#include "latticestructuralelement.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"
#include "sm/CrossSections/latticecrosssection.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "../sm/Materials/structuralmaterial.h"
#endif

namespace oofem {
REGISTER_Element(Lattice3d);

Lattice3d :: Lattice3d(int n, Domain *aDomain) : LatticeStructuralElement(n, aDomain)
{
    numberOfDofMans = 2;
    geometryFlag = 0;
    this->eccS = 0.;
    this->eccT = 0.;
}

Lattice3d :: ~Lattice3d()
{}


void
Lattice3d :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the strain matrix of the receiver.
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    //Assemble Bmatrix (used to compute strains and rotations}
    answer.resize(6, 12);
    answer.zero();

    //Normal displacement jump in x-direction
    //First node
    answer.at(1, 1) = -1.;
    answer.at(1, 2) = 0.;
    answer.at(1, 3) = 0.;
    answer.at(1, 4) = 0.;
    answer.at(1, 5) = -this->eccT;
    answer.at(1, 6) = this->eccS;
    //Second node
    answer.at(1, 7) = 1.;
    answer.at(1, 8) = 0.;
    answer.at(1, 9) = 0.;
    answer.at(1, 10) = 0.;
    answer.at(1, 11) = this->eccT;
    answer.at(1, 12) = -this->eccS;

    //Shear displacement jump in y-plane
    //first node
    answer.at(2, 1) = 0.;
    answer.at(2, 2) = -1.;
    answer.at(2, 3) =  0.;
    answer.at(2, 4) = this->eccT;
    answer.at(2, 5) = 0;
    answer.at(2, 6) = -this->length / 2.;
    //Second node
    answer.at(2, 7) = 0.;
    answer.at(2, 8) = 1.;
    answer.at(2, 9) =  0.;
    answer.at(2, 10) = -this->eccT;
    answer.at(2, 11) = 0;
    answer.at(2, 12) = -this->length / 2.;

    //Shear displacement jump in z-plane
    //first node
    answer.at(3, 1) = 0.;
    answer.at(3, 2) = 0.;
    answer.at(3, 3) = -1.;
    answer.at(3, 4) = -this->eccS;
    answer.at(3, 5) = this->length / 2.;
    answer.at(3, 6) = 0.;
    //Second node
    answer.at(3, 7) = 0.;
    answer.at(3, 8) = 0.;
    answer.at(3, 9) =  1.;
    answer.at(3, 10) = this->eccS;
    answer.at(3, 11) = this->length / 2.;
    answer.at(3, 12) = 0.;

    //Rotation around x-axis
    //First node
    answer.at(4, 1) = 0.;
    answer.at(4, 2) = 0;
    answer.at(4, 3) = 0.;
    answer.at(4, 4) = -1.;
    //answer.at(4, 4) = -sqrt(Ip / this->area);
    answer.at(4, 5) = 0.;
    answer.at(4, 6) = 0.;
    //Second node
    answer.at(4, 7) = 0.;
    answer.at(4, 8) = 0.;
    answer.at(4, 9) = 0.;
    answer.at(4, 10) = 1.;
    //    answer.at(4, 10) = sqrt(Ip / this->area);
    answer.at(4, 11) = 0.;
    answer.at(4, 12) = 0.;

    //Rotation around y-axis
    //First node
    answer.at(5, 1) = 0.;
    answer.at(5, 2) = 0.;
    answer.at(5, 3) = 0.;
    answer.at(5, 4) = 0.;
    answer.at(5, 5) = -1.;
    //    answer.at(5, 5) = -sqrt(I1 / this->area);
    answer.at(5, 6) = 0.;
    //Second node
    answer.at(5, 7) = 0.;
    answer.at(5, 8) = 0.;
    answer.at(5, 9) =  0.;
    answer.at(5, 10) = 0.;
    answer.at(5, 11) = 1.;
    //answer.at(5, 11) = sqrt(I1 / this->area);
    answer.at(5, 12) = 0.;

    //Rotation around z-axis
    //First node
    answer.at(6, 1) = 0.;
    answer.at(6, 2) = 0.;
    answer.at(6, 3) = 0.;
    answer.at(6, 4) = 0.;
    answer.at(6, 5) = 0.;
    answer.at(6, 6) = -1.;
    //    answer.at(6, 6) = -sqrt(I2 / this->area);
    //Second node
    answer.at(6, 7) = 0.;
    answer.at(6, 8) = 0.;
    answer.at(6, 9) =  0.;
    answer.at(6, 10) = 0.;
    answer.at(6, 11) = 0.;
    answer.at(6, 12) = 1.;
    //    answer.at(6, 12) = sqrt(I2 / this->area);
    //    answer.times(1. / this->length);
    return;
}


    void
    Lattice3d::giveInternalForcesVector(FloatArray &answer,
                                             TimeStep *tStep, int useUpdatedGpRecord)
    {
        FloatMatrix b, bt;
        FloatArray u, stress, strain;

        this->length = giveLength();

        this->computeVectorOf(VM_Total, tStep, u);

        if ( initialDisplacements ) {
            u.subtract(* initialDisplacements);
        }

        // zero answer will resize accordingly when adding first contribution
        answer.clear();

        this->computeBmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), b);

        bt.beTranspositionOf(b);

        if ( useUpdatedGpRecord == 1 ) {
            LatticeMaterialStatus *lmatStat = dynamic_cast < LatticeMaterialStatus * > ( integrationRulesArray [ 0 ]->getIntegrationPoint(0)->giveMaterialStatus() );
            stress = lmatStat->giveLatticeStress();
        } else {
            if ( !this->isActivated(tStep) ) {
                strain.zero();
            }
            strain.beProductOf(b, u);
            strain.times(1. / this->length);
            this->computeStressVector(stress, strain, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);
        }

        answer.beProductOf(bt, stress);
        if ( !this->isActivated(tStep) ) {
            answer.zero();
            return;
        }
    }


void
Lattice3d :: giveGPCoordinates(FloatArray &coords)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }
    coords.resize(3);
    coords = this->globalCentroid;
    return;
}

 
double
Lattice3d :: giveLength()
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    return this->length;
}


int
Lattice3d :: giveCrackFlag()
{
    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    LatticeMaterialStatus *status = static_cast< LatticeMaterialStatus * >( gp->giveMaterialStatus() );
    int crackFlag = 0;
    crackFlag = status->giveCrackFlag();

    return crackFlag;
}


double
Lattice3d :: giveCrackWidth()
{
    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    LatticeMaterialStatus *status = static_cast< LatticeMaterialStatus * >( gp->giveMaterialStatus() );
    double crackWidth = 0;
    crackWidth = status->giveCrackWidth();

    return crackWidth;
}

void
Lattice3d :: givePlasticStrain(FloatArray &plasticStrain)
{
    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    LatticeMaterialStatus *status = static_cast< LatticeMaterialStatus * >( this->giveMaterial()->giveStatus(gp) );
    plasticStrain = status->givePlasticLatticeStrain();
    return;
}

void
Lattice3d :: giveOldPlasticStrain(FloatArray &plasticStrain)
{
    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    LatticeMaterialStatus *status = static_cast< LatticeMaterialStatus * >( this->giveMaterial()->giveStatus(gp) );
    plasticStrain = status->giveOldPlasticLatticeStrain();
    return;
}

void
Lattice3d :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer = static_cast< LatticeCrossSection * >( this->giveCrossSection() )->give3dStiffnessMatrix(rMode, gp, tStep);
}

void
Lattice3d :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    answer = static_cast< LatticeCrossSection * >( this->giveCrossSection() )->giveLatticeStress3d(strain, gp, tStep);
}

void
Lattice3d :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                    TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    FloatMatrix d, bi, bj, bjt, dbj, dij;

    answer.resize(12, 12);
    answer.zero();
    this->computeBmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), bj);
    this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

    
    //    double volume = this->computeVolumeAround(integrationRulesArray [ 0 ]->getIntegrationPoint(0) );

    /* for ( int i = 1; i <= 6; i++ ) { */
    /*     d.at(i, i) *= volume; */
    /* } */

    dbj.beProductOf(d, bj);
    bjt.beTranspositionOf(bj);
    answer.beProductOf(bjt, dbj);
    answer.times(1./this->giveLength());

    return;
}

void Lattice3d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    integrationRulesArray.resize(1);
    integrationRulesArray [ 0 ].reset(new GaussIntegrationRule(1, this, 1, 3) );
    integrationRulesArray [ 0 ]->SetUpPointsOnLine(1, _3dLattice);
}



double Lattice3d :: giveArea() {
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    return this->area;
}


bool
Lattice3d :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    FloatMatrix lcs;
    int i, j;

    answer.resize(12, 12);
    answer.zero();

    this->giveLocalCoordinateSystem(lcs);
    for ( i = 1; i <= 3; i++ ) {
        for ( j = 1; j <= 3; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(i + 3, j + 3) = lcs.at(i, j);
            answer.at(i + 6, j + 6) = lcs.at(i, j);
            answer.at(i + 9, j + 9) = lcs.at(i, j);
        }
    }

    return 1;
}


double
Lattice3d :: giveNormalStress()
{
    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    LatticeMaterialStatus *status = static_cast< LatticeMaterialStatus * >( gp->giveMaterialStatus() );
    double normalStress = 0;
    normalStress = status->giveNormalLatticeStress();

    return normalStress;
}


int
Lattice3d :: giveLocalCoordinateSystem(FloatMatrix &answer)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer = this->localCoordinateSystem;

    return 1;
}


double
Lattice3d :: computeVolumeAround(GaussPoint *aGaussPoint)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    return this->area * this->length;
}

void
Lattice3d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, D_w, R_u, R_v, R_w
    };
}

void
Lattice3d :: initializeFrom(InputRecord &ir)
{
    LatticeStructuralElement :: initializeFrom(ir);
    
    minLength = 1.e-20;
    IR_GIVE_OPTIONAL_FIELD(ir, minLength, _IFT_Lattice3d_mlength);

    polygonCoords.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, polygonCoords, _IFT_Lattice3d_polycoords);
    numberOfPolygonVertices = polygonCoords.giveSize() / 3.;

    couplingFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, couplingFlag, _IFT_Lattice3d_couplingflag);

    IR_GIVE_OPTIONAL_FIELD(ir, couplingNumbers, _IFT_Lattice3d_couplingnumber);

    pressures.resize(numberOfPolygonVertices);
    pressures.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, pressures, _IFT_Lattice3d_pressures);
}


int
Lattice3d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer.resize(3);
    answer = this->globalCentroid;

    return 1;
}


void
Lattice3d :: computeGeometryProperties()
{
    //coordinates of the two nodes
    Node *nodeA, *nodeB;
    FloatArray coordsA(3), coordsB(3);

    nodeA  = this->giveNode(1);
    nodeB  = this->giveNode(2);

    for ( int i = 0; i < 3; i++ ) {
        coordsA.at(i + 1) =  nodeA->giveCoordinate(i + 1);
        coordsB.at(i + 1) =  nodeB->giveCoordinate(i + 1);
    }

    //Construct an initial temporary local coordinate system
    FloatArray s(3), t(3);

    //Calculate normal vector
    this->normal.resize(3);
    for ( int i = 0; i < 3; i++ ) {
        this->normal.at(i + 1) = coordsB.at(i + 1) - coordsA.at(i + 1);
    }

    this->length  = sqrt(pow(normal.at(1), 2.) + pow(normal.at(2), 2.) + pow(normal.at(3), 2.) );

    // Compute midpoint
    this->midPoint.resize(3);
    for ( int i = 0; i < 3; i++ ) {
        this->midPoint.at(i + 1) = 0.5 * ( coordsB.at(i + 1) + coordsA.at(i + 1) );
    }

    for ( int i = 0; i < 3; i++ ) {
        this->normal.at(i + 1) /= length;
    }

    computeCrossSectionProperties();

    this->geometryFlag = 1;
    
    return;
}

 
void
  Lattice3d :: computeCrossSectionProperties() {
   //@Todo: Extend this function so that it can deal with beams and shells.
   //Currently, polygon vertices are read in which are then used to compute cross-sectional properties.
   //What we now want is to distinguish three cases: 1) Frame element. No polygon vertices are provided. All properties are read from cross-section, 2) Shell element. Two polygon vertices are provided in one direction is read from two nodes. Thickness is provided by cross-section, 3) Full 3D. Three or more polygon vertices are provided. 
  if ( this->numberOfPolygonVertices > 0 ) { //Either Full 3D or shell
     
     if(this->numberOfPolygonVertices == 1){
       OOFEM_ERROR("Only one cross-section is provided. Check meshing approach.\n");
     }
     else if (this->numberOfPolygonVertices == 2){ //Shell case. Exactly two vertices and thickness are provided.
      FloatArray th(1);
      this->giveCrossSection()->give(CS_Thickness, th, this);
      double t = th.at(1);      
      if(t <= 0.){
	OOFEM_ERROR("Two cross-section vertices are provided but no valid shell thickness\n");
      }
      
      // Preserve the two given mid-line points (global coords) BEFORE resize
      FloatArray m1(3), m2(3);
      m1.at(1) = polygonCoords.at(1);
      m1.at(2) = polygonCoords.at(2);
      m1.at(3) = polygonCoords.at(3);
      m2.at(1) = polygonCoords.at(4);
      m2.at(2) = polygonCoords.at(5);
      m2.at(3) = polygonCoords.at(6);

      // e = element axis (assumed normal to cross-section), already normalized
      //      const FloatArray &e = this->normal;

    // w = direction along the provided mid-line, projected to be perpendicular to e
      FloatArray w = m2; w.subtract(m1);
      double we = w.dotProduct(this->normal);
      FloatArray proj = this->normal;
      proj.times(we);
      w.subtract(proj);

    if ( w.computeNorm() < 1e-12 ) {
      OOFEM_ERROR("Shell cross-section mid-line is (nearly) parallel to element axis; width direction cannot be formed.\n");
    }
    w.normalize();

    // n = thickness direction in the cross-section plane
    FloatArray n(3);
    n.beVectorProductOf(this->normal, w);
    if ( n.computeNorm() < 1e-12 ) {
        OOFEM_ERROR("Thickness direction cannot be formed (degenerate cross product).\n");
    }
    n.normalize();

    // Build 4 corners: (m1±t/2 n, m2±t/2 n)
    const double h = 0.5 * t;
    FloatArray off = n; off.times(h);

    FloatArray c1 = m1;
    c1.add(off);
    FloatArray c2 = m2;
    c2.add(off);
    FloatArray c3 = m2;
    c3.subtract(off);
    FloatArray c4 = m1;
    c4.subtract(off);

    // Overwrite polygonCoordinates with 4 vertices (global coords)
    this->numberOfPolygonVertices = 4;
    this->polygonCoords.resize(12);

    //C1
    for(int i = 1;i<=3;i++){
      this->polygonCoords.at(i) = c1.at(i);
    }

    //C2
    for(int i = 1;i<=3;i++){
      this->polygonCoords.at(3+i) = c2.at(i);
    }
    
    //C3
    for(int i = 1;i<=3;i++){
      this->polygonCoords.at(6+i) = c3.at(i);
    }

    //C4
    for(int i = 1;i<=3;i++){
      this->polygonCoords.at(9+i) = c4.at(i);
    }
        
     }

     //      Debug
    
     // --- Compute section reference point x0 on the element axis from polygonCoords ---
     // coordsA must exist here; if not, rebuild it from nodeA as you do elsewhere.
     Node *nodeA = this->giveNode(1);
     FloatArray coordsA(3);
     for ( int i = 0; i < 3; ++i ) {
       coordsA.at(i+1) = nodeA->giveCoordinate(i+1);
     }

     FloatArray x0 = coordsA;
     double alphaMean = 0.0;

     for ( int k = 0; k < this->numberOfPolygonVertices; ++k ) {
       FloatArray v(3);
       v.at(1) = this->polygonCoords.at(3*k + 1);
       v.at(2) = this->polygonCoords.at(3*k + 2);
       v.at(3) = this->polygonCoords.at(3*k + 3);
       
       v.subtract(coordsA);
       alphaMean += v.dotProduct(this->normal);
     }
     alphaMean /= double(this->numberOfPolygonVertices);

     FloatArray tmp = this->normal;
     tmp.times(alphaMean);
     x0.add(tmp);
// --- end x0 ---
    
    // Compute the cross-section properties from the geometry provided by the vertices.

    //Construct two perpendicular axis so that n is normal to the plane which they create
    //Check, if one of the components of the normal-direction is zero
    FloatArray s(3), t(3);
    if ( this->normal.at(1) == 0 ) {
        s.at(1) = 0.;
        s.at(2) = this->normal.at(3);
        s.at(3) = -this->normal.at(2);
    } else if ( this->normal.at(2) == 0 ) {
        s.at(1) = this->normal.at(3);
        s.at(2) = 0.;
        s.at(3) = -this->normal.at(1);
    } else {
        s.at(1) = this->normal.at(2);
        s.at(2) = -this->normal.at(1);
        s.at(3) = 0.;
    }

    s.normalize();

    t.beVectorProductOf(this->normal, s);
    t.normalize();

    //Set up rotation matrix
    FloatMatrix lcs(3, 3);

    for ( int i = 1; i <= 3; i++ ) {
        lcs.at(1, i) = this->normal.at(i);
        lcs.at(2, i) = s.at(i);
        lcs.at(3, i) = t.at(i);
    }


    //Calculate the local coordinates of the polygon vertices
    FloatArray help(3), test(3);
    FloatArray lpc(3 * this->numberOfPolygonVertices);
    for ( int k = 0; k < this->numberOfPolygonVertices; k++ ) {
        for ( int n = 0; n < 3; n++ ) {
            help(n) = this->polygonCoords(3 * k + n);
        }


	FloatArray rel = help;
	rel.subtract(x0);
	test.beProductOf(lcs, rel);
 
	//        test.beProductOf(lcs, help);
        for ( int n = 0; n < 3; n++ ) {
            lpc(3 * k + n) = test(n);
        }
    }

    this->area = 0.;

    for ( int k = 0; k < numberOfPolygonVertices; k++ ) {
        if ( k < numberOfPolygonVertices - 1 ) {
            this->area += lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 2) - lpc(3 * ( k + 1 ) + 1) * lpc(3 * k + 2);
        } else {   //Back to zero for n+1
            this->area += lpc(3 * k + 1) * lpc(2) - lpc(1) * lpc(3 * k + 2);
        }
    }

    this->area *= 0.5;

    FloatArray tempCoords(3 * numberOfPolygonVertices);
    if ( this->area < 0 ) { //Set area to a positive value and rearrange the coordinate entries
        this->area *= -1.;

        for ( int k = 0; k < numberOfPolygonVertices; k++ ) {
            for ( int m = 0; m < 3; m++ ) {
                tempCoords.at(3 * k + m + 1) = this->polygonCoords.at(3 * ( numberOfPolygonVertices - k - 1 ) + m + 1);
            }
        }

        this->polygonCoords = tempCoords;

        // Calculate again local co-ordinate system for different order
        for ( int k = 0; k < this->numberOfPolygonVertices; k++ ) {
            for ( int n = 0; n < 3; n++ ) {
                help(n) = this->polygonCoords(3 * k + n);
            }

	    FloatArray rel = help;
	    rel.subtract(x0);
	    test.beProductOf(lcs, rel);
	    
	    //            test.beProductOf(lcs, help);
	    for ( int n = 0; n < 3; n++ ) {
                lpc(3 * k + n) = test(n);
            }
        }
    }

    if ( this->area < pow(minLength, 2.) ) {
        this->area = pow(minLength, 2.);
    }

    //Calculate centroids
    centroid.resize(3);
    for ( int k = 0; k < this->numberOfPolygonVertices; k++ ) {
        if ( k < this->numberOfPolygonVertices - 1 ) {
            centroid.at(2) += ( lpc(3 * k + 1) + lpc(3 * ( k + 1 ) + 1) ) * ( lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 2) - lpc(3 * ( k + 1 ) + 1) * lpc(3 * k + 2) );
            centroid.at(3) += ( lpc(3 * k + 2) + lpc(3 * ( k + 1 ) + 2) ) * ( lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 2) - lpc(3 * ( k + 1 ) + 1) * lpc(3 * k + 2) );
        } else { //Back to zero for n+1
            centroid.at(2) += ( lpc(3 * k + 1) + lpc(1) ) * ( lpc(3 * k + 1) * lpc(2) - lpc(1) * lpc(3 * k + 2) );
            centroid.at(3) += ( lpc(3 * k + 2) + lpc(2) ) * ( lpc(3 * k + 1) * lpc(2) - lpc(1) * lpc(3 * k + 2) );
        }
    }

    centroid.times(1. / ( 6. * this->area ) );

    //    centroid.at(1) = lpc.at(1); //The first component of all lpcs should be the same
    centroid.at(1) = 0.0;
    
    //Shift coordinates to centroid
    for ( int k = 0; k < this->numberOfPolygonVertices; k++ ) {
        for ( int l = 0; l < 3; l++ ) {
            lpc(3 * k + l) -= centroid(l);
        }
    }

    //Compute second moments of area.
    //This is for the temporary coordinate system
    double Ixx = 0.;
    double Iyy = 0.;
    double Ixy = 0.;
    double a;

    for ( int k = 0; k < this->numberOfPolygonVertices; k++ ) {
        if ( k < this->numberOfPolygonVertices - 1 ) {
            a = lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 2) - lpc(3 * ( k + 1 ) + 1) * lpc(3 * k + 2);

            Ixx += ( ( pow(lpc(3 * k + 2), 2.) + lpc(3 * k + 2) * lpc(3 * ( k + 1 ) + 2) + pow(lpc(3 * ( k + 1 ) + 2), 2.) ) * a ) / 12.;

            Iyy += ( ( pow(lpc(3 * k + 1), 2.) + lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 1) + pow(lpc(3 * ( k + 1 ) + 1), 2.) ) * a ) / 12.;

            Ixy += ( ( lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 2) + 2. * lpc(3 * k + 1) * lpc(3 * k + 2) +
                       2 * lpc(3 * ( k + 1 ) + 1) * lpc(3 * ( k + 1 ) + 2) + lpc(3 * ( k + 1 ) + 1) * lpc(3 * k + 2) ) * a ) / 24.;
        } else {   //Back to zero for n+1
            a = lpc(3 * k + 1) * lpc(2) - lpc(1) * lpc(3 * k + 2);

            Ixx += ( ( pow(lpc(3 * k + 2), 2.) + lpc(3 * k + 2) * lpc(2) + pow(lpc(2), 2.) ) * a ) / 12.;

            Iyy += ( ( pow(lpc(3 * k + 1), 2.) + lpc(3 * k + 1) * lpc(1) + pow(lpc(1), 2.) ) * a ) / 12.;

            Ixy += ( ( lpc(3 * k + 1) * lpc(2) + 2. * lpc(3 * k + 1) * lpc(3 * k + 2) +
                       2 * lpc(1) * lpc(2) + lpc(1) * lpc(3 * k + 2) ) * a ) / 24.;
        }
    }

    //Compute main axis of the cross-section
    double angleChange = 0.;
    double sum = fabs(Ixx + Iyy);
    double pi = 3.14159265;
    if ( ( fabs(Ixx - Iyy) / sum > 1.e-6 ) && fabs(Ixy) / sum > 1.e-6 ) {
        angleChange = 0.5 * atan(-2 * Ixy / ( Ixx - Iyy ) );
    } else if ( ( fabs(Ixx - Iyy) / sum < 1.e-6 ) && fabs(Ixy) / sum > 1.e-6 ) {
        angleChange = pi / 4.;
    }

    if ( Iyy > Ixx ) {
        angleChange = angleChange + pi / 2.;
    }

    //Moment of inertias saved in the element

    this->I1 = ( Ixx + Iyy ) / 2. + sqrt(pow( ( Ixx - Iyy ) / 2., 2. ) + pow(Ixy, 2.) );
    this->I2 = ( Ixx + Iyy ) / 2. - sqrt(pow( ( Ixx - Iyy ) / 2., 2. ) + pow(Ixy, 2.) );

    this->Ip = I1 + I2;

    //Rotation around normal axis by angleChange

    FloatMatrix rotationChange(3, 3);
    rotationChange.zero();

    rotationChange.at(1, 1) = 1.;
    rotationChange.at(2, 2) = cos(angleChange);
    rotationChange.at(2, 3) = -sin(angleChange);

    rotationChange.at(3, 2) = sin(angleChange);
    rotationChange.at(3, 3) = cos(angleChange);

    this->localCoordinateSystem.beProductOf(rotationChange, lcs);

    //Calculate the polygon vertices in the new coordinate system
    for ( int k = 0; k < this->numberOfPolygonVertices; k++ ) {
        for ( int n = 0; n < 3; n++ ) {
            help(n) = this->polygonCoords(3 * k + n);
        }


	FloatArray rel = help;
	rel.subtract(x0);
	test.beProductOf(this->localCoordinateSystem, rel);
	//        test.beProductOf(this->localCoordinateSystem, help);
        for ( int n = 0; n < 3; n++ ) {
            lpc(3 * k + n) = test(n);
        }
    }

    //Calculate centroid again in local coordinate system
    centroid.zero();
    for ( int k = 0; k < this->numberOfPolygonVertices; k++ ) {
        if ( k < this->numberOfPolygonVertices - 1 ) {
            centroid.at(2) += ( lpc(3 * k + 1) + lpc(3 * ( k + 1 ) + 1) ) * ( lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 2) - lpc(3 * ( k + 1 ) + 1) * lpc(3 * k + 2) );
            centroid.at(3) += ( lpc(3 * k + 2) + lpc(3 * ( k + 1 ) + 2) ) * ( lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 2) - lpc(3 * ( k + 1 ) + 1) * lpc(3 * k + 2) );
        } else {   //Back to zero for n+1
            centroid.at(2) += ( lpc(3 * k + 1) + lpc(1) ) * ( lpc(3 * k + 1) * lpc(2) - lpc(1) * lpc(3 * k + 2) );
            centroid.at(3) += ( lpc(3 * k + 2) + lpc(2) ) * ( lpc(3 * k + 1) * lpc(2) - lpc(1) * lpc(3 * k + 2) );
        }
    }

    centroid.times(1. / ( 6. * this->area ) );

    //    centroid.at(1) = lpc.at(1); //The first component of all lpcs should be the same
    centroid.at(1) = 0.0;


    this->eccS = centroid.at(2);
    this->eccT = centroid.at(3);
    
    /* FloatArray midPointLocal(3); */
    /* midPointLocal.beProductOf(this->localCoordinateSystem, midPoint); */

    /* //eccentricities stored in the element */
    /* this->eccS = centroid.at(2) - midPointLocal.at(2); */
    /* this->eccT = centroid.at(3) - midPointLocal.at(3); */

    FloatMatrix transposeLCS;
    transposeLCS.beTranspositionOf(this->localCoordinateSystem);

    this->globalCentroid.beProductOf(transposeLCS, centroid);
    this->globalCentroid.add(x0);
  } //end of Full 3D case
  else{ //Frame case. No vertices. Check for degenerated meshes (1 vertex).
    if(this->numberOfPolygonVertices == 1){
      OOFEM_ERROR("Too small number of polygon vertices for 3D and shell. Too many for frame element. Check mesh.\n");
      return;
    }
    //Read properties from cross-section
    FloatArray lc(1);

    this->giveCrossSection()->give(CS_Area, lc, this);
    this->area = lc.at(1);
    
    this->giveCrossSection()->give(CS_InertiaMomentY, lc, this);
    this->I1 = lc.at(1);
    
    this->giveCrossSection()->give(CS_InertiaMomentZ, lc, this);
    this->I2 = lc.at(1);

    this->giveCrossSection()->give(CS_TorsionConstantX, lc, this);
    this->Ip = lc.at(1);

    //@TODO. Need to add more. We will need to proide Poisson's ratio if we want to give properties for shear.
    //But this has to come from the material.
  }
  return;
}

double Lattice3d :: giveIy() {
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }
    return this->I1;
}

double Lattice3d :: giveIz() {
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    return this->I2;
}

double Lattice3d :: giveIk() {
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    return this->Ip;
}



void
Lattice3d :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    double density = static_cast< LatticeCrossSection * >( this->giveCrossSection() )->give('d', gp);
    double halfMass = density * computeVolumeAround(gp) / 2.;
    answer.resize(12, 12);
    answer.zero();
    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = halfMass;
    answer.at(7, 7) = answer.at(8, 8) = answer.at(9, 9) = halfMass;
}



void
Lattice3d :: saveContext(DataStream &stream, ContextMode mode)
{
    LatticeStructuralElement :: saveContext(stream, mode);

    contextIOResultType iores;

    if ( ( mode & CM_Definition ) ) {
        if ( ( iores = polygonCoords.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( !stream.write(couplingFlag) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( ( iores = couplingNumbers.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }
}


void
Lattice3d :: restoreContext(DataStream &stream, ContextMode mode)
{
    LatticeStructuralElement :: restoreContext(stream, mode);

    contextIOResultType iores;

    if ( mode & CM_Definition ) {
        if ( ( iores = this->polygonCoords.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( !stream.read(this->couplingFlag) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( ( iores = this->couplingNumbers.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }
}

#ifdef __OOFEG

void
Lattice3d :: drawYourself(oofegGraphicContext &gc, TimeStep *tStep)
{
    OGC_PlotModeType mode = gc.giveIntVarPlotMode();

    if ( mode == OGC_rawGeometry ) {
        this->drawRawGeometry(gc, tStep);
        this->drawRawCrossSections(gc, tStep);
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




void Lattice3d :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;

    WCRec p [ 2 ]; /* points */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor(gc.getElementColor() );
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


void
Lattice3d :: drawRawCrossSections(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;

    //Create as many points as we have polygon vertices
    this->numberOfPolygonVertices = this->polygonCoords.giveSize() / 3.;
    WCRec p [ numberOfPolygonVertices ]; /* poin */

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);

    EASValsSetLayer(OOFEG_RAW_CROSSSECTION_LAYER);

    for ( int i = 0; i < numberOfPolygonVertices; i++ ) {
        p [ i ].x = ( FPNum ) polygonCoords(3 * i);
        p [ i ].y = ( FPNum ) polygonCoords(3 * i + 1);
        p [ i ].z = ( FPNum ) polygonCoords(3 * i + 2);
    }


    WCRec pTemp [ 2 ]; /* points */
    for ( int i = 0; i < numberOfPolygonVertices; i++ ) {
        if ( i < numberOfPolygonVertices - 1 ) {
            pTemp [ 0 ] = p [ i ];
            pTemp [ 1 ] = p [ i + 1 ];
        } else {
            pTemp [ 0 ] = p [ i ];
            pTemp [ 1 ] = p [ 0 ];
        }

        go = CreateLine3D(pTemp);
        EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
        EGAttachObject(go, ( EObjectP ) this);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}



void Lattice3d :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    double defScale = gc.getDefScale();

    WCRec p [ 2 ]; /* points */

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor(gc.getDeformedElementColor() );
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
