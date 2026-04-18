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
 *               Copyright () 1993 - 2019   Borek Patzak
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

    //Stress is now converted to sectional forces
        FloatArray s;
        convertStressToResultants3d(s,stress, integrationRulesArray [ 0 ]->getIntegrationPoint(0));


        answer.beProductOf(bt, s);
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
    LatticeMaterialStatus *status = static_cast< LatticeMaterialStatus * >( this->giveCrossSection()->giveMaterial(gp)->giveStatus(gp) );
    plasticStrain = status->givePlasticLatticeStrain();
    return;
}

void
Lattice3d :: giveOldPlasticStrain(FloatArray &plasticStrain)
{
    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    LatticeMaterialStatus *status = static_cast< LatticeMaterialStatus * >( this->giveCrossSection()->giveMaterial(gp)->giveStatus(gp) );
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
    FloatMatrix d, ds, bi, bj, bjt, dbj, dij;

    answer.resize(12, 12);
    answer.zero();
    this->computeBmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), bj);
    this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

    convertTangentToResultantTangent3d(ds, d, integrationRulesArray [ 0 ]->getIntegrationPoint(0));

    dbj.beProductOf(ds, bj);
    bjt.beTranspositionOf(bj);
    answer.beProductOf(bjt, dbj);
    answer.times(1./this->giveLength());

    return;
}

void Lattice3d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    this->numberOfGaussPoints = 1;

    integrationRulesArray.resize(1);
    integrationRulesArray [ 0 ].reset(new GaussIntegrationRule(1, this, 1, 3) );
    integrationRulesArray [ 0 ]->SetUpPointsOnLine(1, _3dLattice);
}



double Lattice3d :: giveArea(GaussPoint *gp) {
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
    LatticeStructuralElement ::initializeFrom(ir);
    
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

    thickness = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, thickness, _IFT_Lattice3d_thickness);

//Introduce here the geometry calculation
//computeGeometryProperties();

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
Lattice3d::giveGpCoordinates(FloatArray &answer)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }
    answer.resize(3);
    answer = this->globalCentroid;
    return;
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


  auto *cs = dynamic_cast< LatticeCrossSection * >( this->giveCrossSection() );
  if (!cs) {
    OOFEM_ERROR("Expected LatticeCrossSection");
  }
  
  if (cs->giveShape() == 1) {
    double r = cs->giveRadius();
    
    this->area = M_PI * r * r;
    this->I1   = M_PI * pow(r, 4) / 4.0;
    this->I2   = this->I1;
    this->Ip   = M_PI * pow(r, 4) / 2.0;
    
    const double k = 0.9;
    this->shearArea1 = k * this->area;
    this->shearArea2 = k * this->area;
    
    FloatArray s(3), t(3), ref(3);
    ref.resize(3);
    ref.zero();
    ref.at(3) = 1.0;
    
    if (fabs(normal.dotProduct(ref)) > 0.99) {
      ref.zero();
      ref.at(2) = 1.0;
    }
    
    s.beVectorProductOf(ref, normal);
    s.normalize();
    
    t.beVectorProductOf(normal, s);
    t.normalize();
    
    this->localCoordinateSystem.resize(3,3);
    for (int i = 1; i <= 3; ++i) {
      this->localCoordinateSystem.at(1,i) = normal.at(i);
      this->localCoordinateSystem.at(2,i) = s.at(i);
      this->localCoordinateSystem.at(3,i) = t.at(i);
    }
    
    centroid.resize(3);
    centroid.zero();

    FloatArray midPointLocal(3);
    midPointLocal.beProductOf(this->localCoordinateSystem, midPoint);
    centroid = midPointLocal;

    this->eccS = 0.0;
    this->eccT = 0.0;
    this->globalCentroid = midPoint;

    return;
  }
  
  //Shape is 0 (polygon based cross-section
  
  if(this->numberOfPolygonVertices < 3){
    OOFEM_ERROR("Too small number of polygon vertices. Check meshing approach.\n");
  }
  
    //Construct two perpendicular axis so that n is normal to the plane which they create
    //Check, if one of the components of the normal-direction is zero
  FloatArray s(3), t(3);
  FloatArray ref(3);
  ref = {0,0,1};
  
  if (fabs(normal.dotProduct(ref)) > 0.99) {
    ref = {0,1,0};
  }
  
  s.beVectorProductOf(ref, normal);
  s.normalize();
  
  t.beVectorProductOf(normal, s);
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
    FloatArray lpc(3 * numberOfPolygonVertices);
    for ( int k = 1; k <= numberOfPolygonVertices; k++ ) {
        for ( int n = 1; n <= 3; n++ ) {
	  help.at(n) = polygonCoords.at(3 * (k - 1) + n);
        }

	test.beProductOf(lcs, help);

	for ( int n = 1; n <= 3; n++ ) {
	  lpc.at(3 * (k - 1) + n) = test.at(n);
        }
    }

    this->area = 0.;
    for ( int k = 1; k <= numberOfPolygonVertices; k++ ) {
        if ( k < numberOfPolygonVertices ) {
	  this->area += lpc.at(3 * (k - 1) + 2) * lpc.at(3 * ( k ) + 3) - lpc.at(3 * ( k ) + 2) * lpc.at(3 * (k - 1) + 3);
        } else {   //Back to zero for n+1
	  this->area += lpc.at(3 * (k-1) + 2) * lpc.at(3) - lpc.at(2) * lpc.at(3 * (k-1) + 3);
        }
    }

    this->area *= 0.5;

    FloatArray tempCoords(3 * numberOfPolygonVertices);
    if ( this->area < 0 ) { //Set area to a positive value and rearrange the coordinate entries
        this->area *= -1.;

        for ( int k = 0; k < numberOfPolygonVertices; k++ ) {
            for ( int m = 0; m < 3; m++ ) {
                tempCoords.at(3 * k + m + 1) = polygonCoords.at(3 * ( numberOfPolygonVertices - k - 1 ) + m + 1);
            }
        }

        polygonCoords = tempCoords;

        // Calculate again local co-ordinate system for different order
        for ( int k = 1; k <= numberOfPolygonVertices; k++ ) {
            for ( int n = 1; n <= 3; n++ ) {
	      help.at(n) = polygonCoords.at(3 * (k-1) + n );
            }

            test.beProductOf(lcs, help);
            for ( int n = 1; n <= 3; n++ ) {
	      lpc.at(3 * (k-1) + n) = test.at(n);
            }
        }
    }

    if ( this->area < pow(minLength, 2.) ) {
        this->area = pow(minLength, 2.);
    }

    //Calculate centroids
    centroid.resize(3);
    centroid.zero();
    for ( int k = 1; k <= numberOfPolygonVertices; k++ ) {
        if ( k < numberOfPolygonVertices ) {
	  centroid.at(2) += ( lpc.at(3 * (k-1) + 2) + lpc.at(3 * ( k ) + 2) ) * ( lpc.at(3 * (k-1) + 2) * lpc.at(3 * ( k ) + 3) - lpc.at(3 * ( k ) + 2) * lpc.at(3 * (k-1) + 3) );	 
	  centroid.at(3) += ( lpc.at(3 * (k-1) + 3) + lpc.at(3 * ( k ) + 3) ) * ( lpc.at(3 * (k-1) +2) * lpc.at(3 * ( k ) + 3) - lpc.at(3 * ( k ) + 2) * lpc.at(3 * (k-1) + 3) );
        } else { //Back to zero for n+1
	  centroid.at(2) += ( lpc.at(3 * (k-1) + 2) + lpc.at(2) ) * ( lpc.at(3 * (k-1) + 2) * lpc.at(3) - lpc.at(2) * lpc.at(3 * (k-1) + 3) );
	  centroid.at(3) += ( lpc.at(3 * (k-1) + 3) + lpc.at(3) ) * ( lpc.at(3 * (k-1) + 2) * lpc.at(3) - lpc.at(2) * lpc.at(3 * (k-1) + 3) );
        }
    }

    centroid.times(1. / ( 6. * this->area ) );

    centroid.at(1) = lpc.at(1); //The first component of all lpcs should be the same

    //Shift coordinates to centroid
    for ( int k = 1; k <= numberOfPolygonVertices; k++ ) {
        for ( int l = 1; l <= 3; l++ ) {
	  lpc.at(3 * (k-1) + l) -= centroid.at(l);
        }
    }

    //Compute second moments of area.
    //This is for the temporary coordinate system
    double Ixx = 0.;
    double Iyy = 0.;
    double Ixy = 0.;
    double a;

    for ( int k = 1; k <= numberOfPolygonVertices; k++ ) {
        if ( k < numberOfPolygonVertices ) {
	  a = lpc.at(3 * (k-1) + 2) * lpc.at(3 * ( k ) + 3) - lpc.at(3 * ( k ) + 2) * lpc.at(3 * (k-1) + 3);

	  Ixx += ( ( pow(lpc.at(3 * (k-1) + 3), 2.) + lpc.at(3 * (k-1) + 3) * lpc.at(3 * ( k ) + 3) + pow(lpc.at(3 * ( k ) + 3), 2.) ) * a ) / 12.;

	  Iyy += ( ( pow(lpc.at(3 * (k-1) + 2), 2.) + lpc.at(3 * (k-1) + 2) * lpc.at(3 * ( k ) + 2) + pow(lpc.at(3 * ( k ) + 2), 2.) ) * a ) / 12.;

	  Ixy += ( ( lpc.at(3 * (k-1) + 2) * lpc.at(3 * ( k ) + 3) + 2. * lpc.at(3 * (k-1) + 2) * lpc.at(3 * (k-1) + 3) +
		     2 * lpc.at(3 * ( k ) + 2) * lpc.at(3 * ( k ) + 3) + lpc.at(3 * ( k ) + 2) * lpc.at(3 * (k-1) + 3) ) * a ) / 24.;
        } else {   //Back to zero for n+1
	  a = lpc.at(3 * (k-1) + 2) * lpc.at(3) - lpc.at(2) * lpc.at(3 * (k-1) + 3);

	  Ixx += ( ( pow(lpc.at(3 * (k-1) + 3), 2.) + lpc.at(3 * (k-1) + 3) * lpc.at(3) + pow(lpc.at(3), 2.) ) * a ) / 12.;

	  Iyy += ( ( pow(lpc.at(3 * (k-1) + 2), 2.) + lpc.at(3 * (k-1) + 2) * lpc.at(2) + pow(lpc.at(2), 2.) ) * a ) / 12.;

	  Ixy += ( ( lpc.at(3 * (k-1) + 2) * lpc.at(3) + 2. * lpc.at(3 * (k-1) + 2) * lpc.at(3 * (k-1) + 3) +
		     2 * lpc.at(2) * lpc.at(3) + lpc.at(2) * lpc.at(3 * (k-1) + 3) ) * a ) / 24.;
        }
    }

    //Compute main axis of the cross-section
    double angleChange = 0.;
    double tol = 1e-12 * (fabs(Ixx) + fabs(Iyy));

    if ( fabs(Ixx - Iyy) < tol && fabs(Ixy) < tol ) {
      angleChange = 0.0;  // orientation arbitrary
    } else {
      angleChange = 0.5 * atan2(-2.0 * Ixy, (Ixx - Iyy));
    }

    //Moment of inertias saved in the element
    this->I1 = ( Ixx + Iyy ) / 2. + sqrt(pow( ( Ixx - Iyy ) / 2., 2. ) + pow(Ixy, 2.) );
    this->I2 = ( Ixx + Iyy ) / 2. - sqrt(pow( ( Ixx - Iyy ) / 2., 2. ) + pow(Ixy, 2.) );

    this->Ip = I1 + I2;

if (cs->giveShape() == 2) {
  // shell-type rectangle given by polygon vertices
// Use effective rotational/torsional parameter J = b*t^3/6,
// where t is the physical shell thickness and b is the in-plane tributary width.
// This preserves shell-like thickness scaling and avoids the non-physical
// mesh dependence.

  
    if (this->numberOfPolygonVertices != 4) {
        OOFEM_ERROR("Rectangle cross-section must have 4 vertices.");
    }

    auto edgeLength = [&lpc](int k1, int k2) {
        double dy = lpc.at(3 * (k2 - 1) + 2) - lpc.at(3 * (k1 - 1) + 2);
        double dz = lpc.at(3 * (k2 - 1) + 3) - lpc.at(3 * (k1 - 1) + 3);
        return std::sqrt(dy * dy + dz * dz);
    };

    const double e1 = edgeLength(1, 2);
    const double e2 = edgeLength(2, 3);
    const double e3 = edgeLength(3, 4);
    const double e4 = edgeLength(4, 1);

    const double rectTol = 1e-6 * (e1 + e2 + e3 + e4);
    if (std::fabs(e1 - e3) > rectTol || std::fabs(e2 - e4) > rectTol) {
        OOFEM_ERROR("Shape 2 cross-section is not a rectangle.");
    }

    const double d1 = 0.5 * (e1 + e3);
    const double d2 = 0.5 * (e2 + e4);
    const double t  = this->giveThickness();

    const double thickTol = 1e-6 * (d1 + d2);

    double b = 0.0;
    double h = 0.0;

    if (std::fabs(d1 - t) < thickTol) {
        h = d1;
        b = d2;
    } else if (std::fabs(d2 - t) < thickTol) {
        h = d2;
        b = d1;
    } else {
        OOFEM_ERROR("Rectangle cross-section does not match element thickness.");
    }

    this->J = b * h * h * h / 6.0;

      /* if (b <= 0.0 || h <= 0.0) { */
      /* 	OOFEM_ERROR("Invalid rectangle dimensions in cross-section."); */
      /* } */
      
      /* if (h > b) { */
      /*   std::swap(b, h); */
      /* } */
      
      /* const double r = h / b; */
      /* this->J = (b * h * h * h / 3.0) * (1.0 - 0.63 * r + 0.052 * pow(r, 5)); */
      

    } else {
      this->J = this->Ip; // temporary fallback for generic polygons
    }

    

    this->shearArea1 = area;
    this->shearArea2 = area;
    
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
    for ( int k = 1; k <= numberOfPolygonVertices; k++ ) {
        for ( int n = 1; n <= 3; n++ ) {
	  help.at(n) = polygonCoords.at(3 * (k-1) + n);
        }

        test.beProductOf(this->localCoordinateSystem, help);
        for ( int n = 1; n <= 3; n++ ) {
	  lpc.at(3 * (k-1) + n) = test.at(n);
        }
    }

    //Calculate centroid again in local coordinate system
    centroid.zero();
    for ( int k = 1; k <= numberOfPolygonVertices; k++ ) {
        if ( k < numberOfPolygonVertices ) {
	  centroid.at(2) += ( lpc.at(3 * (k-1) + 2) + lpc.at(3 * ( k ) + 2) ) * ( lpc.at(3 * (k-1) + 2) * lpc.at(3 * ( k ) + 3) - lpc.at(3 * ( k ) + 2) * lpc.at(3 * (k-1) + 3) );
	  centroid.at(3) += ( lpc.at(3 * (k-1) + 3) + lpc.at(3 * ( k ) + 3) ) * ( lpc.at(3 * (k-1) + 2) * lpc.at(3 * ( k ) + 3) - lpc.at(3 * ( k ) + 2) * lpc.at(3 * (k-1) + 3) );
        } else {   //Back to zero for n+1
	  centroid.at(2) += ( lpc.at(3 * (k-1) + 2) + lpc.at(2) ) * ( lpc.at(3 * (k-1) + 2) * lpc.at(3) - lpc.at(2) * lpc.at(3 * (k-1) + 3) );
          centroid.at(3) += ( lpc.at(3 * (k-1) + 3) + lpc.at(3) ) * ( lpc.at(3 * (k-1) + 2) * lpc.at(3) - lpc.at(2) * lpc.at(3 * (k-1) + 3) );
        }
    }

    centroid.times(1. / ( 6. * this->area ) );

    centroid.at(1) = lpc.at(1); //The first component of all lpcs should be the same

    FloatArray midPointLocal(3);
    midPointLocal.beProductOf(this->localCoordinateSystem, midPoint);

    //eccentricities stored in the element
    this->eccS = centroid.at(2) - midPointLocal.at(2);
    this->eccT = centroid.at(3) - midPointLocal.at(3);

    FloatMatrix transposeLCS;
    transposeLCS.beTranspositionOf(this->localCoordinateSystem);

    this->globalCentroid.beProductOf(transposeLCS, centroid);
     
  return;
}


bool
Lattice3d::giveRectangularSectionDimensions(double &by, double &bz, GaussPoint *gp) const
{
    if (this->numberOfPolygonVertices != 4) {
        return false;
    }

    FloatArray global(3), local(3);
    std::vector<std::pair<double,double>> yz;
    yz.reserve(4);

    for (int k = 1; k <= 4; ++k) {
        for (int i = 1; i <= 3; ++i) {
            global.at(i) = this->polygonCoords.at(3 * (k - 1) + i);
        }

        local.beProductOf(this->localCoordinateSystem, global);
        yz.emplace_back(local.at(2), local.at(3));
    }

    auto edgeLength = [&yz](int k1, int k2) {
        const double dy = yz[k2].first  - yz[k1].first;
        const double dz = yz[k2].second - yz[k1].second;
        return std::sqrt(dy * dy + dz * dz);
    };

    const double e1 = edgeLength(0, 1);
    const double e2 = edgeLength(1, 2);
    const double e3 = edgeLength(2, 3);
    const double e4 = edgeLength(3, 0);

    const double scale = std::max({e1, e2, e3, e4, 1.0});
    const double tol = 1e-8 * scale;

    if (std::fabs(e1 - e3) > tol || std::fabs(e2 - e4) > tol) {
        return false;
    }

    by = 0.5 * (e1 + e3);
    bz = 0.5 * (e2 + e4);

    return (by > tol && bz > tol);
}


double Lattice3d :: giveI1(GaussPoint *gp) {
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }
    return this->I1;
}

double Lattice3d :: giveI2(GaussPoint *gp) {
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    return this->I2;
}


double Lattice3d :: giveJ(GaussPoint *gp) {
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }
    return this->J;
}



double Lattice3d :: giveShearArea1(GaussPoint *gp) {
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }
    //Temporary assumption. Ideally, shear area should be less than area.
    return this->shearArea1;
}

double Lattice3d :: giveShearArea2(GaussPoint *gp) {
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }
    //Temporary assumption. Ideally, shear area should be less than area.
    return this->shearArea2;
}

double Lattice3d :: giveTributaryWidth(GaussPoint *gp)
{
    if (geometryFlag == 0) {
        computeGeometryProperties();
    }

    if (this->thickness > 0.0) {
        return this->area / this->thickness;
    }

    return 1.0;
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
