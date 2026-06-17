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

#include "vtkxmllatticeexportmodule.h"
#include "element.h"
#include "gausspoint.h"
#include "timestep.h"
#include "engngm.h"
#include "node.h"
#include "materialinterface.h"
#include "mathfem.h"
#include "cltypes.h"
#include "material.h"
#include "classfactory.h"
#include "crosssection.h"
#include "dof.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "dofiditem.h"
#include <cmath>
#include <map>

#ifdef __SM_MODULE
 #include "../sm/Elements/LatticeElements/latticestructuralelement.h"
 #include "../sm/Elements/LatticeElements/lattice3d.h"
 #include "../sm/Elements/LatticeElements/latticelink3d.h"
 #include "../sm/Elements/LatticeElements/latticelink3dboundary.h"
 #include "../sm/CrossSections/latticecrosssection.h"
#endif

namespace {
using namespace oofem;

// R = I + sin(theta)*K + (1-cos(theta))*K^2 where K = skew(spin/theta).
// Degrades to I + skew(spin) for small theta.
void rodriguesRotation(FloatMatrix &R, const FloatArray &spin)
{
    R.resize(3, 3); R.zero();
    R.at(1,1) = R.at(2,2) = R.at(3,3) = 1.0;

    const double theta = spin.computeNorm();
    if ( theta < 1e-12 ) {
        R.at(1,2) = -spin.at(3); R.at(2,1) =  spin.at(3);
        R.at(1,3) =  spin.at(2); R.at(3,1) = -spin.at(2);
        R.at(2,3) = -spin.at(1); R.at(3,2) =  spin.at(1);
        return;
    }
    const double s = std::sin(theta) / theta;
    const double c = ( 1.0 - std::cos(theta) ) / ( theta * theta );
    const double k1 = spin.at(1), k2 = spin.at(2), k3 = spin.at(3);
    R.at(1,1) += c * ( -k2*k2 - k3*k3 );
    R.at(1,2) += -s * k3 + c * ( k1*k2 );
    R.at(1,3) +=  s * k2 + c * ( k1*k3 );
    R.at(2,1) +=  s * k3 + c * ( k1*k2 );
    R.at(2,2) += c * ( -k1*k1 - k3*k3 );
    R.at(2,3) += -s * k1 + c * ( k2*k3 );
    R.at(3,1) += -s * k2 + c * ( k1*k3 );
    R.at(3,2) +=  s * k1 + c * ( k2*k3 );
    R.at(3,3) += c * ( -k1*k1 - k2*k2 );
}

// v_def = node_ref + node_disp + R(node_rot) * (v_ref - node_ref).
void rigidBodyTransform(FloatArray &v_def, const FloatArray &v_ref,
                        const FloatArray &node_ref, const FloatArray &node_disp,
                        const FloatArray &node_rot)
{
    FloatMatrix R;
    rodriguesRotation(R, node_rot);
    FloatArray rel(3);
    for ( int i = 1; i <= 3; ++i ) rel.at(i) = v_ref.at(i) - node_ref.at(i);
    FloatArray rotated(3);
    rotated.beProductOf(R, rel);
    v_def.resize(3);
    for ( int i = 1; i <= 3; ++i ) v_def.at(i) = node_ref.at(i) + node_disp.at(i) + rotated.at(i);
}

// Subdivide a 4-vertex shell polygon into N strips along the shell-normal direction.
// stripVerts is (4*N) x 3; strip k occupies rows 4k+1 .. 4k+4 (CCW order).
// Falls back to copying the original polygon if N==1 or the polygon isn't 4-vertex.
void subdivideShellPolygon(const FloatArray &polyRef, int polyNV,
                            const FloatArray &normal, int N,
                            FloatMatrix &stripVerts)
{
    if ( N <= 1 || polyNV != 4 || normal.giveSize() != 3 ) {
        stripVerts.resize(polyNV, 3);
        for ( int k = 0; k < polyNV; ++k ) {
            for ( int j = 1; j <= 3; ++j ) stripVerts.at(k+1, j) = polyRef.at(3*k + j);
        }
        return;
    }
    FloatArray centroid(3); centroid.zero();
    for ( int k = 0; k < 4; ++k ) for ( int j = 1; j <= 3; ++j ) centroid.at(j) += polyRef.at(3*k + j) / 4.0;
    FloatArray proj(4);
    for ( int k = 0; k < 4; ++k ) {
        double p = 0;
        for ( int j = 1; j <= 3; ++j ) p += ( polyRef.at(3*k + j) - centroid.at(j) ) * normal.at(j);
        proj.at(k+1) = p;
    }
    double meanP = 0; for ( int k = 1; k <= 4; ++k ) meanP += proj.at(k) / 4.0;
    int top[2] = { -1, -1 }, bot[2] = { -1, -1 }, nt = 0, nb = 0;
    for ( int k = 0; k < 4; ++k ) {
        if ( proj.at(k+1) > meanP ) { if ( nt < 2 ) top[nt++] = k; }
        else { if ( nb < 2 ) bot[nb++] = k; }
    }
    if ( nt != 2 || nb != 2 ) {  // degenerate; fall back
        stripVerts.resize(4, 3);
        for ( int k = 0; k < 4; ++k ) for ( int j = 1; j <= 3; ++j ) stripVerts.at(k+1, j) = polyRef.at(3*k + j);
        return;
    }
    auto vdist2 = [&](int a, int b) {
        double d = 0;
        for ( int j = 1; j <= 3; ++j ) {
            double diff = polyRef.at(3*a + j) - polyRef.at(3*b + j);
            d += diff * diff;
        }
        return d;
    };
    // Pair top[0] with whichever bot is closest in 3D (i.e. same in-plane position).
    int botPair0 = ( vdist2(top[0], bot[0]) < vdist2(top[0], bot[1]) ) ? bot[0] : bot[1];
    int botPair1 = ( botPair0 == bot[0] ) ? bot[1] : bot[0];
    int topPair0 = top[0], topPair1 = top[1];

    stripVerts.resize(4 * N, 3);
    for ( int k = 0; k < N; ++k ) {
        double tb = static_cast< double >( k ) / N;
        double tt = static_cast< double >( k + 1 ) / N;
        for ( int j = 1; j <= 3; ++j ) {
            double b0 = polyRef.at(3*botPair0 + j), t0 = polyRef.at(3*topPair0 + j);
            double b1 = polyRef.at(3*botPair1 + j), t1 = polyRef.at(3*topPair1 + j);
            stripVerts.at(4*k + 1, j) = b0 + tb * ( t0 - b0 );
            stripVerts.at(4*k + 2, j) = b1 + tb * ( t1 - b1 );
            stripVerts.at(4*k + 3, j) = b1 + tt * ( t1 - b1 );
            stripVerts.at(4*k + 4, j) = b0 + tt * ( t0 - b0 );
        }
    }
}

#ifdef __SM_MODULE
// Returns nLayers for hybrid shell (>1), else 1. Also fills shellNormal if applicable.
int getElementStripCount(Element *el, FloatArray &shellNormal)
{
    shellNormal.resize(0);
    auto *lat3d = dynamic_cast< Lattice3d * >( el );
    if ( !lat3d || !lat3d->isHybridShell() ) return 1;
    shellNormal = lat3d->giveShellNormal();
    auto *lcs = dynamic_cast< LatticeCrossSection * >( el->giveCrossSection() );
    if ( !lcs ) return 1;
    return lcs->giveNLayers();
}
#else
int getElementStripCount(Element *, FloatArray &shellNormal)
{
    shellNormal.resize(0);
    return 1;
}
#endif

// Build a frame (e_x = given axis, e_y, e_z). Pick world Z as up reference; world Y if vertical.
void buildFrameFromAxis(const FloatArray &axis, FloatArray &ex, FloatArray &ey, FloatArray &ez)
{
    ex.resize(3);
    double len = axis.computeNorm();
    if ( len < 1e-12 ) { ex.at(1) = 1; ex.at(2) = 0; ex.at(3) = 0; }
    else for ( int i = 1; i <= 3; ++i ) ex.at(i) = axis.at(i) / len;
    FloatArray up(3); up.at(1) = 0; up.at(2) = 0; up.at(3) = 1;
    if ( std::fabs(ex.dotProduct(up)) > 0.99 ) { up.at(1) = 0; up.at(2) = 1; up.at(3) = 0; }
    double exDotUp = ex.dotProduct(up);
    ey.resize(3);
    for ( int i = 1; i <= 3; ++i ) ey.at(i) = up.at(i) - exDotUp * ex.at(i);
    ey.times(1.0 / ey.computeNorm());
    ez.resize(3);
    ez.beVectorProductOf(ex, ey);
}

// Build a local frame (e_x, e_y, e_z) given the element axis. Pick world Z as up reference
// unless axis is nearly vertical, in which case use world Y.
void buildElementFrame(const FloatArray &nodeA, const FloatArray &nodeB,
                        FloatArray &ex, FloatArray &ey, FloatArray &ez)
{
    ex.resize(3);
    for ( int i = 1; i <= 3; ++i ) ex.at(i) = nodeB.at(i) - nodeA.at(i);
    double len = ex.computeNorm();
    if ( len < 1e-12 ) {
        ex.at(1) = 1; ex.at(2) = 0; ex.at(3) = 0;
        ey.resize(3); ey.at(1) = 0; ey.at(2) = 1; ey.at(3) = 0;
        ez.resize(3); ez.at(1) = 0; ez.at(2) = 0; ez.at(3) = 1;
        return;
    }
    ex.times(1.0 / len);
    FloatArray up(3); up.at(1) = 0; up.at(2) = 0; up.at(3) = 1;
    if ( std::fabs(ex.dotProduct(up)) > 0.99 ) { up.at(1) = 0; up.at(2) = 1; up.at(3) = 0; }
    double exDotUp = ex.dotProduct(up);
    ey.resize(3);
    for ( int i = 1; i <= 3; ++i ) ey.at(i) = up.at(i) - exDotUp * ex.at(i);
    ey.times(1.0 / ey.computeNorm());
    ez.resize(3);
    ez.beVectorProductOf(ex, ey);
}

// Synthesise a rectangular polygon (4 vertices CCW) of dimensions b × h, centred at midpoint,
// in the (e_y, e_z) plane perpendicular to the element axis.
void synthesiseRectangle(const FloatArray &mid, const FloatArray &ey, const FloatArray &ez,
                          double b, double h, FloatArray &polyCoords)
{
    polyCoords.resize(12);
    double signsY[4] = { +1, -1, -1, +1 };
    double signsZ[4] = { +1, +1, -1, -1 };
    for ( int k = 0; k < 4; ++k ) {
        for ( int j = 1; j <= 3; ++j ) {
            polyCoords.at(3*k + j) = mid.at(j)
                                    + 0.5 * b * signsY[k] * ey.at(j)
                                    + 0.5 * h * signsZ[k] * ez.at(j);
        }
    }
}

// Synthesise an N-gon (N vertices CCW) of given radius, centred at midpoint, in (e_y, e_z) plane.
void synthesisePolygonNgon(const FloatArray &mid, const FloatArray &ey, const FloatArray &ez,
                            double radius, int N, FloatArray &polyCoords)
{
    polyCoords.resize(3 * N);
    for ( int k = 0; k < N; ++k ) {
        double angle = 2.0 * M_PI * k / N;
        double cy = radius * std::cos(angle);
        double cz = radius * std::sin(angle);
        for ( int j = 1; j <= 3; ++j ) {
            polyCoords.at(3*k + j) = mid.at(j) + cy * ey.at(j) + cz * ez.at(j);
        }
    }
}

// Try to synthesise a cross-section polygon for an element without polygonCoords.
// Returns # vertices placed in polyCoords (0 if synthesis failed → treat as line for 2-node, else skip).
#ifdef __SM_MODULE
int synthesiseElementPolygon(Element *el, FloatArray &polyCoords)
{
    if ( el->giveNumberOfDofManagers() < 2 ) return 0;
    // Link/bond elements: render as line, not as a synthesised cross-section.
    if ( dynamic_cast< LatticeLink3d * >( el ) || dynamic_cast< LatticeLink3dBoundary * >( el ) ) return 0;
    const FloatArray &cA = el->giveNode(1)->giveCoordinates();
    const FloatArray &cB = el->giveNode(2)->giveCoordinates();
    FloatArray mid(3);
    for ( int i = 1; i <= 3; ++i ) mid.at(i) = 0.5 * ( cA.at(i) + cB.at(i) );
    FloatArray ex, ey, ez;
    buildElementFrame(cA, cB, ex, ey, ez);

    auto *lcs = dynamic_cast< LatticeCrossSection * >( el->giveCrossSection() );
    if ( lcs && lcs->giveShape() == 1 && lcs->giveRadius() > 0 ) {
        // Circular: 8-gon.
        synthesisePolygonNgon(mid, ey, ez, lcs->giveRadius(), 8, polyCoords);
        return 8;
    }

    auto *lse = dynamic_cast< LatticeStructuralElement * >( el );
    if ( !lse ) return 0;
    IntegrationRule *iRule = el->giveDefaultIntegrationRulePtr();
    if ( !iRule || iRule->giveNumberOfIntegrationPoints() == 0 ) return 0;
    GaussPoint *gp = iRule->getIntegrationPoint(0);
    double A = lse->giveArea(gp);
    if ( A <= 0 ) return 0;
    double I1 = lse->giveI1(gp);
    double I2 = lse->giveI2(gp);
    double b, h;
    if ( I1 > 0 ) {
        b = std::sqrt(12.0 * I1 / A);
        h = A / b;
    } else if ( I2 > 0 ) {
        h = std::sqrt(12.0 * I2 / A);
        b = A / h;
    } else {
        b = h = std::sqrt(A);
    }
    synthesiseRectangle(mid, ey, ez, b, h, polyCoords);
    return 4;
}
#else
int synthesiseElementPolygon(Element *, FloatArray &polyCoords) {
    polyCoords.resize(0);
    return 0;
}
#endif

// Read total nodal displacement (D_u, D_v, D_w) and rotation (R_u, R_v, R_w).
void giveNodeKinematics(Node *node, TimeStep *tStep, FloatArray &disp, FloatArray &rot)
{
    disp.resize(3); disp.zero();
    rot.resize(3); rot.zero();
    if ( !node || !tStep ) return;
    IntArray dispIds = { D_u, D_v, D_w };
    IntArray rotIds  = { R_u, R_v, R_w };
    FloatArray vec;
    node->giveUnknownVector(vec, dispIds, VM_Total, tStep, true);
    for ( int i = 1; i <= std::min(3, vec.giveSize()); ++i ) disp.at(i) = vec.at(i);
    node->giveUnknownVector(vec, rotIds, VM_Total, tStep, true);
    for ( int i = 1; i <= std::min(3, vec.giveSize()); ++i ) rot.at(i) = vec.at(i);
}

} // anonymous namespace

#ifdef __SM_MODULE
 #include "../sm/Elements/LatticeElements/latticestructuralelement.h"
 #include "../sm/Elements/LatticeElements/lattice2dboundary.h"
 #include "../sm/Elements/LatticeElements/lattice3dboundary.h"
 #include "../sm/Elements/LatticeElements/latticelink3dboundary.h"
 #include "../sm/Elements/LatticeElements/latticelink3d.h"
 #include "../sm/Elements/LatticeElements/latticeframe3d.h"
 #include "../sm/Elements/LatticeElements/latticeframe3dnl.h"
 #include "../sm/Elements/LatticeElements/lattice3d.h"
 #include "../sm/CrossSections/latticecrosssection.h"
#endif

#ifdef __TM_MODULE
 #include "../tm/Elements/LatticeElements/latticetransportelement.h"
#endif


namespace oofem {
REGISTER_ExportModule(VTKXMLLatticeExportModule)

VTKXMLLatticeExportModule::VTKXMLLatticeExportModule(int n, EngngModel *e) : VTKXMLExportModule(n, e)
{
    //Find out the element type and elemNodes
}


VTKXMLLatticeExportModule::~VTKXMLLatticeExportModule()
{}


void
VTKXMLLatticeExportModule::initializeFrom(InputRecord &ir)
{
    VTKXMLExportModule::initializeFrom(ir);
    this->crossSectionExportFlag = false;
    IR_GIVE_OPTIONAL_FIELD(ir, this->crossSectionExportFlag, _IFT_VTKXMLLatticeExportModule_cross);
    this->crossOutputStep = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, this->crossOutputStep, _IFT_VTKXMLLatticeExportModule_crossstep);
    if ( this->crossOutputStep < 1 ) this->crossOutputStep = 1;
    this->crossOutputCounter = 0;
    this->pvdBufferLine.clear();
    this->pvdBufferCross.clear();
}

std::string
VTKXMLLatticeExportModule::giveOutputFileNameCross(TimeStep *tStep)
{
    return this->giveOutputBaseFileName(tStep) + ".cross.vtu";
}

std::ofstream
VTKXMLLatticeExportModule::giveOutputStreamCross(TimeStep *tStep)
{
    std::string fileName = giveOutputFileNameCross(tStep);
    std::ofstream streamC;

    if ( pythonExport ) {
        streamC = std::ofstream(NULL_DEVICE);//do not write anything
    } else {
        streamC = std::ofstream(fileName);
    }

    if ( !streamC.good() ) {
        OOFEM_ERROR("failed to open file %s", fileName.c_str() );
    }

    streamC.fill('0');//zero padding
    return streamC;
}

void
VTKXMLLatticeExportModule::giveSwitches(IntArray &answer, int location) {
    answer.resize(3);
    answer.zero();
    int counter = 1;
    for ( int x = -1; x <  2; x++ ) {
        for ( int y = -1; y <  2; y++ ) {
            for ( int z = -1; z <  2; z++ ) {
                if ( !( z == 0 && y == 0 && x == 0 ) ) {
                    if ( counter == location ) {
                        answer(0) = x;
                        answer(1) = y;
                        answer(2) = z;
                    }
                    counter++;
                }
            }
        }
    }

    return;
}

void
VTKXMLLatticeExportModule::setupVTKPiece(ExportRegion &vtkPiece, TimeStep *tStep, Set &region)
{
    // Stores all neccessary data (of a region) in a VTKPiece so it can be exported later.

    Domain *d  = emodel->giveDomain(1);
    Element *elem;

    int nnodes = d->giveNumberOfDofManagers();

    // Assemble local->global and global->local region map and get number of
    // single cells to process, the composite cells exported individually.
    this->initRegionNodeNumbering(vtkPiece, d, tStep, region);
    const IntArray& mapG2L = vtkPiece.getMapG2L();
    const IntArray& mapL2G = vtkPiece.getMapL2G();
    const int numNodes = vtkPiece.giveNumberOfNodes();
    const int numRegionEl = vtkPiece.giveNumberOfCells();



    if ( numNodes > 0 && numRegionEl > 0 ) {
        // Export nodes as vtk vertices
        vtkPiece.setNumberOfNodes(numNodes);
        for ( int inode = 1; inode <= numNodes; inode++ ) {
            if ( mapL2G.at(inode) <= nnodes && mapL2G.at(inode) != 0 ) { //DofManagers in domain (input file)
                const auto &coords = d->giveNode(mapL2G.at(inode) )->giveCoordinates();
                vtkPiece.setNodeCoords(inode, coords);
            } else if ( mapL2G.at(inode) > nnodes && mapL2G.at(inode) <= numNodes && mapL2G.at(inode) != 0 ) { //extra image nodes
                FloatArray helpArray(3);
                helpArray.at(1) = uniqueNodeTable.at(mapL2G.at(inode), 1);
                helpArray.at(2) = uniqueNodeTable.at(mapL2G.at(inode), 2);
                helpArray.at(3) = uniqueNodeTable.at(mapL2G.at(inode), 3);
                vtkPiece.setNodeCoords(inode, helpArray);
            } else {
                FloatArray helpArray(3);
                helpArray.zero();
                vtkPiece.setNodeCoords(inode, helpArray);
            }
        }


        //-------------------------------------------
        // Export all the cell data for the piece
        //-------------------------------------------
        IntArray cellNodes;
        vtkPiece.setNumberOfCells(numRegionEl);

        int offset = 0;
        int cellNum = 0;
        IntArray elems = region.giveElementList();
        int helpCounter = 0;
        for ( int ei = 1; ei <= elems.giveSize(); ei++ ) {
            int elNum = elems.at(ei);
            elem = d->giveElement(elNum);

            // Skip elements that:
            // are inactivated or of composite type ( these are exported individually later)
            if ( this->isElementComposite(elem) || !elem->isActivated(tStep) ) {
                continue;
            }

            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

            cellNum++;

            // Set the connectivity
#ifdef __SM_MODULE
            if ( dynamic_cast< Lattice3dBoundary * >( elem ) || dynamic_cast< LatticeLink3dBoundary * >( elem ) || dynamic_cast< Lattice2dBoundary * >( elem ) ) {
                cellNodes.resize(2);
                IntArray loc = elem->giveLocation();
                for ( int ielnode = 1; ielnode <= 2; ielnode++ ) {
                    if ( loc.at(ielnode) != 0 ) {
                        helpCounter++;
                        cellNodes.at(ielnode) = regionToUniqueMap.at(nnodes + helpCounter);
                    } else {
                        cellNodes.at(ielnode) = elem->giveNode(ielnode)->giveNumber();
                    }
                }
            } else {  //Standard case
#endif
            this->giveElementCell(cellNodes, elem);
#ifdef __SM_MODULE
        }
#endif
            // Map from global to local node numbers for the current piece
            int numElNodes = cellNodes.giveSize();


            IntArray connectivity(numElNodes);
            for ( int i = 1; i <= numElNodes; i++ ) {
                connectivity.at(i) = mapG2L.at(cellNodes.at(i) );
            }

            vtkPiece.setConnectivity(cellNum, connectivity);

            vtkPiece.setCellType(cellNum, this->giveCellType(elem) );  // VTK cell type

            offset += numElNodes;
            vtkPiece.setOffset(cellNum, offset);
        }
    } // end of default piece for simple geometry elements
}


void
VTKXMLLatticeExportModule::setupVTKPieceCross(ExportRegion &vtkPieceCross, TimeStep *tStep, Set& region)
{
    // Per lattice element with a polygon: emit 3 cells -- A-side, B-side, midline.
    // A and B copies move with their owner node's rigid-body kinematics (total disp + Rodrigues rot).
    // Midline = average of the two deformed copies; carries the edge-level cell data.
    // Elements without a polygon (frame, transport-only) fall back to a single VTK_VERTEX.

    Domain *domain = emodel->giveDomain( 1 );
    IntArray elements    = region.giveElementList();
    int numberOfElements = elements.giveSize();

    auto getPolygonNodeCount = [&](int ielem) -> int {
        Element *el = domain->giveElement(ielem);
        if ( auto *le = dynamic_cast<LatticeStructuralElement *>(el) ) {
            return le->giveNumberOfCrossSectionNodes();
        } else if ( auto *le = dynamic_cast<LatticeTransportElement *>(el) ) {
            return le->giveNumberOfCrossSectionNodes();
        }
        return 0;
    };

    auto getPolygonCoords = [&](int ielem, FloatArray &coords) {
        Element *el = domain->giveElement(ielem);
        if ( auto *le = dynamic_cast<LatticeStructuralElement *>(el) ) {
            le->giveCrossSectionCoordinates(coords);
        } else if ( auto *le = dynamic_cast<LatticeTransportElement *>(el) ) {
            le->giveCrossSectionCoordinates(coords);
        }
    };

    auto getGpCoords = [&](int ielem, FloatArray &coords) {
        Element *el = domain->giveElement(ielem);
        if ( auto *le = dynamic_cast<LatticeStructuralElement *>(el) ) {
            le->giveGpCoordinates(coords);
        } else if ( auto *le = dynamic_cast<LatticeTransportElement *>(el) ) {
            le->giveGpCoordinates(coords);
        }
    };

    // Per-element kind: 0 = real polygon, 1 = synthesised polygon, 2 = line (bond), 3 = skip.
    IntArray elemKind(numberOfElements);
    IntArray vertsPerCell(numberOfElements);   // # vertices per emitted cell
    IntArray nStripsPerElem(numberOfElements); // # strips per copy (only meaningful for kind 0)
    IntArray cellsPerElem(numberOfElements);   // # cells per element

    // First pass: classify each element; for kind-1, just count incidence per node (for cap detection).
    std::map< int, int > nodeIncidence;
    FloatArray scratchPoly;
    for ( int ie = 1; ie <= numberOfElements; ie++ ) {
        Element *el = domain->giveElement(elements.at(ie));
        int npoly = getPolygonNodeCount(elements.at(ie));
        if ( npoly >= 3 ) {
            elemKind.at(ie) = 0;
            continue;
        }
        int synVerts = synthesiseElementPolygon(el, scratchPoly);
        if ( synVerts >= 3 ) {
            elemKind.at(ie) = 1;
            vertsPerCell.at(ie) = synVerts;
            ++nodeIncidence[el->giveDofManagerNumber(1)];
            ++nodeIncidence[el->giveDofManagerNumber(2)];
            continue;
        }
        if ( el->giveNumberOfDofManagers() == 2 ) {
            elemKind.at(ie) = 2;
            continue;
        }
        elemKind.at(ie) = 3;
    }

    // Stage 4b (fan approach): for each shell element with a 4-vertex polygon, emit 4 closure
    // triangles per element (A_top, B_top, A_bot, B_bot). Each triangle has 3 fresh vertices,
    // all under one node's rigid-body kinematics. No watertightness, no sorting, no dedup.
    // Counting happens in the second pass below.

    // Detect whether a shell element sits at the specimen boundary. The converter places one
    // pair of polygon TOP-BOT vertices at the edge midpoint (xm ± b·snorm) when there is no
    // adjacent triangle on that side — i.e. the component of (vertex - xm) orthogonal to snorm
    // collapses to zero. Returns true iff exactly one TOP-BOT pair has this property → emit
    // 4 side tris. Using the shell-normal-orthogonal magnitude avoids dependence on the
    // converter's polygon lateral axis (which is not snorm × elemAxis on curved shells).
    auto shellHasBoundarySide = [&](Element *el) -> bool {
#ifdef __SM_MODULE
        auto *l3d = dynamic_cast< Lattice3d * >( el );
        if ( !l3d || !l3d->isShellElement() ) return false;
        FloatArray poly;
        l3d->giveCrossSectionCoordinates(poly);
        if ( poly.giveSize() != 12 ) return false;
        const FloatArray &snorm = l3d->giveShellNormal();
        const FloatArray &cA = l3d->giveNode(1)->giveCoordinates();
        const FloatArray &cB = l3d->giveNode(2)->giveCoordinates();
        FloatArray xm(3), e(3);
        for ( int j = 1; j <= 3; ++j ) {
            xm.at(j) = 0.5 * ( cA.at(j) + cB.at(j) );
            e.at(j)  = cB.at(j) - cA.at(j);
        }
        double elemLen = e.computeNorm();
        if ( elemLen <= 0 ) return false;
        const double tol = 1e-3 * elemLen;
        int nBndTop = 0, nBndBot = 0;
        for ( int k = 0; k < 4; ++k ) {
            double dn = 0;
            for ( int j = 1; j <= 3; ++j ) dn += ( poly.at(3*k + j) - xm.at(j) ) * snorm.at(j);
            double latSq = 0;
            for ( int j = 1; j <= 3; ++j ) {
                double dp = ( poly.at(3*k + j) - xm.at(j) ) - dn * snorm.at(j);
                latSq += dp * dp;
            }
            if ( std::sqrt(latSq) < tol ) { ( dn > 0 ) ? ++nBndTop : ++nBndBot; }
        }
        return ( nBndTop == 1 && nBndBot == 1 );
#else
        (void) el;
        return false;
#endif
    };

    // Second pass: count cells (now we know terminal/interior status per node).
    int totalNodes = 0, totalCells = 0;
    for ( int ie = 1; ie <= numberOfElements; ie++ ) {
        Element *el = domain->giveElement(elements.at(ie));
        const int kind = elemKind.at(ie);
        if ( kind == 0 ) {
            int npoly = getPolygonNodeCount(elements.at(ie));
            FloatArray dummyNormal;
            int nStrips = getElementStripCount(el, dummyNormal);
            int vpc = ( nStrips > 1 && npoly == 4 ) ? 4 : npoly;
            vertsPerCell.at(ie) = vpc;
            nStripsPerElem.at(ie) = nStrips;
            int nClosure = 0;
#ifdef __SM_MODULE
            // 4 fan triangles per shell element with a 4-vertex polygon (A_top, B_top, A_bot, B_bot).
            auto *l3dCnt = dynamic_cast< Lattice3d * >( el );
            if ( l3dCnt && l3dCnt->isShellElement() && npoly == 4 ) {
                nClosure = 4;
                // +4 fan triangles for the boundary side (A_top→poly_top→poly_bot→A_bot and B-mirror).
                if ( shellHasBoundarySide(el) ) nClosure += 4;
            }
#endif
            cellsPerElem.at(ie) = 3 * nStrips + nClosure;
            totalNodes += 3 * nStrips * vpc + 3 * nClosure;
            totalCells += 3 * nStrips + nClosure;
        } else if ( kind == 1 ) {
            int N = vertsPerCell.at(ie);
            int nA = el->giveDofManagerNumber(1);
            int nB = el->giveDofManagerNumber(2);
            int capA = ( nodeIncidence[nA] == 1 ) ? 1 : 0;
            int capB = ( nodeIncidence[nB] == 1 ) ? 1 : 0;
            nStripsPerElem.at(ie) = 1;
            // Cells: 3 midpoint + 2N lateral quads + capA + capB
            int nCells = 3 + 2 * N + capA + capB;
            // Vertices: 3N midpoint copies + 2N joint vertices (caps share joint vertices)
            int nV = 5 * N;
            cellsPerElem.at(ie) = nCells;
            totalNodes += nV;
            totalCells += nCells;
        } else if ( kind == 2 ) {
            vertsPerCell.at(ie) = 2;
            nStripsPerElem.at(ie) = 1;
            cellsPerElem.at(ie) = 1;
            totalNodes += 2;
            totalCells += 1;
        }
    }
    if ( totalNodes == 0 || totalCells == 0 ) return;

    // Stash per-element layout into members for exportCellVarsCross.
    elemKindCross = elemKind;
    cellsPerElemCross = cellsPerElem;
    vertsPerCellCross = vertsPerCell;
    nStripsPerElemCross = nStripsPerElem;
    capCountCross.resize(numberOfElements);
    for ( int ie = 1; ie <= numberOfElements; ie++ ) {
        if ( elemKind.at(ie) != 1 ) { capCountCross.at(ie) = 0; continue; }
        Element *el = domain->giveElement(elements.at(ie));
        int nA = el->giveDofManagerNumber(1);
        int nB = el->giveDofManagerNumber(2);
        int c = 0;
        if ( nodeIncidence[nA] == 1 ) ++c;
        if ( nodeIncidence[nB] == 1 ) ++c;
        capCountCross.at(ie) = c;
    }

    // (Per-element extrusion uses the element's own axis directly — no per-node joint axis needed.)

    vtkPieceCross.setNumberOfNodes(totalNodes);
    vtkPieceCross.setNumberOfCells(totalCells);

    // Pre-size the per-cell polygonRole tag (0 = A-side, 1 = B-side, 2 = midline, -1 = fallback point).
    polygonRoleCross.resize(totalCells);

    // Per-vertex displacement (so ParaView Warp by Vector can show the deformed bodies).
    displacementCross.resize(totalNodes, 3); displacementCross.zero();

    FloatArray polyRef, vRef(3), vDef(3);
    FloatArray nodeRefA(3), nodeRefB(3), dispA, dispB, rotA, rotB;
    IntArray connectivity;
    FloatArray coords(3);

    syntheticCross.resize(totalCells); syntheticCross.zero();

    int nodeOffset = 0, cellIdx = 0, vtkOffset = 0;
    for ( int ie = 1; ie <= numberOfElements; ie++ ) {
        Element *el = domain->giveElement(elements.at(ie));
        const int vpc = vertsPerCell.at(ie);
        const int nStrips = nStripsPerElem.at(ie);
        const int kind = elemKind.at(ie);

        if ( kind == 3 ) continue;  // skipped

        if ( kind == 0 ) {
            getPolygonCoords(elements.at(ie), polyRef);
            const int polyNV = polyRef.giveSize() / 3;
            FloatArray shellNormal;
            getElementStripCount(el, shellNormal);

            // Subdivide polygon (reference frame) into nStrips quads when hybrid.
            FloatMatrix stripVerts;
            subdivideShellPolygon(polyRef, polyNV, shellNormal, nStrips, stripVerts);
            const int totalStripV = nStrips * vpc;

            Node *nA = el->giveNode(1);
            Node *nB = el->giveNumberOfDofManagers() >= 2 ? el->giveNode(2) : nA;
            const FloatArray &cA = nA->giveCoordinates();
            const FloatArray &cB = nB->giveCoordinates();
            for ( int i = 1; i <= 3; ++i ) {
                nodeRefA.at(i) = cA.at(i);
                nodeRefB.at(i) = cB.at(i);
            }
            giveNodeKinematics(nA, tStep, dispA, rotA);
            giveNodeKinematics(nB, tStep, dispB, rotB);

            // Deformed strip vertices under A's and B's rigid-body kinematics.
            FloatMatrix vA(totalStripV, 3), vB(totalStripV, 3);
            for ( int k = 0; k < totalStripV; ++k ) {
                for ( int j = 1; j <= 3; ++j ) vRef.at(j) = stripVerts.at(k+1, j);
                rigidBodyTransform(vDef, vRef, nodeRefA, dispA, rotA);
                for ( int j = 1; j <= 3; ++j ) vA.at(k+1, j) = vDef.at(j);
                rigidBodyTransform(vDef, vRef, nodeRefB, dispB, rotB);
                for ( int j = 1; j <= 3; ++j ) vB.at(k+1, j) = vDef.at(j);
            }

            // Emit 3 * nStrips cells: copy ∈ {0,1,2} × strip ∈ 0..nStrips-1.
            for ( int copy = 0; copy < 3; ++copy ) {
                for ( int s = 0; s < nStrips; ++s ) {
                    connectivity.resize(vpc);
                    for ( int k = 1; k <= vpc; ++k ) {
                        const int stripVertRow = s * vpc + k;  // 1-based
                        int nodeIdx = nodeOffset + k;
                        for ( int j = 1; j <= 3; ++j ) coords.at(j) = stripVerts.at(stripVertRow, j);
                        vtkPieceCross.setNodeCoords(nodeIdx, coords);
                        for ( int j = 1; j <= 3; ++j ) {
                            double def = ( copy == 0 ) ? vA.at(stripVertRow, j)
                                       : ( copy == 1 ) ? vB.at(stripVertRow, j)
                                       : 0.5 * ( vA.at(stripVertRow, j) + vB.at(stripVertRow, j) );
                            displacementCross.at(nodeIdx, j) = def - stripVerts.at(stripVertRow, j);
                        }
                        connectivity.at(k) = nodeIdx;
                    }
                    ++cellIdx;
                    vtkPieceCross.setConnectivity(cellIdx, connectivity);
                    vtkPieceCross.setCellType(cellIdx, 7);  // VTK_POLYGON
                    vtkOffset += vpc;
                    vtkPieceCross.setOffset(cellIdx, vtkOffset);
                    polygonRoleCross.at(cellIdx) = copy;
                    nodeOffset += vpc;
                }
            }

            // Stage 4b fan closures: per shell element with 4-vertex polygon, emit 4 triangles:
            // (A_top, polyTop1, polyTop2), (B_top, ...), (A_bot, polyBot1, polyBot2), (B_bot, ...).
            // Each triangle's 3 vertices use one node's rigid-body kinematics. No watertightness.
#ifdef __SM_MODULE
            auto *l3dShell = dynamic_cast< Lattice3d * >( el );
            if ( l3dShell && l3dShell->isShellElement() && polyNV == 4 ) {
                const FloatArray &snorm = l3dShell->giveShellNormal();
                // Identify TOP / BOT polygon vertex indices by shell-normal projection.
                FloatArray polyCntr(3); polyCntr.zero();
                for ( int k = 0; k < polyNV; ++k ) for ( int j = 1; j <= 3; ++j ) polyCntr.at(j) += polyRef.at(3*k + j) / polyNV;
                int topIdx[2] = { -1, -1 };
                int botIdx[2] = { -1, -1 };
                int nT = 0, nB = 0;
                double maxProj = 0, minProj = 0;
                for ( int k = 0; k < polyNV; ++k ) {
                    double p = 0;
                    for ( int j = 1; j <= 3; ++j ) p += ( polyRef.at(3*k + j) - polyCntr.at(j) ) * snorm.at(j);
                    if ( p > 0 ) {
                        if ( nT < 2 ) topIdx[nT++] = k;
                        if ( p > maxProj ) maxProj = p;
                    } else {
                        if ( nB < 2 ) botIdx[nB++] = k;
                        if ( p < minProj ) minProj = p;
                    }
                }
                if ( nT == 2 && nB == 2 ) {
                    const double halfH_top =  maxProj;
                    const double halfH_bot = -minProj;
                    auto emitFanTriangle = [&](int nodeOwner, int surfId) {
                        const FloatArray &nodeRefOwn = ( nodeOwner == 0 ) ? nodeRefA : nodeRefB;
                        const FloatArray &dispOwn    = ( nodeOwner == 0 ) ? dispA    : dispB;
                        const FloatArray &rotOwn     = ( nodeOwner == 0 ) ? rotA     : rotB;
                        const int *idx = ( surfId == 0 ) ? topIdx : botIdx;
                        const double off = ( surfId == 0 ) ? halfH_top : -halfH_bot;
                        // Three reference positions: node_top/bot + 2 polygon vertices.
                        FloatArray nodeCap(3);
                        for ( int j = 1; j <= 3; ++j ) nodeCap.at(j) = nodeRefOwn.at(j) + off * snorm.at(j);
                        FloatArray triRef[3];
                        triRef[0] = nodeCap;
                        for ( int t = 0; t < 2; ++t ) {
                            triRef[t+1].resize(3);
                            for ( int j = 1; j <= 3; ++j ) triRef[t+1].at(j) = polyRef.at(3 * idx[t] + j);
                        }
                        connectivity.resize(3);
                        for ( int t = 0; t < 3; ++t ) {
                            int nodeIdx = nodeOffset + t + 1;
                            for ( int j = 1; j <= 3; ++j ) coords.at(j) = triRef[t].at(j);
                            vtkPieceCross.setNodeCoords(nodeIdx, coords);
                            FloatArray vDefFan;
                            rigidBodyTransform(vDefFan, triRef[t], nodeRefOwn, dispOwn, rotOwn);
                            for ( int j = 1; j <= 3; ++j ) displacementCross.at(nodeIdx, j) = vDefFan.at(j) - triRef[t].at(j);
                            connectivity.at(t + 1) = nodeIdx;
                        }
                        ++cellIdx;
                        vtkPieceCross.setConnectivity(cellIdx, connectivity);
                        vtkPieceCross.setCellType(cellIdx, 7);  // VTK_POLYGON (3-vertex)
                        vtkOffset += 3;
                        vtkPieceCross.setOffset(cellIdx, vtkOffset);
                        polygonRoleCross.at(cellIdx) = ( surfId == 0 ) ? 5 : 6;
                        syntheticCross.at(cellIdx) = 1;
                        nodeOffset += 3;
                    };
                    emitFanTriangle(0, 0);  // A-top
                    emitFanTriangle(1, 0);  // B-top
                    emitFanTriangle(0, 1);  // A-bot
                    emitFanTriangle(1, 1);  // B-bot

                    // Side closure: locate the TOP-BOT polygon pair at the edge midpoint
                    // (boundary side) and emit 4 triangles tying it to A_top/A_bot and B_top/B_bot.
                    // Boundary criterion: shell-normal-orthogonal component of (vertex - xm)
                    // is ~0 — independent of the converter's lateral basis (which deviates from
                    // snorm × elemAxis on curved shells).
                    {
                        FloatArray e(3);
                        for ( int j = 1; j <= 3; ++j ) e.at(j) = nodeRefB.at(j) - nodeRefA.at(j);
                        double elemLen = e.computeNorm();
                        FloatArray xm(3);
                        for ( int j = 1; j <= 3; ++j ) xm.at(j) = 0.5 * ( nodeRefA.at(j) + nodeRefB.at(j) );
                        if ( elemLen > 0 ) {
                            const double tol = 1e-3 * elemLen;
                            auto vertLat = [&](int idx) {
                                double dn = 0;
                                for ( int j = 1; j <= 3; ++j ) dn += ( polyRef.at(3 * idx + j) - xm.at(j) ) * snorm.at(j);
                                double latSq = 0;
                                for ( int j = 1; j <= 3; ++j ) {
                                    double dp = ( polyRef.at(3 * idx + j) - xm.at(j) ) - dn * snorm.at(j);
                                    latSq += dp * dp;
                                }
                                return std::sqrt(latSq);
                            };
                            int bndTop = -1, bndBot = -1;
                            int nBndTop = 0, nBndBot = 0;
                            for ( int t = 0; t < 2; ++t ) if ( vertLat(topIdx[t]) < tol ) { bndTop = topIdx[t]; ++nBndTop; }
                            for ( int b = 0; b < 2; ++b ) if ( vertLat(botIdx[b]) < tol ) { bndBot = botIdx[b]; ++nBndBot; }
                            if ( nBndTop == 1 && nBndBot == 1 ) {
                                FloatArray polyBndTop(3), polyBndBot(3);
                                for ( int j = 1; j <= 3; ++j ) {
                                    polyBndTop.at(j) = polyRef.at(3 * bndTop + j);
                                    polyBndBot.at(j) = polyRef.at(3 * bndBot + j);
                                }
                                auto emitSideTriangle = [&](int nodeOwner, const FloatArray &triR0, const FloatArray &triR1, const FloatArray &triR2) {
                                    const FloatArray &nodeRefOwn = ( nodeOwner == 0 ) ? nodeRefA : nodeRefB;
                                    const FloatArray &dispOwn    = ( nodeOwner == 0 ) ? dispA    : dispB;
                                    const FloatArray &rotOwn     = ( nodeOwner == 0 ) ? rotA     : rotB;
                                    const FloatArray *triRef[3] = { &triR0, &triR1, &triR2 };
                                    connectivity.resize(3);
                                    for ( int t = 0; t < 3; ++t ) {
                                        int nodeIdx = nodeOffset + t + 1;
                                        for ( int j = 1; j <= 3; ++j ) coords.at(j) = triRef[t]->at(j);
                                        vtkPieceCross.setNodeCoords(nodeIdx, coords);
                                        FloatArray vDefSide;
                                        rigidBodyTransform(vDefSide, *triRef[t], nodeRefOwn, dispOwn, rotOwn);
                                        for ( int j = 1; j <= 3; ++j ) displacementCross.at(nodeIdx, j) = vDefSide.at(j) - triRef[t]->at(j);
                                        connectivity.at(t + 1) = nodeIdx;
                                    }
                                    ++cellIdx;
                                    vtkPieceCross.setConnectivity(cellIdx, connectivity);
                                    vtkPieceCross.setCellType(cellIdx, 7);  // VTK_POLYGON (3-vertex)
                                    vtkOffset += 3;
                                    vtkPieceCross.setOffset(cellIdx, vtkOffset);
                                    polygonRoleCross.at(cellIdx) = 7;
                                    syntheticCross.at(cellIdx) = 1;
                                    nodeOffset += 3;
                                };
                                FloatArray aTop(3), aBot(3), bTop(3), bBot(3);
                                for ( int j = 1; j <= 3; ++j ) {
                                    aTop.at(j) = nodeRefA.at(j) + halfH_top * snorm.at(j);
                                    aBot.at(j) = nodeRefA.at(j) - halfH_bot * snorm.at(j);
                                    bTop.at(j) = nodeRefB.at(j) + halfH_top * snorm.at(j);
                                    bBot.at(j) = nodeRefB.at(j) - halfH_bot * snorm.at(j);
                                }
                                // A-side quad split into two triangles, both under A's kinematics.
                                emitSideTriangle(0, aTop, polyBndTop, polyBndBot);
                                emitSideTriangle(0, aTop, polyBndBot, aBot);
                                // B-side quad split into two triangles, both under B's kinematics.
                                emitSideTriangle(1, bTop, polyBndTop, polyBndBot);
                                emitSideTriangle(1, bTop, polyBndBot, bBot);
                            }
                        }
                    }
                }
            }
#endif
        } else if ( kind == 1 ) {
            // Synthesised polygon (frame/rebar): 3 midpoint cells + 2N lateral quads + 0–2 caps.
            // Vertex layout (per element): 0..N-1 = midpoint A copy, N..2N-1 = B copy, 2N..3N-1 = midline copy,
            // 3N..4N-1 = joint at A, 4N..5N-1 = joint at B.
            synthesiseElementPolygon(el, polyRef);
            const int N = vpc;
            Node *nA = el->giveNode(1);
            Node *nB = el->giveNode(2);
            const FloatArray &cA = nA->giveCoordinates();
            const FloatArray &cB = nB->giveCoordinates();
            for ( int i = 1; i <= 3; ++i ) {
                nodeRefA.at(i) = cA.at(i);
                nodeRefB.at(i) = cB.at(i);
            }
            giveNodeKinematics(nA, tStep, dispA, rotA);
            giveNodeKinematics(nB, tStep, dispB, rotB);

            // Polygon radius: rotational distance from polygon centroid to its first vertex.
            FloatArray polyCentroid(3); polyCentroid.zero();
            for ( int k = 0; k < N; ++k ) for ( int j = 1; j <= 3; ++j ) polyCentroid.at(j) += polyRef.at(3*k + j) / N;
            // Vertex 0 distance to centroid → polygon "radius" (for octagon = R; for rect we just use vertex0).
            FloatArray vec0(3);
            for ( int j = 1; j <= 3; ++j ) vec0.at(j) = polyRef.at(j) - polyCentroid.at(j);
            const double polyR = vec0.computeNorm();
            // Phase angle of vertex 0 in the element's (e_y, e_z) plane.
            FloatArray elemEx, elemEy, elemEz;
            buildElementFrame(nodeRefA, nodeRefB, elemEx, elemEy, elemEz);
            double phase0 = std::atan2( vec0.dotProduct(elemEz), vec0.dotProduct(elemEy) );

            // Per-element extrusion: end-polygons are the midpoint polygon translated to each node,
            // using the element's own (e_y, e_z) frame. This automatically allows any incidence count
            // (no per-node joint axis averaging), at the cost of non-watertight angled junctions.
            auto buildEndVerts = [&](const FloatArray &nodePos, FloatMatrix &out) {
                out.resize(N, 3);
                for ( int k = 0; k < N; ++k ) {
                    double angle = phase0 + 2.0 * M_PI * k / N;
                    double cy = polyR * std::cos(angle);
                    double cz = polyR * std::sin(angle);
                    for ( int j = 1; j <= 3; ++j ) out.at(k+1, j) = nodePos.at(j) + cy * elemEy.at(j) + cz * elemEz.at(j);
                }
            };

            FloatMatrix jointRefA, jointRefB;
            buildEndVerts(nodeRefA, jointRefA);
            buildEndVerts(nodeRefB, jointRefB);

            FloatMatrix vA(N, 3), vB(N, 3), vJointA(N, 3), vJointB(N, 3);
            for ( int k = 0; k < N; ++k ) {
                for ( int j = 1; j <= 3; ++j ) vRef.at(j) = polyRef.at(3*k + j);
                rigidBodyTransform(vDef, vRef, nodeRefA, dispA, rotA);
                for ( int j = 1; j <= 3; ++j ) vA.at(k+1, j) = vDef.at(j);
                rigidBodyTransform(vDef, vRef, nodeRefB, dispB, rotB);
                for ( int j = 1; j <= 3; ++j ) vB.at(k+1, j) = vDef.at(j);
                for ( int j = 1; j <= 3; ++j ) vRef.at(j) = jointRefA.at(k+1, j);
                rigidBodyTransform(vDef, vRef, nodeRefA, dispA, rotA);
                for ( int j = 1; j <= 3; ++j ) vJointA.at(k+1, j) = vDef.at(j);
                for ( int j = 1; j <= 3; ++j ) vRef.at(j) = jointRefB.at(k+1, j);
                rigidBodyTransform(vDef, vRef, nodeRefB, dispB, rotB);
                for ( int j = 1; j <= 3; ++j ) vJointB.at(k+1, j) = vDef.at(j);
            }

            // Emit 3 midpoint copies (A, B, midline) — vertex indices 1..N, N+1..2N, 2N+1..3N.
            const int baseMidA   = nodeOffset;
            const int baseMidB   = nodeOffset + N;
            const int baseMidM   = nodeOffset + 2*N;
            const int baseJointA = nodeOffset + 3*N;
            const int baseJointB = nodeOffset + 4*N;
            for ( int copy = 0; copy < 3; ++copy ) {
                const int base = nodeOffset + copy * N;
                connectivity.resize(N);
                for ( int k = 1; k <= N; ++k ) {
                    int nodeIdx = base + k;
                    for ( int j = 1; j <= 3; ++j ) coords.at(j) = polyRef.at(3*(k-1) + j);
                    vtkPieceCross.setNodeCoords(nodeIdx, coords);
                    for ( int j = 1; j <= 3; ++j ) {
                        double def = ( copy == 0 ) ? vA.at(k, j)
                                   : ( copy == 1 ) ? vB.at(k, j)
                                   : 0.5 * ( vA.at(k, j) + vB.at(k, j) );
                        displacementCross.at(nodeIdx, j) = def - polyRef.at(3*(k-1) + j);
                    }
                    connectivity.at(k) = nodeIdx;
                }
                ++cellIdx;
                vtkPieceCross.setConnectivity(cellIdx, connectivity);
                vtkPieceCross.setCellType(cellIdx, 7);  // VTK_POLYGON
                vtkOffset += N;
                vtkPieceCross.setOffset(cellIdx, vtkOffset);
                polygonRoleCross.at(cellIdx) = copy;
                syntheticCross.at(cellIdx) = 1;
            }

            // Joint vertex sets at A and B (indices 3N+1..4N, 4N+1..5N).
            for ( int k = 1; k <= N; ++k ) {
                int idxA = baseJointA + k;
                for ( int j = 1; j <= 3; ++j ) coords.at(j) = jointRefA.at(k, j);
                vtkPieceCross.setNodeCoords(idxA, coords);
                for ( int j = 1; j <= 3; ++j ) displacementCross.at(idxA, j) = vJointA.at(k, j) - jointRefA.at(k, j);
                int idxB = baseJointB + k;
                for ( int j = 1; j <= 3; ++j ) coords.at(j) = jointRefB.at(k, j);
                vtkPieceCross.setNodeCoords(idxB, coords);
                for ( int j = 1; j <= 3; ++j ) displacementCross.at(idxB, j) = vJointB.at(k, j) - jointRefB.at(k, j);
            }
            nodeOffset += 5 * N;

            // Lateral quads on A-side: midpoint A-vertex k → joint A-vertex k. Quad k spans k..k+1.
            for ( int k = 0; k < N; ++k ) {
                int kn = ( k + 1 ) % N;
                connectivity.resize(4);
                connectivity.at(1) = baseMidA   + k  + 1;
                connectivity.at(2) = baseMidA   + kn + 1;
                connectivity.at(3) = baseJointA + kn + 1;
                connectivity.at(4) = baseJointA + k  + 1;
                ++cellIdx;
                vtkPieceCross.setConnectivity(cellIdx, connectivity);
                vtkPieceCross.setCellType(cellIdx, 7);
                vtkOffset += 4;
                vtkPieceCross.setOffset(cellIdx, vtkOffset);
                polygonRoleCross.at(cellIdx) = 3;  // lateral
                syntheticCross.at(cellIdx) = 1;
            }
            // Lateral quads on B-side.
            for ( int k = 0; k < N; ++k ) {
                int kn = ( k + 1 ) % N;
                connectivity.resize(4);
                connectivity.at(1) = baseMidB   + k  + 1;
                connectivity.at(2) = baseMidB   + kn + 1;
                connectivity.at(3) = baseJointB + kn + 1;
                connectivity.at(4) = baseJointB + k  + 1;
                ++cellIdx;
                vtkPieceCross.setConnectivity(cellIdx, connectivity);
                vtkPieceCross.setCellType(cellIdx, 7);
                vtkOffset += 4;
                vtkPieceCross.setOffset(cellIdx, vtkOffset);
                polygonRoleCross.at(cellIdx) = 3;
                syntheticCross.at(cellIdx) = 1;
            }
            // Caps at terminal nodes.
            if ( nodeIncidence[el->giveDofManagerNumber(1)] == 1 ) {
                connectivity.resize(N);
                for ( int k = 1; k <= N; ++k ) connectivity.at(k) = baseJointA + k;
                ++cellIdx;
                vtkPieceCross.setConnectivity(cellIdx, connectivity);
                vtkPieceCross.setCellType(cellIdx, 7);
                vtkOffset += N;
                vtkPieceCross.setOffset(cellIdx, vtkOffset);
                polygonRoleCross.at(cellIdx) = 4;  // cap
                syntheticCross.at(cellIdx) = 1;
            }
            if ( nodeIncidence[el->giveDofManagerNumber(2)] == 1 ) {
                connectivity.resize(N);
                for ( int k = 1; k <= N; ++k ) connectivity.at(k) = baseJointB + k;
                ++cellIdx;
                vtkPieceCross.setConnectivity(cellIdx, connectivity);
                vtkPieceCross.setCellType(cellIdx, 7);
                vtkOffset += N;
                vtkPieceCross.setOffset(cellIdx, vtkOffset);
                polygonRoleCross.at(cellIdx) = 4;
                syntheticCross.at(cellIdx) = 1;
            }
        } else if ( kind == 2 ) {
            // Line (bond link): VTK_LINE between the two nodes, each endpoint with its own kinematics.
            Node *nA = el->giveNode(1);
            Node *nB = el->giveNode(2);
            const FloatArray &cA = nA->giveCoordinates();
            const FloatArray &cB = nB->giveCoordinates();
            giveNodeKinematics(nA, tStep, dispA, rotA);
            giveNodeKinematics(nB, tStep, dispB, rotB);

            for ( int j = 1; j <= 3; ++j ) coords.at(j) = cA.at(j);
            vtkPieceCross.setNodeCoords(nodeOffset + 1, coords);
            for ( int j = 1; j <= 3; ++j ) displacementCross.at(nodeOffset + 1, j) = dispA.at(j);

            for ( int j = 1; j <= 3; ++j ) coords.at(j) = cB.at(j);
            vtkPieceCross.setNodeCoords(nodeOffset + 2, coords);
            for ( int j = 1; j <= 3; ++j ) displacementCross.at(nodeOffset + 2, j) = dispB.at(j);

            connectivity.resize(2);
            connectivity.at(1) = nodeOffset + 1;
            connectivity.at(2) = nodeOffset + 2;
            ++cellIdx;
            vtkPieceCross.setConnectivity(cellIdx, connectivity);
            vtkPieceCross.setCellType(cellIdx, 3);  // VTK_LINE
            vtkOffset += 2;
            vtkPieceCross.setOffset(cellIdx, vtkOffset);
            polygonRoleCross.at(cellIdx) = -2;
            syntheticCross.at(cellIdx) = 1;
            nodeOffset += 2;
        }
    }

    // Cell variables: edge data only on the midline cell (role 2); closures get zeros.
    this->exportCellVarsCross(vtkPieceCross, region, cellVarsToExport, tStep);
}


void
VTKXMLLatticeExportModule::exportPrimaryVarsCross(ExportRegion &vtkPiece, Set &region, IntArray& primaryVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    FloatArray valueArray;
    smoother.clear(); // Makes sure primary smoother is up-to-date with potentially new mesh.


    //This is a way of plotting cell variables in the form of points in the current state.
    //Displacements are at nodes, but points are at midpoints of elements. Interpolate the displacements so that a fictious displacement at midpoint is obtained. This does not show th jump of rotations.
    
    //const IntArray& mapG2L = vtkPiece.getMapG2L();
    //    const IntArray& mapL2G = vtkPiece.getMapL2G();

    //    vtkPiece.setNumberOfPrimaryVarsToExport(primaryVarsToExport, mapL2G.giveSize() );
    vtkPiece.setNumberOfPrimaryVarsToExport(primaryVarsToExport, d->giveNumberOfElements() );

    FloatArray saveAverage;
    for ( int i = 1, n = primaryVarsToExport.giveSize(); i <= n; i++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);
        int nodeCounter = 0;
        for ( int ielem = 1; ielem <= d->giveNumberOfElements(); ielem++ ) {



            saveAverage.zero();
	  Element *elem = d->giveElement(ielem);

          if (!(dynamic_cast<LatticeFrame3d *>(elem) || dynamic_cast<LatticeFrame3dNL *>(elem))) {
              continue;
          }

	  for (int inode = 1; inode<=elem->giveNumberOfNodes();inode++){
	    DofManager *dman = elem->giveNode(inode);
	    this->getNodalVariableFromPrimaryField(valueArray, dman, tStep, type, region, smoother);
	    saveAverage += valueArray;
	  }
          saveAverage.times(1.0 / elem->giveNumberOfNodes());

          //This is the problem. I need to store it in the cross-section node and not in the element.
          //For the special case, there are
          //Find the cross-section nodes and then loop over
          //Debug

         int numberOfCrossSectionNodes = 0;
          if (  dynamic_cast< LatticeStructuralElement * >(elem) ) {
            numberOfCrossSectionNodes =  ( static_cast< LatticeStructuralElement * >( elem) )->giveNumberOfCrossSectionNodes();
          } else if ( dynamic_cast< LatticeTransportElement * >( elem ) ) {
              numberOfCrossSectionNodes = ( static_cast<LatticeTransportElement *>( elem ) )->giveNumberOfCrossSectionNodes();
          }
          for(int k = 1; k<=numberOfCrossSectionNodes; k++){
              FloatArray average;
              average = saveAverage;
	    nodeCounter++;
            printf("nodeCounter = %d\n", nodeCounter);
	    vtkPiece.setPrimaryVarInNode(type, nodeCounter, std::move(average) );
          }
        }
    }
}



void
VTKXMLLatticeExportModule::doOutput(TimeStep *tStep, bool forcedOutput)
{
    // Did the main (line) VTU actually get written? It writes iff testTimeStepOutput passes.
    const bool lineWillWrite = ( testTimeStepOutput(tStep) || forcedOutput );

    this->doOutputNormal(tStep, forcedOutput);

    if ( crossSectionExportFlag ) {
        // Cross VTU is gated by both the line-output gate AND the crossOutputStep multiplier.
        // The counter ticks once per line-output trigger, so crossstep=N means "1 in N line outputs".
        bool crossWillWrite = false;
        if ( lineWillWrite ) {
            if ( forcedOutput || ( this->crossOutputCounter % this->crossOutputStep == 0 ) ) {
                crossWillWrite = true;
            }
            ++this->crossOutputCounter;
        }
        if ( crossWillWrite ) {
            this->doOutputCross(tStep, forcedOutput);
        }
    }
    this->defaultVTKPiece.clear();
}


void
VTKXMLLatticeExportModule::doOutputCross(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }

    this->fileStreamCross = this->giveOutputStreamCross(tStep);
    struct tm *current;
    time_t now;
    time(& now);
    current = localtime(& now);

    this->fileStreamCross << "<!-- TimeStep " << tStep->giveTargetTime() * timeScale << " Computed " << current->tm_year + 1900 << "-" << setw(2) << current->tm_mon + 1 << "-" << setw(2) << current->tm_mday << " at " << current->tm_hour << ":" << current->tm_min << ":" << setw(2) << current->tm_sec << " -->\n";
    this->fileStreamCross << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    this->fileStreamCross << "<UnstructuredGrid>\n";

    int nPiecesToExport = this->giveNumberOfRegions(); //old name: region, meaning: sets
    int anyPieceNonEmpty = 0;

    for ( int pieceNum = 1; pieceNum <= nPiecesToExport; pieceNum++ ) {
        Set* region = this->giveRegionSet(pieceNum);
        // Fills a data struct (VTKPiece) with all the necessary data.
        this->setupVTKPieceCross(this->defaultVTKPieceCross, tStep, *region);

        // Write the VTK piece to file.
        anyPieceNonEmpty += this->writeVTKPieceCross(this->defaultVTKPieceCross, tStep);
        this->defaultVTKPieceCross.clear();
    }


    if ( anyPieceNonEmpty == 0 ) {
        // write empty piece, Otherwise ParaView complains if the whole vtu file is without <Piece></Piece>
        this->fileStreamCross << "<Piece NumberOfPoints=\"0\" NumberOfCells=\"0\">\n";
        this->fileStreamCross << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"> </DataArray>\n</Cells>\n";
        this->fileStreamCross << "</Piece>\n";
    }

    this->fileStreamCross << "</UnstructuredGrid>\n</VTKFile>";
    this->fileStreamCross.close();

    this->appendPvdEntryCross(tStep, this->giveOutputFileNameCross(tStep));
}


void
VTKXMLLatticeExportModule::doOutputNormal(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }

    this->fileStream = this->giveOutputStream(tStep);
    struct tm *current;
    time_t now;
    time(& now);
    current = localtime(& now);

    this->fileStream << "<!-- TimeStep " << tStep->giveTargetTime() * timeScale << " Computed " << current->tm_year + 1900 << "-" << setw(2) << current->tm_mon + 1 << "-" << setw(2) << current->tm_mday << " at " << current->tm_hour << ":" << current->tm_min << ":" << setw(2) << current->tm_sec << " -->\n";

    this->fileStream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    this->fileStream << "<UnstructuredGrid>\n";

    this->giveSmoother(); // make sure smoother is created, Necessary? If it doesn't exist it is created /JB

    int nPiecesToExport = this->giveNumberOfRegions();     //old name: region, meaning: sets
    int anyPieceNonEmpty = 0;

    NodalRecoveryModel *smoother = giveSmoother();
    NodalRecoveryModel *primVarSmoother = givePrimVarSmoother();

    for ( int pieceNum = 1; pieceNum <= nPiecesToExport; pieceNum++ ) {
        // Fills a data struct (VTKPiece) with all the necessary data.
        Set *region = this->giveRegionSet(pieceNum);
        this->setupVTKPiece(this->defaultVTKPiece, tStep, *region);
        this->writeVTKPieceProlog(this->defaultVTKPiece, tStep);  
        // Export primary, internal and XFEM variables as nodal quantities
        this->exportPrimaryVars(this->defaultVTKPiece, *region, primaryVarsToExport, *primVarSmoother, tStep);
        this->exportIntVars(this->defaultVTKPiece, *region, internalVarsToExport, *smoother, tStep);
        //const IntArray &elements = region.giveElementList();
        this->exportCellVars(this->defaultVTKPiece, *region, cellVarsToExport, tStep);
        // Write the VTK piece to file.
        anyPieceNonEmpty += this->writeVTKPieceVariables(this->defaultVTKPiece, tStep);
        this->writeVTKPieceEpilog(this->defaultVTKPiece, tStep);  
    }

    /*
     * Output all composite elements - one piece per composite element
     * Each element is responsible of setting up a VTKPiece which can then be exported
     */
    Domain *d = emodel->giveDomain(1);

    for ( int pieceNum = 1; pieceNum <= nPiecesToExport; pieceNum++ ) {
        const IntArray &elements = this->giveRegionSet(pieceNum)->giveElementList();
        for ( int i = 1; i <= elements.giveSize(); i++ ) {
            Element *el = d->giveElement(elements.at(i) );
            if ( this->isElementComposite(el) ) {
                if ( el->giveParallelMode() != Element_local ) {
                    continue;
                }

                //this->exportCompositeElement(this->defaultVTKPiece, el, tStep);
                this->exportCompositeElement(this->defaultVTKPieces, el, tStep);

                for ( int j = 0; j < ( int ) this->defaultVTKPieces.size(); j++ ) {
                    this->writeVTKPieceProlog(this->defaultVTKPieces[j], tStep);  
                    anyPieceNonEmpty += this->writeVTKPieceVariables(this->defaultVTKPieces [ j ], tStep);
                    this->writeVTKPieceEpilog(this->defaultVTKPieces[j], tStep);  
                }
            }
        }
    }     // end loop over composite elements

    if ( anyPieceNonEmpty == 0 ) {
        // write empty piece, Otherwise ParaView complains if the whole vtu file is without <Piece></Piece>
        this->fileStream << "<Piece NumberOfPoints=\"0\" NumberOfCells=\"0\">\n";
        this->fileStream << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"> </DataArray>\n</Cells>\n";
        this->fileStream << "</Piece>\n";
    }

    this->fileStream << "</UnstructuredGrid>\n</VTKFile>";
    this->fileStream.close();

    this->appendPvdEntryLine(tStep, this->giveOutputFileName(tStep));
}


void
VTKXMLLatticeExportModule::appendPvdEntryLine(TimeStep *tStep, const std::string &vtuFilename)
{
    std::ostringstream entry;
    entry << "<DataSet timestep=\"" << tStep->giveTargetTime() * this->timeScale
          << "\" group=\"\" part=\"\" file=\"" << vtuFilename << "\"/>";
    this->pvdBufferLine.push_back(entry.str() );
    const std::string pvdName = this->emodel->giveOutputBaseFileName() + ".m" + std::to_string(this->number) + ".pvd";
    this->writePvdCollection(pvdName, this->pvdBufferLine);
}


void
VTKXMLLatticeExportModule::appendPvdEntryCross(TimeStep *tStep, const std::string &vtuFilename)
{
    std::ostringstream entry;
    entry << "<DataSet timestep=\"" << tStep->giveTargetTime() * this->timeScale
          << "\" group=\"\" part=\"\" file=\"" << vtuFilename << "\"/>";
    this->pvdBufferCross.push_back(entry.str() );
    const std::string pvdName = this->emodel->giveOutputBaseFileName() + ".m" + std::to_string(this->number) + ".cross.pvd";
    this->writePvdCollection(pvdName, this->pvdBufferCross);
}


void
VTKXMLLatticeExportModule::writePvdCollection(const std::string &pvdFilename, const std::vector< std::string > &buffer)
{
    if ( this->pythonExport ) {
        return; // python harness suppresses on-disk writes.
    }
    std::ofstream stream(pvdFilename.c_str() );
    if ( !stream.good() ) {
        OOFEM_ERROR("failed to open file %s", pvdFilename.c_str() );
    }
    stream << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n";
    for ( const auto &entry : buffer ) {
        stream << entry << "\n";
    }
    stream << "</Collection>\n</VTKFile>";
    stream.close();
}


int
VTKXMLLatticeExportModule::initRegionNodeNumbering(ExportRegion& vtkPiece, Domain *domain, TimeStep *tStep, Set& region)
{
    int nnodes = domain->giveNumberOfDofManagers();
    int elementNode, node;
    int currOffset = 1;
    Element *element;

    int  regionDofMans = 0;
    int regionSingleCells = 0;

    IntArray elements = region.giveElementList();

    int extraNodes = 0.;
    for ( int ie = 1; ie <= elements.giveSize(); ie++ ) {
        element = domain->giveElement(elements.at(ie) );

#ifdef __SM_MODULE
        if ( dynamic_cast< Lattice3dBoundary * >( element ) || dynamic_cast< LatticeLink3dBoundary * >( element ) || dynamic_cast< Lattice2dBoundary * >( element ) ) {
            IntArray loc = element->giveLocation();
            for ( int ielnode = 1; ielnode <= 2; ielnode++ ) {
                if ( loc.at(ielnode) != 0 ) {
                    extraNodes++;
                }
            }
        }
#endif
    }
    IntArray& regionG2LNodalNumbers = vtkPiece.getMapG2L();
    IntArray& regionL2GNodalNumbers = vtkPiece.getMapL2G();

    regionG2LNodalNumbers.resize(nnodes + extraNodes);
    regionG2LNodalNumbers.zero();

    periodicMap.resize(nnodes + extraNodes);
    regionToUniqueMap.resize(nnodes + extraNodes);
    locationMap.resize(nnodes + extraNodes);

    uniqueNodeTable.resize(nnodes + extraNodes, 3);
    uniqueNodeTable.zero();

    int totalNodes = 0;
    int uniqueNodes = 0;

    for ( int ie = 1; ie <= elements.giveSize(); ie++ ) {
        int ielem = elements.at(ie);
        element = domain->giveElement(ielem);

        if ( this->isElementComposite(element) ) {
            continue;                                    // composite cells exported individually
        }

        if ( !element->isActivated(tStep) ) {                    //skip inactivated elements
            continue;
        }

        if ( element->giveParallelMode() != Element_local ) {
            continue;
        }

        regionSingleCells++;
        elemNodes = element->giveNumberOfNodes();

#ifdef __SM_MODULE
        if ( dynamic_cast< Lattice3dBoundary * >( element ) || dynamic_cast< LatticeLink3dBoundary * >( element ) || dynamic_cast< Lattice2dBoundary * >( element ) ) {
            elemNodes = 2;
        }
#endif

        for ( elementNode = 1; elementNode <= elemNodes; elementNode++ ) {
#ifdef __SM_MODULE
            if ( dynamic_cast< Lattice3dBoundary * >( element ) || dynamic_cast< LatticeLink3dBoundary * >( element ) || dynamic_cast< Lattice2dBoundary * >( element ) ) { //beam elements - only unique nodes
                IntArray loc = element->giveLocation();
                FloatArray nodeCoords;
                if ( loc.at(elementNode) != 0 ) { //boundary element with mirrored node
                    totalNodes++;
                    regionDofMans++;
                    regionG2LNodalNumbers.at(nnodes + totalNodes) = 1;
                    node = element->giveNode(elementNode)->giveNumber();
                    if ( regionG2LNodalNumbers.at(node) == 0 ) { // assign new number
                        /* mark for assignment. This is done later, as it allows to preserve
                         * natural node numbering.
                         */
                        regionG2LNodalNumbers.at(node) = 1;
                        regionDofMans++;
                    }
                    if ( regionG2LNodalNumbers.at(nnodes) == 0  ) {
                        regionG2LNodalNumbers.at(nnodes) = 1;
                        regionDofMans++;
                    }

                    periodicMap.at(nnodes + totalNodes) = node;
                    locationMap.at(nnodes + totalNodes) = loc.at(elementNode);

                    element->recalculateCoordinates(elementNode, nodeCoords);

                    //get only the unique nodes
                    int repeatFlag = 0;
                    for ( int j = nnodes + 1; j <= nnodes + uniqueNodes; j++ ) {
                        double dx = fabs(uniqueNodeTable.at(j, 1) - nodeCoords.at(1) );
                        double dy = fabs(uniqueNodeTable.at(j, 2) - nodeCoords.at(2) );
                        double dz = fabs(uniqueNodeTable.at(j, 3) - nodeCoords.at(3) );
                        if ( dx < 1e-9 && dy < 1e-9 && dz < 1e-9 ) {//node already present
                            repeatFlag++;
                            regionToUniqueMap.at(nnodes + totalNodes) = j;
                        }
                    }

                    if ( repeatFlag == 0 ) {//new unique node
                        uniqueNodes++;
                        uniqueNodeTable.at(nnodes + uniqueNodes, 1) = nodeCoords.at(1);
                        uniqueNodeTable.at(nnodes + uniqueNodes, 2) = nodeCoords.at(2);
                        uniqueNodeTable.at(nnodes + uniqueNodes, 3) = nodeCoords.at(3);
                        regionToUniqueMap.at(nnodes + totalNodes) = nnodes + uniqueNodes;
                    }
                }
            } else {   //regular element
#endif
            node = element->giveNode(elementNode)->giveNumber();
            if ( regionG2LNodalNumbers.at(node) == 0 ) {     // assign new number
                /* mark for assignment. This is done later, as it allows to preserve
                 * natural node numbering.
                 */
                regionG2LNodalNumbers.at(node) = 1;
                regionDofMans++;
            }
#ifdef __SM_MODULE
        }
#endif
        }
    }

    uniqueNodeTable.resizeWithData(nnodes + uniqueNodes, 3);
    regionDofMans = nnodes + uniqueNodes;
    vtkPiece.setNumberOfNodes(regionDofMans);   
    vtkPiece.setNumberOfCells(regionSingleCells);

    regionG2LNodalNumbers.resizeWithValues(regionDofMans);
    regionL2GNodalNumbers.resize(regionDofMans);

    for ( int i = 1; i <= regionDofMans; i++ ) {
        if ( regionG2LNodalNumbers.at(i) ) {
            regionG2LNodalNumbers.at(i) = currOffset++;
            regionL2GNodalNumbers.at(regionG2LNodalNumbers.at(i) ) = i;
        }
    }
    return 1;
}


void
VTKXMLLatticeExportModule::exportPrimaryVars(ExportRegion &vtkPiece, Set& region, IntArray& primaryVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    int nnodes = d->giveNumberOfDofManagers();
    FloatArray valueArray;
    smoother.clear(); // Makes sure primary smoother is up-to-date with potentially new mesh.

    //const IntArray& mapG2L = vtkPiece.getMapG2L();
    const IntArray& mapL2G = vtkPiece.getMapL2G();

    vtkPiece.setNumberOfPrimaryVarsToExport(primaryVarsToExport, mapL2G.giveSize() );

    //Get the macroscopic field (deformation gradients, curvatures etc.)
    DofManager *controlNode = d->giveNode(nnodes);   //assuming the control node is last
    IntArray dofIdArray;
    controlNode->giveCompleteMasterDofIDArray(dofIdArray);

    FloatArray macroField(controlNode->giveNumberOfDofs() );
    for ( int j = 1; j <= controlNode->giveNumberOfDofs(); j++ ) {
        macroField.at(j) = controlNode->giveDofWithID(dofIdArray.at(j) )->giveUnknown(VM_Total, tStep);
    }

    //Get unit cell size
    const auto unitCellSize = controlNode->giveCoordinates();

    for ( int i = 1, n = primaryVarsToExport.giveSize(); i <= n; i++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);

        for ( int inode = 1; inode <= mapL2G.giveSize(); inode++ ) {
            if ( inode <= nnodes && mapL2G.at(inode) <= nnodes && mapL2G.at(inode) != 0 ) { //no special treatment for master nodes
                DofManager *dman = d->giveNode(mapL2G.at(inode) );

                this->getNodalVariableFromPrimaryField(valueArray, dman, tStep, type, region, smoother);
                vtkPiece.setPrimaryVarInNode(type, inode, std::move(valueArray) );
            } else { //special treatment for image nodes
                //find the periodic node, enough to find the first occurrence
                int pos = 0;
                if ( mapL2G.at(inode) != 0 ) {
                    pos = regionToUniqueMap.findFirstIndexOf(mapL2G.at(inode) );
                }
                if ( pos ) {
                    DofManager *dman = d->giveNode(periodicMap.at(pos) );
                    IntArray switches;
                    giveSwitches(switches, locationMap.at(pos) );
                    //get the master unknown
                    FloatArray helpArray;
                    this->getNodalVariableFromPrimaryField(helpArray, dman, tStep, type, region, smoother);
                    //recalculate the image unknown
                    if ( type == DisplacementVector ) {
                        if ( dofIdArray.giveSize() == 3 && dofIdArray.at(1) == E_xx && dofIdArray.at(2) == E_yy && dofIdArray.at(3) == G_xy ) { //Macroscale: 2D solid using Voigt notation
                            valueArray.resize(helpArray.giveSize() );
                            valueArray.at(1) = helpArray.at(1) + unitCellSize.at(1) * switches.at(1) * macroField.at(1);
                            valueArray.at(2) = helpArray.at(2) + unitCellSize.at(2) * switches.at(2) * macroField.at(2) + unitCellSize.at(1) * switches.at(1) * macroField.at(3);
                            valueArray.at(3) = helpArray.at(3);
                        } else if ( dofIdArray.giveSize() == 1 && dofIdArray.at(1) == E_xx ) { //Macroscale: 3D truss. 1 DOF: EXX
                            valueArray.resize(helpArray.giveSize() );
                            valueArray.at(1) = helpArray.at(1) + unitCellSize.at(1) * switches.at(1) * macroField.at(1);
                            valueArray.at(2) = helpArray.at(2);
                            valueArray.at(3) = helpArray.at(3);
                        } else if ( dofIdArray.giveSize() == 6 && dofIdArray.at(1) == E_xx && dofIdArray.at(2) == E_yy && dofIdArray.at(3) == E_zz && dofIdArray.at(4) == G_yz && dofIdArray.at(5) == G_xz && dofIdArray.at(6) == G_xy ) { //Macroscale: 3D solid using Voigt notation
                            valueArray.resize(helpArray.giveSize() );
                            valueArray.at(1) = helpArray.at(1) + unitCellSize.at(1) * switches.at(1) * macroField.at(1) +
                                               unitCellSize.at(3) * switches.at(3) * macroField.at(5) + unitCellSize.at(2) * switches.at(2) * macroField.at(6);
                            valueArray.at(2) = helpArray.at(2) + unitCellSize.at(2) * switches.at(2) * macroField.at(2) +
                                               unitCellSize.at(3) * switches.at(3) * macroField.at(4);
                            valueArray.at(3) = helpArray.at(3) + unitCellSize.at(3) * switches.at(3) * macroField.at(3);
                        } else {
                            OOFEM_ERROR("Unknown periodic element type\n");
                        }
                    }
                } else {
                    valueArray.resize(3);
                }
                vtkPiece.setPrimaryVarInNode(type, inode, std::move(valueArray) );
            }
        }
    }
}


void
VTKXMLLatticeExportModule::exportIntVars(ExportRegion &vtkPiece, Set& region, IntArray& internalVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    int nnodes = d->giveNumberOfDofManagers();
    InternalStateType isType;
    FloatArray answer;

    smoother.clear(); // Makes sure smoother is up-to-date with potentially new mesh.
    //const IntArray& mapG2L = vtkPiece.getMapG2L();
    const IntArray& mapL2G = vtkPiece.getMapL2G();
    // Export of Internal State Type fields
    vtkPiece.setNumberOfInternalVarsToExport(internalVarsToExport, mapL2G.giveSize() );
    for ( int field = 1; field <= internalVarsToExport.giveSize(); field++ ) {
        isType = ( InternalStateType ) internalVarsToExport.at(field);

        for ( int nodeNum = 1; nodeNum <= mapL2G.giveSize(); nodeNum++ ) {
            if ( nodeNum <= nnodes && mapL2G.at(nodeNum) <= nnodes && mapL2G.at(nodeNum) != 0 ) { //no special treatment for master nodes
                Node *node = d->giveNode(mapL2G.at(nodeNum) );
                this->getNodalVariableFromIS(answer, node, tStep, isType, region, smoother);
                vtkPiece.setInternalVarInNode(isType, nodeNum, answer);
            } else { //special treatment for image nodes
                //find the periodic node, enough to find the first occurrence
                int pos = 0;
                if ( mapL2G.at(nodeNum) != 0 ) {
                    pos = regionToUniqueMap.findFirstIndexOf(mapL2G.at(nodeNum) );
                }

                if ( pos ) {
                    Node *node = d->giveNode(periodicMap.at(pos) );
                    this->getNodalVariableFromIS(answer, node, tStep, isType, region, smoother);
                    vtkPiece.setInternalVarInNode(isType, nodeNum, answer);
                } else { //fill with zeroes
                    InternalStateValueType valType = giveInternalStateValueType(isType);
                    int ncomponents = giveInternalStateTypeSize(valType);
                    answer.resize(ncomponents);

                    if ( isType == IST_BeamForceMomentTensor ) {
                        answer.resize(6);
                    }

                    answer.zero();
                    vtkPiece.setInternalVarInNode(isType, nodeNum, answer);
                }
            }
        }
    }
}


bool
VTKXMLLatticeExportModule::writeVTKPieceCross(ExportRegion &vtkPieceCross, TimeStep *tStep)
{
    if ( !vtkPieceCross.giveNumberOfCells() ) {
        return false;
    }

    // Write output: node coords
    int numNodes = vtkPieceCross.giveNumberOfNodes();
    int numEl = vtkPieceCross.giveNumberOfCells();
    FloatArray coords;

    this->fileStreamCross << "<Piece NumberOfPoints=\"" << numNodes << "\" NumberOfCells=\"" << numEl << "\">\n";
    this->fileStreamCross << "<Points>\n <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"> ";

    for ( int inode = 1; inode <= numNodes; inode++ ) {
        coords = vtkPieceCross.giveNodeCoords(inode);
        ///@todo move this below into setNodeCoords since it should alwas be 3 components anyway
        for ( int i = 1; i <= coords.giveSize(); i++ ) {
            this->fileStreamCross << scientific << coords.at(i) << " ";
        }

        for ( int i = coords.giveSize() + 1; i <= 3; i++ ) {
            this->fileStreamCross << scientific << 0.0 << " ";
        }
    }

    this->fileStreamCross << "</DataArray>\n</Points>\n";
    this->fileStreamCross << "<Cells>\n";
    this->fileStreamCross << " <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"> ";

    IntArray cellNodes;
    for ( int ielem = 1; ielem <= numEl; ielem++ ) {
        cellNodes = vtkPieceCross.giveCellConnectivity(ielem);

        for ( int i = 1; i <= cellNodes.giveSize(); i++ ) {
            this->fileStreamCross << cellNodes.at(i) - 1 << " ";
        }
        this->fileStreamCross << " ";
    }

    this->fileStreamCross << "</DataArray>\n";

    // output the offsets (index of individual element data in connectivity array)
    this->fileStreamCross << " <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"> ";

    for ( int ielem = 1; ielem <= numEl; ielem++ ) {
        this->fileStreamCross << vtkPieceCross.giveCellOffset(ielem) << " ";
    }

    this->fileStreamCross << "</DataArray>\n";


    // output cell (element) types
    this->fileStreamCross << " <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"> ";
    for ( int ielem = 1; ielem <= numEl; ielem++ ) {
        this->fileStreamCross << vtkPieceCross.giveCellType(ielem) << " ";
    }

    this->fileStreamCross << "</DataArray>\n";
    this->fileStreamCross << "</Cells>\n";


    ///@todo giveDataHeaders is currently not updated wrt the new structure -> no file names in headers /JB
    std::string pointHeader, cellHeader;
    this->giveDataHeaders(pointHeader, cellHeader);

    this->fileStreamCross << pointHeader.c_str();

    // Displacement field so ParaView's Warp by Vector shows the deformed bodies.
    if ( displacementCross.giveNumberOfRows() == numNodes ) {
        this->fileStreamCross << " <DataArray type=\"Float64\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\"> ";
        for ( int inode = 1; inode <= numNodes; ++inode ) {
            this->fileStreamCross << scientific << displacementCross.at(inode, 1) << " "
                                  << scientific << displacementCross.at(inode, 2) << " "
                                  << scientific << displacementCross.at(inode, 3) << " ";
        }
        this->fileStreamCross << "</DataArray>\n";
    }

    this->fileStreamCross << "</PointData>\n";
    this->fileStreamCross << cellHeader.c_str();

    this->writeCellVarsCross(vtkPieceCross);

    this->fileStreamCross << "</CellData>\n";
    this->fileStreamCross << "</Piece>\n";

    vtkPieceCross.clear();
    return true;
}


void
VTKXMLLatticeExportModule::writePrimaryVarsCross(ExportRegion &vtkPiece)
{
    for ( int i = 1; i <= primaryVarsToExport.giveSize(); i++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);
        InternalStateValueType valType = giveInternalStateValueType(type);
        int ncomponents = giveInternalStateTypeSize(valType);
        ( void ) ncomponents; //silence the warning
        int numNodes = vtkPiece.giveNumberOfNodes();
        const char *name = __UnknownTypeToString(type);
        ( void ) name; //silence the warning

        if (!this->fileStreamCross.is_open()) {
            OOFEM_ERROR("fileStreamCross is not open");
        }

        this->fileStreamCross << " <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"" << ncomponents << "\" format=\"ascii\"> ";
        for ( int inode = 1; inode <= numNodes; inode++ ) {
            FloatArray &valueArray = vtkPiece.givePrimaryVarInNode(type, inode);

            if (valueArray.giveSize() == 0) {
                OOFEM_WARNING("Empty valueArray at node %d", inode);
                continue;
            }
	    for ( int k = 1; k <= valueArray.giveSize(); k++ ) {
	      this->fileStreamCross << scientific << valueArray.at(k) << " ";
	    }
        }
        this->fileStreamCross << "</DataArray>\n";
    }
}


void
VTKXMLLatticeExportModule::exportCellVarsCross(ExportRegion &vtkPiece, Set &region, IntArray &cellVarsToExport, TimeStep *tStep)
{
    // Layout per element (from setupVTKPieceCross):
    //   * Polygon: 3 copies (A, B, midline) × nStrips cells, all 3*nStrips cells contiguous.
    //   * Fallback: 1 cell.
    // Per copy, strips are emitted s=0..nStrips-1 in IP order: for hybrid shells, midline
    // strip s carries layer-IP s's data; A and B strips get zeros.
    Domain *d = emodel->giveDomain(1);
    const IntArray &elems = region.giveElementList();
    const int nElem = elems.giveSize();
    const int nCells = vtkPiece.giveNumberOfCells();

    vtkPiece.setNumberOfCellVarsToExport(cellVarsToExport, nCells);

    for ( int field = 1; field <= cellVarsToExport.giveSize(); ++field ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(field);
        InternalStateValueType valType = giveInternalStateValueType(type);
        int ncomponents = giveInternalStateTypeSize(valType);
        FloatArray zeroVec(ncomponents);
        zeroVec.zero();
        FloatArray valueArray;

        int cellIdx = 0;
        for ( int ie = 1; ie <= nElem; ++ie ) {
            Element *el = d->giveElement(elems.at(ie));
            if ( cellIdx >= nCells ) break;
            const int kind = ( ie <= elemKindCross.giveSize() ) ? elemKindCross.at(ie) : 0;
            const int ncells_this = ( ie <= cellsPerElemCross.giveSize() ) ? cellsPerElemCross.at(ie) : 0;
            const int nStrips = ( ie <= nStripsPerElemCross.giveSize() ) ? nStripsPerElemCross.at(ie) : 1;
            const int N = ( ie <= vertsPerCellCross.giveSize() ) ? vertsPerCellCross.at(ie) : 0;
            const int capCount = ( ie <= capCountCross.giveSize() ) ? capCountCross.at(ie) : 0;
            if ( ncells_this == 0 ) continue;

            if ( el->giveParallelMode() != Element_local ) {
                for ( int k = 1; k <= ncells_this; ++k ) vtkPiece.setCellVar(type, ++cellIdx, zeroVec);
                continue;
            }

            // Kind 2 (line bond): single cell with element data.
            if ( kind == 2 ) {
                this->getCellVariableFromIS(valueArray, el, type, tStep);
                vtkPiece.setCellVar(type, ++cellIdx, valueArray);
                continue;
            }

            // Kind 0 + kind 1 share the "A zeros / B zeros / midline data" pattern for the first 3 cells.
            // Kind 0 has nStrips midline cells (per-IP for hybrid shell).
            // Kind 1 has 1 midline cell + 2N lateral quads (element data) + capCount caps (element data).

            // A copies: zeros.
            for ( int s = 0; s < nStrips; ++s ) vtkPiece.setCellVar(type, ++cellIdx, zeroVec);
            // B copies: zeros.
            for ( int s = 0; s < nStrips; ++s ) vtkPiece.setCellVar(type, ++cellIdx, zeroVec);
            // Midline copies: per-IP value when nStrips>1 (kind 0 hybrid), else element-level.
            if ( nStrips == 1 ) {
                this->getCellVariableFromIS(valueArray, el, type, tStep);
                vtkPiece.setCellVar(type, ++cellIdx, valueArray);
            } else {
                IntegrationRule *iRule = el->giveDefaultIntegrationRulePtr();
                const int nGp = iRule ? iRule->giveNumberOfIntegrationPoints() : 0;
                for ( int s = 0; s < nStrips; ++s ) {
                    if ( s < nGp ) {
                        GaussPoint *gp = iRule->getIntegrationPoint(s);
                        valueArray.resize(ncomponents); valueArray.zero();
                        int ok = el->giveIPValue(valueArray, gp, type, tStep);
                        if ( !ok || valueArray.giveSize() != ncomponents ) valueArray = zeroVec;
                        vtkPiece.setCellVar(type, ++cellIdx, valueArray);
                    } else {
                        vtkPiece.setCellVar(type, ++cellIdx, zeroVec);
                    }
                }
            }

            // Kind 0 closure cells (TOP/BOT + SIDE for boundary shells): zero data.
            // These were emitted per shell element in setupVTKPieceCross; cellsPerElemCross
            // already counts them, so skip the remainder of ncells_this here, otherwise the
            // NEXT element's data would land in the current element's closure slots.
            if ( kind == 0 ) {
                const int nClosure = ncells_this - 3 * nStrips;
                for ( int k = 0; k < nClosure; ++k ) vtkPiece.setCellVar(type, ++cellIdx, zeroVec);
            }

            // Kind 1: lateral quads (2N) + caps (0..2), all carrying element-level data.
            if ( kind == 1 ) {
                this->getCellVariableFromIS(valueArray, el, type, tStep);
                for ( int k = 0; k < 2 * N + capCount; ++k ) {
                    vtkPiece.setCellVar(type, ++cellIdx, valueArray);
                }
            }
        }

        // Trailing safety: any leftover cells (none expected) get zero data.
        while ( cellIdx < nCells ) vtkPiece.setCellVar(type, ++cellIdx, zeroVec);
    }
}


void
VTKXMLLatticeExportModule::writeCellVarsCross(ExportRegion &vtkPiece)
{
    FloatArray valueArray;
    int numCells = vtkPiece.giveNumberOfCells();
    for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);
        InternalStateValueType valType = giveInternalStateValueType(type);
        int ncomponents = giveInternalStateTypeSize(valType);
        const char *name = __InternalStateTypeToString(type);
        ( void ) name; //silence the warning

        this->fileStreamCross << " <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"" << ncomponents << "\" format=\"ascii\"> ";
        valueArray.resize(ncomponents);
        for ( int ielem = 1; ielem <= numCells; ielem++ ) {
            valueArray = vtkPiece.giveCellVar(type, ielem);

            if (valueArray.giveSize() == 0){
                OOFEM_WARNING("Empty valueArray at node %d", ielem);
                continue;
            }

            for ( int k = 1; k <= valueArray.giveSize(); k++ ) {
                this->fileStreamCross << valueArray.at(k) << " ";
            }
        }
        this->fileStreamCross << "</DataArray>\n";
    }

    // polygonRole: 0 = A-side, 1 = B-side, 2 = midline, -2 = line (bond), -1 = point fallback.
    if ( polygonRoleCross.giveSize() == numCells ) {
        this->fileStreamCross << " <DataArray type=\"Int32\" Name=\"polygonRole\" NumberOfComponents=\"1\" format=\"ascii\"> ";
        for ( int ielem = 1; ielem <= numCells; ielem++ ) {
            this->fileStreamCross << polygonRoleCross.at(ielem) << " ";
        }
        this->fileStreamCross << "</DataArray>\n";
    }
    // syntheticGeometry: 1 = polygon synthesised from sectional properties or line for bond, 0 = real polygon.
    if ( syntheticCross.giveSize() == numCells ) {
        this->fileStreamCross << " <DataArray type=\"Int32\" Name=\"syntheticGeometry\" NumberOfComponents=\"1\" format=\"ascii\"> ";
        for ( int ielem = 1; ielem <= numCells; ielem++ ) {
            this->fileStreamCross << syntheticCross.at(ielem) << " ";
        }
        this->fileStreamCross << "</DataArray>\n";
    }
}
} // end namespace oofem
