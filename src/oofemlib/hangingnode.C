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

#include "hangingnode.h"
#include "slavedof.h"
#include "floatarray.h"
#include "intarray.h"
#include "element.h"
#include "feinterpol.h"
#include "spatiallocalizer.h"
#include "classfactory.h"

#include "array"
#include "vector"
#include <cmath>

namespace oofem {
REGISTER_DofManager(HangingNode);

HangingNode :: HangingNode(int n, Domain *aDomain) : Node(n, aDomain)
{
#ifdef __OOFEG
    initialized = false;
#endif
}

void HangingNode :: initializeFrom(InputRecord &ir)
{
    Node :: initializeFrom(ir);
    this->masterElement = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, this->masterElement, _IFT_HangingNode_masterElement);
    this->masterElementRot = -1;
    IR_GIVE_OPTIONAL_FIELD( ir, this->masterElementRot, _IFT_HangingNode_masterElement );
    this->masterRegion = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->masterRegion, _IFT_HangingNode_masterRegion);
}

int HangingNode :: checkConsistency()
{
    if (this->masterElement !=-1){
        return true;
    }
        
    Element *e = this->domain->giveGlobalElement(this->masterElement);
    int result = Node :: checkConsistency();

#if 0
    // Check if master is in same mode
    if ( parallel_mode != DofManager_local ) {
        for ( int i = 1; i <= countOfMasterNodes; i++ ) {
            if ( e->giveNode(i)->giveParallelMode() != parallel_mode ) {
                OOFEM_WARNING("Mismatch in parallel mode of HangingNode and master");
                return false;
            }
        }
    }
#endif

    // Check local coordinate systems
    for ( int i = 1; i <= e->giveNumberOfNodes(); ++i ) {
        if ( !this->hasSameLCS( e->giveNode(i) ) ) {
            OOFEM_WARNING("Different lcs for master/slave nodes.");
            result = false;
        }
    }
    return result;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//2d
std::tuple<int, int, double> HangingNode::findMasterNodesWithAlignedCoordinates( const IntArray &masterNodes, FloatArray localLcoords, Element *e ) const
{
    FEInterpolation *fei = e->giveInterpolation();
    FloatArray localLcoords1, localLcoords2;
    if ( masterNodes.giveSize() < 2 ) {
        OOFEM_ERROR( "Not enough master nodes for interpolation." );
    }

    // Hanging node global coordinates
    FloatArray hangingGlobalCoords = this->giveCoordinates();
    double xh                      = hangingGlobalCoords.at( 1 );
    double yh                      = hangingGlobalCoords.at( 2 );

    int node1 = -1, node2 = -1;
    bool alignedWithX = false; // True if aligned along x-axis, false if y-axis

    // Identify two master nodes aligned with the hanging node (in global coordinates)
    for ( int i = 1; i <= masterNodes.giveSize(); i++ ) {
        FloatArray masterCoords = this->domain->giveNode( masterNodes.at( i ) )->giveCoordinates();
        double xm               = masterCoords.at( 1 );
        double ym               = masterCoords.at( 2 );

        if ( fabs( xm - xh ) < 1e-12 ) { // Same x-coordinate (aligned vertically)
            if ( node1 == -1 ) {
                node1 = masterNodes.at( i );
            } else {
                node2        = masterNodes.at( i );
                alignedWithX = false; // Aligned along y-direction
                break;
            }
        } else if ( fabs( ym - yh ) < 1e-12 ) { // Same y-coordinate (aligned horizontally)
            if ( node1 == -1 ) {
                node1 = masterNodes.at( i );
            } else {
                node2        = masterNodes.at( i );
                alignedWithX = true; // Aligned along x-direction
                break;
            }
        }
    }

    if ( node1 == -1 || node2 == -1 ) {
        OOFEM_ERROR( "Could not find two suitable master nodes for interpolation." );
    }

    // Convert master node global coordinates to local (isoparametric) coordinates
    FloatArray masterLcoords1, masterLcoords2;
    FloatArray masterGlobalCoords1 = this->domain->giveNode( node1 )->giveCoordinates();
    FloatArray masterGlobalCoords2 = this->domain->giveNode( node2 )->giveCoordinates();

    fei->global2local( localLcoords1, masterGlobalCoords1, FEIElementGeometryWrapper( e ) );
    fei->global2local( localLcoords2, masterGlobalCoords2, FEIElementGeometryWrapper( e ) );

    //double ksi1 = localLcoords1.at( 1 );
    //double eta1 = localLcoords1.at( 2 );
    //double ksi2 = localLcoords2.at( 1 );
    //double eta2 = localLcoords2.at( 2 );

    //// Compute dl using local coordinates
    //double dl = alignedWithX ? fabs( ksi2 - ksi1 ) : fabs( eta2 - eta1 );
    double dl = alignedWithX ? fabs( masterGlobalCoords1.at( 1 ) - hangingGlobalCoords.at( 1 ) ) : fabs( masterGlobalCoords1.at( 2 ) - hangingGlobalCoords.at( 2 ) ); 
    //double dl = alignedWithX ? fabs( ksi1 - localLcoords.at( 1 ) ) : fabs( eta1 - localLcoords.at( 2 ) ); 
    return { node1, node2, dl };
}

FloatArray HangingNode::computeMasterContributionForRv( const IntArray &masterNodes, const FloatArray &translationContribution, const FloatArray &lcoords, Element *e ) const
{
    FloatArray localLcoords = lcoords;
    // Find master nodes and parametric distance (dl)
    auto [masterNode1, masterNode2, dl] = findMasterNodesWithAlignedCoordinates( masterNodes, localLcoords, e );

    if ( fabs( dl ) < 1e-12 ) {
        OOFEM_ERROR( "Division by zero detected in rotation contribution computation." );
    }

    int masterSize = masterNodes.giveSize();
    FloatArray contribution;
    contribution.resize( masterSize );
    contribution.zero();

    // Compute rotation contributions using translation contributions in parametric space
    //double theta = (translationContribution.at( masterNode1 ) - translationContribution.at( masterNode2 )) / dl;
    //double theta = ( 1.0 - translationContribution.at( masterNode1 ) ) / dl;
    //double theta2 = translationContribution.at( masterNode2 ) / dl;

    // Assign contributions to the correct master nodes
    for ( int i = 1; i <= masterSize; ++i ) {
        if ( masterNodes.at( i ) == masterNode1 ) {
            contribution.at( i ) = ( 1.0 - translationContribution.at( i ) ) / dl;
            //contribution.at( i ) = theta; // Contribution from masterNode1
        } else if ( masterNodes.at( i ) == masterNode2 ) {
            //contribution.at( i ) = -( 1.0 - translationContribution.at( i ) ) / dl;
            //contribution.at( i ) = -theta; // Contribution from masterNode2
            contribution.at( i ) =  - translationContribution.at( i )  / dl;
        } else {
            contribution.at( i ) = 0.0; // No contribution from other nodes
        }
    }

    return contribution;
}
/////////////////////////////////////////////////////////////////////////////////////////
// 3d
// Compute the Euclidean distance between two 3D points
double HangingNode::computeDistance( const FloatArray &masterGlobalCoords1, const FloatArray &masterGlobalCoords2 )
{
    //FloatArray masterGlobalCoords1 = n1.giveCoordinates();
    //FloatArray masterGlobalCoords2 = n2.giveCoordinates();
    return std::sqrt( std::pow( masterGlobalCoords1.at( 1 ) - masterGlobalCoords2.at( 1 ), 2 ) + std::pow( masterGlobalCoords1.at( 2 ) - masterGlobalCoords2.at( 2 ), 2 ) + std::pow( masterGlobalCoords1.at( 3 ) - masterGlobalCoords2.at( 3 ), 2 ) );
}

double HangingNode::distance( const std::array<double, 3> &a, const std::array<double, 3> &b )
{
    return std::sqrt( std::pow( b[0] - a[0], 2 ) + std::pow( b[1] - a[1], 2 ) + std::pow( b[2] - a[2], 2 ) );
}

int HangingNode::findNormalDirection(const FloatArray& localLcoords1, const FloatArray& localLcoords2, const FloatArray& localLcoords3) {
    // Tolerance to handle floating-point precision
    const double EPSILON = 1e-6;

    // Check if all three points have the same zeta (ζ) coordinate
    if ( std::fabs( localLcoords1[2] - localLcoords2[2] ) < EPSILON && std::fabs( localLcoords2[2] - localLcoords3[2] ) < EPSILON ) {
        return 2; // Normal along ζ → Points lie in ξ-η plane
    }

    // Check if all three points have the same eta (η) coordinate
    if ( std::fabs( localLcoords1[1] - localLcoords2[1] ) < EPSILON && std::fabs( localLcoords2[1] - localLcoords3[1] ) < EPSILON ) {
        return 1; // Normal along η → Points lie in ξ-ζ plane
    }

    // Check if all three points have the same ksi (ξ) coordinate
    if ( std::fabs( localLcoords1[0] - localLcoords2[0] ) < EPSILON && std::fabs( localLcoords2[0] - localLcoords3[0] ) < EPSILON ) {
        return 0; // Normal along ξ → Points lie in η-ζ plane
    }

    return -1; // Not aligned in a single coordinate plane
}


// Function to find the three closest nodes to a hanging node
std::vector<int> HangingNode::findThreeClosestNodes( FloatArray &hangingNodeCoords, const IntArray &hexaNodes )
{
    std::vector<std::pair<double, int> > distances;
    std::array<double, 3> hangCoords = { hangingNodeCoords.at( 1 ), hangingNodeCoords.at( 2 ), hangingNodeCoords.at( 3 ) };
    // Compute distances and store index positions
    for ( size_t i = 0; i < hexaNodes.giveSize(); ++i ) {
        const auto &masternode = hexaNodes[i];
        FloatArray masterGlobalCoords = this->domain->giveNode( masternode )->giveCoordinates();
        //masternode.giveCoordinates(); 
        // Convert FloatArray to std::array<double, 3> if needed
        std::array<double, 3> masterCoords = { masterGlobalCoords.at( 1 ), masterGlobalCoords.at( 2 ), masterGlobalCoords.at( 3 ) };
        double dist                        = distance( hangCoords, masterCoords );        
        distances.emplace_back( dist, i ); // Store distance and index
    }

    // Find the 3 smallest distances using nth_element
    if ( distances.size() < 3 ) {
        throw std::runtime_error( "Not enough nodes to find three closest ones." );
    }

    std::nth_element( distances.begin(), distances.begin() + 3, distances.end() );

    // Return the indices of the three closest nodes
    return { distances[0].second, distances[1].second, distances[2].second };
}

//// Function to create the transformation matrix of the master triangle
//double [3][3] createTransformationTABC( const std::vector<double> &coords )
//{
//    double x1 = coords[0], y1 = coords[1], z1 = coords[2];
//    double x2 = coords[3], y2 = coords[4], z2 = coords[5];
//    double exx = coords[6], exy = coords[7], exz = coords[8];
//
//    double elemLength = std::sqrt( std::pow( x2 - x1, 2 ) + std::pow( y2 - y1, 2 ) + std::pow( z2 - z1, 2 ) );
//
//    double L1 = ( x2 - x1 ) / elemLength;
//    double L2 = ( y2 - y1 ) / elemLength;
//    double L3 = ( z2 - z1 ) / elemLength;
//
//    double DXM = exy * ( z2 - z1 ) - exz * ( y2 - y1 );
//    double DYM = -exx * ( z2 - z1 ) + exz * ( x2 - x1 );
//    double DZM = exx * ( y2 - y1 ) - exy * ( x2 - x1 );
//    double DM  = std::sqrt( DXM * DXM + DYM * DYM + DZM * DZM );
//
//    double M1 = ( exx != 0 ) ? exx / std::abs( exx ) : 0;
//    double M2 = ( exy != 0 ) ? exy / std::abs( exy ) : 0;
//    double M3 = ( exz != 0 ) ? exz / std::abs( exz ) : 0;
//
//    double N1 = (DM != 0) ? DXM / DM : 0;
//    double N2 = ( DM != 0 ) ? DYM / DM : 0 ;
//    double N3 = ( DM != 0 ) ? DZM / DM : 0 ;
//
//    double 
//
//    TABC[0][0] = L1;
//    TABC[0][1] = L2;
//    TABC[0][2] = L3;
//    TABC[1][0] = M1;
//    TABC[1][1] = M2;
//    TABC[1][2] = M3;
//    TABC[2][0] = N1;
//    TABC[2][1] = N2;
//    TABC[2][2] = N3; // Placeholder for N values if needed
//}

FloatArray HangingNode::computeTriangleRotations( DofIDItem id, FloatArray &hangingNodeCoords, const IntArray &hexaNodes, const FloatArray &translationContribution, Element *e )
{
    int numCols = hexaNodes.giveSize(); // Assuming T has 9 columns, adjust if needed
    FloatArray localLcoords1, localLcoords2, localLcoords3, globalLcoords1, globalLcoords2, globalLcoords3;
    FEInterpolation *fei = e->giveInterpolation();

    // Find the three closest nodes to the hanging node
    std::vector<int> nodeSurface = findThreeClosestNodes( hangingNodeCoords, hexaNodes );

    // Extract the triangle nodes based on the closest nodes
    FloatArray masterGlobalCoords1 = this->domain->giveNode( hexaNodes[nodeSurface[0]] )->giveCoordinates();
    FloatArray masterGlobalCoords2 = this->domain->giveNode( hexaNodes[nodeSurface[1]] )->giveCoordinates();
    FloatArray masterGlobalCoords3 = this->domain->giveNode( hexaNodes[nodeSurface[2]] )->giveCoordinates();

    // Compute the edge lengths
    double La1    = computeDistance( masterGlobalCoords1, masterGlobalCoords2 );
    double Lb1    = computeDistance( masterGlobalCoords2, masterGlobalCoords3 );
    double Lc1    = computeDistance( masterGlobalCoords1, masterGlobalCoords3 );
    double S1     = ( La1 + Lb1 + Lc1 ) / 2.0;
    double Omega1 = std::sqrt( S1 * ( S1 - La1 ) * ( S1 - Lb1 ) * ( S1 - Lc1 ) );

    //std::array<double, 6> Array_Axes = { 0, 1, 2, 0, 1, 2 };
    fei->global2local( localLcoords1, masterGlobalCoords1, FEIElementGeometryWrapper( e ) );
    fei->global2local( localLcoords2, masterGlobalCoords2, FEIElementGeometryWrapper( e ) );
    fei->global2local( localLcoords3, masterGlobalCoords3, FEIElementGeometryWrapper( e ) );

    int normal_dir = findNormalDirection( localLcoords1, localLcoords2, localLcoords3 );
 //    FloatArray A = { localLcoords1[Array_Axes[normal_dir]], localLcoords1[Array_Axes[normal_dir+1]], localLcoords1[Array_Axes[normal_dir+2]] };
 //    FloatArray B = { localLcoords2[Array_Axes[normal_dir]], localLcoords2[Array_Axes[normal_dir+1]], localLcoords2[Array_Axes[normal_dir+2]] };
 //    FloatArray C = { localLcoords3[Array_Axes[normal_dir]], localLcoords3[Array_Axes[normal_dir+1]], localLcoords3[Array_Axes[normal_dir+2]] };

 ///*   FloatArray A = { Local_X1[0], Local_Y1[0], Local_Z1[0] };
 //   FloatArray B = { Local_X1[1], Local_Y1[1], Local_Z1[1] };
 //   FloatArray C = { Local_X1[2], Local_Y1[2], Local_Z1[2] };*/
 //   fei->local2global( globalLcoords1, A, FEIElementGeometryWrapper( e ) );
 //   fei->local2global( globalLcoords2, B, FEIElementGeometryWrapper( e ) );
 //   fei->local2global( globalLcoords3, C, FEIElementGeometryWrapper( e ) );
    std::array<double, 3> X1 = { masterGlobalCoords3[0] - masterGlobalCoords2[0], masterGlobalCoords1[0] - masterGlobalCoords3[0], masterGlobalCoords2[0] - masterGlobalCoords1[0] };
    std::array<double, 3> Y1 = { masterGlobalCoords3[1] - masterGlobalCoords2[1], masterGlobalCoords1[1] - masterGlobalCoords3[1], masterGlobalCoords2[1] - masterGlobalCoords1[1] };
    std::array<double, 3> Z1 = { masterGlobalCoords3[2] - masterGlobalCoords2[2], masterGlobalCoords1[2] - masterGlobalCoords3[2], masterGlobalCoords2[2] - masterGlobalCoords1[2] };
    
    std::array<int, 6> Array_NodesContr = { nodeSurface[0], nodeSurface[1], nodeSurface[2], nodeSurface[0], nodeSurface[1], nodeSurface[2] };
    FloatArray T(numCols);
    if ( !std::isnan( Omega1 ) && Omega1 != 0 ) {    
            for ( int i = 0; i < 3; i++ ) {
                if ( normal_dir == 0 ) // Normal in ksi (y) -> Plane is eta-zeta (x-z)
                {
                    if ( id == R_w ) { 
                        T.at( hexaNodes[Array_NodesContr[i]] ) += X1[i] / ( 2.0 * Omega1 ) * ( 1.0 - translationContribution.at( hexaNodes[Array_NodesContr[i]] ) );
                        T.at( hexaNodes[Array_NodesContr[i+1]] ) += X1[i] / ( 2.0 * Omega1 ) * ( - translationContribution.at( hexaNodes[Array_NodesContr[i+1]] ) );
                        T.at( hexaNodes[Array_NodesContr[i + 2]] ) += X1[i] / ( 2.0 * Omega1 ) * ( -translationContribution.at( hexaNodes[Array_NodesContr[i + 2]] ) );
                    } else if ( id == R_v ) {
                        T.at( hexaNodes[Array_NodesContr[i]] ) += Z1[i] / ( 2.0 * Omega1 ) * ( 1.0 - translationContribution.at( hexaNodes[Array_NodesContr[i]] ) );
                        T.at( hexaNodes[Array_NodesContr[i+1]] ) += Z1[i] / ( 2.0 * Omega1 ) * (  - translationContribution.at( hexaNodes[Array_NodesContr[i+1]] ) );
                        T.at( hexaNodes[Array_NodesContr[i + 2]] ) += Z1[i] / ( 2.0 * Omega1 ) * ( -translationContribution.at( hexaNodes[Array_NodesContr[i + 2]] ) );
                    } else if ( id == R_u ) {
                        T.at( hexaNodes[Array_NodesContr[i]] ) += ( - Z1[i] / ( 2.0 * Omega1 ) - X1[i] / ( 2.0 * Omega1 )) * ( 1.0 - translationContribution.at( hexaNodes[Array_NodesContr[i]] ) );
                        T.at( hexaNodes[Array_NodesContr[i+1]] ) += ( -Z1[i] / ( 2.0 * Omega1 ) - X1[i] / ( 2.0 * Omega1 ) ) * (  - translationContribution.at( hexaNodes[Array_NodesContr[i+1]] ) );
                        T.at( hexaNodes[Array_NodesContr[i+2]] )     += ( -Z1[i] / ( 2.0 * Omega1 ) - X1[i] / ( 2.0 * Omega1 ) ) * (  - translationContribution.at( hexaNodes[Array_NodesContr[i+2]] ) );
                    }
                } else if ( normal_dir == 1 ) // Normal in eta (x) -> Plane is ksi-zeta (y-z)
                {
                    if ( id == R_w ) { 
                        T.at( hexaNodes[Array_NodesContr[i]] )     += Z1[i] / ( 2.0 * Omega1 ) * ( 1.0 - translationContribution.at( hexaNodes[Array_NodesContr[i]] ) );
                        T.at( hexaNodes[Array_NodesContr[i + 1]] ) += Z1[i] / ( 2.0 * Omega1 ) * ( -translationContribution.at( hexaNodes[Array_NodesContr[i + 1]] ) );
                        T.at( hexaNodes[Array_NodesContr[i + 2]] ) += Z1[i] / ( 2.0 * Omega1 ) * ( -translationContribution.at( hexaNodes[Array_NodesContr[i + 2]] ) );
                    } else if ( id == R_v ) {
                        T.at( hexaNodes[Array_NodesContr[i]] ) += Y1[i] / ( 2.0 * Omega1 ) * ( 1.0 - translationContribution.at( hexaNodes[Array_NodesContr[i]] ) );
                        T.at( hexaNodes[Array_NodesContr[i+1]] ) += Y1[i] / ( 2.0 * Omega1 ) * (  - translationContribution.at( hexaNodes[Array_NodesContr[i+1]] ) );
                        T.at( hexaNodes[Array_NodesContr[i + 2]] ) += Y1[i] / ( 2.0 * Omega1 ) * ( -translationContribution.at( hexaNodes[Array_NodesContr[i + 2]] ) );
                    } else if ( id == R_u ) {
                        T.at( hexaNodes[Array_NodesContr[i]] ) = ( -Z1[i] / ( 2.0 * Omega1 ) - Y1[i] / ( 2.0 * Omega1 ) ) * ( 1.0 - translationContribution.at( hexaNodes[Array_NodesContr[i]] ) );
                        T.at( hexaNodes[Array_NodesContr[i+1]] ) = ( -Z1[i] / ( 2.0 * Omega1 ) - Y1[i] / ( 2.0 * Omega1 ) ) * ( - translationContribution.at( hexaNodes[Array_NodesContr[i+1]] ) );
                        T.at( hexaNodes[Array_NodesContr[i + 2]] ) = ( -Z1[i] / ( 2.0 * Omega1 ) - Y1[i] / ( 2.0 * Omega1 ) ) * ( -translationContribution.at( hexaNodes[Array_NodesContr[i + 2]] ) );
                    }
                } else if ( normal_dir == 2 ) // Normal in zeta (z) -> Plane is ksi-eta (y-x)
                {
                    if ( id == R_w ) { 
                        T.at( hexaNodes[Array_NodesContr[i]] ) += Y1[i] / ( 2.0 * Omega1 ) * ( 1.0 - translationContribution.at( hexaNodes[Array_NodesContr[i]] ) );
                        T.at( hexaNodes[Array_NodesContr[i + 1]] ) += Y1[i] / ( 2.0 * Omega1 ) * ( -translationContribution.at( hexaNodes[Array_NodesContr[i + 1]] ) );
                        T.at( hexaNodes[Array_NodesContr[i + 2]] ) += Y1[i] / ( 2.0 * Omega1 ) * ( -translationContribution.at( hexaNodes[Array_NodesContr[i + 2]] ) );
                    } else if ( id == R_v ) {
                        T.at( hexaNodes[Array_NodesContr[i]] ) += X1[i] / ( 2.0 * Omega1 ) * ( 1.0 - translationContribution.at( hexaNodes[Array_NodesContr[i]] ) );
                        T.at( hexaNodes[Array_NodesContr[i + 1]] ) += X1[i] / ( 2.0 * Omega1 ) * ( -translationContribution.at( hexaNodes[Array_NodesContr[i + 1]] ) );
                        T.at( hexaNodes[Array_NodesContr[i + 2]] ) += X1[i] / ( 2.0 * Omega1 ) * ( -translationContribution.at( hexaNodes[Array_NodesContr[i + 2]] ) );
                    } else if ( id == R_u ) {
                        T.at( hexaNodes[Array_NodesContr[i]] ) = ( -X1[i] / ( 2.0 * Omega1 ) - Y1[i] / ( 2.0 * Omega1 ) ) * ( 1.0 - translationContribution.at( hexaNodes[Array_NodesContr[i]] ) );
                        T.at( hexaNodes[Array_NodesContr[i+1]] ) = ( -X1[i] / ( 2.0 * Omega1 ) - Y1[i] / ( 2.0 * Omega1 ) ) * ( 1.0 - translationContribution.at( hexaNodes[Array_NodesContr[i+1]] ) );
                        T.at( hexaNodes[Array_NodesContr[i+2]] )     = ( -X1[i] / ( 2.0 * Omega1 ) - Y1[i] / ( 2.0 * Omega1 ) ) * ( 1.0 - translationContribution.at( hexaNodes[Array_NodesContr[i+2]] ) );
                    }
                }
            }
    }
            //multiplyMatrices( TABC1.transpose(), T.block( 3, 0, 3, 3 ) );
            //multiplyMatrices( T.block( 3, 0, 3, 3 ), TABC1 );
    return T;
}


void HangingNode :: postInitialize()
{
    Node :: postInitialize();

    Element *e;
    FEInterpolation *fei;
    FloatArray lcoords, masterContribution;
    FloatArray globalCoords1, globalCoords3; // Store global coordinates of master nodes

#ifdef __OOFEG
    if ( initialized ) {
        return;
    }
    initialized = true;
#endif

    // First check element and interpolation
    if ( masterElement == -1 ) { // Then we find it by taking the closest (probably containing element)
        FloatArray closest;
        SpatialLocalizer *sp = this->domain->giveSpatialLocalizer();
        sp->init();
        // Closest point or containing point? It should be contained, but with numerical errors it might be slightly outside
        // so the closest point is more robust.
        if ( !( e = sp->giveElementClosestToPoint(lcoords, closest, coordinates, this->masterRegion) ) ) {
            OOFEM_ERROR("Couldn't find closest element (automatically).");
        }
        this->masterElement = e->giveNumber();
        #  ifdef DEBUG
        OOFEM_LOG_INFO("Found element %d for hanging node %d\n", e->giveGlobalNumber(), this->giveNumber() );
        #endif
    } else if ( !( e = this->giveDomain()->giveGlobalElement(this->masterElement) ) ) {
        OOFEM_ERROR("Requested element %d doesn't exist.", this->masterElement);
    }
    if ( !( fei = e->giveInterpolation() ) ) {
        OOFEM_ERROR("Requested element %d doesn't have a interpolator.", this->masterElement);
    }

    if ( lcoords.giveSize() == 0 ) { // we don't need to do this again if the spatial localizer was used.
        fei->global2local( lcoords, coordinates, FEIElementGeometryWrapper(e) );
    }

    // Initialize slave dofs (inside check of consistency of receiver and master dof)
    const IntArray &masterNodes = e->giveDofManArray();
    for ( Dof *dof: *this ) {
        SlaveDof *sdof = dynamic_cast< SlaveDof * >(dof);
        if ( sdof ) {
            DofIDItem id = sdof->giveDofID();
            fei          = e->giveInterpolation( id );
            if ( !fei ) {
                OOFEM_ERROR( "Requested interpolation for dof id %d doesn't exist in element %d.",
                    id, this->masterElement );
            }
#if 0 // This won't work (yet), as it requires some more general FEI classes, or something similar.
            if ( fei->hasMultiField() ) {
                FloatMatrix multiContribution;
                IntArray masterDofIDs, masterNodesDup, dofids;
                fei->evalMultiN(multiContribution, dofids, lcoords, FEIElementGeometryWrapper(e), 0.0);
                masterContribution.flatten(multiContribution);
                masterDofIDs.clear();
                for ( int i = 0; i <= multiContribution.giveNumberOfColumns(); ++i ) {
                    masterDofIDs.followedBy(dofids);
                    masterNodesDup.followedBy(masterNodes);
                }
                sdof->initialize(masterNodesDup, & masterDofIDs, masterContribution);
            } else { }
#else
            // Note: There can be more masterNodes than masterContributions, since all the
            // FEI classes are based on that the first nodes correspond to the simpler/linear interpolation.
            // If this assumption is changed in FEIElementGeometryWrapper + friends,
            // masterNode will also need to be modified for each dof accordingly.
            fei->evalN( masterContribution, lcoords, FEIElementGeometryWrapper( e ) );
            // Special handling for RV
            if ( (id == R_v || id == R_w || id == R_u) && (this->domain-> giveDomainType() == _3dMode)) {
                // Compute master contributions for Rv
                 FloatArray translationContribution = masterContribution;
                /*masterContribution                 = this->computeMasterContributionForRv( masterNodes, translationContribution, lcoords, e );*/
                 masterContribution = this->computeTriangleRotations( id, coordinates, masterNodes, translationContribution, e );
            } else if ( ( id == R_v || id == R_w || id == R_u ) && ( this->domain->giveDomainType() == _2dPlaneStressMode || this->domain->giveDomainType() == _PlaneStrainMode ) ) {
                // Compute master contributions for Rv
                FloatArray translationContribution = masterContribution;
                masterContribution                 = this->computeMasterContributionForRv( masterNodes, translationContribution, lcoords, e );
            }
            sdof->initialize( masterNodes, IntArray(), masterContribution );
#endif
     /*   } else if ( (dof->giveDofID() == R_v) && (masterElementRot > 0) ) {
            fei->evalN( masterContribution, lcoords, FEIElementGeometryWrapper( e ) );
            FloatArray translationContribution = masterContribution;
            masterContribution                 = this->computeMasterContributionForRv( masterNodes, translationContribution, lcoords, e );    */        
        }
    }
}
} // end namespace oofem
