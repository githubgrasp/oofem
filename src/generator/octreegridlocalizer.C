/* $Header: /home/cvs/bp/oofem/oofemlib/src/octreelocalizer.C,v 1.16.4.1 2004/04/05 15:19:43 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

#include "generatorerror.h"
#include "octreegridlocalizer.h"
#include "element.h"
#include "grid.h"
#include "integrationrule.h"
#include "vertex.h"
#include "node.h"
#include "mathfem.h"
#ifndef __MAKEDEPEND
 #include <time.h>
 #include <iostream>
#endif

oofemOctantRec :: oofemOctantRec(OctreeGridLocalizer *loc, oofemOctantRec *parent, oofem::FloatArray &origin, double size)
{
    int i, j, k;

    this->localizer = loc;
    this->parent = parent;
    this->origin = origin;
    this->size   = size;
    this->nodeList    = NULL;
    this->elementList = NULL;

    for ( i = 0; i <= 1; i++ ) {
        for ( j = 0; j <= 1; j++ ) {
            for ( k = 0; k <= 1; k++ ) {
                this->child [ i ] [ j ] [ k ] = NULL;
            }
        }
    }
}

oofemOctantRec :: ~oofemOctantRec()
{
    // destructor - delete all receivers childs recursivelly
    // delete childs first, then itself
    int i, j, k;
    for ( i = 0; i <= 1; i++ ) {
        for ( j = 0; j <= 1; j++ ) {
            for ( k = 0; k <= 1; k++ ) {
                if ( this->child [ i ] [ j ] [ k ] ) {
                    delete this->child [ i ] [ j ] [ k ];
                }
            }
        }
    }

    if ( nodeList ) {
        delete nodeList;
    }

    if ( elementList ) {
        delete elementList;
    }
}

std :: list< int > *
oofemOctantRec :: giveNodeList() {
    // std::cout<<"kab00m "<< '\n';
    if ( nodeList == NULL ) {
        nodeList = new std :: list< int >;
    }

    return nodeList;
}

oofemOctantRec *
oofemOctantRec :: giveChild(int xi, int yi, int zi)
{
    if ( ( xi >= 0 ) && ( xi < 2 ) && ( yi >= 0 ) && ( yi < 2 ) && ( zi >= 0 ) && ( zi < 2 ) ) {
        return this->child [ xi ] [ yi ] [ zi ];
    } else {
      generator::error4("oofemOctantRec::giveChild invalid child index (%d,%d,%d)", xi, yi, zi);
    }

    return NULL;
}


int
oofemOctantRec :: containsPoint(const oofem::FloatArray &coords)
{
    int i;

    for ( i = 1; i <= coords.giveSize(); i++ ) {
        if ( localizer->giveOctreeMaskValue(i) ) {
            if ( coords.at(i) < this->origin.at(i) ) {
                return 0;
            }

            if ( coords.at(i) > ( this->origin.at(i) + this->size ) ) {
                return 0;
            }
        }
    }

    return 1;
}


int
oofemOctantRec :: giveChildContainingPoint(oofemOctantRec **child, const oofem::FloatArray &coords)
{
    int i;
    oofem::IntArray ind(3);

    if ( this->containsPoint(coords) ) {
        if ( this->isTerminalOctant() ) {
            * child = NULL;
            return -1;
        }

        for ( i = 1; i <= coords.giveSize(); i++ ) {
            if ( localizer->giveOctreeMaskValue(i) && ( coords.at(i) > ( this->origin.at(i) + this->size / 2. ) ) ) {
                ind.at(i) = 1;
            } else {
                ind.at(i) = 0;
            }
        }

        * child = this->child [ ind.at(1) ] [ ind.at(2) ] [ ind.at(3) ];
        return 1;
    } else {
        * child = NULL;
        return -2;
    }
}


int
oofemOctantRec :: isTerminalOctant() {
    // This implemetation is relying on fact, that child [0][0][0]
    // is created even for all degenerated trees
    if ( this->child [ 0 ] [ 0 ] [ 0 ] ) {
        return 0;
    }

    return 1;
}

int
oofemOctantRec :: divideLocally(int level, const oofem::IntArray &octantMask)
{
    int i, j, k, result = 1;

    if ( this->isTerminalOctant() ) {
        // create corresponding child octants
        int i, j, k;
	oofem::FloatArray childOrigin(3);

        for ( i = 0; i <= octantMask.at(1); i++ ) {
            for ( j = 0; j <= octantMask.at(2); j++ ) {
                for ( k = 0; k <= octantMask.at(3); k++ ) {
                    childOrigin.at(1) = this->origin.at(1) + i * ( this->size / 2. );
                    childOrigin.at(2) = this->origin.at(2) + j * ( this->size / 2. );
                    childOrigin.at(3) = this->origin.at(3) + k * ( this->size / 2. );

                    this->child [ i ] [ j ] [ k ] = new oofemOctantRec(localizer, this, childOrigin, this->size / 2.0);
                }
            }
        }
    }

    int newLevel = level - 1;
    if ( newLevel > 0 ) {
        // propagate message to childs recursivelly with level param decreased
        for ( i = 0; i <= octantMask.at(1); i++ ) {
            for ( j = 0; j <= octantMask.at(2); j++ ) {
                for ( k = 0; k <= octantMask.at(3); k++ ) {
                    if ( this->child [ i ] [ j ] [ k ] ) {
                        result &= this->child [ i ] [ j ] [ k ]->divideLocally(newLevel, octantMask);
                    }
                }
            }
        }
    }

    return result;
}


oofemOctantRec :: boundingBoxStatus
oofemOctantRec :: testBoundingBox(const oofem::FloatArray &coords, double radius)
{
    int i, test = 1, nsd, size = coords.giveSize();
    double dist;

    nsd = localizer->giveOctreeMaskValue(1) + localizer->giveOctreeMaskValue(2) + localizer->giveOctreeMaskValue(3);

    oofem::FloatArray cellCenter = this->origin;
    // first simple test to exclude too far bboxes, but detect preciselly the corner bboxes
    double cellRadius = sqrt( nsd * ( this->size / 2.0 ) * ( this->size / 2. ) );
    for ( i = 1; i < 4; i++ ) {
        cellCenter.at(i) += this->size / 2.0;
    }

    for ( i = 1, dist = 0.0; i <= size; i++ ) {
        if ( localizer->giveOctreeMaskValue(i) ) {
            dist += ( cellCenter.at(i) - coords.at(i) ) * ( cellCenter.at(i) - coords.at(i) );
        }
    }

    dist = sqrt(dist);
    if ( dist > ( cellRadius + radius ) ) {
        return BBS_OUTSIDECELL;
    }

    int centerInside = this->containsPoint(coords);
    if ( centerInside ) { // test if whole bbox inside
        for ( i = 1; i <= size; i++ ) {
            if ( localizer->giveOctreeMaskValue(i) ) {
                if ( ( this->origin.at(i) > ( coords.at(i) - radius ) )  ||  ( ( this->origin.at(i) + this->size ) < ( coords.at(i) + radius ) ) ) {
                    test = 0;
                }
            }
        }

        if ( test ) {
            return BBS_INSIDECELL;
        } else {
            return BBS_CONTAINSCELL;
        }

    } else { // BBox center not inside cell, but may hit cell
        // test if bounding sphere hits boundary surface
        int inBounds = 1;
        for ( i = 1; i <= size; i++ ) {
            if ( localizer->giveOctreeMaskValue(i) && ( fabs( cellCenter.at(i) - coords.at(i) ) > ( this->size / 2. + radius ) ) ) {
                inBounds = 0;
            }
        }

        if ( inBounds ) {
            return BBS_CONTAINSCELL;
        }
    }

    return BBS_OUTSIDECELL;
}


oofemOctantRec *
OctreeGridLocalizer :: findTerminalContaining(oofemOctantRec *startCell, const oofem::FloatArray &coords) {
    int result;
    oofemOctantRec *currCell = startCell;
    if ( startCell->containsPoint(coords) ) {
        // found terminal octant containing node
        while ( !currCell->isTerminalOctant() ) {
            result = currCell->giveChildContainingPoint(& currCell, coords);
            if ( result == -2 ) {
	      generator::error("findTerminalContaining: internal error - octree inconsistency");
            }
        }

        ;

        return currCell;
    } else {
        return NULL;
    }
}


int
OctreeGridLocalizer :: buildOctreeDataStructure()
{
    int i, j, init = 1, nnode = this->grid->giveNumberOfInputVertices();
    double rootSize, resolutionLimit;
    oofem::FloatArray minc(3), maxc(3), * coords;
    Vertex *dman;

    // test if tree already built
    if ( rootCell ) {
        return 1;
    }

    // measure time consumed by octree build phase
    clock_t sc = :: clock();
    clock_t ec;

    // first determine grid extends (bounding box), and check for degenerated grid type
    for ( i = 1; i <= nnode; i++ ) {
        dman = grid->giveInputVertex(i);
        coords = ( ( ( Vertex * ) dman )->giveCoordinates() );
        if ( init ) {
            init = 0;
            for ( j = 1; j <= coords->giveSize(); j++ ) {
                minc.at(j) = maxc.at(j) = coords->at(j);
            }
        } else {
            for ( j = 1; j <= coords->giveSize(); j++ ) {
                if ( coords->at(j) < minc.at(j) ) {
                    minc.at(j) = coords->at(j);
                }

                if ( coords->at(j) > maxc.at(j) ) {
                    maxc.at(j) = coords->at(j);
                }
            }
        }
        //        }
    } // end loop over input nodes

    //Extend cell in all three directions by plus/minus its size

    oofem::FloatArray distances(3);
    for ( int i = 1; i <= 3; i++ ) {
        distances.at(i) = maxc.at(i) - minc.at(i);
        maxc.at(i) += distances.at(i);
        minc.at(i) -= distances.at(i);
    }

    // determine root size
    rootSize = 0.0;
    for ( i = 1; i <= 3; i++ ) {
      rootSize = 1.000001 * std::max( rootSize, maxc.at(i) - minc.at(i) );
    }

    // check for degenerated grid
    resolutionLimit = std::min(1.e-3, rootSize / 1.e6);
    for ( i = 1; i <= 3; i++ ) {
        if ( ( maxc.at(i) - minc.at(i) ) > resolutionLimit ) {
            this->octreeMask.at(i) = 1;
        } else {
            this->octreeMask.at(i) = 0;
        }
    }

    // Create root Octant
    this->rootCell = new oofemOctantRec(this, NULL, minc, rootSize);

    // Build octree tree
    if ( nnode > OCTREE_MAX_NODES_LIMIT ) {
        this->rootCell->divideLocally(1, this->octreeMask);
    }

    ec = :: clock();

    // compute max. tree depth
    int treeDepth = 0;
    this->giveMaxTreeDepthFrom(this->rootCell, treeDepth);
    // compute processor time used by the program
    long nsec = ( ec - sc ) / CLOCKS_PER_SEC;
    printf("Octree init [depth %d in %lds]\n", treeDepth, nsec);

    return 1;
}


void
OctreeGridLocalizer :: insertSequentialNode(int nodeNum, const oofem::FloatArray &coords)
{
    this->insertNodeIntoOctree(this->rootCell, nodeNum, coords);

    return;
}



int
OctreeGridLocalizer :: insertNodeIntoOctree(oofemOctantRec *rootCell, int nodeNum, const oofem::FloatArray &coords)
{
    int nCellItems, cellDepth, result;
    oofemOctantRec *currCell;
    nodeContainerType *cellNodeList;
    nodeContainerType :: const_iterator pos;
    Vertex *dman;
    oofem::FloatArray *nodeCoords;

    // found terminal octant containing node
    currCell = this->findTerminalContaining(rootCell, coords);
    // request cell node list
    cellNodeList = currCell->giveNodeList();
    nCellItems = cellNodeList->size();
    cellDepth  = this->giveCellDepth(currCell);
    // check for refinement criteria
    // should also include max refinement level criteria
    //
    if ( ( nCellItems > OCTREE_MAX_NODES_LIMIT ) && ( cellDepth <= OCTREE_MAX_DEPTH ) ) {
        // refine tree one level
        currCell->divideLocally(1, this->octreeMask);
        // propagate all nodes already assigned to currCell to childs
        for ( pos = cellNodeList->begin(); pos != cellNodeList->end(); ++pos ) {
            dman = grid->giveVertex(* pos);
            nodeCoords = ( ( ( Vertex * ) dman )->giveCoordinates() );
            this->insertNodeIntoOctree(currCell, * pos, * nodeCoords);
        }

        // remove node list at relative root cell
        currCell->deleteNodeList();
        // now the new node record is saved simply to one of created children.
        // Generally, one has to check if child is candidate for further refinement
        // (case, when all parent nodes are belonging to particular child), refine it
        // and check childs again. This can be implemented by recursive procedure, but
        // probability of this case is relatively small.
        // Current implementation simply insert node to child created in first refinement,
        // which can lead ne node violation of refinement criteria.
        // But implementation is simple and effective.
        // Please note, that the next insertion in particular cell cause the refinement.

        // find child containing new node
        result = currCell->giveChildContainingPoint(& currCell, coords);
        if ( result != 1 ) {
	  generator::error("insertNodeIntoOctree: internal error - octree inconsistency");
        }
    }

    currCell->addNode(nodeNum);
    return 1;
}

int
OctreeGridLocalizer :: checkNodesWithinBox(const oofem::FloatArray &coords, const double radius)
{
    nodeContainerType nodeSet;
    giveAllNodesWithinBox(nodeSet, coords, radius);
    if ( nodeSet.empty() ) {
        return 0;
    } else   {
        return 1;
    }
}


void
OctreeGridLocalizer::giveAllNodesWithinBox(nodeContainerType& nodeSet,
                                           const oofem::FloatArray& coords,
                                           const double radius)
{
    this->init();

    // Start from the terminal cell containing the center (or root if outside)
    oofemOctantRec* curr = this->findTerminalContaining(rootCell, coords);
    if (!curr) curr = rootCell;

    // Ascend until the bounding sphere is fully inside the cell
    if (curr != rootCell) {
        while (curr->testBoundingBox(coords, radius) != oofemOctantRec::BBS_INSIDECELL) {
            curr = curr->giveParent();
            if (curr == rootCell) break;
        }
    }

    // Now traverse from 'curr' downward and collect nodes
    this->giveNodesWithinBox(nodeSet, curr, coords, radius);
}




void
OctreeGridLocalizer::giveNodesWithinBox(nodeContainerType& out,
                                        oofemOctantRec* startCell,
                                        const oofem::FloatArray& coords,
                                        const double radius)
{
    // Precompute radius^2
    const double r2 = radius * radius;

    // Simple explicit stack (vector is way faster than recursion here)
    struct Frame { oofemOctantRec* cell; };
    std::vector<Frame> stack;
    stack.reserve(64);
    stack.push_back({ startCell });

    while (!stack.empty()) {
        oofemOctantRec* cell = stack.back().cell;
        stack.pop_back();

        // Cheap prune using your existing test (keeps behavior identical)
        auto hit = cell->testBoundingBox(coords, radius);
        if (hit == oofemOctantRec::BBS_OUTSIDECELL) continue;

        if (cell->isTerminalOctant()) {
            // Terminal: test nodes
            nodeContainerType* cellNodes = cell->giveNodeList();
            if (cellNodes && !cellNodes->empty()) {
                // Compare squared distances; avoid sqrt() in FloatArray::distance
                for (int nid : *cellNodes) {
                    const oofem::FloatArray* nc = grid->giveVertex(nid)->giveCoordinates();

                    double d2 = 0.0;
                    const int D = nc->giveSize(); // usually 3
                    for (int i = 1; i <= D; ++i) {
                        if (!this->giveOctreeMaskValue(i)) continue; // skip degenerate axes
                        const double d = nc->at(i) - coords.at(i);
                        d2 += d * d;
                    }

                    if (d2 <= r2) {
                        out.push_back(nid);
                    }
                }
            }
        } else {
            // Push existing children
            for (int i = 0; i <= octreeMask.at(1); ++i)
            for (int j = 0; j <= octreeMask.at(2); ++j)
            for (int k = 0; k <= octreeMask.at(3); ++k) {
                if (auto* ch = cell->giveChild(i, j, k)) {
                    stack.push_back({ ch });
                }
            }
        }
    }
}

int
OctreeGridLocalizer :: giveCellDepth(oofemOctantRec *cell)
{
    return ( int ) ( log( this->rootCell->giveSize() / cell->giveSize() ) / M_LN2 );
}

void
OctreeGridLocalizer :: giveMaxTreeDepthFrom(oofemOctantRec *root, int &maxDepth)
{
    int i, j, k, depth = this->giveCellDepth(root);
    maxDepth = std::max(maxDepth, depth);

    for ( i = 0; i <= octreeMask.at(1); i++ ) {
        for ( j = 0; j <= octreeMask.at(2); j++ ) {
            for ( k = 0; k <= octreeMask.at(3); k++ ) {
                if ( root->giveChild(i, j, k) ) {
                    this->giveMaxTreeDepthFrom(root->giveChild(i, j, k), maxDepth);
                }
            }
        }
    }
}


void
OctreeGridLocalizer :: giveListOfTerminalCellsInBoundingBox(std :: list< oofemOctantRec * > &cellList, const oofem::FloatArray &coords,
                                                            const double radius, oofemOctantRec *currentCell)
{
    int i, j, k;
    oofemOctantRec :: boundingBoxStatus BBStatus;

    BBStatus = currentCell->testBoundingBox(coords, radius);
    if ( ( BBStatus == oofemOctantRec :: BBS_INSIDECELL ) || ( BBStatus == oofemOctantRec :: BBS_CONTAINSCELL ) ) {
        if ( currentCell->isTerminalOctant() ) {
            cellList.push_back(currentCell);
        } else {
            for ( i = 0; i <= octreeMask.at(1); i++ ) {
                for ( j = 0; j <= octreeMask.at(2); j++ ) {
                    for ( k = 0; k <= octreeMask.at(3); k++ ) {
                        if ( currentCell->giveChild(i, j, k) ) {
                            this->giveListOfTerminalCellsInBoundingBox( cellList, coords, radius, currentCell->giveChild(i, j, k) );
                        }
                    }
                }
            }
        }
    }
}

int
OctreeGridLocalizer :: init(bool force)
{
    if ( force ) {
        if ( rootCell ) {
            delete rootCell;
        }
        rootCell = NULL;
        elementIPListsInitialized = 0;
    }
    if ( !rootCell ) {
        return this->buildOctreeDataStructure();
    } else {
        return 0;
    }
}
