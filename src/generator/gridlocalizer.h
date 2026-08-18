/* $Header: /home/cvs/bp/oofem/oofemlib/src/spatiallocalizer.h,v 1.9.4.1 2004/04/05 15:19:43 bp Exp $ */
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


//   *******************************
//   *** CLASS SPATIAL LOCALIZER ***
//   *******************************

#ifndef gridlocalizer_h
#define gridlocalizer_h

#include "gridcomponent.h"
//#include "compiler.h"

#include "interface.h"
#include "logger.h"
#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <set>
 #include <list>
#endif

class Grid;

class oofemOctantRec;

/**
 * The base class for all grid localizers.
 * The basic task is to provide grid information and localization for domain, to which receiver is associated.
 * Typical services include searching the closes node to give position, serching of an element containing given point, etc.
 * If special element algorithms required, these should be included using interface concept.
 */
class GridLocalizer : public GridComponent
{
protected:

public:

    /**
     * Typedefs to introduce the container type for element numbers, returned by some services
     */
    typedef std::set< int >elementContainerType;
    /**
     * Typedefs to introduce the container type for nodal numbers, returned by some services
     */
    typedef std::list< int >nodeContainerType;


    /// Constructor.
    GridLocalizer(int n, Grid *g) : GridComponent(n, g) { }

    /// Destructor.
    virtual ~GridLocalizer() = default;

    /**
     * Returns the list of grid nodes within an axis-aligned box.
     * @param nodeList output container populated with the matching node ids
     * @param coords centre of the query box
     * @param radius half-width of the box (bounding sphere radius in the
     *        default implementation)
     */
    virtual void giveAllNodesWithinBox(nodeContainerType &nodeList, const oofem::FloatArray &coords, const double radius) = 0;

    /// Returns nonzero if any grid node lies within the given box.
    virtual int checkNodesWithinBox(const oofem::FloatArray &coords, const double radius) = 0;

    /// Insert a new node with the given number and coordinates into the
    /// localiser in sequence (no rebalancing).
    virtual void insertSequentialNode(int nodeNum, const oofem::FloatArray &coords) = 0;

    /**
     * Initialise the localiser data structure if not done previously.
     * If `force` is true, the initialisation is enforced (useful when the
     * grid geometry has changed). Returns nonzero on success.
     */
    virtual int init(bool force = false) { return 1; }

protected:
};

#endif // gridlocalizer_h
