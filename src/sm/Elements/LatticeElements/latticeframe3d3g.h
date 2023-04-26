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

#ifndef latticeframe3d3g_h
#define latticeframe3d3g_h

#include "latticeframe3d.h"

///@name Input fields for LatticeFrame3d3g
//@{
#define _IFT_LatticeFrame3d3g_Name "latticeframe3d3g"

//@}

namespace oofem {
/**
 * This class implements a geometric nonlinear 3-dimensional frame element based on rigid body spring theory presented in Toi 1991 and Toi 1993. It is an extension of a geometric linear 3-dimensional frame element. It belongs to the group of lattice models in OOFEM.
 * Authors: Gumaa Abdelrhim and Peter Grassl
 */

class LatticeFrame3d3g : public LatticeFrame3d
{
protected:

public:
    LatticeFrame3d3g(int n, Domain *);
    virtual ~LatticeFrame3d3g();
    double computeCurrentLength();

    int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    const char *giveInputRecordName() const override { return _IFT_LatticeFrame3d3g_Name; }
    const char *giveClassName() const override { return "latticeframe3d3g"; }

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override;

protected:
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    virtual void  computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;
};
} // end namespace oofem
#endif
