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

#ifndef latticeframe3dnl2_h
#define latticeframe3dnl2_h

#include "latticeframe3d.h"

///@name Input fields for Latticeframe3dnl2
//@{
#define _IFT_LatticeFrame3dNL2_Name "latticeframe3dnl2"

//@}

namespace oofem {
/**
 * This class implements a geometric nonlinear 3-dimensional frame element based on rigid body spring models. It is an extension of a geometric linear 3-dimensional frame element based on rigid body spring theory called LatticeFrame3D presented in Toi 1991 and Toi 1993. It belongs to the group of lattice models in the OOFEM structure, but can be used as a standard 3D beam element.
 * References:
 * Toi, Y. (1991). Shifted integration technique in one‐dimensional plastic collapse analysis using linear and cubic finite elements. International Journal for Numerical Methods in Engineering, 31(8), 1537-1552.
 * Toi, Y., & Isobe, D. (1993). Adaptively shifted integration technique for finite element collapse analysis of framed structures. International Journal for Numerical Methods in Engineering, 36(14), 2323-2339.
 * Authors: Gumaa Abdelrhim and Peter Grassl, 2024
 */

class LatticeFrame3dNL2 : public LatticeFrame3d
{
protected:

public:
    LatticeFrame3dNL2(int n, Domain *);
    virtual ~LatticeFrame3dNL2();

    const char *giveInputRecordName() const override { return _IFT_LatticeFrame3dNL2_Name; }
    const char *giveClassName() const override { return "latticeframe3dnl2"; }

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override;

protected:

  bool computeGtoLStrainRotationMatrix(FloatMatrix &answer);
  bool computeGtoLRotationMatrix(FloatMatrix &) override;
  virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override;
  void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
  virtual void  computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;
};
} // End namespace oofem
#endif
