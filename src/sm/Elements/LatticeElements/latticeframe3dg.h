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

#ifndef latticeframe3dg_h
#define latticeframe3dg_h

#include "latticeframe3d.h"

///@name Input fields for LatticeFrame3dg
//@{
#define _IFT_LatticeFrame3dg_Name "latticeframe3dg"
#define _IFT_LatticeFrame3dg_refnode "refnode"
#define _IFT_LatticeFrame3dg_refangle "refangle"
#define _IFT_LatticeFrame3dg_zaxis "zaxis"
#define _IFT_LatticeFrame3dg_s "s"
//@}

namespace oofem {
/**
 * This class implements a 3-dimensional frame element based on rigid body spring theory presented in Toi 1991 and Toi 1993. It belongs to the group of lattice models in OOFEM.
 * Authors: Gumaa Abdelrhim and Peter Grassl
 */

class LatticeFrame3dg : public LatticeFrame3d
{
protected:
    int referenceNode;
    FloatArray zaxis;
    double referenceAngle = 0;
    double kappa;
    double length = 0.;
    double iy, iz, ik;
    double area, shearareay, shearareaz;
    double s;
    FloatMatrix localCoordinateSystem;

    FloatArray midPoint, globalCentroid, normal;

public:
    LatticeFrame3dg(int n, Domain *);
    virtual ~LatticeFrame3dg();
    int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    const char *giveInputRecordName() const override { return _IFT_LatticeFrame3dg_Name; }
    const char *giveClassName() const override { return "latticeframe3dg"; }


protected:
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override;
  //  virtual void computeBmatrixAt( GaussPoint *aGaussPoint, FloatMatrix &answer) override;
    virtual void  computeStrainVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override;
    bool computeGtoLRotationMatrix(FloatMatrix &) override;

  void computeBFmatrixAt( GaussPoint *aGaussPoint, FloatMatrix &answer);

};
} // end namespace oofem
#endif
