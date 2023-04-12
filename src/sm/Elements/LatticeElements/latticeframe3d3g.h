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

#ifndef latticeframe3d3g_h
#define latticeframe3d3g_h

#include "latticeframe3d.h"

///@name Input fields for LatticeFrame3d3g
//@{
#define _IFT_LatticeFrame3d3g_Name "latticeframe3d3g"

//@}

namespace oofem {
/**
* This class implements a 3-dimensional frame element based on rigid body spring theory presented in Toi 1991 and Toi 1993. It belongs to the group of lattice models in OOFEM.
* Authors: Gumaa Abdelrhim and Peter Grassl
*/

class LatticeFrame3d3g : public LatticeFrame3d
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
   LatticeFrame3d3g(int n, Domain *);
   virtual ~LatticeFrame3d3g();
   int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;

   double giveArea() override;

   double giveIy() override;

   double giveIz() override;

   double giveIk() override;

   double giveShearAreaY() override;

   double giveShearAreaZ() override;
   double computeLength() override;
  double computeCurrentLength();
   void initializeFrom(InputRecord &ir) override;
   void giveDofManDofIDMask(int inode, IntArray &) const override;

   int giveLocalCoordinateSystem(FloatMatrix &answer) override;

   const char *giveInputRecordName() const override { return _IFT_LatticeFrame3d3g_Name; }
   const char *giveClassName() const override { return "latticeframe3d3g"; }

   void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override;

   Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }



   protected:
   void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
   void computeBDmatrixAt(GaussPoint *, FloatMatrix &);
   //  virtual void computeBmatrixAt( GaussPoint *aGaussPoint, FloatMatrix &answer) override;
   virtual void  computeStrainVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;
   bool computeGtoLRotationMatrix(FloatMatrix &) override;
   void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
   void computeBFmatrixAt( GaussPoint *aGaussPoint, FloatMatrix &answer);

   void computeLumpedMassMatrix( FloatMatrix &answer, TimeStep *tStep )override;
};
} // end namespace oofem
#endif
