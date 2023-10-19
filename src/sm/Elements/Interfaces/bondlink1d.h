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

#ifndef bondlink1dtruss_h
#define bondlink1dtruss_h

#include "../structuralelement.h"

///@name Input fields for BondLink1d
//@{
#define _IFT_BondLink1d_Name "bondlink1d"
#define _IFT_BondLink1d_length "length"
#define _IFT_BondLink1d_diameter "diameter"
#define _IFT_BondLink1d_dirvector "dirvector"
//@}

namespace oofem {
/**
 * This class implements a bond link for connecting truss and continuum elements with coinciding nodes.
 *
 * @author: Peter Grassl
 */

class BondLink1d : public StructuralElement
{
protected:
    double bondLength = 0.;

    FloatMatrix localCoordinateSystem;
    double bondDiameter = 0.;
    FloatArray directionVector;
    int geometryFlag = 0;
    FloatArray rigid;
    FloatArray globalCentroid;

public:
    BondLink1d(int n, Domain *);

    double computeVolumeAround(GaussPoint *aGaussPoint) override;

    double giveLength();

    MaterialMode giveMaterialMode() override { return _1dInterface; }

    virtual double giveBondLength();

    virtual double giveBondDiameter();

    int computeNumberOfDofs() override { return 2; }

    void giveDofManDofIDMask(int inode, IntArray &) const override;


    void giveInternalForcesVector(FloatArray &answer,
                                  TimeStep *tStep, int useUpdatedGpRecord) override;

    const char *giveInputRecordName() const override { return _IFT_BondLink1d_Name; }
    const char *giveClassName()  const override { return "BondLink1d"; }
    void initializeFrom(InputRecord &ir) override;

    Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }

    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;

#ifdef __OOFEG
    void drawYourself(oofegGraphicContext &context, TimeStep *tStep) override;
    void drawRawGeometry(oofegGraphicContext &, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &, TimeStep *tStep, UnknownType) override;
#endif


protected:
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override;

    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;

    /**
     * This function is different from the standard computeGlobalCordinates
     * function as it returns the global coordinates of the gausspoint
     * independent to the value of the lcoords.
     */
    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;

  
    void computeGaussPoints() override;
    integrationDomain  giveIntegrationDomain() const override { return _Line; }
};
} // end namespace oofem
#endif
