#ifndef interfacedisk_h
#define interfacedisk_h

#include "inclusion.h"
#include "grid.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <sstream>
#endif


/// 2D circular inclusion — the 2D analog of `InterfaceSphere`. Reads
/// `#@interfacedisk <id> centre 2 cx cy radius r [refine <r>] [itz <t>]`.
/// Places points on the circle plus a sister ring at `radius + itz` so
/// the converter can form an aggregate / ITZ / matrix triple of element
/// materials.
class InterfaceDisk : public Inclusion
{
protected:
    /// `[cx cy]`.
    oofem::FloatArray centre;

    /// Inclusion diameter (used to default `itzThickness` when absent).
    double diameter;

    /// Perpendicular thickness of the ITZ halo.
    double itzThickness;

public:

    /**
     * Constructor. Creates a 2D disk inclusion.
     * @param n inclusion number in the grid
     * @param aGrid grid to which the inclusion belongs
     */
    InterfaceDisk(int n, Grid *aGrid);

    /// Destructor.
    ~InterfaceDisk();

    /// Returns the inclusion radius.
    double giveRadius() { return this->radius; }

    /// Returns the disk centre `[cx cy]`.
    const oofem::FloatArray &giveCentre() const { return this->centre; }

    /// Returns the thickness of the ITZ halo around the inclusion.
    double giveITZThickness() { return this->itzThickness; }

    /// Place a ring of points on the circle plus a sister ring offset
    /// outward by `itzThickness`.
    int generatePoints() override;

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "InterfaceDisk"; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@interfacedisk <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);
};

#endif // interfacedisk_h
