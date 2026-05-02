#ifndef disk_h
#define disk_h


#include "grid.h"
#include "region.h"

#include "floatarray.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <sstream>
#endif


/// Solid 2D disk region — the 2D analog of `Sphere`. Reads
/// `#@disk <id> centre 2 cx cy radius r [refine <r>]`. Coordinates are
/// stored internally as 3D with `z = 0`.
class Disk : public Region
{
protected:
    /// `[cx cy]` (z is implied 0).
    oofem::FloatArray centre;

public:

    /**
     * Constructor. Creates a 2D disk region.
     * @param n region number in the grid
     * @param aGrid grid to which the region belongs
     */
    Disk(int n, Grid *aGrid);

    /// Destructor.
    ~Disk();

    /// Returns the disk radius.
    double giveRadius() { return this->radius; }

    /// Generate points: centre, circumference ring, interior random fill
    /// with mirror across the circle.
    int generatePoints() override;

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Disk"; }

    /// Parse keyword/value tokens from an open istringstream positioned
    /// after the `#@disk <num>` prefix.
    void initializeFromTokens(std::istringstream &iss);
};


#endif // disk_h
