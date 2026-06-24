#ifndef converter_holedisk_h
#define converter_holedisk_h

#include "floatarray.h"

class Grid;


/// A circular hole (void) in the converter mesh — the counterpart of the
/// generator's `HoleDisk` inclusion. The generator places mechanical nodes ON
/// the hole/matrix rim; this class owns the converter-side hole geometry:
/// classifying points inside / on the rim, and dropping the dual cross-section
/// (Voronoi/transport) vertices of the rim edges onto the circle so transport
/// nodes lie on the boundary with no tangential edge. Created by `#@holedisk`
/// or `#@diskinclusion ... delete`. Mirrors the "region owns its geometry"
/// design (cf. `Disk::modifyVoronoiCrossSection`), for a void rather than a
/// domain region.
class HoleDisk
{
protected:
    /// Hole id (the `<num>` of the directive); used to reference the hole, e.g.
    /// from a future `#@coupling hole <id>`.
    int number = 0;
    /// Circle centre.
    double cx = 0., cy = 0.;
    /// Circle radius.
    double radius = 0.;

public:
    HoleDisk(int n, double cx, double cy, double radius);

    int giveNumber() const { return number; }
    double giveCentreX() const { return cx; }
    double giveCentreY() const { return cy; }
    double giveRadius() const { return radius; }

    /// True iff (x,y) lies strictly inside the hole (closer to the centre than
    /// `radius`). Rim points (at the radius) are not "inside" — they are kept.
    bool isStrictlyInside(double x, double y) const;

    /// True iff (x,y) lies on the rim (within `tol` of the circle).
    bool onRim(double x, double y, double tol) const;

    /// Project the dual cross-section Voronoi vertices of this hole's rim
    /// mechanical edges radially onto the circle (and flag them on-boundary),
    /// so the transport nodes land on the rim. The 2D-inclusion analog of
    /// `Disk::modifyVoronoiCrossSection`, for one hole.
    void projectVoronoiCrossSection(Grid *grid, double tol) const;

    /// If the Delaunay chord (dA,dB) is a rim edge of this hole (both endpoints
    /// on the circle), set `m` to its midpoint projected radially onto the
    /// circle and return true; otherwise return false and leave `m` untouched.
    bool boundaryNodeCoord(const oofem::FloatArray &dA, const oofem::FloatArray &dB,
                           double tol, oofem::FloatArray &m) const;
};


#endif // converter_holedisk_h
