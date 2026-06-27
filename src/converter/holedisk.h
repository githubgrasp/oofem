#ifndef converter_holedisk_h
#define converter_holedisk_h

#include "floatarray.h"

#include <vector>

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
    /// Transport-only hole: deleted from the TM mesh but kept in the SM mesh
    /// (the SM keeps the inclusion+ITZ; only the transport domain is voided here,
    /// with its boundary on this circle = the ITZ midline). Created by `#@tmhole`.
    /// Regular holes (`#@holedisk`) are voids in both domains (tmOnly = false).
    bool tmOnly = false;

public:
    HoleDisk(int n, double cx, double cy, double radius);

    int giveNumber() const { return number; }
    double giveCentreX() const { return cx; }
    double giveCentreY() const { return cy; }
    double giveRadius() const { return radius; }
    bool isTmOnly() const { return tmOnly; }
    void setTmOnly(bool t) { tmOnly = t; }

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

    /// One rim mechanical node receiving the hydro-mechanical coupling load.
    struct RimCouplingEntry {
        int delaunayVertex = 0;          ///< global Delaunay-vertex id (mechanical node)
        double dirX = 0., dirY = 0.;     ///< outward radial unit direction at the node
        double tributary = 0.;           ///< tributary boundary length (sum of incident half-edges)
        std::vector< int >neighbourVertices; ///< other endpoint of each incident rim edge (for TM-node lookup)
    };

    /// Identify the mechanical nodes on this hole's rim and fill `entries` with
    /// each node's outward radial direction and tributary boundary length (half
    /// of every incident rim Delaunay edge). Pure geometry on the hole circle —
    /// the grid turns the result into coupling boundary-condition records.
    void computeRimCoupling(Grid *grid, std::vector< RimCouplingEntry > &entries) const;
};


#endif // converter_holedisk_h
