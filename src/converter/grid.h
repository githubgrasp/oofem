//   ********************
//   *** CLASS GRID ***
//   ********************

#ifndef grid_h
#define grid_h

#include "converterlistutils.h"
#include "convertererror.h"
#include "floatarray.h"
#include "intarray.h"

#include <unordered_map>
#include <unordered_set>

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <time.h>
 #include <map>
 #include <list>
 #include <set>
#endif




class Vertex;
class Curve;
class Surface;
class Region;
class Line;
class Inclusion;
class GridLocalizer;
class Fibre;
class BoundarySphere;
class BoundaryCylinder;


/// Raw T3D mesh node. `entType`/`entID` tag the geometric entity (vertex,
/// curve, surface, region) the node belongs to, as emitted by T3D.
struct Node {
    int id;
    double x, y, z;

    int entType = 0;   ///< T3D geometric entity type
    int entID   = 0;   ///< T3D geometric entity ID
};

/// T3D surface triangle. `entType` is 3/5/6 (curve/surface/region).
/// `bndCurveId`/`bndCurveProp` per triangle edge are only populated when
/// T3D was run with `outType & 8` (boundary-curve info requested).
struct Tri {
    int id;
    int n1, n2, n3;
    int entType = 0;  ///< 3=curve, 5=surface, 6=region
    int entID   = 0;
    int entProp = 0;

    /// Edge→curve mapping; only meaningful when (outType & 8).
    /// Indices: 0=(n1-n2), 1=(n2-n3), 2=(n3-n1).
    int bndCurveId[ 3 ]   = { 0, 0, 0 };
    int bndCurveProp[ 3 ] = { 0, 0, 0 };
};

/// T3D tetrahedron with per-face geometric-entity tags for boundary lookup.
struct Tet {
    int id;
    int n1, n2, n3, n4;

    int entType;
    int entID;
    int entProp;

    int faceEntID[ 4 ];
    int faceEntType[ 4 ];
    int faceEntProp[ 4 ];
};

/// Undirected mesh edge linking two nodes. Adjacent triangle/tetra indices
/// are -1 until populated by `buildEdgeAdjacency3D`.
struct Edge {
    int n1, n2;     ///< node ids (sorted ascending)
    int tri1 = -1;  ///< adjacent triangle index (-1 if none)
    int tri2 = -1;
    int tet1 = -1;  ///< adjacent tetra index (-1 if none)
    int tet2 = -1;
};


/// Curve segment derived from a boundary triangle edge: the two endpoints
/// and the T3D curve id the segment belongs to. Produced by
/// `buildCurveSegsFromTris`.
struct CurveSeg {
    int n1, n2;     ///< sorted node ids
    int curveID;
};



/**
 * Description of class grid
 */
class Grid
{
private:

    /** Grid type. This determines the type of input to generate
     */
    enum GridType { _3dSM, _3dTM, _2dSM, _2dTM };

    GridType gridType;

    /// Out-of-plane thickness for 2D lattice elements (`thick` field of the
    /// `lattice2D` directive). Set by `#@thickness <t>`; defaults to 1.0.
    double latticeThickness = 1.0;

    /// Spatial dimension parsed from the `mesh.nodes` header (2 or 3).
    int spatialDim = 3;

    /// Per-Voronoi-vertex flag: true if the vertex was strictly inside
    /// any delete-mode `#@notch` box before the projection step. Indexed
    /// 1..nVoronoiVertices; entry 0 is unused (qhull at-infinity sentinel).
    std::vector< char > voronoiOrigInsideDeletingNotch;

    std::vector< std::vector< int > >edgeToTets;
    std::vector< std::vector< int > >edgeToBoundaryTris;

    std::string controlFileName;

    std::map< int, std::map< int, std::vector< int > > >entityNodes;

    std::map< int, std::map< int, std::vector< int > > >entityTris;

    std::vector< CurveSeg >curveSegs;
    std::unordered_map< int, std::vector< int > >curveToSegIdx;

    bool use3DFrameSection = false;

    double frameRadius = 0.0;


    struct BCRequest {
        int entType;           // 1=vertex, 2=curve, 5=patch
        int entID;
        std::vector< int >extraVertices; // endpoints for curves

        std::vector< int >dofs;
        std::vector< double >values;
        int setID = -1; // NEW
    };

    std::vector< BCRequest >bcRequests;

    struct SetDef {
        int setID;
        std::vector< int >nodeIDs; // sorted
    };

    std::vector< SetDef >generatedNodeSets;


    double defaultThickness = 1.0;

    std::map< int, std::map< int, double > >entityThickness;
    // entityThickness[entType][entID] = thickness

    /// Grid number
    int number;

    // Determines the type of macroscopic element, so that a corresponding number of DOFs is created for the control node
    const char *boundElemName;
    const char *boundBeamElemName;

    double diameter, density;

    double liveLoad_q = 0.0;   // N/m^2

    double shellWidthScale = -1.0; // < 0 means: use area-based scaling


    bool liveLoadEnabled = false;

    /// When true, giveOutput() also writes the auxiliary POV-Ray rendering files
    /// (`*.vor.line.pov`, `*.del.cross.pov`, etc.). Default off — opt in with the
    /// `#@pov` directive in the control file.
    bool emitPovOutput = false;

    /// When true, giveOutput() also writes VTK files for ParaView visualisation
    /// (`*.voronoielement.vtu`, `*.delaunayelement.*.vtu`, etc.). Default off —
    /// opt in with the `#@vtk` directive in the control file.
    bool emitVtkOutput = false;

    oofem::FloatArray liveDir; // size 3

    /// Per-entity element override for the T3D writer. Populated by the
    /// `#@element <entityKind> <entityID> <name> <crossSect> <mat>` directive
    /// in the control file. Empty by default — writers fall back to their
    /// hardcoded default element name. Priority during resolution:
    /// curve (lowest entType) > surface > tetra/region.
    struct EdgeSpec {
        std::string elementName;
        int crossSect = 1;
        int material  = 1;
    };
    std::map< std::pair< int, int >, EdgeSpec >elementSpecsByEntity;

    /// A notch is an axis-aligned box. By default (`material <m>`) matrix
    /// elements whose midpoint falls inside have their material reassigned to
    /// `<m>`. With the `delete` keyword the element is omitted entirely —
    /// useful for sharp notches and small voids (the dual mesh effectively
    /// gains a hole). Populated by either of:
    ///   `#@notch <id> box 6 xmin ymin zmin xmax ymax zmax material <m>`
    ///   `#@notch <id> box 6 xmin ymin zmin xmax ymax zmax delete`
    struct NotchSpec {
        double xmin, ymin, zmin, xmax, ymax, zmax;
        int material = 3;
        bool deleteFlag = false;
    };
    std::vector< NotchSpec >notchSpecs;

    /// A sphere inclusion: elements whose endpoints both fall inside get
    /// `inside` material; elements that straddle the sphere (one endpoint
    /// inside, one outside — the ITZ zone) get `interface` material.
    /// Populated by the `#@sphereinclusion <id> centre 3 x y z radius r
    /// itz t inside <mi> interface <mif>` directive.
    struct SphereInclusionSpec {
        double cx = 0., cy = 0., cz = 0.;
        double radius = 0., itz = 0.;
        int inside    = 2;
        int interface_ = 3;
    };
    std::vector< SphereInclusionSpec >sphereInclusionSpecs;

    /// A straight-axis cylindrical inclusion with ITZ halo. The inclusion
    /// classification uses perpendicular distance from each endpoint to the
    /// (infinite) axis through the two line points. Populated by the
    /// `#@cylinderinclusion <id> line 6 x1 y1 z1 x2 y2 z2 radius r itz t
    ///   inside <mi> interface <mif>` directive.
    struct CylinderInclusionSpec {
        double x1 = 0., y1 = 0., z1 = 0.;
        double x2 = 0., y2 = 0., z2 = 0.;
        double radius = 0., itz = 0.;
        int inside    = 2;
        int interface_ = 3;
    };
    std::vector< CylinderInclusionSpec >cylinderInclusionSpecs;

    /// Per-material bodyload map: material id → bc id. Populated by the
    /// `#@bodyload <mat> <bc_id>` directive. When an element is emitted
    /// with crossSect/mat equal to <mat>, the writer appends
    /// " bodyloads 1 <bc_id>" to its record. Decouples bodyload placement
    /// from inclusion directives: Wong sets bodyload on the matrix
    /// material only; the corrosion cylinder sets it on the interface
    /// material only.
    std::map< int, int >bodyloadByMaterial;

    /// Emit any Delaunay vertex on a specified face of the region's bounding
    /// box as a `rigidarmnode` slaved to a named control vertex. Populated by
    /// `#@rigidarm <master_ctl_id> face <axis> <side> mastermask 6 m1..m6
    /// doftype 6 d1..d6`. axis is 1|2|3 (x|y|z); side is min|max. The master
    /// vertex itself is kept as a regular node. Use-case: cantilever ends
    /// clamped as rigid bodies attached to support/load control nodes.
    struct RigidArmSpec {
        int masterCtlId = 0;
        int axis = 1;       // 1=x, 2=y, 3=z
        bool sideMax = false;
        oofem::IntArray mastermask;   // 6 ints
        oofem::IntArray doftype;      // 6 ints
    };
    std::vector< RigidArmSpec >rigidArmSpecs;

    /// Slave the chosen DOFs of every Delaunay vertex on the named face of
    /// the rect/prism box to a control vertex via `DT_simpleSlave` — the
    /// slave's DOF value is identical to the master's, with no rigid-arm
    /// rotation kinematics. Populated by `#@slaveside <master_ctl_id> face
    /// <axis> <side> dofs <list>`. Use-case: pulling one face of a non-
    /// periodic specimen uniformly without coupling rotation.
    struct SlaveSideSpec {
        int masterCtlId = 0;
        int axis = 1;                  // 1=x, 2=y, 3=z
        bool sideMax = false;
        oofem::IntArray slavedDofs;    // subset of {D_u=1, D_v=2, D_w=3, R_u=4, R_v=5, R_w=6}
    };
    std::vector< SlaveSideSpec >slaveSideSpecs;

    /// Tag every node lying on the chosen face plane of region 1's bounding
    /// box with `bc <n> <id1> <id2> …` so OOFEM applies the listed
    /// `BoundaryCondition` records to those nodes. Populated by
    /// `#@nodebc <bc_id> face <axis> <min|max>`. Multiple specs may name the
    /// same face (each contributes one BC id) or different faces.
    /// Use-case: prescribing saturation on top/bottom of a transport mesh
    /// without a post-process script.
    struct NodeBCSpec {
        int bcId = 0;
        int axis = 1;       // 1=x, 2=y, 3=z
        bool sideMax = false;
    };
    std::vector< NodeBCSpec >nodeBCSpecs;

    /// When true, TM element emitters append `lumpedcapacity 1` to every
    /// `latticemt2D` / `latticemt3D` line. Enabled by `#@lumpedcapacity 1`.
    /// Switches the element capacity matrix from the consistent (coupled)
    /// form to a TPFA-monotone diagonal form — needed for nonlinear c(p)
    /// like Richards-style unsaturated flow.
    bool emitLumpedCapacity = false;

    /// Any Delaunay line with either endpoint equal to the named control
    /// vertex gets the specified material. Populated by `#@material_around
    /// <ctl_id> material <m>`. Precedence: material_around < notch <
    /// inclusion — notch/inclusion directives override. Use-case: cantilever
    /// elastic elements attached to the support/load points.
    struct MaterialAroundSpec {
        int ctlId = 0;
        int material = 1;
    };
    std::vector< MaterialAroundSpec >materialAroundSpecs;

    /// When true, emitters append `couplingflag 1 couplingnumber N <ids>`
    /// to each element record. For an SM lattice3D element (at a Delaunay
    /// line), the ids are the Voronoi-line cross-section elements; for a TM
    /// latticemt3D element (at a Voronoi line), they are the Delaunay-line
    /// cross-section elements. Image-side (flag==1) cross-section elements
    /// get swapped for their periodic partner via `givePeriodicElement()`.
    /// Enabled by the `#@couplingflag` directive.
    bool emitCouplingFlag = false;

    /// Control-vertex definitions from `#@controlvertex <id> coords 3 x y z`.
    /// Each entry declares a specific mesh coordinate whose nearest Delaunay
    /// vertex id is exposed at write time via the `#@CTL<id>` inline placeholder.
    std::vector< std::pair< int, oofem::FloatArray > >controlVertexDefinitions;
    /// Resolved mapping: directive id → raw Delaunay vertex index (1-based).
    std::map< int, int >controlNodeIds;

    int t3dOutType = 0;

    int dimension;

    double TOL;

    oofem::IntArray periodicityFlag;

    int maxIter, randomInteger;

    int periodicNodeCounter;

    //for t3d shells
    double shellThickness = 1.0;

    double beamWidth = -1.0;
    double beamThickness = -1.0;

    std::vector< Vertex * >delaunayVertexList;

    std::vector< Vertex * >voronoiVertexList;

    std::vector< Vertex * >reinforcementNodeList;

    std::vector< Line * >delaunayLineList;

    std::vector< Line * >voronoiLineList;

    std::vector< Fibre * >fibreList;

    std::vector< Line * >latticeBeamList;

    std::vector< Line * >latticeLinkList;

    std::vector< Vertex * >interNodeList;

    std::vector< Vertex * >vertexList;

    std::vector< Curve * >curveList;

    std::vector< Surface * >surfaceList;

    std::vector< Region * >regionList;

    std::vector< Inclusion * >inclusionList;

    typedef std::list< int >nodeContainerType;

    std::vector< Node >nodes;

    std::vector< Tri >tris;

    std::vector< Tet >tets;

    std::vector< Edge >edges;

    std::unordered_map< int, int >nodeIndex;

    std::vector< double >edgeWidth;

    std::vector< int >loadNodeSetID;  // size = nodes.size()
    std::vector< double >loadFx, loadFy, loadFz;


public:

    /// Octree localiser over Delaunay vertices. Owned by Grid, used to answer
    /// spatial queries (nearest-vertex, in-box) during control-file parsing
    /// and during geometry→mesh resolution.
    GridLocalizer *delaunayLocalizer;
    /// Octree localiser over Voronoi vertices (dual of the Delaunay mesh).
    GridLocalizer *voronoiLocalizer;
    /// Octree localiser over reinforcement (fibre) nodes.
    GridLocalizer *reinforcementLocalizer;

    /// Constructor. `n` is the grid number (in practice always 1).
    Grid(int n);

    /// Destructor — deletes owned components and localisers.
    ~Grid();


    /// Returns the coordinates of node `nid` from the raw T3D node list.
    oofem::FloatArray getX(int nid) const;

    /// Returns the outward-pointing unit normal of the triangle at index
    /// `triIndex` in the T3D triangle list.
    oofem::FloatArray triNormal(int triIndex) const;

    /// Populate `edgeWidth` (size = edges.size()) with the width used to
    /// derive 3D lattice cross-section areas. `thickness` is the default
    /// shell thickness; per-entity overrides from `#@thickness` apply first.
    void computeEdgeWidths(double thickness);

    /// Target cross-section area for edge `edgeIndex` given its length `Le`,
    /// based on the stored widths and local geometry.
    double computeTargetArea(int edgeIndex, double Le) const;

    /// Connect boundary-triangle edges to T3D curves so that the T3D curve
    /// numbering can be referenced from control-file directives.
    void buildCurveSegsFromTris();

    /// Connect surface triangles to T3D patches so that patch numbers can
    /// reference triangles from control-file directives.
    void buildBoundaryTrisFromTets();

    /// Populate `edges` with unique edges across tris/tets, recording the
    /// adjacent triangle and tetra indices on each edge (Edge::tri1/2, tet1/2).
    void buildEdgeAdjacency3D();

    /// Barycentre (coord average of 4 vertices) of tetrahedron `tetIndex`.
    oofem::FloatArray tetBarycentre(int tetIndex) const;
    /// Barycentre (coord average of 3 vertices) of triangle `triIndex`.
    oofem::FloatArray faceBarycentre(int triIndex) const;

    /// Rebuild the `entityTris` index after tri list changes (e.g. after
    /// filtering). Maps (entType, entID) → list of triangle indices.
    void rebuildEntityTris();

    /// Signed volume of tetrahedron `tetIndex`; sign depends on node order.
    double tetVolume(int tetIndex) const;

    /// Area of the polygon `polycoords` projected onto the plane through
    /// `xm` with orthonormal basis (r, s). Used to compute lattice cross-
    /// section areas from 3D Voronoi polygons.
    double computePolygonAreaProjected(const oofem::FloatArray &polycoords,
                                       const oofem::FloatArray &xm,
                                       const oofem::FloatArray &r,
                                       const oofem::FloatArray &s) const;


    /// Assemble the Voronoi polygon associated with edge `edgeIndex` and
    /// return its vertices concatenated in `polycoords` (x1 y1 z1 x2 y2 z2 …).
    void buildEdgePolygon3D(int edgeIndex, oofem::FloatArray &polycoords) const;
    /// Write one OOFEM element record for a 3D lattice edge, resolving its
    /// cross-section, material and polygon. Increments `eid`.
    void write3DEdgeSection(std::ostream &out, int &eid, const Edge &e, int edgeIndex);


    /// Euclidean distance between raw mesh nodes `n1` and `n2`.
    double segLength(int n1, int n2) const;

    /// Fill `L` (size = nodes.size()) with the half-length each node
    /// contributes on curve `curveID`. Used to convert per-unit-length
    /// loads into per-node loads for `set`-based load emission.
    void computeNodalLengthsOnCurve(int curveID, std::vector< double > &L) const;

    /// Barycentre of triangle `triIndex` (alias of `faceBarycentre`).
    oofem::FloatArray triBarycentre(int triIndex) const;


    /// Read the control file (`controlFileName`) and dispatch `#@…`
    /// directives into Grid state. Shared by both the T3D and qhull
    /// pipelines — each directive is consumed by whichever writers care
    /// about it; the disjoint directive sets mean the two pipelines never
    /// interfere with one another's state.
    void readControlRecords();

    /// Read a packing-format file (as produced by `src/aggregate/`) and
    /// populate inclusion specs. Each `sphere` line contributes a
    /// `SphereInclusionSpec` whose centre and radius come from the line
    /// while `itz` and the inside/interface material ids come from the
    /// `#@inclusionfile` directive. Each `fibre` line is appended to
    /// `fibreList` (renumbered to avoid collision with `#@fibre` ids).
    /// Ellipsoid lines are unsupported by the converter's material-
    /// classification logic and trigger a warning.
    void readInclusionFile(const std::string &path,
                           double itz,
                           int insideMaterial,
                           int interfaceMaterial);

    /// Resolve the node-id set targeted by a BC request (vertex → single
    /// node; curve → nodes on curve; patch → nodes on surface patch).
    std::set< int >collectBCNodes(const BCRequest &bc) const;

    /// Euclidean length of an Edge entry using the node coordinates.
    double edgeLength(const Edge &e) const;

    /// Area of triangle `triIndex` in the T3D triangle list.
    double triArea(int triIndex) const;

    /// Build per-load sets from `loadRequests`: one set per load entity,
    /// populated with the nodes the load should distribute over.
    void prepareLiveLoadSets();

    /// Thickness for a T3D (entType, entID) pair, falling back to
    /// `defaultThickness` if no `#@thickness` directive targeted it.
    double giveThicknessForEntity(int entType, int entID) const;

    /// True if `t` starts with `#@` (the converter directive prefix).
    bool isConverterDirective(const std::string &t) const
    {
        return t.rfind("#@", 0) == 0;
    }

    /// Debug helper: sum shell-triangle areas against the sum of per-edge
    /// contributions and abort if they disagree beyond tolerance.
    void checkAreaConservation() const;


    /// A live-load request parsed from a `#@LOAD` directive. Targets a
    /// geometric entity (vertex/curve/surface/shell) and carries the load
    /// intensity `q` plus the load-time-function id `ltf`.
    struct LoadRequest {
        int entType = -1; ///< 1=vertex, 2=curve, 3=surface, 5=patch, 6=shell
        int entID   = -1;
        double q    = 0.0;///< N/m² for surfaces, N/m for curves, N for vertices
        int setID   = -1;
        int ltf     = 1;  ///< default load-time function id
    };
    std::vector< LoadRequest >loadRequests;

    /// Map a string entity keyword (e.g. "vertex", "curve", "surface",
    /// "patch", "shell") to its integer T3D entType.
    int entityTypeFromString(const std::string &s) const;


    /// Fill `A` (size = nodes.size()) with the area each node contributes
    /// on a triangular entity (surface/patch/shell), distributing each
    /// triangle's area 1/3 to each vertex. Used to convert surface loads
    /// into per-node loads.
    void computeNodalAreasOnTriEntity(int entType, int entID, std::vector< double > &A) const;

    /// Grid number (in practice always 1).
    int giveNumber() { return this->number; }

    /// Emit NodalLoad records for all live-load sets, using ids starting
    /// at the caller-owned counter `bcID` (incremented per emission).
    void writeLiveLoads(std::ostream &out, int &bcID);

    /// Materialise the sets referenced by `bcRequests` — populates
    /// `bcRequests[i].setID` and appends the resolved set into
    /// `generatedNodeSets`.
    void prepareBCSets();

    /// Emit BoundaryCondition records (Dirichlet) for all BC requests,
    /// referencing the sets previously built by `prepareBCSets`.
    void writeBCRecords(std::ostream &out, int &bcID) const;

    /// For a given edge, return the (entType, entID) of the most-specific T3D
    /// geometric entity it belongs to, with priority curve > surface > region.
    /// Returns {0, 0} if no entity is associated.
    std::pair< int, int >entityForEdge(const Edge &e) const;

    /// Resolve the EdgeSpec for an edge, falling back to `defaultSpec` if no
    /// `#@element` directive matches the edge's entity.
    EdgeSpec resolveEdgeSpec(const Edge &e, const EdgeSpec &defaultSpec) const;

    /// For a matrix lattice line with endpoint coords A/B, return the notch
    /// material if the midpoint is inside any `#@notch ... material <m>` box,
    /// else `defaultMat`. Delete-mode notches are ignored here — see
    /// `notchDeletes` for that path.
    int resolveNotchMaterial(const oofem::FloatArray &A, const oofem::FloatArray &B,
                             int defaultMat) const;

    /// True iff the midpoint of segment A-B lies inside any delete-mode
    /// `#@notch` box. The caller skips emission of the corresponding element.
    bool notchDeletes(const oofem::FloatArray &A, const oofem::FloatArray &B) const;

    /// For a matrix lattice line with endpoint coords A/B, resolve the sphere-
    /// inclusion material. Returns the `inside` material if both endpoints fall
    /// within `radius + itz/2` of a sphere centre; the `interface` material
    /// if exactly one endpoint does; otherwise returns `defaultMat`.
    int resolveInclusionMaterial(const oofem::FloatArray &A, const oofem::FloatArray &B,
                                 int defaultMat) const;

    /// Match each `#@controlvertex` declaration to the nearest Delaunay vertex
    /// and populate `controlNodeIds`. Called after the Delaunay vertex list is loaded.
    void resolveControlVertices();

    void writeGeneratedSets(std::ostream &out) const;

    /// Geometric tolerance (typically 1e-6 of bounding-box extent).
    double giveTol() { return this->TOL; }

    /// Override the grid number. Rarely used (grids are always grid 1).
    void setNumber(int nn) { this->number = nn; }

    /// Maximum iterations for the iterative mesh-generation algorithms.
    int giveMaximumIterations() { return this->maxIter; }

    /// Nominal particle diameter (used by discrete-particle mesh generation).
    double giveDiameter() { return this->diameter; }

    /// Copy the 3-entry periodicity flag into `answer` (per-axis: 0 none,
    /// 1 full, 2 partial/surface).
    void givePeriodicityFlag(oofem::IntArray &answer) { answer = this->periodicityFlag; }

    /// Seed for the random-number generator used during mesh generation.
    int giveRandomInteger() { return this->randomInteger; }

    /// Write Voronoi lines as POV-Ray cylinders to `filename`.
    void giveVoronoiPOVOutput(const std::string &filename);

    /// Write Delaunay lines as POV-Ray cylinders to `filename`.
    void giveDelaunayPOVOutput(const std::string &filename);

    /// Write Voronoi cell geometry as VTK unstructured grid.
    void giveVoronoiCellVTKOutput(FILE *outputStream);

    /// Write Delaunay elements (lattice lines) as VTK unstructured grid.
    void giveDelaunayElementVTKOutput(FILE *outputStream);

    /// Variant of `giveDelaunayElementVTKOutput` restricted to elements
    /// whose material id equals `nb_mtx` — useful for per-material views.
    void giveDelaunayElementVTKOutput2(FILE *outputStream, int nb_mtx);

    /// Write Delaunay cross-section polygons as VTK unstructured grid.
    void giveDelaunayCrossSectionVTKOutput(FILE *outputStream);

    /// Write Voronoi cross-section polygons as VTK unstructured grid.
    void giveVoronoiCrossSectionVTKOutput(FILE *outputStream);

    /// Write Voronoi elements (dual lattice lines) as VTK unstructured grid.
    void giveVoronoiElementVTKOutput(FILE *outputStream);

    /// Set `gridType` from a control-file keyword ("3dSM", "3dTM"). Aborts
    /// on an unrecognised name.
    void resolveGridType(const std::string &name);

    /// Returns the `n`-th Vertex (1-based) or nullptr if out of range.
    Vertex *giveVertex(int n);

    /// Returns the `n`-th Curve (1-based) or nullptr if out of range.
    Curve *giveCurve(int n);

    /// Returns the `n`-th Surface (1-based) or nullptr if out of range.
    Surface *giveSurface(int n);

    /// Returns the `n`-th Region (1-based) or nullptr if out of range.
    Region *giveRegion(int n);

    /// Returns the `n`-th Inclusion (1-based) or nullptr if out of range.
    Inclusion *giveInclusion(int n);
    /// Returns the `n`-th Fibre (1-based) or nullptr if out of range.
    Fibre *giveFibre(int n);


    /// Returns the `n`-th Delaunay vertex (1-based).
    Vertex *giveDelaunayVertex(int n);

    /// Returns the `n`-th Voronoi vertex (1-based).
    Vertex *giveVoronoiVertex(int n);

    /// Returns the `n`-th Delaunay line (1-based). Delaunay lines are the
    /// lattice "beam" elements connecting neighbouring particles.
    Line *giveDelaunayLine(int n);



    /// Returns the `n`-th Voronoi line (1-based). Voronoi lines carry the
    /// dual (transport) mesh.
    Line *giveVoronoiLine(int n);
    /// Returns the `n`-th fibre-coupling link (1-based).
    Line *giveLatticeLink(int n);
    /// Returns the `n`-th lattice beam segment of a fibre (1-based).
    Line * giveLatticeBeam(int n);
    /// Returns the `n`-th reinforcement (fibre) node (1-based).
    Vertex *giveReinforcementNode(int n);
    /// Returns the `n`-th intermediate node (fibre↔matrix coupling helper).
    Vertex *giveInterNode(int n);

    /// Main entry point used by the T3D pipeline. Runs geometry resolution
    /// and writes the requested output files.
    int generateOutput();

    /// In-place sort of `randomNumbers` into ascending `sortedRandomNumbers`.
    void sortRandomNumbers(oofem::FloatArray &sortedRandomNumbers, oofem::FloatArray &randomNumbers);

    /// Build `rankVector[i]` = rank of `randomNumbers[i]` (1-based).
    void createRankTable(oofem::IntArray &rankVector, oofem::FloatArray &randomNumbers);


    /** Order orientation of cross-section vertices to be counterclockwise.
     *  This should only be needed for cross-sections of Voronoi elements with more than 3 vertices.
     *  (Special case of boundary elements) */
    void orderDelaunayCrossSectionVertices(int elementNumber);


    /// Qhull pipeline entry point. Reads the qhull `mesh.nodes` +
    /// `mesh.voronoi` pair, parses the control template, and emits the
    /// OOFEM input. Returns 0 on success.
    int instanciateYourselfFromQhull(const std::string &controlFile, const char *nodeFileName, const char *voronoiFileName);


    /// Read a T3D mesh from file `fn` into `nodes`/`tris`/`tets`. Returns
    /// false on parse error.
    bool readT3d(const std::string &fn,
                 std::vector< Node > &nodes,
                 std::vector< Tri > &tris,
                 std::vector< Tet > &tets);

    /// Derive the unique undirected edge list for `tris` and `tets` into
    /// `edges`, with adjacency indices filled in.
    void buildEdges(const std::vector< Tri > &tris,
                    const std::vector< Tet > &tets,
                    std::vector< Edge > &edges);


    /// Write one OOFEM `node` record per T3D node to `out`.
    void writeT3dNodesOofem(std::ostream &out);

    /// Write one OOFEM element record per T3D edge to `out`.
    void writeT3dElemsOofem(std::ostream &out);

    /// T3D pipeline entry point. Reads `t3dFileName`, parses control
    /// directives from `controlFileName`, and emits the OOFEM input.
    /// Returns 0 on success.
    int instanciateYourselfFromT3d(const std::string &t3dFileName, const std::string &controlFileName);

    /// Numerical Recipes ran1 uniform random generator. `idum` is a
    /// caller-owned seed state.
    double ran1(long *idum);

    /// Returns number of dof managers in grid.
    int                giveNumberOfVertices() { return converter::size1(vertexList); }
    /// Returns number of lines in grid.
    int                giveNumberOfCurves() { return converter::size1(curveList); }
    /// Returns number of material models in grid
    int                giveNumberOfSurfaces() { return converter::size1(surfaceList); }
    /// Returns number of cross section models in grid
    int                giveNumberOfRegions() { return converter::size1(regionList); }

    int                giveNumberOfInclusions() { return converter::size1(inclusionList); }

    int                giveNumberOfDelaunayVertices() { return converter::size1(delaunayVertexList); }

    int                giveNumberOfVoronoiVertices() { return converter::size1(voronoiVertexList); }
    int                giveNumberOfDelaunayLines() { return converter::size1(delaunayLineList); }

    int                giveNumberOfVoronoiLines() { return converter::size1(voronoiLineList); }

    int                giveNumberOfLatticeBeams() { return converter::size1(latticeBeamList); }

    int                giveNumberOfLatticeLinks() { return converter::size1(latticeLinkList); }

    int                giveNumberOfReinforcementNode() { return converter::size1(reinforcementNodeList); }
    int                giveNumberOfInterNodes() { return converter::size1(interNodeList); }

    int giveNumberOfFibres() { return converter::size1(fibreList); }


    /// Resize the Delaunay vertex list to `_newSize` entries.
    void resizeDelaunayVertices(int _newSize);
    /// Resize the Voronoi vertex list to `_newSize` entries.
    void resizeVoronoiVertices(int _newSize);
    /// Resize the Delaunay line list to `_newSize` entries.
    void resizeDelaunayLines(int _newSize);
    /// Resize the Voronoi line list to `_newSize` entries.
    void resizeVoronoiLines(int _newSize);
    //@{
    /// Resizes the internal data structure to accomodate space for _newSize dofManagers
    void resizeVertices(int _newSize);
    /// Resizes the internal data structure to accomodate space for _newSize lines
    void resizeCurves(int _newSize);
    /// Resizes the internal data structure to accomodate space for _newSize cross section models
    void resizeSurfaces(int _newSize);
    /// Resizes the internal data structure to accomodate space for _newSize materials
    void resizeRegions(int _newSize);
    /// Resizes the internal data structure to accomodate space for _newSize materials
    void resizeInclusions(int _newSize);


    /// Install `obj` as the `i`-th Delaunay vertex (1-based). Grid takes ownership.
    void setDelaunayVertex(int i, Vertex *obj);
    /// Install `obj` as the `i`-th Voronoi vertex (1-based). Grid takes ownership.
    void setVoronoiVertex(int i, Vertex *obj);
    /// Install `obj` as the `i`-th Delaunay line (1-based). Grid takes ownership.
    void setDelaunayLine(int i, Line *obj);
    /// Install `obj` as the `i`-th Voronoi line (1-based). Grid takes ownership.
    void setVoronoiLine(int i, Line *obj);


    /// Sets i-th componet. The component will be futher managed and maintained by grid object.
    void setVertex(int i, Vertex *obj);
    /// Sets i-th componet. The component will be futher managed and maintained by grid object.
    void setCurve(int i, Curve *obj);
    /// Sets i-th componet. The component will be futher managed and maintained by grid object.
    void setSurface(int i, Surface *obj);
    /// Sets i-th componet. The component will be futher managed and maintained by grid object.
    void setRegion(int i, Region *obj);

    /// Sets i-th componet. The component will be futher managed and maintained by grid object.
    void setInclusion(int i, Inclusion *obj);

    /// Install `obj` as the `i`-th fibre (1-based). Grid takes ownership.
    void setFibre(int i, Fibre *obj);
    /// Install `obj` as the `i`-th lattice link (1-based). Grid takes ownership.
    void setLatticeLink(int i, Line *obj);
    /// Install `obj` as the `i`-th lattice beam (1-based). Grid takes ownership.
    void setLatticeBeam(int i, Line *obj);
    /// Install `obj` as the `i`-th reinforcement node (1-based). Grid takes ownership.
    void setReinforcementNode(int i, Vertex *obj);
    /// Install `obj` as the `i`-th intermediate node (1-based). Grid takes ownership.
    void setInterNode(int i, Vertex *obj);


    /// Class name ("Grid") for debug/log output.
    const char *giveClassName() const { return "Grid"; }

    /// Main dispatch for the qhull pipeline — switches on `gridType` and
    /// calls the matching writer (`give3DSMOutput` / `give3DTMOutput`).
    void giveOutput(const std::string &fileName);
    /// T3D pipeline dispatch — writes the OOFEM input generated from T3D
    /// meshes to `fileName`.
    void giveOutputT3d(const std::string &fileName);

    /// Write the full OOFEM `.in` input to `fileName` (header, counts,
    /// nodes, elements, cross-sections, materials, BCs, load-time
    /// functions, sets, extractor block).
    void giveOofemOutput(const std::string &fileName);

    /// Write the Voronoi / Delaunay / cross-section VTK files alongside
    /// `fileName` (base name for the .vtu outputs).
    void giveVtkOutput(const std::string &fileName);
    /// Alternative VTK writer that emits one set of VTK files per material
    /// (0 … nb_of_mt-1) for per-material visualisation.
    void giveVtkOutput2(const std::string &fileName, int nb_of_mt);

    /// Write the POV-Ray rendering files alongside `fileName`.
    void givePOVOutput(const std::string &fileName);

    /// Qhull writer for 3D structural-mechanics analyses — emits lattice3D
    /// lines (Delaunay edges) with cross-section polygons and materials
    /// resolved from notch/inclusion/material_around directives. Also
    /// handles rigidarm nodes and couplingflag emission.
    void give3DSMOutput(const std::string &fileName);

    /// 2D structural-mechanics writer for the qhull pipeline. Emits one
    /// `lattice2D` line per Delaunay edge whose midpoint is inside the
    /// rectangle (Stage 6 MVP — periodicity is not yet implemented). Each
    /// element is parameterised with `gpCoords 2 gx gy` (segment midpoint),
    /// `width <w>` (length of the dual Voronoi edge), and
    /// `thick <t>` (the `#@thickness` directive value). Notch deletion and
    /// `#@diskinclusion` material resolution are applied via the same
    /// `notchDeletes` / `resolveNotchMaterial` / `resolveInclusionMaterial`
    /// pipeline as the 3D writer.
    void give2DSMOutput(const std::string &fileName);

    /// 2D mass-transport writer for the qhull pipeline. Emits one
    /// `latticemt2D` line per Voronoi edge whose endpoints are both inside
    /// the rectangle. The "nodes" of a transport element are Voronoi
    /// vertices (the dual mesh); the cross-section width is the length of
    /// the Delaunay edge that the Voronoi edge is dual to. gpCoords is the
    /// midpoint of the Voronoi-edge segment. Same `#@thickness` directive
    /// drives the out-of-plane thickness.
    void give2DTMOutput(const std::string &fileName);

    /// Snap Voronoi vertices that fall outside the 2D `#@rect` region to
    /// the nearest rectangle face by clamping coordinates, and snap
    /// Voronoi vertices that fall strictly inside any 2D `#@notch` box to
    /// the nearest notch face. Mirrors the 3D `Prism::modifyVoronoiCross`
    /// projection step so transport nodes (Voronoi vertices) sit on the
    /// physical boundaries. Idempotent — call once after qhull parsing.
    void project2DVoronoiVerticesToBoundaries();

    /// 3D notch counterpart of the 2D projection: Voronoi vertices that
    /// fall strictly inside any 3D `#@notch` box AND bound a crossing
    /// Voronoi edge (one endpoint inside the notch, one outside) get
    /// snapped to the nearest notch face. The 3D *outer* boundary
    /// projection is handled separately by `Prism::modifyVoronoiCross`;
    /// this method only adds the inner-boundary (notch) projection so 3D
    /// matches the behaviour now active in 2D. Idempotent.
    void project3DVoronoiVerticesToNotches();

    /// True iff Voronoi vertex `id` was strictly inside any delete-mode
    /// `#@notch` box BEFORE the projection step ran. Used by the TM
    /// writers to delete elements whose two Voronoi endpoints were both
    /// originally inside a void — distinct from "midpoint inside notch"
    /// because crossing edges have one endpoint snapped to the surface
    /// after projection, yet should be kept (snapped) rather than
    /// deleted. Populated alongside the projection passes.
    bool wasVoronoiOriginallyInsideDeletingNotch(int id) const;

    /// Qhull writer for 3D transport-mechanics analyses — emits
    /// latticemt3D elements on the Voronoi dual of the Delaunay mesh,
    /// with the same directive-driven material/inclusion handling.
    void give3DTMOutput(const std::string &fileName);

    /// Walk fibreList: discretise each fibre, build reinforcement nodes,
    /// lattice beams (fibre segments), and lattice links (fibre↔matrix coupling).
    void discretizeFibres();


    /// Evaluate a degree-n polynomial with coefficients `a[0..n]` at `x`.
    double  dpolyValue(int n, double a[], double x);

    /// Numerical Recipes indexx: fill `indx[1..n]` with the 1-based indices
    /// of `arr` in ascending-value order.
    void indexx(int n, double arr[], int indx[]);

    /// Numerical Recipes integer vector allocator, range [nl..nh].
    int *ivector(int nl, int nh);

    /// Numerical Recipes paired deallocator for `ivector`.
    void free_ivector(int *v, int nl, int nh);

    /// Create a reinforcement (fibre) node at `coordR` and register it
    /// with the grid. Used by Fibre during discretisation.
    Vertex *createReinfNode(oofem::FloatArray coordR);
    /// Create an intermediate (fibre↔matrix coupling) node at `coordS`
    /// and register it with the grid.
    Vertex *createInterNode(oofem::FloatArray coordS);

    /// VTK writer for fibre lattice-beam segments.
    void giveBeamElementVTKOutput(FILE *outputStream);
    /// VTK writer for fibre↔matrix lattice-link elements.
    void giveLinkElementVTKOutput(FILE *outputStream);

    /// Localiser helper: return the 1-based ids of Delaunay vertices
    /// inside the axis-aligned tolerance box centred at `coord`.
    oofem::IntArray findDelaunayNodesWithinBox(oofem::FloatArray coord, double TOL);

    /// Numerical Recipes Gaussian elimination solver (`n`×`n` system with
    /// `rhs_num` right-hand sides, in-place on `a`).
    int r8mat_solve(int n, int rhs_num, double a[]);
};

#endif // grid_h
