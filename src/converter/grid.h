//   ********************
//   *** CLASS GRID ***
//   ********************

#ifndef grid_h
#define grid_h

#include "converterlistutils.h"
#include "converterdatareader.h"
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
#define _IFT_Grid_type "grid"
#define _IFT_Grid_macrotype "macrotype"
#define _IFT_Grid_diam "diam"
#define _IFT_Grid_ranint "ranint"
#define _IFT_Grid_regflag "regflag"
#define _IFT_Grid_perflag "perflag"

#define _IFT_Grid_mechmean "mechmean"
#define _IFT_Grid_mechcov "mechcov"
#define _IFT_Grid_mechmax "mechmax"
#define _IFT_Grid_mechmin "mechmin"

#define _IFT_Grid_young "young"
#define _IFT_Grid_gamma1 "gamma1"
#define _IFT_Grid_gamma2 "gamma2"
#define _IFT_Grid_conf "conf"
#define _IFT_Grid_tanbeta "tanbeta"
#define _IFT_Grid_tanphi "tanphi"
#define _IFT_Grid_deltarad "deltarad"

#define _IFT_Grid_nvertex "nvertex"
#define _IFT_Grid_ncurve "ncurve"
#define _IFT_Grid_nsurface "nsurface"
#define _IFT_Grid_nregion "nregion"
#define _IFT_Grid_ninclusion "ninclusion"
#define _IFT_Grid_nfiber "nfibre"

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


struct Node {
    int id;
    double x, y, z;

    int entType = 0;   // T3D geometric entity type
    int entID   = 0;   // T3D geometric entity ID
};

struct Tri {
    int id;
    int n1, n2, n3;
    int entType = 0;  // 3,5,6
    int entID   = 0;
    int entProp = 0; // <-- add this


    // only meaningful when (outType & 8)
    int bndCurveId[ 3 ]   = { 0, 0, 0 };   // edge (n1-n2), (n2-n3), (n3-n1)
    int bndCurveProp[ 3 ] = { 0, 0, 0 };
};

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

struct Edge {
    int n1, n2;     // node ids (sorted)
    int tri1 = -1;  // adjacent triangle index
    int tri2 = -1;
    int tet1 = -1;  // adjacent tetra index
    int tet2 = -1;
};


struct CurveSeg {
    int n1, n2;     // sorted node ids
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
    enum GridType { _3dSM, _3dTM, _3dSMTM, _3dPerSM, _3dPerTM, _3dPerSMTM, _3dWong, _3dCantSM, _3dCantTM, _3dCantExtraTM, _3dCantSMTM, _3dKupfer, _3dImran, _3dNotch };

    GridType gridType;

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

    /// A notch is an axis-aligned box within which matrix elements whose midpoint
    /// falls inside get a softened/alternative material. Populated by the
    /// `#@notch <id> box 6 xmin ymin zmin xmax ymax zmax material <m>` directive.
    struct NotchSpec {
        double xmin, ymin, zmin, xmax, ymax, zmax;
        int material = 3;
    };
    std::vector< NotchSpec >notchSpecs;

    /// A sphere inclusion: elements whose endpoints both fall inside get
    /// `inside` material; elements that straddle the sphere (one endpoint
    /// inside, one outside — the ITZ zone) get `interface` material and an
    /// optional `bodyload` id. Populated by the
    /// `#@sphereinclusion <id> centre 3 x y z radius r itz t inside <mi>
    ///   interface <mif> [bodyload <b>]` directive.
    struct SphereInclusionSpec {
        double cx = 0., cy = 0., cz = 0.;
        double radius = 0., itz = 0.;
        int inside    = 2;
        int interface_ = 3;
        int bodyload = -1;
    };
    std::vector< SphereInclusionSpec >sphereInclusionSpecs;

    /// A straight-axis cylindrical inclusion with ITZ halo. The inclusion
    /// classification uses perpendicular distance from each endpoint to the
    /// (infinite) axis through the two line points. Populated by the
    /// `#@cylinderinclusion <id> line 6 x1 y1 z1 x2 y2 z2 radius r itz t
    ///   inside <mi> interface <mif> [bodyload <b>]` directive.
    struct CylinderInclusionSpec {
        double x1 = 0., y1 = 0., z1 = 0.;
        double x2 = 0., y2 = 0., z2 = 0.;
        double radius = 0., itz = 0.;
        int inside    = 2;
        int interface_ = 3;
        int bodyload = -1;
    };
    std::vector< CylinderInclusionSpec >cylinderInclusionSpecs;

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

    GridLocalizer *delaunayLocalizer;
    GridLocalizer *voronoiLocalizer;
    GridLocalizer *reinforcementLocalizer;

    /// constructor
    Grid(int n);

    ///  Destructor
    ~Grid();


    oofem::FloatArray getX(int nid) const;

    oofem::FloatArray triNormal(int triIndex) const;

    void computeEdgeWidths(double thickness);

    double computeTargetArea(int edgeIndex, double Le) const;

    //This function connects edges to t3d curves, so that later the curve numbering of t3d can be used to control the input.
    void buildCurveSegsFromTris();

    //This function connects triangles to t3d patches, so that patch number can be used to connect patch to triangles.
    void buildBoundaryTrisFromTets();

    void buildEdgeAdjacency3D();

    oofem::FloatArray tetBarycentre(int tetIndex) const;
    oofem::FloatArray faceBarycentre(int triIndex) const;

    void rebuildEntityTris();

    double tetVolume(int tetIndex) const;

    double computePolygonAreaProjected(const oofem::FloatArray &polycoords,
                                       const oofem::FloatArray &xm,
                                       const oofem::FloatArray &r,
                                       const oofem::FloatArray &s) const;


    void buildEdgePolygon3D(int edgeIndex, oofem::FloatArray &polycoords) const;
    void write3DEdgeSection(std::ostream &out, int &eid, const Edge &e, int edgeIndex);


    double segLength(int n1, int n2) const;

    void computeNodalLengthsOnCurve(int curveID, std::vector< double > &L) const;

    oofem::FloatArray triBarycentre(int triIndex) const;


    // function is used to read input from control file to convert mesh from mesh generation output to oofem input
    void readControlRecords();

    std::set< int >collectBCNodes(const BCRequest &bc) const;

    double edgeLength(const Edge &e) const;

    double triArea(int triIndex) const;

    void prepareLiveLoadSets();

    double giveThicknessForEntity(int entType, int entID) const;

    bool isConverterDirective(const std::string &t) const
    {
        return t.rfind("#@", 0) == 0;
    }

    void checkAreaConservation() const;


    struct LoadRequest {
        int entType = -1; // 1 vertex, 2 curve, 3 surface, 5 patch, 6 shell
        int entID   = -1;
        double q    = 0.0;// N/m^2 for tri-entities, N/m for curves, N for vertexint setID   = -1;
        int setID   = -1;
        int ltf     = 1; // default load time function
    };
    std::vector< LoadRequest >loadRequests;

    int entityTypeFromString(const std::string &s) const;


    void computeNodalAreasOnTriEntity(int entType, int entID, std::vector< double > &A) const;
    //  void computeNodalAreas(std::vector<double> &A) const;

    int giveNumber() { return this->number; }

    void writeLiveLoads(std::ostream &out, int &bcID);

    void prepareBCSets();

    void writeBCRecords(std::ostream &out, int &bcID) const;

    /// For a given edge, return the (entType, entID) of the most-specific T3D
    /// geometric entity it belongs to, with priority curve > surface > region.
    /// Returns {0, 0} if no entity is associated.
    std::pair< int, int >entityForEdge(const Edge &e) const;

    /// Resolve the EdgeSpec for an edge, falling back to `defaultSpec` if no
    /// `#@element` directive matches the edge's entity.
    EdgeSpec resolveEdgeSpec(const Edge &e, const EdgeSpec &defaultSpec) const;

    /// For a matrix lattice line with endpoint coords A/B, return the notch
    /// material if the midpoint is inside any `#@notch` box, else `defaultMat`.
    int resolveNotchMaterial(const oofem::FloatArray &A, const oofem::FloatArray &B,
                             int defaultMat) const;

    /// For a matrix lattice line with endpoint coords A/B, resolve the sphere-
    /// inclusion material. Returns the `inside` material if both endpoints fall
    /// within `radius + itz/2` of a sphere centre; the `interface` material
    /// (and populates `bodyloadOut` if set, otherwise -1) if exactly one endpoint
    /// does; otherwise returns `defaultMat`.
    int resolveInclusionMaterial(const oofem::FloatArray &A, const oofem::FloatArray &B,
                                 int defaultMat, int &bodyloadOut) const;

    /// Match each `#@controlvertex` declaration to the nearest Delaunay vertex
    /// and populate `controlNodeIds`. Called after the Delaunay vertex list is loaded.
    void resolveControlVertices();

    void writeGeneratedSets(std::ostream &out) const;

    double giveTol() { return this->TOL; }

    void setNumber(int nn) { this->number = nn; }

    int giveMaximumIterations() { return this->maxIter; }

    double giveDiameter() { return this->diameter; }

    void givePeriodicityFlag(oofem::IntArray &answer) { answer = this->periodicityFlag; }

    int giveRandomInteger() { return this->randomInteger; }

    void giveVoronoiPOVOutput(const std::string &filename);

    void giveDelaunayPOVOutput(const std::string &filename);

    void giveVoronoiCellVTKOutput(FILE *outputStream);

    void giveDelaunayElementVTKOutput(FILE *outputStream);

    void giveDelaunayElementVTKOutput2(FILE *outputStream, int nb_mtx);// added to visualize only a selected material

    void giveDelaunayCrossSectionVTKOutput(FILE *outputStream);

    void giveVoronoiCrossSectionVTKOutput(FILE *outputStream);

    void giveVoronoiElementVTKOutput(FILE *outputStream);

    void resolveGridType(const std::string &name);

    Vertex *giveVertex(int n);

    Curve *giveCurve(int n);

    Surface *giveSurface(int n);

    Region *giveRegion(int n);

    Inclusion *giveInclusion(int n);
    Fibre *giveFibre(int n);


    Vertex *giveDelaunayVertex(int n);

    Vertex *giveVoronoiVertex(int n);

    Line *giveDelaunayLine(int n);



    Line *giveVoronoiLine(int n);
    Line *giveLatticeLink(int n);
    Line * giveLatticeBeam(int n);
    Vertex *giveReinforcementNode(int n);
    Vertex *giveInterNode(int n);

    int generateOutput();

    void sortRandomNumbers(oofem::FloatArray &sortedRandomNumbers, oofem::FloatArray &randomNumbers);

    void createRankTable(oofem::IntArray &rankVector, oofem::FloatArray &randomNumbers);


    /** Order orientation of cross-section vertices to be counterclockwise.
     *  This should only be needed for cross-sections of Voronoi elements with more than 3 vertices.
     *  (Special case of boundary elements) */
    void orderDelaunayCrossSectionVertices(int elementNumber);


    int instanciateYourselfFromQhull(const std::string &controlFile, const char *nodeFileName, const char *voronoiFileName);

    void readQhullControlRecords(const std::string &controlFile);


    bool readT3d(const std::string &fn,
                 std::vector< Node > &nodes,
                 std::vector< Tri > &tris,
                 std::vector< Tet > &tets);

    void buildEdges(const std::vector< Tri > &tris,
                    const std::vector< Tet > &tets,
                    std::vector< Edge > &edges);


    void writeT3dNodesOofem(std::ostream &out);

    void writeT3dElemsOofem(std::ostream &out);

    //Main functions which converts T3D meshes into lattices
    int instanciateYourselfFromT3d(const std::string &t3dFileName, const std::string &controlFileName);

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


    void resizeDelaunayVertices(int _newSize);
    void resizeVoronoiVertices(int _newSize);
    void resizeDelaunayLines(int _newSize);
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


    void setDelaunayVertex(int i, Vertex *obj);
    void setVoronoiVertex(int i, Vertex *obj);
    void setDelaunayLine(int i, Line *obj);
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

    void setFibre(int i, Fibre *obj);
    void setLatticeLink(int i, Line *obj);
    void setLatticeBeam(int i, Line *obj);
    void setReinforcementNode(int i, Vertex *obj);
    void setInterNode(int i, Vertex *obj);


    const char *giveClassName() const { return "Grid"; }

    void giveOutput(const std::string &fileName);
    void giveOutputT3d(const std::string &fileName);

    void giveOofemOutput(const std::string &fileName);

    void giveVtkOutput(const std::string &fileName);
    void giveVtkOutput2(const std::string &fileName, int nb_of_mt); // alternative function to obtain VTK files for each material

    void givePOVOutput(const std::string &fileName);

    void give3DSMOutput(const std::string &fileName);

    void give3DSMTMOutput(const std::string &fileName);

    void give3DPeriodicSMTMOutput(const std::string &fileName);

    void give3DTMOutput(const std::string &fileName);

    void give3DPeriodicSMOutput(const std::string &fileName);

    void give3DPeriodicTMOutput(const std::string &fileName);

    /// Walk fibreList: discretise each fibre, build reinforcement nodes,
    /// lattice beams (fibre segments), and lattice links (fibre↔matrix coupling).
    /// Shared by both ConverterTXTDataReader and qhull-template paths.
    void discretizeFibres();


    //Analysis for Concreep paper
    void give3DWongOutput(const std::string &fileName);

    //SM mesh for Cantilver benchmark of 3D Lattice paper
    void give3DCantileverSMOutput(const std::string &fileName);

    //TM mesh for Cantilver benchmark of 3D Lattice paper for CMAME
    void give3DCantileverTMOutput(const std::string &fileName);

    //TM mesh for extra Cantilver benchmark of 3D Lattice paper for CMAME
    void give3DCantileverTMExtraOutput(const std::string &fileName);

    //SMTM mesh for Cantilver benchmark of 3D Lattice paper for CMAME
    void give3DCantileverSMTMOutput(const std::string &fileName);


    //Random SM mesh for Kupfer experiment
    void give3DKupferOutput(const std::string &fileName);

    //Random SM mesh for Imran experiment
    void give3DImranOutput(const std::string &fileName);

    //Random SM mesh for Notch conference
    void give3DNotchOutput(const std::string &fileName);

    double  dpolyValue(int n, double a[], double x);

    void indexx(int n, double arr[], int indx[]);

    int *ivector(int nl, int nh);

    void free_ivector(int *v, int nl, int nh);

    // tool functions to enable fibre (subobjects) to create Nodes at the scale of the grid
    Vertex *createReinfNode(oofem::FloatArray coordR);
    Vertex *createInterNode(oofem::FloatArray coordS);

    // VTK outputs for fibres
    void giveBeamElementVTKOutput(FILE *outputStream);
    void giveLinkElementVTKOutput(FILE *outputStream);

    oofem::IntArray findDelaunayNodesWithinBox(oofem::FloatArray coord, double TOL); // function created to allow other objects to use the localizer

    int r8mat_solve(int n, int rhs_num, double a[]);
};

#endif // grid_h
