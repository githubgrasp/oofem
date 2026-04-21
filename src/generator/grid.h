#ifndef grid_h
#define grid_h

#include "generatorlistutils.h"

#include "floatarray.h"

#include "gridtype.h"
#include "logger.h"
#include <iostream>

#include "domain.h"
#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <time.h>
 #include <map>
#endif

class Vertex;
class Curve;
class Surface;
class Region;
class Inclusion;
class GridLocalizer;
class Refinement;
/**
 * Grid provides services for reading components description from
 * input stream and instantiating corresponding components acordingly. The basic Grid task are following
 * <UL>
 * <LI>
 * Reading its description from input and creating corresponding objects.</LI>
 * <LI>
 * Provides services for accessing its particular components.</LI>
 * <LI>
 * Checking yourself.</LI>
 * </UL>
 */
class Grid
{
private:

    /** Grid type. This determines the dimension of the problem
     */
    int dimension;

    int vertexCount = 0;

    /// Grid number
    int number;

    double density;

    //The type of grid used ( if BCC regType = 1, else if FCC regType = 2)
    int regType;

    int maxIter, randomInteger;

    int regularFlag;

    oofem::IntArray periodicityFlag;

    int randomFlag;

    /// Opt in to writing points.vtk alongside nodes.dat (default off).
    int vtkFlag = 0;

    /// Output filename for the generated point list. Populated either by
    /// the first non-directive line (legacy path) or by the `#@output`
    /// directive (new `#@` path).
    std::string outputFileName;

    GridLocalizer *gridLocalizer;

public:

    /// Number of cells along each axis used by the regular-grid point
    /// generators (BCC / FCC). Populated from `#@xyzedges` when present
    /// and `#@regflag 1`; unused in random-placement mode.
    oofem::IntArray xyzEdges;

    /// Geometric tolerance (typically `1e-6 * diameter`) used in all
    /// float-coordinate comparisons.
    double TOL;

    /// Target nominal spacing between generated points. Set by `#@diam`.
    double diameter;

    /// Final, generated vertex list (the points the generator writes out).
    std::vector< Vertex * >vertexList;

    /// Input vertices read from `#@vertex` directives (bounding-box
    /// corners, explicit seed points).
    std::vector< Vertex * >inputVertexList;

    /// Control vertices read from `#@controlvertex` — named node locations
    /// the converter references for BCs / loads.
    std::vector< Vertex * >controlVertexList;

    /// Curves read from `#@curve`.
    std::vector< Curve * >curveList;

    /// Surfaces read from `#@surface`.
    std::vector< Surface * >surfaceList;

    /// Regions read from `#@prism` / `#@cylinder` / `#@sphere`.
    std::vector< Region * >regionList;

    /// Inclusions read from `#@intersphere` / `#@interfacecylinder` /
    /// `#@interfaceplane` / `#@interfacesurface`.
    std::vector< Inclusion * >inclusionList;

    /// Local refinement boxes read from `#@refineprism`.
    std::vector< Refinement * >refinementList;


    /// Constructor. Creates an empty `n`-th grid.
    Grid(int n);

    /// Destructor. Deletes owned components.
    ~Grid();

    /// Returns the grid number (1-based; currently always 1).
    int giveNumber() { return this->number; }

    /// Overrides the grid number.
    void setNumber(int nn) { this->number = nn; }

    /// Returns the current RNG seed (always negative once the grid is
    /// fully initialised; see `#@ranint` in the user manual).
    int giveRandomInteger() { return this->randomInteger; }

    /// Maximum number of random-placement attempts before aborting.
    int giveMaximumIterations() { return this->maxIter; }

    /// Append a point to `vertexList` at the given coordinates. Used
    /// internally by the per-region / per-surface / per-curve point
    /// generators.
    void addVertex(const oofem::FloatArray &coords);

    /// Dump the generated vertex list as an ASCII VTK PolyData file at
    /// `path`. Triggered by the `#@vtk` directive.
    void exportVTK(const std::string &path);

    /// Returns the `n`-th (1-based) generated vertex, or nullptr.
    Vertex *giveVertex(int n);

    /// Returns the `n`-th (1-based) input vertex (from `#@vertex`).
    Vertex *giveInputVertex(int n);

    /// Returns the `n`-th (1-based) control vertex (from `#@controlvertex`).
    Vertex *giveControlVertex(int n);

    /// Returns the `n`-th (1-based) curve.
    Curve *giveCurve(int n);

    /// Returns the `n`-th (1-based) surface.
    Surface *giveSurface(int n);

    /// Returns the `n`-th (1-based) region.
    Region *giveRegion(int n);

    /// Returns the `n`-th (1-based) inclusion.
    Inclusion *giveInclusion(int n);

    /// Returns the `n`-th (1-based) refinement entry.
    Refinement *giveRefinement(int n);

    /// Returns the receiver's associated spatial localiser (octree).
    GridLocalizer *giveGridLocalizer();

    /// Seed `vertexList` with the `inputVertexList` entries (corner /
    /// bounding-box points). Called before region point generation.
    void generateInputPoints();

    /// Seed `vertexList` with the `controlVertexList` entries (named
    /// support / load points).
    void generateControlPoints();

    /// Emit periodic-image partners of the control vertices that lie on
    /// a periodic boundary. No-op when the grid is non-periodic.
    void generatePeriodicControlPoints();

    /// Drive the per-region periodic point generators. Returns 1 on
    /// success.
    int generatePeriodicPoints();

    /// Drive the per-region regular-grid point generators (BCC / FCC,
    /// selected by `regType`). Returns 1 on success.
    int generateRegularPoints();

    /// Drive the per-region random-placement point generators. Returns
    /// 1 on success.
    int generateRandomPoints();

    /// Drive the per-region mixed-periodic point generators (some axes
    /// periodic, others not).
    int generateMixedPoints();

    /// Top-level point-generation entry point. Runs input / control /
    /// region / inclusion / refinement passes in the right order for
    /// the active periodicity and randomness flags.
    int generatePoints();


    /// Parse a `#@` directive control file into Grid state. Populates
    /// `outputFileName` from a `#@output` directive.
    int readControlRecords(const std::string &controlFile);

    /// Read a packing-format file (as produced by `src/aggregate/`) and
    /// instantiate one `InterfaceSphere` per `sphere` line. `itz` and
    /// `refinement` are applied uniformly to every inclusion read.
    /// Ellipsoid lines are unsupported by the generator and trigger a
    /// warning; fibre lines are silently ignored (they belong to the
    /// converter, not the generator).
    void readInclusionFile(const std::string &path, double itz, double refinement);

    /// The output path the generator should write `nodes.dat` to.
    const std::string &giveOutputFileName() const { return outputFileName; }

    /// Numerical Recipes `ran1` uniform RNG. `idum` is a caller-owned
    /// seed state; negative values reinitialise the generator state.
    double ran1(int *idum);

    /// Returns number of vertices.
    int                giveNumberOfVertices() const { return vertexCount; }

    /// Returns number of input vertices.
    int                giveNumberOfInputVertices() { return generator::size1(inputVertexList); }

    /// Returns number of input vertices.
    int                giveNumberOfControlVertices() { return generator::size1(controlVertexList); }

    /// Returns number of curves.
    int                giveNumberOfCurves() { return generator::size1(curveList); }

    /// Returns number of surfaces
    int                giveNumberOfSurfaces() { return generator::size1(surfaceList); }

    /// Returns number of regions
    int                giveNumberOfRegions() { return generator::size1(regionList); }

    /// Returns number of cross section models in grid
    int                giveNumberOfInclusions() { return generator::size1(inclusionList); }

    /// Returns number of refinements in grid
    int                giveNumberOfRefinements() { return generator::size1(refinementList); }


    /// Returns 1 iff mesh generation should use the regular-grid path
    /// (`#@regflag 1`), 0 for random placement.
    int giveRegularFlag() { return this->regularFlag; }

    /// Copies the 3-entry periodicity flag into `answer` (one entry per
    /// axis: 0 non-periodic, 1 periodic).
    void givePeriodicityFlag(oofem::IntArray &answer) { answer = this->periodicityFlag; }

    /// Returns the random-placement strategy: 0 Bolander, 1 Grassl
    /// (points on boundary), 2 Grassl with periodic partners.
    int giveRandomFlag() { return this->randomFlag; }
    /// Returns 1 iff `#@vtk` was set (opt in to `points.vtk` output).
    int giveVtkFlag() { return this->vtkFlag; }


    /// Map a direction index into the corresponding coordinate index.
    /// Used by the periodic-mesh generators to walk axes in a fixed
    /// order.
    int giveCorrespondingCoordinateIndex(int);

    /**
     *@name Advanced grid manipulation methods.
     */
    //@{
    /// Resizes the internal data structure to accomodate space for _newSize dofManagers
    void resizeVertices(int _newSize);

    /// Resizes the internal data structure to accomodate space for _newSize dofManagers
    void resizeInputVertices(int _newSize);

    /// Resizes the internal data structure to accomodate space for _newSize dofManagers
    void resizeControlVertices(int _newSize);

    /// Resizes the internal data structure to accomodate space for _newSize elements
    void resizeCurves(int _newSize);

    /// Resizes the internal data structure to accomodate space for _newSize cross section models
    void resizeSurfaces(int _newSize);

    /// Resizes the internal data structure to accomodate space for _newSize materials
    void resizeRegions(int _newSize);

    /// Resizes the internal data structure to accomodate space for _newSize materials
    void resizeInclusions(int _newSize);

    /// Resizes the internal data structure to accomodate space for _newSize materials
    void resizeRefinements(int _newSize);



    /// Sets i-th componet. The component will be futher managed and maintained by grid object.
    void setVertex(int i, Vertex *obj);

    /// Sets i-th componet. The component will be futher managed and maintained by grid object.
    void setInputVertex(int i, Vertex *obj);

    /// Sets i-th componet. The component will be futher managed and maintained by grid object.
    void setControlVertex(int i, Vertex *obj);

    /// Sets i-th componet. The component will be futher managed and maintained by grid object.
    void setCurve(int i, Curve *obj);

    /// Sets i-th componet. The component will be futher managed and maintained by grid object.
    void setSurface(int i, Surface *obj);

    /// Sets i-th componet. The component will be futher managed and maintained by grid object.
    void setRegion(int i, Region *obj);

    /// Sets i-th componet. The component will be futher managed and maintained by grid object.
    void setInclusion(int i, Inclusion *obj);

    /// Sets i-th componet. The component will be futher managed and maintained by grid object.
    void setRefinement(int i, Refinement *obj);


    /// Returns the grid's spatial dimension (stored in `dimension`, set
    /// via `#@domain`).
    int giveGridType() { return this->dimension; }

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Grid"; }

    /// Fill `boundaries` with the global bounding box of all regions in
    /// `[xmin xmax ymin ymax zmin zmax]` order.
    void defineBoundaries(oofem::FloatArray &boundaries);

    /// Returns the local target diameter at point `coords`, taking any
    /// active `#@refineprism` entry into account (else returns the
    /// grid-wide `diameter`).
    double giveDiameter(oofem::FloatArray &coords);

    /// Write one `x y z` coordinate triple per generated vertex to
    /// `outputStream`, prefixed with a `3` (dimension) and the vertex
    /// count — the format consumed by `qvoronoi p Fv`.
    void giveOutput(FILE *outputStream);
};

#endif // grid_h
