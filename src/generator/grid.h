#ifndef grid_h
#define grid_h

#include "generatorlistutils.h"

#include "floatarray.h"

//#include "alist.h"
#include "generatordatareader.h"
#include "gridtype.h"
//#include "statecountertype.h"
#include "logger.h"
#include <iostream>

#include "domain.h"
#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <time.h>
 #include <map>
#endif

#define _IFT_Grid_type "type"
#define _IFT_Grid_diam            "diam"
#define _IFT_Grid_maxiter         "maxiter"
#define _IFT_Grid_ranint          "ranint"
#define _IFT_Grid_perflag         "perflag"
#define _IFT_Grid_regflag         "regflag"
#define _IFT_Grid_regtype         "regtype"
#define _IFT_Grid_xyzedges        "xyzedges"
#define _IFT_Grid_aggflag         "aggflag"
#define _IFT_Grid_target          "target"
#define _IFT_Grid_nvertex         "nvertex"
#define _IFT_Grid_ncontrolvertex  "ncontrolvertex"
#define _IFT_Grid_ncurve          "ncurve"
#define _IFT_Grid_nsurface        "nsurface"
#define _IFT_Grid_nregion         "nregion"
#define _IFT_Grid_ninclusion      "ninclusion"
#define _IFT_Grid_nrefinement     "nrefinement"
#define _IFT_Grid_ranflag     "ranflag"

class Vertex;
class Curve;
class Surface;
class Region;
class Inclusion;
class GridLocalizer;
class TimeStep;
class Refinement;

#ifdef __PARALLEL_MODE
class ProcessCommunicator;
class LoadBalancer;
#endif
/**
 * Class and object Grid. Grid contains grid description, or if program runs in parrallel then it contains
 * description of grid associated to particular processor or thread of execution. Generally, it contain and
 * manages lists of Dof managers, elements, boundary conditions, cross sections and materials - these describe
 * the geometry of problem, its constitutive properties and applied boundary conditions. Services for accessing
 * these objects are provided. Grid is attribute of engineering model - which represent type of analysis
 * which should be performed.
 *
 * Grid also provides services for reading its description from
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
    gridType dType;

    /**
     * Associated Engineering model. An abstraction for type of analysis which will be prformed.
     * @see EngngModel
     */
    // EngngModel *engineeringModel;

    /// Grid number
    int number;

    double density;
    //These are the number of edges that are going to be changes into an array
    // xedges, yedges, zedges;
    // Moved to public in order not to retype them in the input files

    //The type of grid used ( if BCC regType = 1, else if FCC regType = 2)
    int regType;

    int maxIter, randomInteger;

    int aggregateFlag;
    int targetAggregates;

    int regularFlag;

  oofem::IntArray periodicityFlag;

    int couplingFlag;

    int randomFlag;

    GridLocalizer *gridLocalizer;
  
public:

  oofem::IntArray xyzEdges;
    double TOL;

    double diameter;


    std::vector<Vertex*> vertexList;
    std::vector<Vertex*> inputVertexList;
    std::vector<Vertex*> controlVertexList;
    std::vector<Curve*> curveList;
    std::vector<Surface*> surfaceList;
    std::vector<Region*> regionList;
    std::vector<Inclusion*> inclusionList;
    std::vector<Refinement*> refinementList;

  
    // /// Vertex list
    // AList< Vertex > *vertexList;

    // /// Input Vertex list
    // AList< Vertex > *inputVertexList;

    // /// Control Vertex list
    // AList< Vertex > *controlVertexList;

    // /// Curve list
    // AList< Curve > *curveList;


    // /// Surface list
    // AList< Surface > *surfaceList;

    // /// Region list
    // AList< Region > *regionList;


    // /// Inclusion list
    // AList< Inclusion > *inclusionList;


    // /// Inclusion list
    // AList< Refinement > *refinementList;


    /**
     * Constructor. Creates empty n-th grid belonging to given ProblemManager
     */
    Grid(int n);  // constructors
    ///  Destructor
    ~Grid();                            // destructor

    /// Returns grid number
    int               giveNumber() { return this->number; }
    /// Returns grid number
    void               setNumber(int nn) { this->number = nn; }

    int Iterations() { return this->maxIter; }


    int giveRandomInteger() { return this->randomInteger; }

    int giveMaximumIterations() { return this->maxIter; }


    Vertex *giveVertex(int n);

    Vertex *giveInputVertex(int n);

    Vertex *giveControlVertex(int n);
    
    Curve *giveCurve(int n);

    Surface *giveSurface(int n);

    Region *giveRegion(int n);

    Inclusion *giveInclusion(int n);

    Refinement *giveRefinement(int n);

    /**
     * Returns receiver's associated spatial localizer.
     */
    GridLocalizer *giveGridLocalizer();

    void generateInputPoints();
      
    void generateControlPoints();

    void generatePeriodicControlPoints();
    
    int generatePeriodicPoints();

    int generateRegularPoints();

    int generateRandomPoints();

    //some directions are with mirroring others with shift. Ideally, this should be the only function later on. 
    int generateMixedPoints();
    
    int generatePoints();


    int instanciateYourself(GeneratorDataReader *dr);

    double ran1(int *idum);

    /// Returns number of vertices.
  int                giveNumberOfVertices() { return generator::size1(vertexList); }

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
    

    /// Returns flag to indicate if mesh generation should be regular
    int                giveRegularFlag() { return this->regularFlag; }


    /// Returns flag to indicate if mesh generation should be regular
  void                givePeriodicityFlag(oofem::IntArray& answer) { answer = this->periodicityFlag;}

    
        /// Returns flag to indicate if mesh generation should be regular
    int                giveRandomFlag() { return this->randomFlag; }



    int                giveCorrespondingCoordinateIndex(int);
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
    

    gridType         giveGridType()        { return dType; }
    /// Sets grid type
    void               setGridType(gridType _dType)         { this->dType = _dType; }

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Grid"; }

  void defineBoundaries(oofem::FloatArray &boundaries);

  double giveDiameter(oofem::FloatArray &coords);
    
  void giveOutput(FILE *outputStream);

};

#endif // grid_h
