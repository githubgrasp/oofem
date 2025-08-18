//   ********************
//   *** CLASS GRID ***
//   ********************

#ifndef grid_h
#define grid_h

#include "converterlistutils.h"
#include "converterdatareader.h"
#include "convertererror.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#include <time.h>
#include <map>
#include <list>
#endif

#define _IFT_Grid_diam "diam"
#define _IFT_Grid_ranint "ranint"
#define _IFT_Grid_regflag "regflag"
#define _IFT_Grid_perflag "perflag"
#define _IFT_Grid_macrotype "macrotype"
#define _IFT_Grid_coupflag "coupflag"
#define _IFT_Grid_meshtype "meshtype"
#define _IFT_Grid_poremean "poremean"
#define _IFT_Grid_porecov "porecov"
#define _IFT_Grid_poremax "poremax"
#define _IFT_Grid_poremin "poremin"
#define _IFT_Grid_throatmean "throatmean"
#define _IFT_Grid_throatcov "throatcov"
#define _IFT_Grid_throatmax "throatmax"
#define _IFT_Grid_throatmin "throatmin"

#define _IFT_Grid_contactangle "contactangle"
#define _IFT_Grid_surfacetension "surfacetension"

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
class Tetra;
class Inclusion;
class GridLocalizer;
class Fibre;
class BoundarySphere;
class BoundaryCylinder;

/**
 * Description of class grid
 */
class Grid
{
private:

    /** Grid type. This determines the type of input to generate
     */
  enum GridType { _3dSM, _3dTM, _3dSMTM, _3dPerSM, _3dPerTM, _3dPerSMTM, _3dFPZ, _3dFPZFibre, _3dFibreBenchmark, _3dWong, _3dPerPoreTM, _3dPerPoreSM, _3dPerPoreSMTM, _3dCantSM, _3dCantTM, _3dCantExtraTM, _3dCantSMTM, _3dBentoniteSM, _3dBentoniteTM, _3dBentoniteSMTM, _3dSphere, _3dCylinder, _3dTetraSM, _3dPerTetraSM, _3dRCPerSM, _3dRCPer2SM, _3dRCSM, _3dTension, _3dGopSha, _3dKupfer, _3dImran, _3dNotch };

    GridType gridType;

    /// Grid number
    int number;

    // Determines the type of macroscopic element, so that a corresponding number of DOFs is created for the control node
    enum MacroType { _Truss, _Membrane, _Beam, _Plate, _3dVoigt, _3d };
    MacroType macroType;
    const char* boundElemName;
    const char* boundBeamElemName;

    double diameter, density;

    int dimension;

    double TOL;

    oofem::IntArray periodicityFlag;

    int meshType;
    
    int regularFlag;

    int couplingFlag;

    int maxIter, randomInteger;

    double throatMean, throatCOV, throatMax, throatMin;

    double poreMean, poreCOV, poreMax, poreMin;

    double mechMean, mechCOV, mechMax, mechMin;

    //elastic
    double youngModulus, gammaOne, gammaTwo;

    //plastic
    double tanBeta, tanPhi;

    //plastic
    double confinement;

    double deltarad;

    double contactAngle, surfaceTension;

    int periodicNodeCounter;

  //AList< Vertex > *delaunayVertexList;
  std::vector< Vertex * >delaunayVertexList;

  //  AList< Vertex > *voronoiVertexList;
    std::vector< Vertex * >voronoiVertexList;

  //  AList< Vertex > *reinforcementNodeList;
  std::vector< Vertex * >reinforcementNodeList;

  //  AList< Line > *delaunayLineList;
    std::vector< Line * >delaunayLineList;

  //  AList< Line > *voronoiLineList;
    std::vector< Line * >voronoiLineList;
    
  //  AList< Fibre > *fibreList;
    std::vector< Line * >fibreList;

  // AList< Tetra > *delaunayTetraList;
    std::vector< Tetra * >delaunayTetraList;
    
  //    AList< Line > *latticeBeamList;
    std::vector< Line * >latticeBeamList;

  //  AList< Line > *latticeLinkList;
    std::vector< Line * >latticeLinkList;
    
  //    AList< Vertex > *interNodeList;
    std::vector< Vertex * >interNodeList;

    /// Line list
  //    AList< Vertex > *vertexList;
    std::vector< Vertex * >vertexList;
  
    /// Curve list
  //    AList< Curve > *curveList;
    std::vector< Curve * >curveList;

    /// Surface list
  //    AList< Surface > *surfaceList;
    std::vector< Surface * >surfaceList;

    /// Cross section list
  //    AList< Region > *regionList;
    std::vector< Region * >regionList;

  //    AList< Inclusion > *inclusionList;
    std::vector< Inclusion * >inclusionList;

    typedef std :: list< int >nodeContainerType;

public:

    GridLocalizer *delaunayLocalizer;
    GridLocalizer *voronoiLocalizer;
    GridLocalizer *reinforcementLocalizer;


    /**
     * Constructor. Creates empty n-th grid belonging to given ProblemManager
     */
    Grid(int n);  // constructors
    ///  Destructor
    ~Grid();
    // destructor

    /// Returns grid number
    int               giveNumber() { return this->number; }

    double giveTol(){return this->TOL;}

    /// Returns grid number
    void               setNumber(int nn) { this->number = nn; }

    int giveMaximumIterations() { return this->maxIter; }
    double giveDiameter() { return this->diameter; }

  void givePeriodicityFlag(oofem::IntArray& answer) {answer = this->periodicityFlag; }

    int giveRegularFlag(){return this->regularFlag;}

    int giveMeshType(){return this->meshType;}
    
    int giveRandomInteger() { return this->randomInteger; }

    void giveVoronoiPOVOutput(char *filename);

    void giveDelaunayPOVOutput(char *filename);

    void giveVoronoiCellVTKOutput(FILE *outputStream);

    void giveDelaunayElementVTKOutput(FILE *outputStream);

    void giveDelaunayElementVTKOutput2(FILE *outputStream, int nb_mtx);// added to visualize only a selected material

    void giveDelaunayCrossSectionVTKOutput(FILE *outputStream);

    void giveVoronoiCrossSectionVTKOutput(FILE *outputStream);

    void giveTetraElementVTKOutput(FILE *outputStream);

    void giveVoronoiElementVTKOutput(FILE *outputStream);

  void resolveGridType(const std::string &name);

  
    void resolveMacroType(const std::string &name);

    Vertex *giveVertex(int n);

    Curve *giveCurve(int n);

    Surface *giveSurface(int n);

    Region *giveRegion(int n);

    Inclusion *giveInclusion(int n);
    Fibre *giveFibre(int n);


    Vertex *giveDelaunayVertex(int n);

    Vertex *giveVoronoiVertex(int n);

    Tetra *giveDelaunayTetra(int n);
    
    Line *giveDelaunayLine(int n);



    Line *giveVoronoiLine(int n);
    Line *givelatticeLink(int n);
    Line * givelatticeBeam(int n);
    Vertex *giveReinforcementNode(int n);
    Vertex *giveInterNode(int n);
    
    int generateOutput();

    void sortRandomNumbers(oofem::FloatArray &sortedRandomNumbers, oofem::FloatArray &randomNumbers);

    void createRankTable(oofem::IntArray &rankVector, oofem::FloatArray &randomNumbers);


    /** Order orientation of cross-section vertices to be counterclockwise.
     *  This should only be needed for cross-sections of Voronoi elements with more than 3 vertices.
     *  (Special case of boundary elements) */
    void orderDelaunayCrossSectionVertices(int elementNumber);


    int instanciateYourself(ConverterDataReader *dr, char nodeFileName[], char delaunayFileName[], char voronoiFileName[]);

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
    
  int                giveNumberOfFibres() { return converter::size1(fibreList); }

  int                giveNumberOfDelaunayTetras() {  return converter::size1(delaunayTetraList); }
    

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

    void resizeDelaunayTetras(int _newSize);


    void setDelaunayVertex(int i, Vertex *obj);
    void setVoronoiVertex(int i, Vertex *obj);
    void setDelaunayLine(int i, Line *obj);
    void setVoronoiLine(int i, Line *obj);

    void setDelaunayTetra(int i, Tetra *obj);

    
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
    void setlatticeLink(int i,Line *obj);
    void setlatticeBeam(int i,Line *obj);
    void setreinforcementNode(int i,Vertex *obj);
    void setinterNode(int i,Vertex *obj);


    const char *giveClassName() const { return "Grid"; }

    void giveOutput(char *fileName);

    void giveOofemOutput(char *fileName);

    void giveVtkOutput(char *fileName);
    void giveVtkOutput2(char *fileName, int nb_of_mt );// alternative function to obtain VTK files for each material
    
    void giveVtkOutputTetra(char *fileName, int nb_of_mt);
    
    void givePOVOutput(char *fileName);


    void give3DSMOutput(char *fileName);

    void give3DSMTMOutput(char *fileName);

    void give3DPeriodicSMTMOutput(char *fileName);

    void give3DTMOutput(char *fileName);

    void give3DPeriodicSMOutput(char *fileName);

    void give3DPeriodicTMOutput(char *fileName);

    void give3DPeriodicPoreTMOutput(char *fileName);

    void give3DPeriodicPoreSMOutput(char *fileName);

    void give3DPeriodicPoreSMTMOutput(char *fileName);

    void give3DBentoniteSMOutput(char *fileName);

    void give3DBentoniteTMOutput(char *fileName);

    void give3DBentoniteCoupledOutput(char *fileName);

    void give3DFPZOutput(char *fileName);

    void give3DFPZFibreOutput(char *fileName);

    //Analysis to test if single fibre solution can be reproduced
    void give3DFibreBenchmarkOutput(char *fileName);
    
    //Analysis for Concreep paper
    void give3DWongOutput(char *fileName);

    //SM mesh for Cantilver benchmark of 3D Lattice paper
    void give3DCantileverSMOutput(char *fileName);

    //TM mesh for Cantilver benchmark of 3D Lattice paper for CMAME
    void give3DCantileverTMOutput(char *fileName);

    //TM mesh for extra Cantilver benchmark of 3D Lattice paper for CMAME
    void give3DCantileverTMExtraOutput(char *fileName);

    //SMTM mesh for Cantilver benchmark of 3D Lattice paper for CMAME
    void give3DCantileverSMTMOutput(char *fileName);

    void give3DSphereOutput(char *fileName);

    void give3DCylinderOutput(char *fileName);

    //Periodic SM mesh for tetrahedra
    void give3DPeriodicTetraSMOutput(char *fileName);

    //Random SM mesh for tetrahedra
    void  give3DTetraSMOutput(char *fileName);

    //Periodic SM mesh for reinforced concrete (concrete with tetrahedra, reinforcement with beams, bond with link elements)
    void give3DRCPeriodicSMOutput(char *fileName);

    //Alternative periodic SM mesh for reinforced concrete (concrete with tetrahedra, reinforcement with beams, bond with link elements). Hanging nodes are used.
    void give3DRCPeriodicSMOutput2(char *fileName);

    //Random SM mesh for reinforced concrete
    void give3DRCSMOutput(char *fileName);

    //SM mesh direct tension
    void give3DTensionOutput(char *fileName);

    //Random SM mesh for GopSha experiment
    void give3DGopShaOutput(char *fileName);

    //Random SM mesh for Kupfer experiment
    void give3DKupferOutput(char *fileName);

    //Random SM mesh for Imran experiment
    void give3DImranOutput(char *fileName);

    //Random SM mesh for Notch conference
    void give3DNotchOutput(char *fileName);
  
    double normalCdfInverse(double cdf, double a, double b);

    double  normal01CdfInverse(double p);

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
    
    oofem::IntArray findDelaunayNodesWithinBox(oofem::FloatArray coord,double TOL); // function created to allow other objects to use the localizer

    void giveTetrahedronBarycentres(oofem::FloatArray &centres, oofem::FloatArray &tetraCoords, oofem::FloatArray &pointCoords);
    
    double* tetrahedron_barycentric ( double tetra[3*4], double p[3] );

    int r8mat_solve ( int n, int rhs_num, double a[] );
    
};

#endif // grid_h
