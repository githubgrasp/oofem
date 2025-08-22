#include "converterlistutils.h"
#include "grid.h"
#include "ellipsoid.h"
#include "vertex.h"
#include "curve.h"
#include "surface.h"
#include "region.h"
#include "interfacesphere.h"
#include "boundarysphere.h"
#include "interfacecylinder.h"
#include "cylinder.h"
#include "inclusion.h"
#include "fibre.h"
#include "line.h"
#include "tetra.h"
#include "converterdatareader.h"
#include "octreegridlocalizer.h"
#include "error.h"
#include "floatarray.h"
#include "prism.h"
#include "convertererror.h"
#include <sstream>

#ifndef __MAKEDEPEND
 #include <string.h>
 #include <stdarg.h>
 #ifdef HAVE_STRINGS_H
 #endif
 #include <strings.h>
 #include <math.h>
 #include <stdio.h>
 #include <cstdlib>
 #include <ctype.h>
 #include <fstream>
 #include <iostream>
#endif

#define IA 16807
#define IM 2147483647
#define AM ( 1.0 / IM )
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV ( 1 + ( IM - 1 ) / NTAB )
#define EPS 1.2e-7
#define RNMX ( 1.0 - EPS )


Grid::Grid(int i)
{
    delaunayLocalizer      = NULL;
    voronoiLocalizer       = NULL;
    reinforcementLocalizer = NULL;
}

Grid::~Grid()
{
    for (auto v : delaunayLineList) {
        delete v;
    }

    for (auto v : voronoiLineList) {
        delete v;
    }

    for (auto v : latticeBeamList) {
        delete v;
    }

    for (auto v : latticeLinkList) {
        delete v;
    }

    for (auto v : delaunayVertexList) {
        delete v;
    }

    for (auto v : voronoiVertexList) {
        delete v;
    }

    for (auto v : reinforcementNodeList) {
        delete v;
    }

    for (auto v : vertexList) {
        delete v;
    }

    for (auto v : curveList) {
        delete v;
    }

    for (auto v : surfaceList) {
        delete v;
    }
    for (auto v : regionList) {
        delete v;
    }

    for (auto v : inclusionList) {
        delete v;
    }

    for (auto v : fibreList) {
        delete v;
    }

    for (auto v : interNodeList) {
        delete v;
    }

    for (auto v : delaunayTetraList) {
        delete v;
    }

    //Localizers
    delete delaunayLocalizer;
    delete voronoiLocalizer;
    delete reinforcementLocalizer;
}



void Grid::resolveGridType(const std::string &name)
{
    //Find the right grid type
    if ( !strncasecmp(name.c_str(), "3dsm", 4) ) {
        gridType = _3dSM;
    } else if ( !strncasecmp(name.c_str(), "3dtm", 4) ) {
        gridType = _3dTM;
    } else if ( !strncasecmp(name.c_str(), "3dcoupledsmtm", 13) ) {
        gridType = _3dSMTM;
    } else if ( !strncasecmp(name.c_str(), "3dpersm", 7) ) {
        gridType = _3dPerSM;
    } else if ( !strncasecmp(name.c_str(), "3dpertm", 7) ) {
        gridType = _3dPerTM;
    } else if ( !strncasecmp(name.c_str(), "3dpercoupledsmtm", 15) ) {
        gridType = _3dPerSMTM;
    } else if ( !strncasecmp(name.c_str(), "3dfpz", 5) ) {
        gridType = _3dFPZ;
    } else if ( !strncasecmp(name.c_str(), "3dfibrefpz", 10) ) {
        gridType = _3dFPZFibre;
    } else if ( !strncasecmp(name.c_str(), "3dfibrebenchmark", 10) ) {
        gridType = _3dFibreBenchmark;
    } else if ( !strncasecmp(name.c_str(), "3dwong", 6) ) {
        gridType = _3dWong;
    } else if ( !strncasecmp(name.c_str(), "3dperporetm", 11) ) {
        gridType = _3dPerPoreTM;
    } else if ( !strncasecmp(name.c_str(), "3dperporesmtm", 13) ) {
        gridType = _3dPerPoreSMTM;
    } else if ( !strncasecmp(name.c_str(), "3dperporesm", 11) ) {
        gridType = _3dPerPoreSM;
    } else if ( !strncasecmp(name.c_str(), "3dbentonitesm", 13) ) {
        gridType = _3dBentoniteSM;
    } else if ( !strncasecmp(name.c_str(), "3dbentonitetm", 13) ) {
        gridType = _3dBentoniteTM;
    } else if ( !strncasecmp(name.c_str(), "3dbentonitecoupled", 18) ) {
        gridType = _3dBentoniteSMTM;
    } else if ( !strncasecmp(name.c_str(), "3dcantSM", 8) ) {
        gridType = _3dCantSM;
    } else if ( !strncasecmp(name.c_str(), "3dcantTM", 8) ) {
        gridType = _3dCantTM;
    } else if ( !strncasecmp(name.c_str(), "3dcantextratM", 13) ) {
        gridType = _3dCantExtraTM;
    } else if ( !strncasecmp(name.c_str(), "3dcantcoupledsmtm", 13) ) {
        gridType = _3dCantSMTM;
    } else if ( !strncasecmp(name.c_str(), "3dsphere", 8) ) {
        gridType = _3dSphere;
    } else if ( !strncasecmp(name.c_str(), "3dcylinder", 10) ) {
        gridType = _3dCylinder;
    } else if ( !strncasecmp(name.c_str(), "3dpertetrasm", 12) ) {
        gridType = _3dPerTetraSM;
    } else if ( !strncasecmp(name.c_str(), "3dtetrasm", 9) ) {
        gridType = _3dTetraSM;
    } else if ( !strncasecmp(name.c_str(), "3drcpersm", 9) ) {
        gridType = _3dRCPerSM;
    } else if ( !strncasecmp(name.c_str(), "3drcper2sm", 9) ) {
        gridType = _3dRCPer2SM;
    } else if ( !strncasecmp(name.c_str(), "3drcsm", 6) ) {
        gridType = _3dRCSM;
    } else if ( !strncasecmp(name.c_str(), "3dtension", 9) ) {
        gridType = _3dTension;
    } else if ( !strncasecmp(name.c_str(), "3dgopsha", 8) ) {
        gridType = _3dGopSha;
    } else if ( !strncasecmp(name.c_str(), "3dkupfer", 8) ) {
        gridType = _3dKupfer;
    } else if ( !strncasecmp(name.c_str(), "3dimran", 8) ) {
        gridType = _3dImran;
    } else    {
        converter::errorf("Unknown grid type %s\n", name.c_str() );
    }
    return;
}

void Grid::resolveMacroType(const std::string &name)
{
    //Find the right macro element type
    if ( !strncasecmp(name.c_str(), "truss", 5) ) {
        macroType = _Truss;
        boundElemName = "ltrspaceboundarytruss";
        boundBeamElemName = "libeam3dboundarytruss";
    } else if ( !strncasecmp(name.c_str(), "membrane", 8) ) {
        macroType = _Membrane;
        boundElemName = "ltrspaceboundarymembrane";
        boundBeamElemName = "libeam3dboundarymembrane";
    } else if ( !strncasecmp(name.c_str(), "beam", 4) ) {
        macroType = _Beam;
        boundElemName = "ltrspaceboundarybeam";
        boundBeamElemName = "libeam3dboundarybeam";
    } else if ( !strncasecmp(name.c_str(), "plate", 5) ) {
        macroType = _Plate;
        boundElemName = "ltrspaceboundaryplate";
        boundBeamElemName = "libeam3dboundaryplate";
    } else if ( !strncasecmp(name.c_str(), "3dVoigt", 7) ) {
        macroType = _3dVoigt;
        boundElemName = "ltrspaceboundaryvoigt";
        boundBeamElemName = "libeam3dboundaryvoigt";
    } else if ( !strncasecmp(name.c_str(), "3d", 2) ) {
        macroType = _3d;
        boundElemName = "ltrspaceboundary";
        boundBeamElemName = "libeam3dboundary";
    } else {
        macroType = _3d;
        boundElemName = "ltrspaceboundary";
        boundBeamElemName = "libeam3dboundary";
    }
    return;
}

int Grid::instanciateYourself(ConverterDataReader *dr, const char nodeFileName[], const char delaunayFileName[], const char voronoiFileName[])
{
    int i, num;

    Vertex *vertex;
    Curve *curve;
    Surface *surface;
    Region *region;
    InterfaceSphere *interfacesphere;
    Prism *prism;
    BoundarySphere *boundarysphere;
    Fibre *fibre;
    Ellipsoid *ellipsoid;
    InterfaceCylinder *interfacecylinder;
    Cylinder *cylinder;

    Vertex *delaunayVertex;
    Vertex *voronoiVertex;

    Line *delaunayLine;
    Line *voronoiLine;

    Line *linkLine;
    Line *beamLine;

    Tetra *delaunayTetra;

    auto &irDomainRec = dr->giveInputRecord(ConverterDataReader::CIR_domainRec, 1);

    std::string gridTypeName;
    irDomainRec.giveField(gridTypeName, _IFT_Grid_type);
    irDomainRec.finish();

    for (char &c : gridTypeName) {
        c = std::tolower(static_cast < unsigned char > ( c ) );
    }
    resolveGridType(gridTypeName);

    // read control data
    auto &irControlRec = dr->giveInputRecord(ConverterDataReader::CIR_controlRec, 1);


    IR_GIVE_FIELD(irControlRec, diameter, _IFT_Grid_diam); // Macro
    TOL = 1.e-6 * diameter;


    randomInteger = 0;
    IR_GIVE_OPTIONAL_FIELD(irControlRec, randomInteger, _IFT_Grid_ranint); // Macro
    if ( randomInteger >= 0 ) {
        randomInteger = -time(NULL);
    }

    //when regularFlag=0, generator should have placed nodes irregularly
    //when regularFlag=1 regular generator is used.
    //when regularFlag=2, generator with nodes on surfaces and edges was used.
    regularFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(irControlRec, regularFlag, _IFT_Grid_regflag); // Macro
    if ( regularFlag < 0 || regularFlag > 2 ) {
        converter::error("Error: Unknown regular flag. Should be 0, 1 or 2.\n");
    }

    //when periodicityFlag=1 periodic cell generator is used.
    //Periodicity can be enforced in certain directions only.
    periodicityFlag.zero();
    IR_GIVE_OPTIONAL_FIELD(irControlRec, periodicityFlag, _IFT_Grid_perflag); // Macro
    if ( periodicityFlag.giveSize() != 3 ) {
        converter::error("Error: Unknown size of periodicity flag. It needs to be of size 3. If periodicity in all directions is desired then enter \"perflag 3 1 1 1\"\n");
    }

    std::string macroType;

    try {
        if ( irControlRec.hasField(_IFT_Grid_macrotype) ) {
            irControlRec.giveField(macroType, _IFT_Grid_macrotype);
            // normalize to lowercase for case-insensitive matching
            for (char &c : macroType) {
                c = static_cast < char > ( std::tolower(static_cast < unsigned char > ( c ) ) );
            }
            resolveMacroType(macroType);
        } else {
            resolveMacroType("");          // default
        }
    }
    catch(const InputException &e) {
        resolveMacroType("");
    }

    //when couplingFlag=1 coupling is implemented
    couplingFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(irControlRec, couplingFlag, _IFT_Grid_coupflag); // Macro
    if ( couplingFlag < 0 || couplingFlag > 1 ) {
        converter::error("Error: Unknown coupling flag. Should be 0 or 1.\n");
    }


    //meshType: 0 for lattice, 1 for continuum (tetrahedra)
    meshType = 0;
    IR_GIVE_OPTIONAL_FIELD(irControlRec, meshType, _IFT_Grid_meshtype); // Macro
    if ( meshType < 0 || meshType > 1 ) {
        converter::error("Error: Unknown meshtype. Should be\n0 for lattice\n or\n 1 for tetra.\n");
    }

    //Input for the pores
    IR_GIVE_OPTIONAL_FIELD(irControlRec, poreMean, _IFT_Grid_poremean); // Macro
    IR_GIVE_OPTIONAL_FIELD(irControlRec, poreCOV, _IFT_Grid_porecov); // Macro
    IR_GIVE_OPTIONAL_FIELD(irControlRec, poreMax, _IFT_Grid_poremax); // Macro
    IR_GIVE_OPTIONAL_FIELD(irControlRec, poreMin, _IFT_Grid_poremin); // Macro

    //Input for the throats
    IR_GIVE_OPTIONAL_FIELD(irControlRec, throatMean, _IFT_Grid_throatmean); // Macro
    IR_GIVE_OPTIONAL_FIELD(irControlRec, throatCOV, _IFT_Grid_throatcov); // Macro
    IR_GIVE_OPTIONAL_FIELD(irControlRec, throatMax, _IFT_Grid_throatmax); // Macro
    IR_GIVE_OPTIONAL_FIELD(irControlRec, throatMin, _IFT_Grid_throatmin); // Macro

    //Lammat parameters
    this->contactAngle = 0.;
    IR_GIVE_OPTIONAL_FIELD(irControlRec, contactAngle, _IFT_Grid_contactangle); // Macro
    this->surfaceTension = 0.0;
    IR_GIVE_OPTIONAL_FIELD(irControlRec, surfaceTension, _IFT_Grid_surfacetension); // Macro

    //Input for the mechanical elements cross sections
    IR_GIVE_OPTIONAL_FIELD(irControlRec, mechMean, _IFT_Grid_mechmean); // Macro
    IR_GIVE_OPTIONAL_FIELD(irControlRec, mechCOV, _IFT_Grid_mechcov); // Macro
    IR_GIVE_OPTIONAL_FIELD(irControlRec, mechMax, _IFT_Grid_mechmax); // Macro
    IR_GIVE_OPTIONAL_FIELD(irControlRec, mechMin, _IFT_Grid_mechmin); // Macro

    //Input for elastic material
    IR_GIVE_OPTIONAL_FIELD(irControlRec, youngModulus, _IFT_Grid_young); // Macro
    this->gammaOne = 1.;
    IR_GIVE_OPTIONAL_FIELD(irControlRec, gammaOne, _IFT_Grid_gamma1); // Macro

    this->gammaTwo = 1.;
    IR_GIVE_OPTIONAL_FIELD(irControlRec, gammaTwo, _IFT_Grid_gamma2); // Macro

    this->confinement = 0.;
    IR_GIVE_OPTIONAL_FIELD(irControlRec, confinement, _IFT_Grid_conf); // Macro

    //Inut for plastic variables
    this->tanBeta = 0.3;
    IR_GIVE_OPTIONAL_FIELD(irControlRec, tanBeta, _IFT_Grid_tanbeta); // Macro

    this->tanPhi = 0.3;
    IR_GIVE_OPTIONAL_FIELD(irControlRec, tanPhi, _IFT_Grid_tanphi); // Macro

    //input for throat change rate
    IR_GIVE_OPTIONAL_FIELD(irControlRec, deltarad, _IFT_Grid_deltarad); // Macro

    irControlRec.finish();

    // read grid components
    int nvertex, nelem, ncurve, nsurface, nregion, ninclusion, nfibre;
    auto &irDomainCompRec = dr->giveInputRecord(ConverterDataReader::CIR_domainCompRec, 1);
    IR_GIVE_FIELD(irDomainCompRec, nvertex, _IFT_Grid_nvertex); // Macro
    IR_GIVE_FIELD(irDomainCompRec, ncurve, _IFT_Grid_ncurve); // Macro
    IR_GIVE_FIELD(irDomainCompRec, nsurface, _IFT_Grid_nsurface); // Macro
    IR_GIVE_FIELD(irDomainCompRec, nregion, _IFT_Grid_nregion); // Macro
    IR_GIVE_FIELD(irDomainCompRec, ninclusion, _IFT_Grid_ninclusion); // Macro
    IR_GIVE_OPTIONAL_FIELD(irDomainCompRec, nfibre, _IFT_Grid_nfiber); // Macro

    irDomainCompRec.finish();

    //=========================
    //Instanciate localizers
    //=========================

    if ( delaunayLocalizer == NULL ) {
        delaunayLocalizer = new OctreeGridLocalizer(1, this, 0);
    }

    if ( voronoiLocalizer == NULL ) {
        voronoiLocalizer = new OctreeGridLocalizer(1, this, 1);
    }

    if ( reinforcementLocalizer == NULL ) {
        reinforcementLocalizer = new OctreeGridLocalizer(1, this, 2);
    }


    vertexList.resize(nvertex, nullptr);
    for ( i = 0; i < nvertex; i++ ) {
        auto &irVertexRec = dr->giveInputRecord(ConverterDataReader::CIR_vertexRec, i + 1);

        std::string name;
        int num = 0;
        IR_GIVE_RECORD_KEYWORD_FIELD(irVertexRec, name, num);

        if ( num < 1 || num > nvertex ) {
            std::cerr << "instanciateYourself: Invalid vertex number (num=" << num << ")\n";
            std::exit(EXIT_FAILURE);
        }

        if ( converter::includes1(vertexList, num) ) {
            std::cerr << "instanciateYourself: Curve entry already exists (num=" << num << ")\n";
            std::exit(EXIT_FAILURE);
        }

        Vertex *vertex = new Vertex(num, this);
        vertex->initializeFrom(irVertexRec);
        setVertex(num, vertex);
        irVertexRec.finish();
    }




    curveList.resize(ncurve, nullptr);
    for ( i = 0; i < ncurve; i++ ) {
        auto &irCurveRec = dr->giveInputRecord(ConverterDataReader::CIR_curveRec, i + 1);

        std::string name;
        int num = 0;
        IR_GIVE_RECORD_KEYWORD_FIELD(irCurveRec, name, num);

        if ( num < 1 || num > ncurve ) {
            std::cerr << "instanciateYourself: Invalid curve number (num=" << num << ")\n";
            std::exit(EXIT_FAILURE);
        }

        if ( converter::includes1(curveList, num) ) {
            std::cerr << "instanciateYourself: Curve entry already exists (num=" << num << ")\n";
            std::exit(EXIT_FAILURE);
        }

        Curve *curve = new Curve(num, this);
        curve->initializeFrom(irCurveRec);
        setCurve(num, curve);
        irCurveRec.finish();
    }

    surfaceList.resize(nsurface, nullptr);
    for ( i = 0; i < nsurface; i++ ) {
        auto &irSurfaceRec = dr->giveInputRecord(ConverterDataReader::CIR_surfaceRec, i + 1);
        std::string name;
        int num = 0;
        IR_GIVE_RECORD_KEYWORD_FIELD(irSurfaceRec, name, num);

        if ( num < 1 || num > nsurface ) {
            std::cerr << "instanciateYourself: Invalid curve number (num=" << num << ")\n";
            std::exit(EXIT_FAILURE);
        }

        if ( converter::includes1(surfaceList, num) ) {
            std::cerr << "instanciateYourself: Curve entry already exists (num=" << num << ")\n";
            std::exit(EXIT_FAILURE);
        }

        Surface *surface = new Surface(num, this);
        surface->initializeFrom(irSurfaceRec);
        setSurface(num, surface);
        irSurfaceRec.finish();
    }


    regionList.resize(nregion, nullptr);

    for (int i = 0; i < nregion; ++i) {
        auto &irRegionRec = dr->giveInputRecord(ConverterDataReader::CIR_regionRec, i + 1);

        std::string name;
        int num = 0;
        IR_GIVE_RECORD_KEYWORD_FIELD(irRegionRec, name, num);

        if ( num < 1 || num > nregion ) {
            std::cerr << "instanciateYourself: Invalid region number (num=" << num << ")\n";
            std::exit(EXIT_FAILURE);
        }
        if ( converter::includes1(regionList, num) ) {
            std::cerr << "instanciateYourself: Region entry already exists (num=" << num << ")\n";
            std::exit(EXIT_FAILURE);
        }

        Region *r = nullptr;
        if ( name == "prism" || name == "region" ) {
            r = new Prism(num, this);
        } else if ( name == "cylinder" ) {
            r = new Cylinder(num, this);
        } else if ( name == "bsphere" ) {
            r = new BoundarySphere(num, this);
        } else {
            std::cerr << "instanciateYourself: Unknown region type '" << name << "'\n";
            std::exit(EXIT_FAILURE);
        }

        r->initializeFrom(irRegionRec);
        setRegion(num, r);
        irRegionRec.finish();
    }

    inclusionList.resize(ninclusion, nullptr);

    for (int i = 0; i < ninclusion; ++i) {
        auto &irInclusionRec = dr->giveInputRecord(ConverterDataReader::CIR_inclusionRec, i + 1);

        std::string name;
        int num = 0;
        IR_GIVE_RECORD_KEYWORD_FIELD(irInclusionRec, name, num);

        // 1) validate number
        if ( num < 1 || num > ninclusion ) {
            std::cerr << "instanciateYourself: Invalid inclusion number (num=" << num << ")\n";
            std::exit(EXIT_FAILURE);
        }

        // 2) duplicate check
        if ( converter::includes1(inclusionList, num) ) {
            std::cerr << "instanciateYourself: Inclusion entry already exists (num=" << num << ")\n";
            std::exit(EXIT_FAILURE);
        }

        // 3) create the right type
        Inclusion *inc = nullptr;
        if ( name == "intersphere" ) {
            inc = new InterfaceSphere(num, this);
        } else if ( name == "interfacecylinder" ) {
            inc = new InterfaceCylinder(num, this);
        } else if ( name == "ellipsoid" ) {
            inc = new Ellipsoid(num, this);
        } else {
            std::cerr << "instanciateYourself: Unknown inclusion type '" << name << "'\n";
            std::exit(EXIT_FAILURE);
        }

        // 4) init + store
        inc->initializeFrom(irInclusionRec);
        setInclusion(num, inc);

        irInclusionRec.finish();
    }

    //====================================
    // read delaunayVertices
    //====================================
    std::ifstream vertexField(nodeFileName);

    if ( !vertexField.is_open() ) {
        std::cout << "In grid.C: Unable to open file " << nodeFileName << "\n";
        std::exit(1);
    }


    //Read vertex file
    int junk, nDelaunayVertices, nDelaunayTetras;
    oofem::FloatArray coords(3);
    oofem::FloatArray coordsOne(3), coordsTwo(3);


    vertexField.precision(16);

    vertexField >> junk;
    vertexField >> nDelaunayVertices;


    delaunayVertexList.resize(nDelaunayVertices, nullptr);
    for (int i = 0; i < nDelaunayVertices; ++i) {
        double x, y, z;
        if ( !( vertexField >> x >> y >> z ) ) {
            std::cerr << "instanciateYourself: failed to read coordinates for Delaunay vertex "
                      << ( i + 1 ) << "\n";
            std::exit(EXIT_FAILURE);
        }

        if ( converter::includes1(delaunayVertexList, i + 1) ) {
            std::cerr << "instanciateYourself: DelaunayVertex entry already exists (num="
                      << ( i + 1 ) << ")\n";
            std::exit(EXIT_FAILURE);
        }

        // create and store vertex
        auto *v = new Vertex(i + 1, this);
        oofem::FloatArray coords(3);
        coords.at(1) = x;
        coords.at(2) = y;
        coords.at(3) = z;
        v->setCoordinates(coords);
        setDelaunayVertex(i + 1, v);
    }

    delaunayLocalizer->init(true);

    printf("Finished Delaunay vertices\n");


    //=====================================
    // read Delaunay tetrahedra
    //====================================
    // This part is needed for tetrahedral models but not for lattice meshes

    printf("meshType = %d\n", meshType);

    if ( meshType == 1 ) {
        std::ifstream delaunayField(delaunayFileName);
        if ( !delaunayField.is_open() ) {
            std::cout << "In grid.C: Unable to open file " << delaunayFileName << "\n";
            std::exit(1);
        }


        // First value in delaunay.dat is the number of tets
        int nDelaunayTetras = 0;
        if ( !( delaunayField >> nDelaunayTetras ) ) {
            std::cerr << "Failed to read number of Delaunay tetrahedra.\n";
            std::exit(EXIT_FAILURE);
        }

        // 1‑based storage with nullptr padding
        delaunayTetraList.resize(nDelaunayTetras, nullptr);

        oofem::IntArray delaunayVertices(4);
        for (int i = 0; i < nDelaunayTetras; ++i) {
            int a, b, c, d;
            if ( !( delaunayField >> a >> b >> c >> d ) ) {
                std::cerr << "Failed to read vertex indices for Delaunay tetra " << ( i + 1 ) << ".\n";
                std::exit(EXIT_FAILURE);
            }

            // input is 0‑based; convert to 1‑based for internal use
            delaunayVertices.at(1) = a + 1;
            delaunayVertices.at(2) = b + 1;
            delaunayVertices.at(3) = c + 1;
            delaunayVertices.at(4) = d + 1;

            // create and store the tetra
            auto *t = new Tetra(i + 1, this);
            t->setLocalVertices(delaunayVertices);
            setDelaunayTetra(i + 1, t);

            // register connectivity on each referenced vertex
            for (int k = 1; k <= 4; ++k) {
                const int vid = delaunayVertices.at(k);
                auto *v = giveDelaunayVertex(vid);
                if ( !v ) {
                    std::cerr << "Delaunay tetra " << ( i + 1 )
                              << " references missing Delaunay vertex " << vid << ".\n";
                    std::exit(EXIT_FAILURE);
                }
                v->setLocalTetra(i + 1);
            }
        }
    }

    // =========================================
    // read Voronoi vertices
    // =========================================

    std::ifstream voronoiField(voronoiFileName);
    if ( !voronoiField.is_open() ) {
        converter::errorf("In grid.C: Unable to open file %s", voronoiFileName);
    }
    voronoiField.precision(16);

    junk = 0;
    int nVoronoiVertices = 0;

    // First two header lines in voronoi.dat
    if ( !( voronoiField >> junk ) ) {
        converter::error("Voronoi file: failed to read header line 1");
    }
    if ( !( voronoiField >> nVoronoiVertices ) || nVoronoiVertices < 0 ) {
        converter::error("Voronoi file: invalid vertex count on header line 2");
    }

    voronoiVertexList.resize(nVoronoiVertices, nullptr);

    for (int i = 0; i < nVoronoiVertices; ++i) {
        if ( !( voronoiField >> coords.at(1) >> coords.at(2) >> coords.at(3) ) ) {
            converter::errorf("Voronoi file: unexpected EOF reading vertex %d/%d",
                              i + 1, nVoronoiVertices);
        }

        const int id = i + 1; // 1-based id

        if ( converter::includes1(voronoiVertexList, id) ) {
            converter::errorf("Voronoi vertex %d duplicated", id);
        }

        auto * v = new Vertex(id, this);
        v->setCoordinates(coords);
        setVoronoiVertex(id, v);
    }

    if ( voronoiLocalizer ) {
        voronoiLocalizer->init(true);
    }

    std::printf("Finished Voronoi vertices (%d)\n", nVoronoiVertices);

    //==========================
    //Read Delaunay lines
    //==========================
    int nDelaunayLines;
    voronoiField >> nDelaunayLines; //1st line after the Voronoi nodes Coordinates

    int newSize, outsideFlag, boundaryTypeFlag = 0;
    oofem::FloatArray temp;

    int size;
    oofem::IntArray delaunayNodes(2);
    double boundaryDiameter, radius;
    oofem::FloatArray centreLine;
    int realNumber = 0;
    int edgeFlagCounter = 0;
    double area = 0;


    delaunayLineList.resize(nDelaunayLines, nullptr);

    voronoiLineList.reserve(100 * nDelaunayLines);

    int voronoiLineCounter = 0;

    for (int i = 0; i < nDelaunayLines; ++i) {
        int size = 0;
        voronoiField >> size;

        oofem::IntArray delaunayNodes(2);
        voronoiField >> delaunayNodes.at(1) >> delaunayNodes.at(2);

        delaunayNodes.at(1) += 1;
        delaunayNodes.at(2) += 1;


        auto * delaunayLine = new Line(i + 1, this);
        delaunayLine->setVertices(delaunayNodes);

        // Register this line at both Delaunay vertices
        this->giveDelaunayVertex(delaunayNodes.at(1) )->setLocalLine(i + 1);
        this->giveDelaunayVertex(delaunayNodes.at(2) )->setLocalLine(i + 1);

        // ------------------------
        // Read Voronoi nodes list
        // ------------------------
        // 'size' counts: 2 (delaunay nodes) + NvoronoiNodes
        const int nVorNodes = size - 2;
        oofem::IntArray voronoiNodes(nVorNodes);
        for (int k = 0; k < nVorNodes; ++k) {
            voronoiField >> voronoiNodes.at(k + 1);
            // Note: Voronoi node indices are already 1-based in Qhull, 0 means at infinity.
        }

        // Stash polygon on the Delaunay line for cross-section
        delaunayLine->updateCrossSectionVertices(voronoiNodes);

        // --------------------------------
        // Build Voronoi edges around face
        // --------------------------------
        oofem::IntArray nodesA(2);
        for (int m = 0; m < nVorNodes; ++m) {
            // consecutive pair, wrapping at the end
            nodesA.at(1) = voronoiNodes.at(m + 1);
            nodesA.at(2) = ( m < nVorNodes - 1 ) ? voronoiNodes.at(m + 2) : voronoiNodes.at(1);

            // Collect candidate Voronoi lines touching one of the endpoints (prefer a non-zero node)
            oofem::IntArray localVoronoiLines;
            if ( nodesA.at(1) != 0 ) {
                this->giveVoronoiVertex(nodesA.at(1) )->giveLocalLines(localVoronoiLines);
            } else if ( nodesA.at(2) != 0 ) {
                this->giveVoronoiVertex(nodesA.at(2) )->giveLocalLines(localVoronoiLines);
            } else {
                std::fprintf(stderr, "error: cannot have two zero Voronoi nodes\n");
                std::exit(1);
            }

            // See if this Voronoi edge already exists (order-insensitive match)
            bool exists = false;
            for (int k = 0; k < localVoronoiLines.giveSize(); ++k) {
                const int lid = localVoronoiLines.at(k + 1);
                oofem::IntArray localVertices;
                this->giveVoronoiLine(lid)->giveLocalVertices(localVertices);

                const bool same =
                    ( localVertices.at(1) == nodesA.at(1) && localVertices.at(2) == nodesA.at(2) ) ||
                    ( localVertices.at(1) == nodesA.at(2) && localVertices.at(2) == nodesA.at(1) );

                if ( same ) {
                    exists = true;
                    // Update cross-section coupling both ways
                    this->giveVoronoiLine(lid)->updateCrossSectionVertices(delaunayNodes);
                    this->giveVoronoiLine(lid)->updateCrossSectionElement(i + 1);
                    delaunayLine->updateCrossSectionElement(lid);
                    break;
                }
            }

            // Create new Voronoi line if it doesn't exist
            if ( !exists ) {
                const int newId = ++voronoiLineCounter;
                auto * vorLine = new Line(newId, this);
                vorLine->setVertices(nodesA);

                // Couple to the Delaunay line (cross-section)
                vorLine->updateCrossSectionVertices(delaunayNodes);
                vorLine->updateCrossSectionElement(i + 1);
                delaunayLine->updateCrossSectionElement(newId);

                converter::put1_replace(voronoiLineList, newId, vorLine);

                // Register this Voronoi line at its endpoints (if finite)
                if ( nodesA.at(1) != 0 ) {
                    this->giveVoronoiVertex(nodesA.at(1) )->setLocalLine(newId);
                }
                if ( nodesA.at(2) != 0 ) {
                    this->giveVoronoiVertex(nodesA.at(2) )->setLocalLine(newId);
                }
            }
        }

        // Finally store the Delaunay line (1-based)
        converter::put1_replace(delaunayLineList, i + 1, delaunayLine);
    }


    std::printf("Finished Delaunay and Voronoi lines\n");




    //Go through all Delaunay vertices and set Voronoi nodes belonging to the cell.
    //This is needed later for VTK files. However, polyhedra are still not supported in VTK files
    oofem::IntArray localLines;
    localLines.zero();
    oofem::IntArray crossSectionNodes;
    crossSectionNodes.zero();
    for ( int i = 0; i < nDelaunayVertices; i++ ) {
        this->giveDelaunayVertex(i + 1)->giveLocalLines(localLines);
        for ( int m = 0; m < localLines.giveSize(); m++ ) {
            this->giveDelaunayLine(localLines.at(m + 1) )->giveCrossSectionVertices(crossSectionNodes);
            this->giveDelaunayVertex(i + 1)->updateCellVertices(crossSectionNodes);
        }
    }


    //Go through all Voronoi vertices and set Delaunay elements belonging to the cell.
    //This is needed for some coupling approaches.
    oofem::IntArray crossSectionElements;
    crossSectionElements.zero();
    for ( int i = 0; i < nVoronoiVertices; i++ ) {
        this->giveVoronoiVertex(i + 1)->giveLocalLines(localLines);
        for ( int m = 0; m < localLines.giveSize(); m++ ) {
            this->giveVoronoiLine(localLines.at(m + 1) )->giveCrossSectionElements(crossSectionElements);
            this->giveVoronoiVertex(i + 1)->updateCellElements(crossSectionElements);
        }
    }


    fibreList.resize(nfibre, nullptr);

    std::printf("\n number of fibres detected : %d \n", ( int ) fibreList.size() );

    for (int i = 0; i < nfibre; ++i) {
        auto &irFibreRec = dr->giveInputRecord(ConverterDataReader::CIR_fibreRec, i + 1);

        std::string kw;
        int num = 0;
        irFibreRec.giveRecordKeywordField(kw, num);

        if ( num < 1 || num > nfibre ) {
            converter::errorf("instanciateYourself: Invalid fibre number (num=%d)", num);
        }
        if ( converter::includes1(fibreList, num) ) {
            converter::errorf("instanciateYourself: Fibre entry already exists (num=%d)", num);
        }

        auto *fibre = new Fibre(num, this);
        fibre->initializeFrom(irFibreRec);

        converter::put1_replace(fibreList, num, fibre);
    }

    //==================================
    // generation of the Lines between Reinforcements nodes (beam elements) , and between Reinf and Del vertices (links)
    //=================================

    int numberOfBeams, numberOfLinks, numberOfReinforcementNodes;
    int indexOfBeamElements, indexOfLinkElements, indexOfReinforcementNodes;
    int globalIndex = 0;
    int beamElementCounter = 0, linkElementCounter = 0;
    oofem::IntArray beamNodes(2), linkNodes(2);
    double portionOfFibre;
    oofem::FloatArray coordP1, coordP2;
    double fibre_diameter;
    oofem::FloatArray dir_vector;

    printf("\n generation of beam elements for fibres and link elements in progress... \n ");

    for ( i = 1; i <= nfibre; i++ ) {
        // 0) collect info about fibre which will be added to elements
        fibre_diameter = giveFibre(i)->giveDiameter();
        dir_vector = giveFibre(i)->giveDirVector();

        // 1) discretization and creation of reinforcement nodes
        /**TODO: This has been written for the lattice generation.
         * Therefore, the nodes of the rebars are placed in the centre of the section crossing
         * the Voronoi cell. For meshtype=1 where tetras are used it would be better to place
         * the node at the intersection with Delaunay tetrahedra of the segment crossing the delaunay tetra.
         **/
        giveFibre(i)->discretizeYouself();

        // 2) creation of link and beam elements associated to the fibre

        numberOfReinforcementNodes = giveFibre(i)->NbOfReinfNodes();
        numberOfBeams = ( giveFibre(i)->NbOfReinfNodes() ) - 1;
        numberOfLinks = ( giveFibre(i)->NbOfReinfNodes() );

        indexOfLinkElements = converter::size1(latticeLinkList);
        indexOfBeamElements = converter::size1(latticeBeamList);
        indexOfReinforcementNodes = converter::size1(reinforcementNodeList);

        latticeLinkList.resize(indexOfLinkElements + numberOfLinks, nullptr);
        latticeBeamList.resize(indexOfBeamElements + numberOfBeams, nullptr);


        for (int j = 1; j <= numberOfBeams; j++ ) {
            beamElementCounter++;

            beamLine = ( Line * ) ( Line(beamElementCounter + 1, this).ofType() );
            setLatticeBeam(beamElementCounter, beamLine);

            beamNodes.at(1) = giveFibre(i)->giveNumberReinforcementNode(j);
            beamNodes.at(2) = giveFibre(i)->giveNumberReinforcementNode(j + 1);

            beamLine->setVertices(beamNodes);
            this->giveReinforcementNode(beamNodes.at(1) )->setLocalLine(beamElementCounter);
            this->giveReinforcementNode(beamNodes.at(2) )->setLocalLine(beamElementCounter);

            beamLine->setDiameter(fibre_diameter);
            beamLine->setDirVector(dir_vector);
        }

        for (int j = 1; j <= numberOfReinforcementNodes; j++ ) {
            linkElementCounter++;
            linkLine = ( Line * ) ( Line(linkElementCounter + 1, this).ofType() );
            setLatticeLink(linkElementCounter, linkLine);

            linkNodes.at(1) = giveFibre(i)->giveNumberReinforcementNode(j);
            linkNodes.at(2) = giveFibre(i)->giveNumberDelaunayNode(j);

            linkLine->setVertices(linkNodes);

            //set of the length of fibre associated to the link
            this->giveInterNode(giveFibre(i)->giveNumberIntersectionPoint(j) )->giveCoordinates(coordP1);
            this->giveInterNode(giveFibre(i)->giveNumberIntersectionPoint(j + 1) )->giveCoordinates(coordP2);
            portionOfFibre = Fibre::computedistance(coordP1, coordP2);
            linkLine->setAssociatedLength(portionOfFibre);

            this->giveReinforcementNode(linkNodes.at(1) )->setLocalLink(linkElementCounter);
            this->giveDelaunayVertex(linkNodes.at(2) )->setLocalLink(linkElementCounter);

            linkLine->setDiameter(fibre_diameter);
            linkLine->setDirVector(dir_vector);

            linkLine->setL_end(giveFibre(i)->giveL_end(j) );
        }
    }

    // reinforcementLocalizer->init(true);

    printf("finished initializing\n");

    return 1;
}

Vertex *Grid::createReinfNode(oofem::FloatArray coordR)
// function to create reinforcement nodes, directly with global index in the grid (not only for the fibre...)
{
    int index(this->giveNumberOfReinforcementNode() + 1);
    reinforcementNodeList.resize(index, nullptr);

    Vertex *reinforcementNode;
    reinforcementNode = ( Vertex * ) ( Vertex(index, this).ofType() );
    reinforcementNode->setCoordinates(coordR);
    setReinforcementNode(index, reinforcementNode);

    return reinforcementNode;
};

Vertex *Grid::createInterNode(oofem::FloatArray coordS)
// function to create reinforcement nodes, directly with global index in the grid (not only for the fibre...)
{
    int index(this->giveNumberOfInterNodes() + 1);
    interNodeList.resize(index, nullptr);

    Vertex *interNode;
    interNode = ( Vertex * ) ( Vertex(index, this).ofType() );
    interNode->setCoordinates(coordS);
    setInterNode(index, interNode);
    return interNode;
};



Vertex * Grid::giveVertex(int n) {
    return converter::require_at1(vertexList, n, "giveVertex");
}


Curve * Grid::giveCurve(int n) {
    return converter::require_at1(curveList, n, "giveCurve");
}


Surface * Grid::giveSurface(int n) {
    return converter::require_at1(surfaceList, n, "giveSurface");
}


Region * Grid::giveRegion(int n) {
    return converter::require_at1(regionList, n, "giveRegion");
}


Inclusion * Grid::giveInclusion(int n) {
    return converter::require_at1(inclusionList, n, "giveInclusion");
}


Fibre * Grid::giveFibre(int n) {
    return converter::require_at1(fibreList, n, "giveFibre");
}


Line * Grid::giveDelaunayLine(int n) {
    return converter::require_at1(delaunayLineList, n, "giveDelaunayLine");
}


Tetra * Grid::giveDelaunayTetra(int n) {
    return converter::require_at1(delaunayTetraList, n, "giveDelaunayTetra");
}


Line * Grid::giveVoronoiLine(int n) {
    return converter::require_at1(voronoiLineList, n, "giveVoronoiLine");
}


Vertex * Grid::giveDelaunayVertex(int n) {
    return converter::require_at1(delaunayVertexList, n, "giveDelaunayVertex");
}


Vertex * Grid::giveVoronoiVertex(int n) {
    return converter::require_at1(voronoiVertexList, n, "giveVoronoiVertex");
}

Vertex * Grid::giveReinforcementNode(int n) {
    return converter::require_at1(reinforcementNodeList, n, "giveReinforcementNode");
}

Line * Grid::giveLatticeBeam(int n) {
    return converter::require_at1(latticeBeamList, n, "giveLatticeBeam");
}

Line * Grid::giveLatticeLink(int n) {
    return converter::require_at1(latticeLinkList, n, "giveLatticeLink");
}

Vertex * Grid::giveInterNode(int n) {
    return converter::require_at1(interNodeList, n, "giveInterNode");
}

// --- Resize helpers ---
void Grid::resizeDelaunayLines(int newSize) {
    converter::ensure_size1(delaunayLineList, newSize);
}

void Grid::resizeVoronoiLines(int newSize) {
    converter::ensure_size1(voronoiLineList, newSize);
}

void Grid::resizeDelaunayVertices(int newSize) {
    converter::ensure_size1(delaunayVertexList, newSize);
}

void Grid::resizeVoronoiVertices(int newSize) {
    converter::ensure_size1(voronoiVertexList, newSize);
}

// --- Set helpers ---
void Grid::setDelaunayVertex(int i, Vertex *obj) {
    converter::put1(delaunayVertexList, i, obj);
}

void Grid::setVoronoiVertex(int i, Vertex *obj) {
    converter::put1(voronoiVertexList, i, obj);
}

void Grid::setVoronoiLine(int i, Line *obj) {
    converter::put1(voronoiLineList, i, obj);
}

void Grid::setDelaunayLine(int i, Line *obj) {
    converter::put1(delaunayLineList, i, obj);
}

void Grid::setDelaunayTetra(int i, Tetra *obj) {
    converter::put1(delaunayTetraList, i, obj);
}


void Grid::setVertex(int i, Vertex *obj) {
    converter::put1(vertexList, i, obj);
}

void Grid::setCurve(int i, Curve *obj) {
    converter::put1(curveList, i, obj);
}

void Grid::setSurface(int i, Surface *obj) {
    converter::put1(surfaceList, i, obj);
}

void Grid::setRegion(int i, Region *obj) {
    converter::put1(regionList, i, obj);
}

void Grid::setInclusion(int i, Inclusion *obj) {
    converter::put1(inclusionList, i, obj);
}

void Grid::setFibre(int i, Fibre *obj) {
    converter::put1(fibreList, i, obj);
}

void Grid::setReinforcementNode(int i, Vertex *obj) {
    converter::put1(reinforcementNodeList, i, obj);
}

void Grid::setLatticeBeam(int i, Line *obj) {
    converter::put1(latticeBeamList, i, obj);
}

void Grid::setLatticeLink(int i, Line *obj) {
    converter::put1(latticeLinkList, i, obj);
}

void Grid::setInterNode(int i, Vertex *obj) {
    converter::put1(interNodeList, i, obj);
}


int Grid::generateOutput()
{
    oofem::FloatArray boundaries;
    //Could be extended to multiple boundaries later
    this->giveRegion(1)->defineBoundaries(boundaries);
    this->giveRegion(1)->findOutsiders(boundaries);

    return 1;
}


void
Grid::orderDelaunayCrossSectionVertices(int elementNumber)
{
    oofem::FloatArray coordsA(3), coordsB(3);
    oofem::IntArray vertices, crossSectionVertices;
    this->giveVoronoiLine(elementNumber)->giveLocalVertices(vertices);
    this->giveVoronoiLine(elementNumber)->giveCrossSectionVertices(crossSectionVertices);
    int size = crossSectionVertices.giveSize();

    if ( size <= 3 ) {
        return;
    }

    //Generate local coordinate system. Normal is given by axis of element.

    for ( int i = 0; i < 3; i++ ) {
        coordsA.at(i + 1) =  this->giveVoronoiVertex(vertices.at(1) )->giveCoordinate(i + 1);
        coordsB.at(i + 1) =  this->giveVoronoiVertex(vertices.at(2) )->giveCoordinate(i + 1);
    }

    //Construct an initial temporary local coordinate system
    oofem::FloatArray n(3), s(3), t(3);

    for ( int i = 0; i < 3; i++ ) {
        n.at(i + 1) = coordsB.at(i + 1) - coordsA.at(i + 1);
    }

    // Compute midpoint
    oofem::FloatArray midPoint(3);
    for ( int i = 0; i < 3; i++ ) {
        midPoint.at(i + 1) = 0.5 * ( coordsB.at(i + 1) + coordsA.at(i + 1) );
    }

    double length  = sqrt(pow(n.at(1), 2.) + pow(n.at(2), 2.) + pow(n.at(3), 2.) );

    if ( length < 1.e-20 ) {
        printf("too small length. Cannot fix orientation\n");
        return;
    }

    for ( int i = 0; i < 3; i++ ) {
        n.at(i + 1) /= length;
    }

    //Create t and s
    if ( n.at(1) == 0 ) {
        s.at(1) = 0.;
        s.at(2) = n.at(3);
        s.at(3) = -n.at(2);
    } else if ( n.at(2) == 0 ) {
        s.at(1) = n.at(3);
        s.at(2) = 0.;
        s.at(3) = -n.at(1);
    } else {
        s.at(1) = n.at(2);
        s.at(2) = -n.at(1);
        s.at(3) = 0.;
    }

    s.normalize();

    t.beVectorProductOf(n, s);
    t.normalize();

    //Set up rotation matrix
    oofem::FloatMatrix lcs(3, 3);

    for ( int i = 1; i <= 3; i++ ) {
        lcs.at(1, i) = n.at(i);
        lcs.at(2, i) = s.at(i);
        lcs.at(3, i) = t.at(i);
    }

    //Calculate the local coordinates of the polygon vertices
    oofem::FloatArray help(3), test(3);
    oofem::FloatArray lpc(3 *size);
    for ( int k = 0; k < size; k++ ) {
        for ( int n = 0; n < 3; n++ ) {
            help(n) = this->giveDelaunayVertex(crossSectionVertices.at(k + 1) )->giveCoordinate(n + 1);
        }

        test.beProductOf(lcs, help);
        for ( int n = 0; n < 3; n++ ) {
            lpc(3 * k + n) = test(n);
        }
    }

    //Check now the order

    int j, k, count = 0;
    double z;
    oofem::FloatArray tempCoords(3);
    int temp;

    for ( int i = 0; i < size; i++ ) {
        j = ( i + 1 ) % size;
        k = ( i + 2 ) % size;
        z  = ( lpc.at(3 * j + 2) - lpc.at(3 * i + 2) ) * ( lpc.at(3 * k + 3) - lpc.at(3 * j + 3) );
        z -= ( lpc.at(3 * j + 3) - lpc.at(3 * i + 3) ) * ( lpc.at(3 * k + 2) - lpc.at(3 * j + 2) );
        if ( z < 0 && count <= 0 ) {//clockwise
            count--;
        } else if ( z > 0  && count >= 0 ) {    //counter clockwise
            count++;
        } else if ( (z < 0 && count > 0) || (z > 0 && count < 0) ) {         //detected problem
            //swap data points j and k and start over
            for ( int n = 0; n < 3; n++ ) {
                tempCoords.at(n + 1) = lpc.at(3 * j + n + 1);
                lpc.at(3 * j + n + 1) = lpc.at(3 * k + n + 1);
                lpc.at(3 * k + n + 1) =  tempCoords.at(n + 1);
            }
            temp = crossSectionVertices.at(j + 1);
            crossSectionVertices.at(j + 1) = crossSectionVertices.at(k + 1);
            crossSectionVertices.at(k + 1) = temp;

            i = 0;
            count = 0;
        }
    }

    this->giveVoronoiLine(elementNumber)->setCrossSectionVertices(crossSectionVertices);

    return;
}



void Grid::giveOutput(const std::string &fileName)
{
    printf("\n starting giving ouputs \n");
    giveOofemOutput(fileName);

    if ( gridType == _3dPerTetraSM || gridType == _3dTetraSM || gridType == _3dRCPerSM || gridType == _3dRCPer2SM || gridType == _3dRCSM ) {
        giveVtkOutputTetra(fileName, 3);
    } else   {
        givePOVOutput(fileName);
        giveVtkOutput2(fileName, 3);
    }
}

void Grid::giveOofemOutput(const std::string &fileName)
{
    //Start with oofem output
    printf("starting giving Oofem output... \n");
    if ( gridType == _3dSM ) { //Base implementation
        give3DSMOutput(fileName);
    } else if ( gridType == _3dTM ) { //Base implementation
        give3DTMOutput(fileName);
    } else if ( gridType == _3dSMTM ) { //Base implementation
        give3DSMTMOutput(fileName);
    } else if ( gridType == _3dPerSM ) {  //Base implementation
        give3DPeriodicSMTMOutput(fileName);
    } else if ( gridType == _3dFPZ ) {
        give3DFPZOutput(fileName);
    } else if ( gridType == _3dFPZFibre ) {
        give3DFPZFibreOutput(fileName);
    } else if ( gridType == _3dFibreBenchmark ) {
        give3DFibreBenchmarkOutput(fileName);
    } else if ( gridType == _3dWong )  { //Implementation for microcracking paper
        give3DWongOutput(fileName);
    } else if ( gridType == _3dPerPoreTM ) {  //Implementation for pore scale analysis
        give3DPeriodicPoreTMOutput(fileName);
    } else if ( gridType == _3dPerPoreSMTM ) {  //Implementation for pore scale analysis
        give3DPeriodicPoreSMTMOutput(fileName);
    } else if ( gridType == _3dPerPoreSM ) {  //Implementation for pore scale analysis
        give3DPeriodicPoreSMOutput(fileName);
    } else if ( gridType == _3dBentoniteSM ) {  //Implementation for structural bentonite analysis
        give3DBentoniteSMOutput(fileName);
    } else if ( gridType == _3dBentoniteTM ) {  //Implementation for transport bentonite analysis
        give3DBentoniteTMOutput(fileName);
    } else if ( gridType == _3dBentoniteSMTM ) { //Implementation for coupled bentonite analysis
        give3DBentoniteCoupledOutput(fileName);
    } else if ( gridType == _3dCantSM ) {  //Implementation for 3D coupling paper structural part
        give3DCantileverSMOutput(fileName);
    } else if ( gridType == _3dCantTM ) {  //Implementation for 3D coupling paper Transport part
        give3DCantileverTMOutput(fileName);
    } else if ( gridType == _3dCantExtraTM ) {  //Implementation for 3D coupling paper
        give3DCantileverTMExtraOutput(fileName);
    } else if ( gridType == _3dCantSMTM ) {  //Implementation for 3D coupling paper Coupled
        give3DCantileverSMTMOutput(fileName);
    } else if ( gridType == _3dSphere ) {  //Implementation for 3D sphere (Milan and Domenico)
        give3DSphereOutput(fileName);
    } else if ( gridType == _3dCylinder ) {  //Implementation for 3D cylinder (Milan and Domenico)
        give3DCylinderOutput(fileName);
    } else if ( gridType == _3dTetraSM ) {  //Implementation for Tetrahedra (Adam)
        give3DTetraSMOutput(fileName);
    } else if ( gridType == _3dPerTetraSM ) {  //Implementation for Tetrahedra, periodic (Adam)
        give3DPeriodicTetraSMOutput(fileName);
    } else if ( gridType == _3dRCSM ) { //Implementation for reinforced concrete (Adam)
        give3DRCSMOutput(fileName);
    } else if ( gridType == _3dRCPerSM ) { //Implementation for reinforced concrete, periodic (Adam)
        give3DRCPeriodicSMOutput(fileName);
    } else if ( gridType == _3dRCPer2SM ) { //Alternative implementation for reinforced concrete, periodic (Adam)
        give3DRCPeriodicSMOutput2(fileName);
    } else if ( gridType == _3dTension ) {  //Implementation for 3D prism subjected to tension for calibration of corrosion study (Milan and Petr)
        give3DTensionOutput(fileName);
    } else if ( gridType == _3dGopSha ) {  //Implementation of Gop Sha experiment (Ismail)
        give3DGopShaOutput(fileName);
    } else if ( gridType == _3dKupfer )   { //Implementation of Kupfer experiment (Ismail)
        give3DKupferOutput(fileName);
    }
    /* else if ( gridType == _3dNotch ) {  //Implementation of FPZ subjected to multi-axial stress states (Notch conference) */
    /*   give3DNotchOutput(fileName); */
    /* } */
    else {
        converter::error("Unknown grid type\n");
    }
    return;
};

void Grid::giveVtkOutput(const std::string &fileName)
{
    // Write Voronoi *cell* VTU only if not periodic in any dir
    if ( periodicityFlag.at(1) != 1 &&
         periodicityFlag.at(2) != 1 &&
         periodicityFlag.at(3) != 1 ) {
        const std::string filename1 = fileName + ".voronoicell.vtu";
        FILE *f1 = converter::fopen_or_die(filename1, "w");
        giveVoronoiCellVTKOutput(f1);
        std::fclose(f1);
    }

    // Delaunay elements
    const std::string filename2 = fileName + ".delaunayelement.vtu";
    FILE *f2 = converter::fopen_or_die(filename2, "w");
    giveDelaunayElementVTKOutput(f2);
    std::fclose(f2);

    // Voronoi elements
    const std::string filename3 = fileName + ".voronoielement.vtu";
    FILE *f3 = converter::fopen_or_die(filename3, "w");
    giveVoronoiElementVTKOutput(f3);
    std::fclose(f3);
}

void
Grid::give3DSMOutput(const std::string &fileName)
{
    //Template for irregular nonperiodic mechanical models. Do not change for applications

    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    int materialType = 1;
    oofem::IntArray nodes;
    oofem::IntArray crossSectionNodes;

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfLines++;
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Mechanical 3D model\n");
    fprintf(outputStream, "NonLinearStatic nmsteps 1 nsteps 1 contextOutputStep 1\n");
    fprintf(outputStream, "nsteps 1 rtolv 0.001 reqIterations 100 stiffMode 2 maxiter 2000 controllmode 1 stepLength 1. minsteplength 1.e-10 Psi 0.\n");
    fprintf(outputStream, "domain 3dLattice\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_all element_all\n");
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 2\n", numberOfNodes, numberOfLines);

    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
        }
    }

    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "lattice3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, "\n");
        }
    }

    fprintf(outputStream, "simplecs 1\n");
    fprintf(outputStream, "latticedamage3d 1 talpha 0. d 0. e 1.56 e0 87.5e+6 stype 1 wf 40.e+6\n");
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStream, "NodalLoad 2 loadTimeFunction 1 Components 6 0. 0. 0. 0. 0. 0.\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");

    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#NODE number 3 dof 1 unknown d\n");
    fprintf(outputStream, "#LOADLEVEL\n");
    fprintf(outputStream, "##TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");
    return;
}


void
Grid::give3DTMOutput(const std::string &fileName)
{
    //Template for irregular nonperiodic transport models. Do not change for applications


    FILE *outputStream = converter::fopen_or_die(fileName, "w");


    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    oofem::IntArray nodes;
    oofem::IntArray crossSectionNodes;
    Vertex *delaunayVertex;
    Line *voronoiLine;
    Vertex *voronoiVertex;

    int materialType = 1;

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfLines++;
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Transport 3D model\n");
    fprintf(outputStream, "nltransienttransportproblem nsteps 5 deltat 1.0 rtol 0.001 alpha 1. nsmax 200 contextOutputStep 100 nmodules 0\n");
    fprintf(outputStream, "domain 2dMassLatticeTransport\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_all element_all\n");

    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 1 nic 0 nltf 2\n", numberOfNodes, numberOfLines);


    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
            fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
        }
    }

    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "latticemt3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveDelaunayVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, "\n");
        }
    }

    fprintf(outputStream, "simplecs 1\n");
    fprintf(outputStream, "latticetransmat 1 d 1. k 1. vis 1. thetas 1. thetar 0. contype 0 c 0. \n");
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.\n");
    fprintf(outputStream, "BoundaryCondition 2 loadTimeFunction 1 prescribedvalue -1.\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");

    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#LOADLEVEL\n");
    fprintf(outputStream, "##TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");

    return;
}



void
Grid::give3DSMTMOutput(const std::string &fileName)
{
    //Template for irregular nonperiodic coupled mechanical transport models. Do not change for applications

    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    const std::string fileName1 = fileName + ".sm";
    FILE *outputStreamSM = converter::fopen_or_die(fileName1, "w");

    const std::string fileName2 = fileName + ".tm";
    FILE *outputStreamTM = converter::fopen_or_die(fileName2, "w");

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    int materialType = 1;

    oofem::IntArray nodes;
    oofem::IntArray crossSectionNodes;
    oofem::IntArray crossSectionElements;

    //Write the control file
    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Control file for coupling of 3d lattice models\n");
    fprintf(outputStream, "StaggeredProblem nsteps 5 deltat 1. prob1  \"oofem.in.tm\" prob2 \"oofem.in.sm\" contextoutputstep 100 coupling 3 2 1 0\n");
    fprintf(outputStream, "##%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "##INCLUDE oofem.in.tm\n");
    fprintf(outputStream, "##%%END_CHECK%%\n");


    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfLines++;
        }
    }

    //Start with output file for the SM part
    fprintf(outputStreamSM, "oofem.out\n");
    fprintf(outputStreamSM, "Mechanical 3D model\n");
    fprintf(outputStreamSM, "NonLinearStatic nmsteps 1 nsteps 1 contextOutputStep 1\n");
    fprintf(outputStreamSM, "nsteps 1 rtolv 1.e-3 reqIterations 100 stiffMode 2 maxiter 2000 controllmode 1 stepLength 1. minsteplength 1.e-10 Psi 0.\n");
    fprintf(outputStreamSM, "domain 3dLattice\n");
    fprintf(outputStreamSM, "OutputManager tstep_all dofman_all element_all\n");
    fprintf(outputStreamSM, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 2\n", numberOfNodes, numberOfLines);

    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            fprintf(outputStreamSM, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
        }
    }

    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamSM, "lattice3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStreamSM, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            this->giveDelaunayLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            //I will need here to check if element exist in dual network. Not clear how this should be done with the output.
            //Should I return as many elements as vertices and set the nonexistant elements to zero? I think that this is needed, as otherwise the geometry is not defined correctly. Anyway, it is not clear that crosssection nodes are set correctly for the transport elements.
            fprintf(outputStreamSM, " couplingflag 1 couplingelements %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( this->giveVoronoiLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 0 ) {
                    fprintf(outputStreamSM, "%d ", crossSectionElements.at(m + 1) );
                } else {
                    fprintf(outputStreamSM, "0 ");
                }
            }
            fprintf(outputStreamSM, "\n");
        }
    }
    fprintf(outputStreamSM, "simplecs 1\n");
    fprintf(outputStreamSM, "latticedamage3d 1 talpha 0. d 0. e 1.56 e0 87.5e+6 stype 1 wf 40.e+6\n");
    fprintf(outputStreamSM, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStreamSM, "NodalLoad 2 loadTimeFunction 1 Components 6 0. 0. 0. 0. 0. 0.\n");
    fprintf(outputStreamSM, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStreamSM, "PiecewiseLinFunction 2 nPoints 2 t 2 0. 1. f(t) 2 1. 1.\n");
    fprintf(outputStreamSM, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStreamSM, "#NODE number 3 dof 1 unknown d\n");
    fprintf(outputStreamSM, "#LOADLEVEL\n");
    fprintf(outputStreamSM, "##TIME\n");
    fprintf(outputStreamSM, "#%%END_CHECK%%\n");

    //Write the transport input file
    fprintf(outputStreamTM, "oofem.out.tm\n");
    fprintf(outputStreamTM, "Transport part of 3D model\n");
    fprintf(outputStreamTM, "nltransienttransportproblem nsteps 5 deltat 1.0 rtol 0.001 alpha 1. nsmax 200 contextOutputStep 100 nmodules 0\n");
    fprintf(outputStreamTM, "domain 2dMassLatticeTransport\n");
    fprintf(outputStreamTM, "OutputManager tstep_all dofman_all element_all\n");
    fprintf(outputStreamTM, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 1 nic 0 nltf 2\n", numberOfNodes, numberOfLines);

    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
            fprintf(outputStreamTM, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
        }
    }

    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamTM, "latticemt3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStreamTM, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            this->giveVoronoiLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            fprintf(outputStreamTM, " couplingflag 1 couplingelements %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( this->giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 0 ) {
                    fprintf(outputStreamTM, "%d ", crossSectionElements.at(m + 1) );
                } else {
                    fprintf(outputStreamTM, "0 ");
                }
            }

            fprintf(outputStreamTM, "\n");
        }
    }

    fprintf(outputStreamTM, "simplecs 1\n");
    fprintf(outputStreamTM, "latticetransmat 1 d 1. k 1. vis 1. thetas 1. thetar 0. contype 0 c 0. \n");
    fprintf(outputStreamTM, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.\n");
    fprintf(outputStreamTM, "BoundaryCondition 2 loadTimeFunction 1 prescribedvalue -1.\n");
    fprintf(outputStreamTM, "ConstantFunction 1 f(t) 1.\n");

    fprintf(outputStreamTM, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStreamTM, "#LOADLEVEL\n");
    fprintf(outputStreamTM, "##TIME\n");
    fprintf(outputStreamTM, "#%%END_CHECK%%\n");

    return;
}



void
Grid::give3DPeriodicSMTMOutput(const std::string &fileName)
{
    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    const std::string fileName1 = fileName + ".sm";
    FILE *outputStreamSM = converter::fopen_or_die(fileName1, "w");

    const std::string fileName2 = fileName + ".tm";
    FILE *outputStreamTM = converter::fopen_or_die(fileName2, "w");

    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    int materialType = 1;

    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;
    oofem::IntArray crossSectionElements;

    //Write the control file
    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Control file for coupling of 3d lattice models\n");
    fprintf(outputStream, "StaggeredProblem nsteps 5 deltat 1. prob1  \"oofem.in.tm\" prob2 \"oofem.in.sm\" contextoutputstep 100 coupling 3 2 1 0\n");
    fprintf(outputStream, "##%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "##INCLUDE oofem.in.tm\n");
    fprintf(outputStream, "##%%END_CHECK%%\n");


    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfLines++;
        }
    }

    //Start with output file for the SM part
    fprintf(outputStreamSM, "oofem.out\n");
    fprintf(outputStreamSM, "Mechanical 3D model\n");
    fprintf(outputStreamSM, "NonLinearStatic nmsteps 1 nsteps 1 contextOutputStep 1\n");
    fprintf(outputStreamSM, "nsteps 1 rtolv 1.e-3 reqIterations 100 stiffMode 2 maxiter 2000 controllmode 1 stepLength 1. minsteplength 1.e-10 Psi 0.\n");
    fprintf(outputStreamSM, "domain 3dLattice\n");
    fprintf(outputStreamSM, "OutputManager tstep_all dofman_all element_all\n");
    fprintf(outputStreamSM, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 2\n", numberOfNodes + 1, numberOfLines);

    int firstFlag = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                fprintf(outputStreamSM, "node %d coords 3 %e %e %e bc 6 1 1 1 1 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStreamSM, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    //Periodic control node
    fprintf(outputStreamSM, "node %d coords 3 %e %e %e load 1 2\n", this->giveNumberOfDelaunayVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );


    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamSM, "lattice3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStreamSM, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            this->giveDelaunayLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            fprintf(outputStreamSM, " couplingflag 1 couplingelements %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( this->giveVoronoiLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 0 ||  this->giveVoronoiLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 2 ) {
                    fprintf(outputStreamSM, "%d ", crossSectionElements.at(m + 1) );
                } else {
                    fprintf(outputStreamSM, "%d ", this->giveVoronoiLine(crossSectionElements.at(m + 1) )->givePeriodicElement() );
                }
            }
            fprintf(outputStreamSM, "\n");
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveDelaunayVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
            }
            materialType = 1;
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamSM, "lattice3Dboundary %d nodes 3 %d %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), this->giveNumberOfDelaunayVertices() + 1, materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStreamSM, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }

            this->giveDelaunayLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            fprintf(outputStreamSM, " couplingflag 1 couplingelements %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( this->giveVoronoiLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 0 ||  this->giveVoronoiLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 2 ) {
                    fprintf(outputStreamSM, "%d ", crossSectionElements.at(m + 1) );
                } else {
                    fprintf(outputStreamSM, "%d ", this->giveVoronoiLine(crossSectionElements.at(m + 1) )->givePeriodicElement() );
                }
            }
            fprintf(outputStreamSM, " location 2 %d %d", location.at(1), location.at(2) );

            fprintf(outputStreamSM, "\n");
        }
    }
    fprintf(outputStreamSM, "simplecs 1\n");
    fprintf(outputStreamSM, "latticedamage3d 1 talpha 0. d 0. e 1.56 e0 87.5e+6 stype 1 wf 40.e+6\n");
    fprintf(outputStreamSM, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStreamSM, "NodalLoad 2 loadTimeFunction 1 Components 6 0. 0. 0. 0. 0. 0.\n");
    fprintf(outputStreamSM, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStreamSM, "PiecewiseLinFunction 2 nPoints 2 t 2 0. 1. f(t) 2 1. 1.\n");
    fprintf(outputStreamSM, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStreamSM, "#NODE number 3 dof 1 unknown d\n");
    fprintf(outputStreamSM, "#LOADLEVEL\n");
    fprintf(outputStreamSM, "##TIME\n");
    fprintf(outputStreamSM, "#%%END_CHECK%%\n");


    //Determine the number of Voronoi nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of Voronoi lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfLines++;
        }
    }

    //Write the transport input file
    fprintf(outputStreamTM, "oofem.out.tm\n");
    fprintf(outputStreamTM, "Transport part of 3D model\n");
    fprintf(outputStreamTM, "nltransienttransportproblem nsteps 5 deltat 1.0 rtol 0.001 alpha 1. nsmax 200 contextOutputStep 100 nmodules 0\n");
    fprintf(outputStreamTM, "domain 2dMassLatticeTransport\n");
    fprintf(outputStreamTM, "OutputManager tstep_all dofman_all element_all\n");
    fprintf(outputStreamTM, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 1\n", numberOfNodes + 1, numberOfLines);

    firstFlag = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                fprintf(outputStreamTM, "node %d coords 3 %e %e %e bc 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStreamTM, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    //Periodic control node
    fprintf(outputStreamTM, "node %d coords 3 %e %e %e ndofs 3 dofIDmask 3 1 2 3 bc 3 2 1 1\n", this->giveNumberOfVoronoiVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );


    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamTM, "latticemt3D_Discrete %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveDelaunayVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStreamTM, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            this->giveVoronoiLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            fprintf(outputStreamTM, " couplingflag 1 couplingelements %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( this->giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 0 || this->giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 2 ) {
                    fprintf(outputStreamTM, "%d ", crossSectionElements.at(m + 1) );
                } else {
                    fprintf(outputStreamTM, "%d ", this->giveDelaunayLine(crossSectionElements.at(m + 1) )->givePeriodicElement() );
                }
            }
            fprintf(outputStreamTM, "\n");
        } else if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveVoronoiVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
            }
            materialType = 1;
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamTM, "latticemt3Dboundary_Discrete %d nodes 3 %d %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), this->giveNumberOfVoronoiVertices() + 1, materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveDelaunayVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStreamTM, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }

            this->giveVoronoiLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            fprintf(outputStreamTM, " couplingflag 1 couplingelements %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( this->giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 0 || this->giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 2 ) {
                    fprintf(outputStreamTM, "%d ", crossSectionElements.at(m + 1) );
                } else {
                    fprintf(outputStreamTM, "%d ", this->giveDelaunayLine(crossSectionElements.at(m + 1) )->givePeriodicElement() );
                }
            }
            fprintf(outputStreamTM, " location 2 %d %d", location.at(1), location.at(2) );

            fprintf(outputStreamTM, "\n");
        }
    }

    fprintf(outputStreamTM, "simplecs 1\n");
    fprintf(outputStreamTM, "latticetransmat 1 d 1. k 1. vis 1. thetas 1. thetar 0. contype 0 c 0. \n");
    fprintf(outputStreamTM, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.\n");
    fprintf(outputStreamTM, "BoundaryCondition 2 loadTimeFunction 1 prescribedvalue -1.\n");
    fprintf(outputStreamTM, "ConstantFunction 1 f(t) 1.\n");

    fprintf(outputStreamTM, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStreamTM, "#LOADLEVEL\n");
    fprintf(outputStreamTM, "##TIME\n");
    fprintf(outputStreamTM, "#%%END_CHECK%%\n");

    return;
}

void
Grid::give3DPeriodicSMOutput(const std::string &fileName)
{
    //Template for irregular periodic mechanical models. Do not change for applications

    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    /* const std::string fileName1 = fileName + ".sm"; */
    /* FILE *outputStreamSM = converter::fopen_or_die(fileName1, "w"); */

    /* const std::string fileName2 = fileName + ".tm"; */
    /* FILE *outputStreamTM = converter::fopen_or_die(fileName2, "w"); */


    /* FILE *outputStream; */
    /* if ( ( outputStream = fopen(fileName, "w") ) == NULL ) { */
    /*     converter::errorf("Cannot open output file %s", fileName); */
    /* } */

    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    int materialType = 1;
    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfLines++;
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Mechanical 3D model\n");
    fprintf(outputStream, "NonLinearStatic nmsteps 1 nsteps 1 contextOutputStep 1 nmodules 2 profileopt 1 lstype 3 smtype 7\n");
    fprintf(outputStream, "nsteps 100 rtolv 1.e-3 reqIterations 100 stiffMode 2 maxiter 1000 controllmode 0 stepLength 5.e-6 minsteplength 1.e-7 hpcmode 2 hpc 2 %d 2 hpcw 1 1. lstype 3 smtype 7\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "vtkxml primvars 1 1 tstep_all domain_all\n");
    fprintf(outputStream, "gpexportmodule vars 2 59 90 tstep_all domain_all\n");
    fprintf(outputStream, "domain 3dLattice\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_all element_all\n");
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 2\n", numberOfNodes + 1, numberOfLines);

    int firstFlag = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                fprintf(outputStream, "node %d coords 3 %e %e %e bc 6 1 1 1 1 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    //Periodic control node
    fprintf(outputStream, "node %d coords 3 %e %e %e load 1 2\n", this->giveNumberOfDelaunayVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );


    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "lattice3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, "\n");
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveDelaunayVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
            }
            materialType = 1;
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "lattice3Dboundary %d nodes 3 %d %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), this->giveNumberOfDelaunayVertices() + 1, materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, " location 2 %d %d", location.at(1), location.at(2) );

            fprintf(outputStream, "\n");
        }
    }
    fprintf(outputStream, "simplecs 1\n");

    fprintf(outputStream, "latticedamage3d 1 talpha 0. d 0. e 30.e9 e0 100.e-6 stype 3 wf 40.e-6 randvars 1 800 randgen 1 2\n");
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStream, "NodalLoad 2 loadTimeFunction 1 Components 6 0. 1. 0. 0. 0. 0.\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "InterpolatingFunction 2 name field.dat dim 3\n");
    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "#LOADLEVEL\n");
    fprintf(outputStream, "##TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");
    return;
}






void
Grid::give3DPeriodicTMOutput(const std::string &fileName)
{
    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    int materialType = 1;
    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;

    //Determine the number of Voronoi nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of Voronoi lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfLines++;
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Mechanical 3D model\n");
    fprintf(outputStream, "discretetransportproblem ncycles 4 incrementout 0.1 Nsteps 1 deltat 1. rtol 1.e-6 alpha 1. nsmax 2000 contextOutputStep 100 profileopt 1 nmodules 0\n");
    fprintf(outputStream, "domain 3dMassLatticeTransport\n");
    fprintf(outputStream, "OutputManager tstep_all element_all dofman_output {%d}\n", numberOfNodes + 1);
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 1\n", numberOfNodes + 1, numberOfLines);

    int firstFlag = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                fprintf(outputStream, "node %d coords 3 %e %e %e bc 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    //Periodic control node
    fprintf(outputStream, "node %d coords 3 %e %e %e ndofs 3 dofIDmask 3 1 2 3 bc 3 2 1 1\n", this->giveNumberOfVoronoiVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );


    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "latticemt3D_Discrete %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveDelaunayVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, "\n");
        } else if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveVoronoiVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
            }
            materialType = 1;
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "latticemt3Dboundary_Discrete %d nodes 3 %d %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), this->giveNumberOfVoronoiVertices() + 1, materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveDelaunayVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, " location 2 %d %d", location.at(1), location.at(2) );

            fprintf(outputStream, "\n");
        }
    }
    fprintf(outputStream, "simplecs 1\n");

    fprintf(outputStream, "latticelammat 1 d 1000. vis 0.001002 ca 40 st 0.48\n");
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.\n");
    fprintf(outputStream, "BoundaryCondition 2 loadTimeFunction 1 prescribedvalue 1.\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", this->giveNumberOfVoronoiVertices() + 1);
    fprintf(outputStream, "#LOADLEVEL\n");
    fprintf(outputStream, "##TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");
    return;
}


void
Grid::give3DPeriodicPoreSMTMOutput(const std::string &fileName)
{
    //Template for irregular nonperiodic coupled mechanical transport models. Do not change for applications

    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    const std::string fileName1 = fileName + ".sm";
    FILE *outputStreamSM = converter::fopen_or_die(fileName1, "w");

    const std::string fileName2 = fileName + ".tm";
    FILE *outputStreamTM = converter::fopen_or_die(fileName2, "w");


    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfLines;
    oofem::FloatArray coords;
    int materialType = 1;

    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;
    oofem::IntArray crossSectionElements;
    int lastNode;

    int numberOfPipeDiameters;
    int numberOfPeriodicLines = 0;

    //Write the control file
    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Control file for coupling of 3d lattice pore models\n");
    fprintf(outputStream, "StaggeredProblem nsteps 44 deltat 1. prob1  \"oofem.in.tm\" prob2 \"oofem.in.sm\" contextoutputstep 100 coupling 3 2 1 0\n");
    fprintf(outputStream, "##%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "##INCLUDE oofem.in.tm\n");
    fprintf(outputStream, "##%%END_CHECK%%\n");

    int numberOfNodes;

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;

    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
            if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
                numberOfPeriodicLines++;
            }
            numberOfLines++;
        }
    }


    if ( numberOfPeriodicLines % 2 != 0 ) {
        converter::error("In give3DElasticDiscreteSMOutput: The number of perdiodic lines is not even");
    }

    //Determine the number of the diameter values for the pipes that will have to be generated.
    numberOfPipeDiameters = numberOfLines - numberOfPeriodicLines / 2;

    //Create first random numbers for vertices and lines
    //Create random numbers for vertices


    double randomRadius = 0.;

    //Create random numbers for lines

    oofem::FloatArray delaunayLineRadius(numberOfPipeDiameters);
    long randomIntegerOne = this->randomInteger - 1;

    double gaussianMechMean, gaussianMechSTD, gaussianMechCOV;

    gaussianMechMean = log(this->mechMean / sqrt(1 + pow(this->mechCOV, 2) ) );
    gaussianMechSTD = sqrt(log(1 + pow(this->mechCOV, 2) ) );
    gaussianMechCOV = gaussianMechSTD / gaussianMechMean;

    if ( gaussianMechMean < 0. ) {
        gaussianMechSTD = -gaussianMechSTD;
    }

    //Create random numbers for lines
    for ( int i = 0; i < numberOfPipeDiameters; i++ ) {
        randomRadius = normalCdfInverse(ran1(& randomIntegerOne), gaussianMechMean,  gaussianMechSTD);

        //Apply cut-offs if necessary
        if ( randomRadius > mechMax ) {
            randomRadius = mechMax;
        } else if ( randomRadius < mechMin ) {
            randomRadius = mechMin;
        }
        delaunayLineRadius.at(i + 1) = randomRadius;
    }



    //Apply sorted line radii to the elements

    double helpRadius = 0;
    int counter;
    counter = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ) {
            counter++;
            helpRadius = delaunayLineRadius.at(counter);
            this->giveDelaunayLine(i + 1)->setRadius(helpRadius);
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
            //Here only the periodic lines that have been used for generating the rank vector are allowed to pass.
            if ( this->giveDelaunayLine(i + 1)->givePeriodicElement() > i + 1  ) {
                counter++;
                helpRadius = delaunayLineRadius.at(counter);
                this->giveDelaunayLine(i + 1)->setRadius(helpRadius);
                //Assign an equal diameter value to the periodic friend.
                this->giveDelaunayLine(this->giveDelaunayLine(i + 1)->givePeriodicElement() )->setRadius(helpRadius);
            } else if ( this->giveDelaunayLine(i + 1)->givePeriodicElement() < i + 1 && this->giveDelaunayLine(this->giveDelaunayLine(i + 1)->givePeriodicElement() )->givePeriodicElement() !=  i + 1 ) {
                double help = this->giveDelaunayLine(this->giveDelaunayLine(i + 1)->givePeriodicElement() )->giveRadius();
                //                printf("Help is %e %d\n", help, i + 1);
                this->giveDelaunayLine(i + 1)->setRadius(help);
            }
        }
    }


    //Start with output file for the SM part
    fprintf(outputStreamSM, "oofem.out\n");
    fprintf(outputStreamSM, "Mechanical 3D model\n");
    fprintf(outputStreamSM, "NonLinearStatic nmsteps 1 nsteps 44 contextOutputStep 1\n");
    fprintf(outputStreamSM, "nsteps 44 rtolv 1.e-3 stiffMode 2 controllmode 1 refloadmode 0 minsteplength 1.e-10 \n");
    fprintf(outputStreamSM, "domain 3dLattice\n");
    fprintf(outputStreamSM, "OutputManager tstep_all element_all dofman_output {%d}\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStreamSM, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 3 nic 0 nltf 2\n", numberOfNodes + 1, numberOfLines);

    int firstFlag = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                fprintf(outputStreamSM, "node %d coords 3 %e %e %e bc 6 1 1 1 1 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStreamSM, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    //Periodic control node
    fprintf(outputStreamSM, "node %d coords 3 %e %e %e bc 6 0 0 0 0 0 0\n", this->giveNumberOfDelaunayVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );


    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamSM, "lattice3DDiscrete %d nodes 2 %d %d crossSect 1 mat %d rad %.16e ", i + 1, nodes.at(1), nodes.at(2), materialType, this->giveDelaunayLine(i + 1)->giveRadius() );


            this->giveDelaunayLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            fprintf(outputStreamSM, " couplingflag 1 couplingnumber %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( this->giveVoronoiLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 0 ||  this->giveVoronoiLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 2 ) {
                    fprintf(outputStreamSM, "%d ", crossSectionElements.at(m + 1) );
                } else {
                    fprintf(outputStreamSM, "%d ", this->giveVoronoiLine(crossSectionElements.at(m + 1) )->givePeriodicElement() );
                }
            }

            fprintf(outputStreamSM, "\n");
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveDelaunayVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
            }
            materialType = 1;
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamSM, "lattice3Dboundarydiscrete %d nodes 3 %d %d %d crossSect 1 mat %d rad %.16e ", i + 1, nodes.at(1), nodes.at(2), this->giveNumberOfDelaunayVertices() + 1, materialType, this->giveDelaunayLine(i + 1)->giveRadius() );


            this->giveDelaunayLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            fprintf(outputStreamSM, " couplingflag 1 couplingnumber %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( this->giveVoronoiLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 0 ||  this->giveVoronoiLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 2 ) {
                    fprintf(outputStreamSM, "%d ", crossSectionElements.at(m + 1) );
                } else {
                    fprintf(outputStreamSM, "%d ", this->giveVoronoiLine(crossSectionElements.at(m + 1) )->givePeriodicElement() );
                }
            }
            fprintf(outputStreamSM, " location 2 %d %d bodyloads 1 3", location.at(1), location.at(2) );

            fprintf(outputStreamSM, "\n");
        }
    }
    fprintf(outputStreamSM, "simplecs 1\n");
    fprintf(outputStreamSM, "latticelinearelastic 1 talpha 0. d 0. e %e a1 %e a2 %e \n", youngModulus, gammaOne, gammaTwo);
    fprintf(outputStreamSM, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStreamSM, "BoundaryCondition 2 loadTimeFunction 1 prescribedvalue -2.e-12\n");
    fprintf(outputStreamSM, "StructTemperatureLoad 3 loadTimeFunction 1 Components 1 1.0\n");
    fprintf(outputStreamSM, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStreamSM, "PiecewiseLinFunction 2 nPoints 5 t 5 0. 1. 2. 3. 4. f(t) 5 1. 1. 1. 1. 1.\n");
    fprintf(outputStreamSM, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStreamSM, "#NODE number 3 dof 1 unknown d\n");
    fprintf(outputStreamSM, "#LOADLEVEL\n");
    fprintf(outputStreamSM, "##TIME\n");
    fprintf(outputStreamSM, "#%%END_CHECK%%\n");


    //Determine the number of Voronoi nodes in the domain
    numberOfNodes = 0;
    lastNode = 0;
    numberOfPeriodicLines = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
            lastNode = i + 1;
        }
    }

    //Determine the number of Voronoi lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {
            if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {
                numberOfPeriodicLines++;
            }
            numberOfLines++;
        }
    }

    //Check if the number of the periodic Voronoi Lines is odd.
    //The generation continues if the number is even, otherwise it exits.
    if ( numberOfPeriodicLines % 2 != 0 ) {
        converter::error("In give3DPeriodicPoreSMTMOutput: The number of perdiodic lines is not even");
    }

    //Determine the number of the diameter values for the pipes that will have to be generated.
    numberOfPipeDiameters = numberOfLines - numberOfPeriodicLines / 2;

    //Create first random numbers for vertices and lines
    //Create random numbers for vertices
    oofem::FloatArray voronoiVertexRadius(numberOfNodes);
    randomIntegerOne = this->randomInteger;
    randomRadius = 0.;


    randomIntegerOne = this->randomInteger;

    double gaussianPoreMean, gaussianPoreSTD, gaussianPoreCOV;

    gaussianPoreMean = log(this->poreMean / sqrt(1 + pow(this->poreCOV, 2) ) );
    gaussianPoreSTD = sqrt(log(1 + pow(this->poreCOV, 2) ) );
    gaussianPoreCOV = gaussianPoreSTD / gaussianPoreMean;


    if ( gaussianPoreMean < 0. ) {
        gaussianPoreSTD = -gaussianPoreSTD;
    }

    for ( int i = 0; i < numberOfNodes; i++ ) {
        randomRadius = normalCdfInverse(ran1(& randomIntegerOne), gaussianPoreMean, gaussianPoreSTD);


        //Apply cut-offs if necessary
        if ( randomRadius > poreMax ) {
            randomRadius = poreMax;
            //maxCounter++;
        } else if ( randomRadius < poreMin ) {
            randomRadius = poreMin;
            //minCounter++;
        }
        voronoiVertexRadius.at(i + 1) = randomRadius;
    }
    //Create random numbers for lines
    oofem::FloatArray voronoiLineRadius(numberOfPipeDiameters);
    long randomIntegerTwo = this->randomInteger - 1;


    double gaussianThroatMean, gaussianThroatSTD, gaussianThroatCOV;

    gaussianThroatMean = log(this->throatMean / sqrt(1 + pow(this->throatCOV, 2) ) );
    gaussianThroatSTD = sqrt(log(1 + pow(this->throatCOV, 2) ) );
    gaussianThroatCOV = gaussianThroatSTD / gaussianThroatMean;


    if ( gaussianThroatMean < 0. ) {
        gaussianThroatSTD = -gaussianThroatSTD;
    }

    //Create random numbers for lines
    for ( int i = 0; i < numberOfPipeDiameters; i++ ) {
        randomRadius =  normalCdfInverse(ran1(& randomIntegerTwo),  gaussianThroatMean, gaussianThroatSTD);
        //Apply cut-offs if necessary
        if ( randomRadius > throatMax ) {
            randomRadius = throatMax;
        } else if ( randomRadius < throatMin ) {
            randomRadius = throatMin;
        }
        voronoiLineRadius.at(i + 1) = randomRadius;
    }


    /*Strategy to apply random numbers to vertices and elements
     * 1. Sort the line radii (smallest first)
     * 2. Apply the pores radii to pores randomly. This is different to our original idea, but easier to implement.
     * 3. Go trough all elements and find the maximum pore radius for each element. One element has two pores attached to it.
     * 4. Rank the elements according to their maximum pore radius.
     * 5. Apply the smallest element radius to the element with the smallest maximum pore radius using the rank table.
     * Any geometrical violations are not correct. Hopefully there will not be many and will not influence the results too much.
     */

    /*@todo:
     * 1. Check how many violations of the pipe to sphere radii occur.
     */

    //Sort the line radii according to size
    oofem::FloatArray sortedVoronoiLineRadius(numberOfPipeDiameters);
    sortRandomNumbers(sortedVoronoiLineRadius, voronoiLineRadius);

    //Apply the vertex radii to all vertices.
    //This is random as vertices were not sorted and randomly generated.

    counter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            counter++;
            this->giveVoronoiVertex(i + 1)->setRadius(voronoiVertexRadius.at(counter) );
        }
    }

    //Got through all elements and find largest pore
    oofem::IntArray localVertices(2);
    double radius = 0.;
    oofem::FloatArray minRadius(numberOfPipeDiameters);
    counter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) {
            counter++;
            minRadius.at(counter) = 1000000.;
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            for ( int m = 0; m < 2; m++ ) {
                radius = this->giveVoronoiVertex(nodes.at(m + 1) )->giveRadius();
                if ( radius < minRadius.at(counter) ) {
                    minRadius.at(counter) = radius;
                }
            }
        } else if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {
            //A flag was implemented in line.h.
            //This flag shows whether the periodic friend of the line under consideration has been taken into account for the rank vector.
            //flag = 0 if not passed, 1 if passed.
            //Only one of the two periodic friends have to be taken into account.
            if ( this->giveVoronoiLine(i + 1)->givePeriodicElement() > i + 1 ) {
                counter++;
                minRadius.at(counter) = 1000000.;
                this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
                //Go through nodes and replace the ones outside with periodic nodes
                for ( int m = 0; m < 2; m++ ) {
                    if ( this->giveVoronoiVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                        location.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->giveLocation();
                        nodes.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                    radius = this->giveVoronoiVertex(nodes.at(m + 1) )->giveRadius();
                    if ( radius < minRadius.at(counter) ) {
                        minRadius.at(counter) = radius;
                    }
                }
            }
        }
    }


    //Need to sort it now.
    //Firstly create a rank table
    //Create indexx talble for radii
    oofem::IntArray rankVector(numberOfPipeDiameters);
    createRankTable(rankVector, minRadius);

    //Apply sorted line radii to the elements

    helpRadius = 0;
    counter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) {
            counter++;
            helpRadius = sortedVoronoiLineRadius.at(rankVector.at(counter) );
            this->giveVoronoiLine(i + 1)->setRadius(helpRadius);
        } else if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {
            //Here only the periodic lines that have been used for generating the rank vector are allowed to pass.
            if ( this->giveVoronoiLine(i + 1)->givePeriodicElement() > i + 1 ) {
                counter++;
                helpRadius = sortedVoronoiLineRadius.at(rankVector.at(counter) );
                this->giveVoronoiLine(i + 1)->setRadius(helpRadius);
                //Assign an equal diameter value to the periodic friend.
                this->giveVoronoiLine(this->giveVoronoiLine(i + 1)->givePeriodicElement() )->setRadius(helpRadius);
            }
        }
    }

    //Write the transport input file
    fprintf(outputStreamTM, "oofem.out.tm\n");
    fprintf(outputStreamTM, "Transport part of 3D model\n");
    fprintf(outputStreamTM, "discretetransportproblem Nsteps 44 deltat 1. homogenFlag 1 configurationflag 1  rtol 1.e-3 alpha 1. nsmax 2000 contextOutputStep 100000 profileopt 1 nmodules 0\n");
    fprintf(outputStreamTM, "domain 3dMassLatticeTransport\n");
    fprintf(outputStreamTM, "OutputManager tstep_all dofman_output {%d}\n", this->giveNumberOfVoronoiVertices() + 1);
    fprintf(outputStreamTM, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 3 nic 0 nltf 2\n", numberOfNodes + 1, numberOfLines);

    firstFlag = 0;
    counter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                counter++;
                fprintf(outputStreamTM, "pore %d coords 3 %e %e %e bc 1 1 rad %.16e\n", i + 1, coords.at(1), coords.at(2), coords.at(3), voronoiVertexRadius.at(counter) );
            } else {
                counter++;
                fprintf(outputStreamTM, "pore %d coords 3 %e %e %e rad %.16e\n", i + 1, coords.at(1), coords.at(2), coords.at(3), voronoiVertexRadius.at(counter) );
            }
        }
    }

    //Periodic control node
    fprintf(outputStreamTM, "node %d coords 3 %e %e %e ndofs 3 dofIDmask 3 1 2 3 bc 3 2 1 1\n", this->giveNumberOfVoronoiVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );


    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamTM, "latticemt3D_Discrete %d nodes 2 %d %d crossSect 1 mat %d rad %.16e", i + 1, nodes.at(1), nodes.at(2), materialType, this->giveVoronoiLine(i + 1)->giveRadius() );
            this->giveVoronoiLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            fprintf(outputStreamTM, " couplingflag 1 couplingnumber %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( this->giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 0 || this->giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 2 ) {
                    fprintf(outputStreamTM, "%d ", crossSectionElements.at(m + 1) );
                } else {
                    if ( this->giveDelaunayLine(crossSectionElements.at(m + 1) )->givePeriodicElement() == 0 ) {
                        converter::errorf("Clang the element is an outsider one %d", i + 1);
                    }
                    fprintf(outputStreamTM, "%d ", this->giveDelaunayLine(crossSectionElements.at(m + 1) )->givePeriodicElement() );
                }
            }
            fprintf(outputStreamTM, "\n");
        } else if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveVoronoiVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
            }
            materialType = 1;
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamTM, "latticemt3Dboundary_Discrete %d nodes 3 %d %d %d crossSect 1 mat %d rad %.16e", i + 1, nodes.at(1), nodes.at(2), this->giveNumberOfVoronoiVertices() + 1, materialType, this->giveVoronoiLine(i + 1)->giveRadius() );

            this->giveVoronoiLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            fprintf(outputStreamTM, " couplingflag 1 couplingnumber %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( this->giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 0 || this->giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 2 ) {
                    fprintf(outputStreamTM, "%d ", crossSectionElements.at(m + 1) );
                } else {
                    if ( this->giveDelaunayLine(crossSectionElements.at(m + 1) )->givePeriodicElement() == 0 ) {
                        converter::errorf("Clang the element is an outsider one %d", i + 1);
                    }
                    fprintf(outputStreamTM, "%d ", this->giveDelaunayLine(crossSectionElements.at(m + 1) )->givePeriodicElement() );
                }
            }
            fprintf(outputStreamTM, " location 2 %d %d", location.at(1), location.at(2) );

            fprintf(outputStreamTM, "\n");
        }
    }

    fprintf(outputStreamTM, "simplecs 1\n");
    fprintf(outputStreamTM, "latticelammat 1 d 1000. vis 0.001002 poremean 1.000000e-05 st 1. \n");
    fprintf(outputStreamTM, "BoundaryCondition 1 loadTimeFunction 2 prescribedvalue 1.\n");
    fprintf(outputStreamTM, "BoundaryCondition 2 loadTimeFunction 1 prescribedvalue 1.\n");
    fprintf(outputStreamTM, "BoundaryCondition 3 loadTimeFunction 1 prescribedvalue 0.\n");
    fprintf(outputStreamTM, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStreamTM, "PiecewiseLinFunction 2 nPoints 45 t 45 -1. 0. 1. 2. 3. 4. 5. 6. 7. 8. 9. 10. 11. 12. 13. 14. 15. 16. 17. 18. 19. 20. 21. 22. 23. 24. 25. 26. 27. 28. 29. 30. 31. 32. 33. 34. 35. 36. 37. 38. 39. 40. 41. 42. 43. f(t) 45 0.01 0.01 0.01 1.142857e+6 1.142857e+6 1.343e+6 1.343e+6 1.4857e+6 1.4857e+6 1.642857e+6 1.642857e+6 1.8857e+6 1.8857e+6 2.242857e+6 2.242857e+6 2.82857e+6 2.82857e+6 3.828571e+6 3.828571e+6 5.871e+6 5.871e+6 1.61e+7 1.61e+7 2.e+6 2.e+6 3.67e+5 3.67e+5 3.23e+5 3.23e+5 3.e+5 3.e+5 2.777e+5 2.777e+5 2.4966e+5 2.4966e+5 2.1628e+5 2.1628e+5 1.786e+5 1.786e+5 1.38e+5 1.38e+5 9.4793e+4 9.4793e+4 3.2e+4 3.2e+4\n");
    fprintf(outputStreamTM, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStreamTM, "#LOADLEVEL\n");
    fprintf(outputStreamTM, "##TIME\n");
    fprintf(outputStreamTM, "#%%END_CHECK%%\n");

    return;
}


void
Grid::give3DBentoniteCoupledOutput(const std::string &fileName)
{
    //Template for irregular nonperiodic coupled mechanical transport models. Do not change for applications


    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    const std::string fileName1 = fileName + ".sm";
    FILE *outputStreamSM = converter::fopen_or_die(fileName1, "w");

    const std::string fileName2 = fileName + ".tm";
    FILE *outputStreamTM = converter::fopen_or_die(fileName2, "w");

    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfLines;
    oofem::FloatArray coords;
    int materialType = 1;

    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;
    oofem::IntArray crossSectionElements;
    oofem::IntArray cellElements;
    cellElements.zero();
    int lastNode;

    int numberOfPipeDiameters;
    int numberOfPeriodicLines = 0;

    //Write the control file
    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Control file for coupling of 3d lattice models\n");
    fprintf(outputStream, "StaggeredProblem nsteps 1001 deltat 1. prob1  \"oofem.in.tm\" prob2 \"oofem.in.sm\" contextoutputstep 1000 coupling 3 2 1 0\n");
    fprintf(outputStream, "##%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "##INCLUDE oofem.in.tm\n");
    fprintf(outputStream, "##%%END_CHECK%%\n");

    int numberOfNodes;

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;

    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
            if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
                numberOfPeriodicLines++;
            }
            numberOfLines++;
        }
    }


    if ( numberOfPeriodicLines % 2 != 0 ) {
        converter::error("In give3DElasticDiscreteSMOutput: The number of perdiodic lines is not even");
    }

    //Determine the number of the diameter values for the pipes that will have to be generated.
    numberOfPipeDiameters = numberOfLines - numberOfPeriodicLines / 2;

    //Create first random numbers for vertices and lines
    //Create random numbers for vertices


    double randomRadius = 0.;

    //Create random numbers for lines

    oofem::FloatArray delaunayLineRadius(numberOfPipeDiameters);
    long randomIntegerOne = this->randomInteger - 1;

    double gaussianMechMean, gaussianMechSTD, gaussianMechCOV;

    gaussianMechMean = log(this->mechMean / sqrt(1 + pow(this->mechCOV, 2) ) );
    gaussianMechSTD = sqrt(log(1 + pow(this->mechCOV, 2) ) );
    gaussianMechCOV = gaussianMechSTD / gaussianMechMean;

    if ( gaussianMechMean < 0. ) {
        gaussianMechSTD = -gaussianMechSTD;
    }

    //Create random numbers for lines
    for ( int i = 0; i < numberOfPipeDiameters; i++ ) {
        randomRadius = normalCdfInverse(ran1(& randomIntegerOne), gaussianMechMean,  gaussianMechSTD);

        //Apply cut-offs if necessary
        if ( randomRadius > mechMax ) {
            randomRadius = mechMax;
        } else if ( randomRadius < mechMin ) {
            randomRadius = mechMin;
        }
        delaunayLineRadius.at(i + 1) = randomRadius;
    }

    //Apply sorted line radii to the elements

    double helpRadius = 0;
    int counter;
    counter = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ) {
            counter++;
            helpRadius = delaunayLineRadius.at(counter);
            this->giveDelaunayLine(i + 1)->setRadius(helpRadius);
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
            //Here only the periodic lines that have been used for generating the rank vector are allowed to pass.
            if ( this->giveDelaunayLine(i + 1)->givePeriodicElement() > i + 1  ) {
                counter++;
                helpRadius = delaunayLineRadius.at(counter);
                this->giveDelaunayLine(i + 1)->setRadius(helpRadius);
                //Assign an equal diameter value to the periodic friend.
                this->giveDelaunayLine(this->giveDelaunayLine(i + 1)->givePeriodicElement() )->setRadius(helpRadius);
            } else if ( this->giveDelaunayLine(i + 1)->givePeriodicElement() < i + 1 && this->giveDelaunayLine(this->giveDelaunayLine(i + 1)->givePeriodicElement() )->givePeriodicElement() !=  i + 1 ) {
                double help = this->giveDelaunayLine(this->giveDelaunayLine(i + 1)->givePeriodicElement() )->giveRadius();
                printf("Help is %e %d\n", help, i + 1);
                this->giveDelaunayLine(i + 1)->setRadius(help);
            }
        }
    }

    printf("In output: Finished sorting of structural elements\n");

    //Start with output file for the SM part
    fprintf(outputStreamSM, "oofem.out.sm\n");
    fprintf(outputStreamSM, "Mechanical 3D model\n");
    fprintf(outputStreamSM, "NonLinearStatic nmsteps 1 nsteps 1001 contextOutputStep 1000 nmodules 2 profileopt 1\n");
    fprintf(outputStreamSM, "nsteps 1001 rtolv 1.e-3 stiffMode 2 controllmode 1 refloadmode 0 minsteplength 1.e-2 maxiter 1000\n");
    fprintf(outputStreamSM, "vtkxmldiscrete tstep_all domain_all\n");
    fprintf(outputStreamSM, "gpexportmodule vars 1 27 tstep_all domain_all\n");
    fprintf(outputStreamSM, "domain 3dLattice\n");
    fprintf(outputStreamSM, "OutputManager tstep_all dofman_output {%d}\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStreamSM, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 2\n", numberOfNodes + 1, numberOfLines);

    int firstFlag = 0;
    oofem::IntArray cellVertices;
    cellVertices.zero();
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);

            if ( firstFlag == 0 ) {
                firstFlag = 1;
                fprintf(outputStreamSM, "rigidbody %d coords 3 %e %e %e bc 6 1 1 1 1 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStreamSM, "rigidbody %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    //Periodic control node
    fprintf(outputStreamSM, "rigidbody %d coords 3 %e %e %e bc 6 2 2 2 0 0 0\n", this->giveNumberOfDelaunayVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );


    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamSM, "lattice3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStreamSM, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }

            this->giveDelaunayLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            fprintf(outputStreamSM, " couplingflag 1 couplingnumber %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( this->giveVoronoiLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 0 ||  this->giveVoronoiLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 2 ) {
                    fprintf(outputStreamSM, "%d ", crossSectionElements.at(m + 1) );
                } else {
                    fprintf(outputStreamSM, "%d ", this->giveVoronoiLine(crossSectionElements.at(m + 1) )->givePeriodicElement() );
                }
            }

            fprintf(outputStreamSM, "\n");
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveDelaunayVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
            }
            materialType = 1;
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamSM, "lattice3Dboundary %d nodes 3 %d %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), this->giveNumberOfDelaunayVertices() + 1, materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStreamSM, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }

            this->giveDelaunayLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            fprintf(outputStreamSM, " couplingflag 1 couplingnumber %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( this->giveVoronoiLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 0 ||  this->giveVoronoiLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 2 ) {
                    fprintf(outputStreamSM, "%d ", crossSectionElements.at(m + 1) );
                } else {
                    fprintf(outputStreamSM, "%d ", this->giveVoronoiLine(crossSectionElements.at(m + 1) )->givePeriodicElement() );
                }
            }
            fprintf(outputStreamSM, " location 2 %d %d", location.at(1), location.at(2) );

            fprintf(outputStreamSM, "\n");
        }
    }
    fprintf(outputStreamSM, "simplecs 1\n");
    fprintf(outputStreamSM, "latticebentonite 1 talpha 0. d 0. e 30.e9 ft 0. fc 10e6 nodamage 1 a1 0.33 ahard 1.e9 swell 0.1\n");
    fprintf(outputStreamSM, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStreamSM, "BoundaryCondition 2 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStreamSM, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStreamSM, "PiecewiseLinFunction 2 nPoints 2 t 2 -1. 1000. f(t) 2 1. 1.\n");
    fprintf(outputStreamSM, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStreamSM, "#REACTION number %d dof 1 unknown d\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStreamSM, "#REACTION number %d dof 2 unknown d\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStreamSM, "#REACTION number %d dof 3 unknown d\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStreamSM, "#TIME\n");
    fprintf(outputStreamSM, "#%%END_CHECK%%\n");


    printf("In output: Finished structural part\n");

    //Determine the number of Voronoi nodes in the domain
    numberOfNodes = 0;
    lastNode = 0;
    numberOfPeriodicLines = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
            lastNode = i + 1;
        }
    }

    //Determine the number of Voronoi lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {
            if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {
                numberOfPeriodicLines++;
            }
            numberOfLines++;
        }
    }

    //Check if the number of the periodic Voronoi Lines is odd.
    //The generation continues if the number is even, otherwise it exits.
    if ( numberOfPeriodicLines % 2 != 0 ) {
        converter::error("In give3DPeriodicPoreSMTMOutput: The number of perdiodic lines is not even");
    }

    //Determine the number of the diameter values for the pipes that will have to be generated.
    numberOfPipeDiameters = numberOfLines - numberOfPeriodicLines / 2;

    //Create first random numbers for vertices and lines
    //Create random numbers for vertices
    oofem::FloatArray voronoiVertexRadius(numberOfNodes);
    randomIntegerOne = this->randomInteger;
    randomRadius = 0.;


    randomIntegerOne = this->randomInteger;

    double gaussianPoreMean, gaussianPoreSTD, gaussianPoreCOV;

    gaussianPoreMean = log(this->poreMean / sqrt(1 + pow(this->poreCOV, 2) ) );
    gaussianPoreSTD = sqrt(log(1 + pow(this->poreCOV, 2) ) );
    gaussianPoreCOV = gaussianPoreSTD / gaussianPoreMean;


    if ( gaussianPoreMean < 0. ) {
        gaussianPoreSTD = -gaussianPoreSTD;
    }

    for ( int i = 0; i < numberOfNodes; i++ ) {
        randomRadius = normalCdfInverse(ran1(& randomIntegerOne), gaussianPoreMean, gaussianPoreSTD);


        //Apply cut-offs if necessary
        if ( randomRadius > poreMax ) {
            randomRadius = poreMax;
            //maxCounter++;
        } else if ( randomRadius < poreMin ) {
            randomRadius = poreMin;
            //minCounter++;
        }
        voronoiVertexRadius.at(i + 1) = randomRadius;
    }
    //Create random numbers for lines
    oofem::FloatArray voronoiLineRadius(numberOfPipeDiameters);
    long randomIntegerTwo = this->randomInteger - 1;


    double gaussianThroatMean, gaussianThroatSTD, gaussianThroatCOV;

    gaussianThroatMean = log(this->throatMean / sqrt(1 + pow(this->throatCOV, 2) ) );
    gaussianThroatSTD = sqrt(log(1 + pow(this->throatCOV, 2) ) );
    gaussianThroatCOV = gaussianThroatSTD / gaussianThroatMean;


    if ( gaussianThroatMean < 0. ) {
        gaussianThroatSTD = -gaussianThroatSTD;
    }

    //Create random numbers for lines
    //oofem::FloatArray voronoiLineRadius(numberOfPipeDiameters);
    //long randomIntegerTwo = this->randomInteger - 1;
    for ( int i = 0; i < numberOfPipeDiameters; i++ ) {
        randomRadius =  normalCdfInverse(ran1(& randomIntegerTwo),  gaussianThroatMean, gaussianThroatSTD);
        //Apply cut-offs if necessary
        if ( randomRadius > throatMax ) {
            randomRadius = throatMax;
        } else if ( randomRadius < throatMin ) {
            randomRadius = throatMin;
        }
        voronoiLineRadius.at(i + 1) = randomRadius;
    }


    /*Strategy to apply random numbers to vertices and elements
     * 1. Sort the line radii (smallest first)
     * 2. Apply the pores radii to pores randomly. This is different to our original idea, but easier to implement.
     * 3. Go trough all elements and find the maximum pore radius for each element. One element has two pores attached to it.
     * 4. Rank the elements according to their maximum pore radius.
     * 5. Apply the smallest element radius to the element with the smallest maximum pore radius using the rank table.
     * Any geometrical violations are not correct. Hopefully there will not be many and will not influence the results too much.
     */

    /*@todo:
     * 1. Check how many violations of the pipe to sphere radii occur.
     */

    //Sort the line radii according to size
    oofem::FloatArray sortedVoronoiLineRadius(numberOfPipeDiameters);
    sortRandomNumbers(sortedVoronoiLineRadius, voronoiLineRadius);

    //Apply the vertex radii to all vertices.
    //This is random as vertices were not sorted and randomly generated.

    counter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            counter++;
            this->giveVoronoiVertex(i + 1)->setRadius(voronoiVertexRadius.at(counter) );
        }
    }

    //Got through all elements and find largest pore
    oofem::IntArray localVertices(2);
    double radius = 0.;
    oofem::FloatArray minRadius(numberOfPipeDiameters);
    counter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) {
            counter++;
            minRadius.at(counter) = 1000000.;
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            for ( int m = 0; m < 2; m++ ) {
                radius = this->giveVoronoiVertex(nodes.at(m + 1) )->giveRadius();
                if ( radius < minRadius.at(counter) ) {
                    minRadius.at(counter) = radius;
                }
            }
        } else if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {
            //A flag was implemented in line.h.
            //This flag shows whether the periodic friend of the line under consideration has been taken into account for the rank vector.
            //flag = 0 if not passed, 1 if passed.
            //Only one of the two periodic friends have to be taken into account.
            if ( this->giveVoronoiLine(i + 1)->givePeriodicElement() > i + 1 ) {
                counter++;
                minRadius.at(counter) = 1000000.;
                this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
                //Go through nodes and replace the ones outside with periodic nodes
                for ( int m = 0; m < 2; m++ ) {
                    if ( this->giveVoronoiVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                        location.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->giveLocation();
                        nodes.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                    radius = this->giveVoronoiVertex(nodes.at(m + 1) )->giveRadius();
                    if ( radius < minRadius.at(counter) ) {
                        minRadius.at(counter) = radius;
                    }
                }
            }
        }
    }


    //Need to sort it now.
    //Firstly create a rank table
    //Create indexx talble for radii
    oofem::IntArray rankVector(numberOfPipeDiameters);
    createRankTable(rankVector, minRadius);

    //Apply sorted line radii to the elements

    helpRadius = 0;
    counter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) {
            counter++;
            helpRadius = sortedVoronoiLineRadius.at(rankVector.at(counter) );
            this->giveVoronoiLine(i + 1)->setRadius(helpRadius);
        } else if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {
            //Here only the periodic lines that have been used for generating the rank vector are allowed to pass.
            if ( this->giveVoronoiLine(i + 1)->givePeriodicElement() > i + 1 ) {
                counter++;
                helpRadius = sortedVoronoiLineRadius.at(rankVector.at(counter) );
                this->giveVoronoiLine(i + 1)->setRadius(helpRadius);
                //Assign an equal diameter value to the periodic friend.
                this->giveVoronoiLine(this->giveVoronoiLine(i + 1)->givePeriodicElement() )->setRadius(helpRadius);
            }
        }
    }

    printf("In output: Finished sorting of transport elements\n");

    //Write the transport input file
    fprintf(outputStreamTM, "oofem.out.tm\n");
    fprintf(outputStreamTM, "Transport part of 3D model\n");
    fprintf(outputStreamTM, "discretetransportproblem nsteps 101 deltat 1. homogenFlag 1 configurationflag 1  rtol 1.e-3 alpha 1. nsmax 2000 contextOutputStep 100000 profileopt 1 nmodules 1\n");
    fprintf(outputStreamTM, "vtkxmldiscrete tstep_all element_all domain_all cellvars 1 112\n");
    fprintf(outputStreamTM, "domain 3dMassLatticeTransport\n");
    fprintf(outputStreamTM, "OutputManager tstep_all dofman_output {%d}\n", this->giveNumberOfVoronoiVertices() + 1);
    fprintf(outputStreamTM, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 4 nic 0 nltf 2\n", numberOfNodes + 1, numberOfLines);

    firstFlag = 0;
    counter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
            this->giveVoronoiVertex(i + 1)->giveCellElements(cellElements);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                counter++;
                fprintf(outputStreamTM, "pore %d coords 3 %e %e %e bc 1 1 rad %.16e couplingflag 1 couplingnumber %d ", i + 1, coords.at(1), coords.at(2), coords.at(3), voronoiVertexRadius.at(counter), cellElements.giveSize() );
                for (int m = 0; m < cellElements.giveSize(); m++) {
                    if ( this->giveDelaunayLine(cellElements.at(m + 1) )->giveOutsideFlag() == 0 || this->giveDelaunayLine(cellElements.at(m + 1) )->giveOutsideFlag() == 2 ) {
                        fprintf(outputStreamTM, "%d ", cellElements.at(m + 1) );
                    } else {
                        if ( this->giveDelaunayLine(cellElements.at(m + 1) )->givePeriodicElement() == 0 ) {
                            converter::errorf("Error: the element %d is outside", i + 1);
                        }
                        fprintf(outputStreamTM, "%d ", this->giveDelaunayLine(cellElements.at(m + 1) )->givePeriodicElement() );
                    }
                }
                fprintf(outputStreamTM, "\n");
            } else {
                counter++;
                fprintf(outputStreamTM, "pore %d coords 3 %e %e %e rad %.16e couplingflag 1 couplingnumber %d ", i + 1, coords.at(1), coords.at(2), coords.at(3), voronoiVertexRadius.at(counter), cellElements.giveSize() );
                for (int m = 0; m < cellElements.giveSize(); m++) {
                    if ( this->giveDelaunayLine(cellElements.at(m + 1) )->giveOutsideFlag() == 0 || this->giveDelaunayLine(cellElements.at(m + 1) )->giveOutsideFlag() == 2 ) {
                        fprintf(outputStreamTM, "%d ", cellElements.at(m + 1) );
                    } else {
                        if ( this->giveDelaunayLine(cellElements.at(m + 1) )->givePeriodicElement() == 0 ) {
                            converter::errorf("Error: the element %d is outside", i + 1);
                        }
                        fprintf(outputStreamTM, "%d ", this->giveDelaunayLine(cellElements.at(m + 1) )->givePeriodicElement() );
                    }
                }
                fprintf(outputStreamTM, "\n");
            }
        }
    }
    //Periodic control node
    fprintf(outputStreamTM, "node %d coords 3 %e %e %e ndofs 3 dofIDmask 3 1 2 3 bc 3 2 2 2\n", this->giveNumberOfVoronoiVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );


    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamTM, "latticemt3D_Discrete %d nodes 2 %d %d crossSect 1 mat %d rad %.16e", i + 1, nodes.at(1), nodes.at(2), materialType, this->giveVoronoiLine(i + 1)->giveRadius() );
            this->giveVoronoiLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            fprintf(outputStreamTM, " couplingflag 1 couplingnumber %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( this->giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 0 || this->giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 2 ) {
                    fprintf(outputStreamTM, "%d ", crossSectionElements.at(m + 1) );
                } else {
                    if ( this->giveDelaunayLine(crossSectionElements.at(m + 1) )->givePeriodicElement() == 0 ) {
                        converter::errorf("Error: the element %d is outside", i + 1);
                    }
                    fprintf(outputStreamTM, "%d ", this->giveDelaunayLine(crossSectionElements.at(m + 1) )->givePeriodicElement() );
                }
            }
            fprintf(outputStreamTM, "\n");
        } else if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveVoronoiVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
            }
            materialType = 1;
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamTM, "latticemt3Dboundary_Discrete %d nodes 3 %d %d %d crossSect 1 mat %d rad %.16e", i + 1, nodes.at(1), nodes.at(2), this->giveNumberOfVoronoiVertices() + 1, materialType, this->giveVoronoiLine(i + 1)->giveRadius() );

            this->giveVoronoiLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            fprintf(outputStreamTM, " couplingflag 1 couplingnumber %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( this->giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 0 || this->giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 2 ) {
                    fprintf(outputStreamTM, "%d ", crossSectionElements.at(m + 1) );
                } else {
                    if ( this->giveDelaunayLine(crossSectionElements.at(m + 1) )->givePeriodicElement() == 0 ) {
                        converter::errorf("Error: the element is outside %d", i + 1);
                    }
                    fprintf(outputStreamTM, "%d ", this->giveDelaunayLine(crossSectionElements.at(m + 1) )->givePeriodicElement() );
                }
            }
            fprintf(outputStreamTM, " location 2 %d %d", location.at(1), location.at(2) );

            fprintf(outputStreamTM, "\n");
        }
    }

    fprintf(outputStreamTM, "simplecs 1\n");
    fprintf(outputStreamTM, "latticelammat 1 d 1000. vis 0.001002 pipemean 1.e-6 pipemin 1.e-9 poremean 1.e-4 poremin 1.e-8\n");
    fprintf(outputStreamTM, "BoundaryCondition 1 loadTimeFunction 2 prescribedvalue 1.\n");
    fprintf(outputStreamTM, "BoundaryCondition 2 loadTimeFunction 1 prescribedvalue %e\n", specimenDimension.at(1) );
    fprintf(outputStreamTM, "BoundaryCondition 3 loadTimeFunction 1 prescribedvalue %e\n", specimenDimension.at(2) );
    fprintf(outputStreamTM, "BoundaryCondition 4 loadTimeFunction 1 prescribedvalue %e\n", specimenDimension.at(3) );
    fprintf(outputStreamTM, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStreamTM, "PiecewiseLinFunction 2 datafile \"time.in\"\n");
    fprintf(outputStreamTM, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStreamTM, "#LOADLEVEL\n");
    fprintf(outputStreamTM, "##TIME\n");
    fprintf(outputStreamTM, "#%%END_CHECK%%\n");

    return;
}



void
Grid::give3DPeriodicPoreSMOutput(const std::string &fileName)
{
    FILE *outputStream = converter::fopen_or_die(fileName, "w");


    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    int materialType = 1;
    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;
    int numberOfPeriodicLines;
    int numberOfPipeDiameters;

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    //Determine the total number and the number of the periodic Delaunay lines in the domain
    numberOfLines = 0;
    numberOfPeriodicLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
            if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
                numberOfPeriodicLines++;
            }
            numberOfLines++;
        }
    }

    //Check if the number of the periodic Delaunay Lines is odd.
    //The generation continues if the number is even, otherwise it exits.
    if ( numberOfPeriodicLines % 2 != 0 ) {
        converter::error("In give3DPlasticDiscreteTriaxialSMOutput: The number of perdiodic lines is not even");
    }

    //Determine the number of the diameter values for the pipes that will have to be generated.
    numberOfPipeDiameters = numberOfLines - numberOfPeriodicLines / 2;

    //Create first random numbers for vertices and lines
    //Create random numbers for vertices

    long randomIntegerOne = this->randomInteger;
    double randomRadius = 0.;

    //Create random numbers for lines
    oofem::FloatArray delaunayLineRadius(numberOfPipeDiameters);
    long randomIntegerTwo = this->randomInteger - 1;


    double gaussianMechMean, gaussianMechSTD, gaussianMechCOV;

    gaussianMechMean = log(this->mechMean / sqrt(1 + pow(this->mechCOV, 2) ) );
    gaussianMechSTD = sqrt(log(1 + pow(this->mechCOV, 2) ) );
    gaussianMechCOV = gaussianMechSTD / gaussianMechMean;
    printf("mechmean %e mechcov %e mechmin %e asd asd mechmax %e\n", this->mechMean, this->mechCOV, mechMax, mechMin);

    if ( gaussianMechMean < 0. ) {
        gaussianMechSTD = -gaussianMechSTD;
    }

    //Create random numbers for lines
    for ( int i = 0; i < numberOfPipeDiameters; i++ ) {
        randomRadius = normalCdfInverse(ran1(& randomIntegerOne), gaussianMechMean,  gaussianMechSTD);

        //Apply cut-offs if necessary
        if ( randomRadius > mechMax ) {
            randomRadius = mechMax;
        } else if ( randomRadius < mechMin ) {
            randomRadius = mechMin;
        }
        delaunayLineRadius.at(i + 1) = randomRadius;
    }

    //Apply sorted line radii to the elements

    double helpRadius = 0;
    int counter;
    counter = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ) {
            counter++;
            helpRadius = delaunayLineRadius.at(counter);
            this->giveDelaunayLine(i + 1)->setRadius(helpRadius);
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
            //Here only the periodic lines that have been used for generating the rank vector are allowed to pass.
            if ( this->giveDelaunayLine(i + 1)->givePeriodicElement() > i + 1  ) {
                counter++;
                helpRadius = delaunayLineRadius.at(counter);
                this->giveDelaunayLine(i + 1)->setRadius(helpRadius);
                //Assign an equal diameter value to the periodic friend.
                this->giveDelaunayLine(this->giveDelaunayLine(i + 1)->givePeriodicElement() )->setRadius(helpRadius);
            } else if ( this->giveDelaunayLine(i + 1)->givePeriodicElement() < i + 1 && this->giveDelaunayLine(this->giveDelaunayLine(i + 1)->givePeriodicElement() )->givePeriodicElement() !=  i + 1 ) {
                double help = this->giveDelaunayLine(this->giveDelaunayLine(i + 1)->givePeriodicElement() )->giveRadius();
                //                printf("Help is %e %d\n", help, i + 1);
                this->giveDelaunayLine(i + 1)->setRadius(help);
            }
        }
    }


    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Mechanical 3D model\n");
    fprintf(outputStream, "NonLinearStatic nmsteps 1 nsteps 1 contextOutputStep 1 profileopt 1 nmodules 0\n");
    fprintf(outputStream, "nsteps 200 rtolv 1.e-3 maxiter 1000 controllmode 1 refloadmode 0 minsteplength 1.e-2 ddm 6 %d 1 %d 2 %d 3 ddv 3 -4.0e-3 -4.0e-3 -4.0e-3 ddltf 3 stiffmode 2 \n", this->giveNumberOfDelaunayVertices() + 1, this->giveNumberOfDelaunayVertices() + 1, this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "domain 3dLattice\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_output {%d}\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 4 nic 0 nltf 3\n", numberOfNodes + 1, numberOfLines);


    int firstFlag = 0;
    counter = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            //used for stats
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                counter++;
                fprintf(outputStream, "node %d coords 3 %e %e %e bc 6 1 1 1 1 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                counter++;
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    //statistics for the line lengths

    counter = 0;
    //Periodic control node
    fprintf(outputStream, "node %d coords 3 %e %e %e load 3 2 3 4\n", this->giveNumberOfDelaunayVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );


    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);
            fprintf(outputStream, "lattice3DDiscrete %d nodes 2 %d %d crossSect 1 mat %d rad %.16e ", i + 1, nodes.at(1), nodes.at(2), materialType, this->giveDelaunayLine(i + 1)->giveRadius() );
            fprintf(outputStream, "\n");


            counter++;
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            if ( this->giveDelaunayLine(i + 1)->givePeriodicElement() > i + 1 ) {
                counter++;
            }

            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveDelaunayVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
            }
            materialType = 1;

            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "lattice3DBoundaryDiscrete %d nodes 3 %d %d %d crossSect 1 mat %d rad %.16e ", i + 1, nodes.at(1), nodes.at(2), this->giveNumberOfDelaunayVertices() + 1, materialType, this->giveDelaunayLine(i + 1)->giveRadius() );
            fprintf(outputStream, " location 2 %d %d", location.at(1), location.at(2) );
            fprintf(outputStream, "\n");
        }
    }

    fprintf(outputStream, "simplecs 1\n");

    fprintf(outputStream, "latticeplastic 1 talpha 0. d 0. e %e a1 %e a2 %e beta %e phi %e \n", youngModulus, gammaOne, gammaTwo, this->tanBeta, this->tanPhi);
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.\n");
    fprintf(outputStream, "NodalLoad 2 loadTimeFunction 2 Components 6 %e 0. 0. 0. 0. 0.\n", this->confinement * specimenDimension.at(3) * specimenDimension.at(2) );
    fprintf(outputStream, "NodalLoad 3 loadTimeFunction 2 Components 6 0. %e 0. 0. 0. 0.\n", this->confinement * specimenDimension.at(1) * specimenDimension.at(3) );
    fprintf(outputStream, "NodalLoad 4 loadTimeFunction 2 Components 6 0. 0. %e 0. 0. 0.\n", this->confinement * specimenDimension.at(1) * specimenDimension.at(2) );

    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "PiecewiseLinFunction 2 nPoints 4 t 4  -1. 0.1 0.2 200. f(t) 4 1. 1. 0. 0.\n");
    fprintf(outputStream, "PiecewiseLinFunction 3 nPoints 4 t 4  -1. 0.1 0.2 200. f(t) 4 0. 0. 1. 200.\n");

    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "#LOADLEVEL\n");
    fprintf(outputStream, "##TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");


    return;
}




void
Grid::give3DBentoniteSMOutput(const std::string &fileName)
{
    //Output for 3D Bentonite modelling

    FILE *outputStream = converter::fopen_or_die(fileName, "w");


    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    int materialType = 1;
    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfLines++;
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Mechanical part of 3D bentonite model subjected to confinement\n");
    fprintf(outputStream, "NonLinearStatic nmsteps 1 nsteps 1 contextOutputStep 1000 nmodules 2 profileopt 1 lstype 3 smtype 7\n");
    fprintf(outputStream, "nsteps 200 rtolv 1.e-3 reqIterations 100 stiffMode 2 manrmsteps 10 maxiter 200 controllmode 0 stepLength 5.e-8 minsteplength 5.e-8 maxrestarts 0 hpcmode 2 hpc 2 %d 2 hpcw 1 -1. lstype 3 smtype 7\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "vtkxml primvars 1 1 tstep_all domain_all\n");
    fprintf(outputStream, "gpexportmodule vars 1 27 tstep_all domain_all\n");
    fprintf(outputStream, "domain 3dLattice\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_output {%d}\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 1\n", numberOfNodes + 1, numberOfLines);

    int firstFlag = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                fprintf(outputStream, "node %d coords 3 %e %e %e bc 6 1 1 1 1 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    //Periodic control node
    fprintf(outputStream, "node %d coords 3 %e %e %e load 1 2 bc 6 1 0 1 0 0 0\n", this->giveNumberOfDelaunayVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );

    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "lattice3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }

            fprintf(outputStream, "\n");
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveDelaunayVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
            }
            materialType = 1;
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "lattice3Dboundary %d nodes 3 %d %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), this->giveNumberOfDelaunayVertices() + 1, materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, " location 2 %d %d", location.at(1), location.at(2) );

            fprintf(outputStream, "\n");
        }
    }
    fprintf(outputStream, "simplecs 1\n");

    fprintf(outputStream, "latticebentonite 1 talpha 0. d 0. e 30.e9 ft 0.e6 fc 10e6 nodamage 1 a1 0.33 ahard 1.e9\n");
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStream, "NodalLoad 2 loadTimeFunction 1 Components 6 0. 1. 0. 0. 0. 0.\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#REACTION number %d dof 1 unknown d\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "#REACTION number %d dof 2 unknown d\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "#REACTION number %d dof 3 unknown d\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "#TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");
    return;
}


void
Grid::give3DBentoniteTMOutput(const std::string &fileName)
{
    //Template for irregular periodic mechanical models. Do not change for applications

    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    const std::string fileName1 = fileName + "_Lengths.dat";
    FILE *outputStreamLengths = converter::fopen_or_die(fileName1, "w");

    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    int materialType = 1;
    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;
    int numberOfPeriodicLines;
    int numberOfPipeDiameters;

    //Determine the number of Voronoi nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    //Determine the total number and the number of the periodic Voronoi lines in the domain
    numberOfLines = 0;
    numberOfPeriodicLines = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        oofem::IntArray helpNodes(2), locationArray(2), switches(2);
        oofem::FloatArray helpCoords(3);

        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 3 ) {
            this->giveVoronoiLine(i + 1)->giveLocalVertices(helpNodes);
            for ( int m = 0; m < 2; m++ ) {
                giveVoronoiVertex(helpNodes.at(m + 1) )->giveCoordinates(helpCoords);
                locationArray.at(m + 1) = this->giveRegion(1)->giveSwitches(switches, helpCoords);
            }
        }
        //If they do, one of the switches array from the nodes is used to shift the line

        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 3 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 4 ) {
            if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 || ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 3 && locationArray.at(1) == locationArray.at(2) ) ) {
                numberOfPeriodicLines++;
            }
            numberOfLines++;
        }
    }

    //Check if the number of the periodic Voronoi Lines is odd.
    //The generation continues if the number is even, otherwise it exits.
    if ( numberOfPeriodicLines % 2 != 0 ) {
        converter::error("In give3DPeriodicPoreTMOutput: The number of perdiodic lines is not even");
    }

    //Determine the number of the diameter values for the pipes that will have to be generated.
    numberOfPipeDiameters = numberOfLines - numberOfPeriodicLines / 2;

    //Create first random numbers for vertices and lines
    //Create random numbers for vertices
    oofem::FloatArray voronoiVertexRadius(numberOfNodes);
    long randomIntegerOne = this->randomInteger;
    double randomRadius = 0.;

    double gaussianPoreMean, gaussianPoreSTD, gaussianPoreCOV;

    gaussianPoreMean = log(this->poreMean / sqrt(1 + pow(this->poreCOV, 2) ) );
    gaussianPoreSTD = sqrt(log(1 + pow(this->poreCOV, 2) ) );
    gaussianPoreCOV = gaussianPoreSTD / gaussianPoreMean;


    if ( gaussianPoreMean < 0. ) {
        gaussianPoreSTD = -gaussianPoreSTD;
    }

    for ( int i = 0; i < numberOfNodes; i++ ) {
        randomRadius = normalCdfInverse(ran1(& randomIntegerOne), gaussianPoreMean, gaussianPoreSTD);

        //Apply cut-offs if necessary
        if ( randomRadius > poreMax ) {
            randomRadius = poreMax;
            //maxCounter++;
        } else if ( randomRadius < poreMin ) {
            randomRadius = poreMin;
            //minCounter++;
        }
        voronoiVertexRadius.at(i + 1) = randomRadius;
    }
    //Create random numbers for lines
    oofem::FloatArray voronoiLineRadius(numberOfPipeDiameters);
    long randomIntegerTwo = this->randomInteger - 1;


    double gaussianThroatMean, gaussianThroatSTD, gaussianThroatCOV;

    gaussianThroatMean = log(this->throatMean / sqrt(1 + pow(this->throatCOV, 2) ) );
    gaussianThroatSTD = sqrt(log(1 + pow(this->throatCOV, 2) ) );
    gaussianThroatCOV = gaussianThroatSTD / gaussianThroatMean;


    if ( gaussianThroatMean < 0. ) {
        gaussianThroatSTD = -gaussianThroatSTD;
    }

    //Create random numbers for lines
    //oofem::FloatArray voronoiLineRadius(numberOfPipeDiameters);
    //long randomIntegerTwo = this->randomInteger - 1;
    for ( int i = 0; i < numberOfPipeDiameters; i++ ) {
        randomRadius =  normalCdfInverse(ran1(& randomIntegerTwo),  gaussianThroatMean, gaussianThroatSTD);

        //Apply cut-offs if necessary
        if ( randomRadius > throatMax ) {
            randomRadius = throatMax;
        } else if ( randomRadius < throatMin ) {
            randomRadius = throatMin;
        }
        voronoiLineRadius.at(i + 1) = randomRadius;
    }


    /*Strategy to apply random numbers to vertices and elements
     * 1. Sort the line radii (smallest first)
     * 2. Apply the pores radii to pores randomly. This is different to our original idea, but easier to implement.
     * 3. Go trough all elements and find the maximum pore radius for each element. One element has two pores attached to it.
     * 4. Rank the elements according to their maximum pore radius.
     * 5. Apply the smallest element radius to the element with the smallest maximum pore radius using the rank table.
     * Any geometrical violations are not correct. Hopefully there will not be many and will not influence the results too much.
     */

    /*@todo:
     * 1. Check how many violations of the pipe to sphere radii occur.
     */

    //Sort the line radii according to size
    oofem::FloatArray sortedVoronoiLineRadius(numberOfPipeDiameters);
    sortRandomNumbers(sortedVoronoiLineRadius, voronoiLineRadius);

    //Apply the vertex radii to all vertices.
    //This is random as vertices were not sorted and randomly generated.

    int counter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            counter++;
            this->giveVoronoiVertex(i + 1)->setRadius(voronoiVertexRadius.at(counter) );
        }
    }

    //Got through all elements and find largest pore
    oofem::IntArray localVertices(2);
    double radius = 0.;
    oofem::FloatArray minRadius(numberOfPipeDiameters);
    counter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        oofem::IntArray helpNodes(2), locationArray(2), switches(2);
        oofem::FloatArray helpCoords(3);

        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 3 ) {
            this->giveVoronoiLine(i + 1)->giveLocalVertices(helpNodes);
            for ( int m = 0; m < 2; m++ ) {
                giveVoronoiVertex(helpNodes.at(m + 1) )->giveCoordinates(helpCoords);
                locationArray.at(m + 1) = this->giveRegion(1)->giveSwitches(switches, helpCoords);
            }
        }

        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 4 ) {
            counter++;
            minRadius.at(counter) = 1000000.;
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            for ( int m = 0; m < 2; m++ ) {
                //new
                if ( this->giveVoronoiVertex(nodes.at(m + 1) )->giveOutsideFlag() == 2 ) {
                    nodes.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
                radius = this->giveVoronoiVertex(nodes.at(m + 1) )->giveRadius();
                if ( radius < minRadius.at(counter) ) {
                    minRadius.at(counter) = radius;
                }
            }
        } else if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 || ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 3 && locationArray.at(1) == locationArray.at(2) ) ) {
            //A flag was implemented in line.h.
            //This flag shows whether the periodic friend of the line under consideration has been taken into account for the rank vector.
            //flag = 0 if not passed, 1 if passed.
            //Only one of the two periodic friends have to be taken into account.
            if ( this->giveVoronoiLine(i + 1)->givePeriodicElement() > i + 1 ) {
                counter++;
                minRadius.at(counter) = 1000000.;
                this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
                //Go through nodes and replace the ones outside with periodic nodes
                for ( int m = 0; m < 2; m++ ) {
                    if ( this->giveVoronoiVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 || this->giveVoronoiVertex(nodes.at(m + 1) )->giveOutsideFlag() == 2 ) {
                        location.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->giveLocation();
                        nodes.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                    radius = this->giveVoronoiVertex(nodes.at(m + 1) )->giveRadius();
                    if ( radius < minRadius.at(counter) ) {
                        minRadius.at(counter) = radius;
                    }
                }
            }
        }
    }

    int help;
    help = 1;
    //Need to sort it now.
    //Firstly create a rank table
    //Create indexx talble for radii
    oofem::IntArray rankVector(numberOfPipeDiameters);
    createRankTable(rankVector, minRadius);

    //Apply sorted line radii to the elements

    double helpRadius = 0;
    counter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) {
            counter++;
            helpRadius = sortedVoronoiLineRadius.at(rankVector.at(counter) );
            this->giveVoronoiLine(i + 1)->setRadius(helpRadius);
        } else if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 3 ) {
            //Here only the periodic lines that have been used for generating the rank vector are allowed to pass.
            if ( this->giveVoronoiLine(i + 1)->givePeriodicElement() > i + 1 ) {
                counter++;
                helpRadius = sortedVoronoiLineRadius.at(rankVector.at(counter) );
                this->giveVoronoiLine(i + 1)->setRadius(helpRadius);
                //Assign an equal diameter value to the periodic friend.
                this->giveVoronoiLine(this->giveVoronoiLine(i + 1)->givePeriodicElement() )->setRadius(helpRadius);
            }
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Transport discrete 3D model\n");
    fprintf(outputStream, "discretetransportproblem nsteps 101 deltat 1. homogenFlag 1 configurationflag 1  rtol 1.e-3 alpha 1. nsmax 2000 contextOutputStep 100000 profileopt 1 nmodules 0\n");
    fprintf(outputStream, "domain 3dMassLatticeTransport\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_output {%d}\n", this->giveNumberOfVoronoiVertices() + 1);
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 4 nic 0 nltf 2\n", numberOfNodes + 1, numberOfLines);

    int firstFlag = 0;
    counter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                counter++;
                fprintf(outputStream, "pore %d coords 3 %e %e %e bc 1 1 rad %.16e\n", i + 1, coords.at(1), coords.at(2), coords.at(3), voronoiVertexRadius.at(counter) );
            } else {
                counter++;
                fprintf(outputStream, "pore %d coords 3 %e %e %e rad %.16e\n", i + 1, coords.at(1), coords.at(2), coords.at(3), voronoiVertexRadius.at(counter) );
            }
        }
    }

    //Periodic control node
    fprintf(outputStream, "node %d coords 3 %e %e %e ndofs 3 dofIDmask 3 1 2 3 bc 3 2 3 4\n", this->giveNumberOfVoronoiVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );


    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);
            fprintf(outputStream, "latticemt3D_Discrete %d nodes 2 %d %d crossSect 1 mat %d rad %.16e ", i + 1, nodes.at(1), nodes.at(2), materialType, this->giveVoronoiLine(i + 1)->giveRadius() );
            fprintf(outputStream, "\n");
        } else if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveVoronoiVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
            }
            materialType = 1;
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "latticemt3Dboundary_Discrete %d nodes 3 %d %d %d crossSect 1 mat %d rad %.16e ", i + 1, nodes.at(1), nodes.at(2), this->giveNumberOfVoronoiVertices() + 1, materialType, this->giveVoronoiLine(i + 1)->giveRadius() );
            fprintf(outputStream, " location 2 %d %d", location.at(1), location.at(2) );
            fprintf(outputStream, "\n");
        }
    }

    double maxSuction;
    maxSuction = 2. * 1. / sortedVoronoiLineRadius.at(1);

    fprintf(outputStream, "simplecs 1\n");
    fprintf(outputStream, "latticelammat 1 d 1000. vis 0.001002 rate %e poremean %e pipemin 0.01e-9 poremin 0.1e-9\n", this->deltarad, this->poreMean);
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 2 prescribedvalue 1.\n");
    fprintf(outputStream, "BoundaryCondition 2 loadTimeFunction 1 prescribedvalue %e\n", specimenDimension.at(1) );
    fprintf(outputStream, "BoundaryCondition 3 loadTimeFunction 1 prescribedvalue %e\n", specimenDimension.at(2) );
    fprintf(outputStream, "BoundaryCondition 4 loadTimeFunction 1 prescribedvalue %e\n", specimenDimension.at(3) );
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "PiecewiseLinFunction 2 datafile \"time.in\"\n");
    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", this->giveNumberOfVoronoiVertices() + 1);
    fprintf(outputStream, "#LOADLEVEL\n");
    fprintf(outputStream, "##TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");
    return;
}



void
Grid::give3DCantileverTMOutput(const std::string &fileName)
{
    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    oofem::IntArray nodes;
    oofem::IntArray crossSectionNodes;
    Vertex *delaunayVertex;
    Line *voronoiLine;
    Vertex *voronoiVertex;

    int materialType = 1;

    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) {
            this->giveRegion(1)->modifyVoronoiCrossSection(i + 1);
        }
    }

    //Make sure that outsiders are correctly identified and run find outsiders once more.
    this->giveRegion(1)->findOutsiders(boundaries);

    //Determine the number of Voronoi lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) {
            for ( int n = 0; n < 2; n++ ) {
                this->giveVoronoiVertex(this->giveVoronoiLine(i + 1)->giveLocalVertex(n + 1) )->setPrintFlag(1);
            }
            numberOfLines++;
            orderDelaunayCrossSectionVertices(i + 1);
        }
    }

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->givePrintFlag() == 1 ) {
            numberOfNodes++;
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "3D transport benchmark for 2015 CMAME paper\n");
    fprintf(outputStream, "nltransienttransportproblem nsteps 9 deltat 0.025 rtol 1.e-10 alpha 0.7 nsmax 200 contextOutputStep 100 nmodules 2 profileopt 1\n");
    fprintf(outputStream, "vtkxml primvars 1 5 tstep_all domain_all\n");
    fprintf(outputStream, "dofman primvars 1 5 tstep_all domain_all\n");
    fprintf(outputStream, "domain 3dMassLatticeTransport\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_all element_all\n");
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 1 nic %d nltf 1\n", numberOfNodes, numberOfLines, numberOfNodes);
    int nodeCounter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->givePrintFlag() == 1 ) {
            this->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
            nodeCounter++;
            if ( fabs(coords.at(3) - boundaries.at(5) ) < TOL || fabs(coords.at(3) - boundaries.at(6) ) < TOL ) {
                fprintf(outputStream, "node %d coords 3 %e %e %e ic 1 %d bc 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3), nodeCounter);
            } else {
                fprintf(outputStream, "node %d coords 3 %e %e %e ic 1 %d\n", i + 1, coords.at(1), coords.at(2), coords.at(3), nodeCounter);
            }
        }
    }

    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "latticemt3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveDelaunayVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, " mlength 1.e-8");
            fprintf(outputStream, "\n");
        }
    }

    fprintf(outputStream, "simplecs 1\n");
    fprintf(outputStream, "latticetransmat 1 d 1. k 1. vis 1. thetas 1. thetar 0. contype 0 c 1. \n");
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.\n");
    //Loop
    double initialValue = 0;
    double pi = 3.1415926535897;
    int conditionCounter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->givePrintFlag() == 1 ) {
            conditionCounter++;
            this->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
            initialValue = sin(pi * coords.at(3) / specimenDimension.at(3) );
            fprintf(outputStream, "InitialCondition %d conditions 1 u %e\n", conditionCounter, initialValue);
        }
    }



    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");

    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#LOADLEVEL\n");
    fprintf(outputStream, "##TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");

    return;
}


void
Grid::give3DCantileverTMExtraOutput(const std::string &fileName)
{
    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);


    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) {
            this->giveRegion(1)->modifyVoronoiCrossSection(i + 1);
        }
    }

    //Make sure that outsiders are correctly identified and run find outsiders once more.
    this->giveRegion(1)->findOutsiders(boundaries);

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    oofem::IntArray nodes;
    oofem::IntArray crossSectionNodes;
    Vertex *delaunayVertex;
    Line *voronoiLine;
    Vertex *voronoiVertex;

    int materialType = 1;

    //Determine the number of Delaunay lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) {
            for ( int n = 0; n < 2; n++ ) {
                this->giveVoronoiVertex(this->giveVoronoiLine(i + 1)->giveLocalVertex(n + 1) )->setPrintFlag(1);
            }
            numberOfLines++;
            orderDelaunayCrossSectionVertices(i + 1);
        }
    }

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    int numberOfBC = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->givePrintFlag() == 1 ) {
            numberOfNodes++;
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "3D transport benchmark for 2015 CMAME paper\n");
    fprintf(outputStream, "nltransienttransportproblem nsteps 9 deltat 0.025 rtol 1.e-10 alpha 1. nsmax 200 contextOutputStep 100 nmodules 2 profileopt 1\n");
    fprintf(outputStream, "vtkxml primvars 1 5 tstep_all domain_all\n");
    fprintf(outputStream, "dofman primvars 1 5 tstep_all domain_all\n");
    fprintf(outputStream, "domain 3dMassLatticeTransport\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_all element_all\n");
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 1\n", numberOfNodes, numberOfLines);
    int nodeCounter = 0;
    int bcCounter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->givePrintFlag() == 1 ) {
            this->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
            nodeCounter++;

            if ( fabs(coords.at(3) - boundaries.at(5) ) < TOL ) {
                fprintf(outputStream, "node %d coords 3 %e %e %e bc 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else if ( fabs(coords.at(3) - boundaries.at(6) ) < TOL ) {
                fprintf(outputStream, "node %d coords 3 %e %e %e bc 1 2\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }


    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "latticemt3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveDelaunayVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, "\n");
            //	}
        }
    }

    fprintf(outputStream, "simplecs 1\n");
    fprintf(outputStream, "latticetransmat 1 d 1. k 1. vis 1. thetas 1. thetar 0. contype 0 c 0. \n");


    fprintf(outputStream, "BoundaryCondition 1  loadTimeFunction 1 prescribedvalue 0\n");
    fprintf(outputStream, "BoundaryCondition 2  loadTimeFunction 1 prescribedvalue 1\n");

    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");

    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#LOADLEVEL\n");
    fprintf(outputStream, "##TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");

    return;
}

void
Grid::give3DCantileverSMOutput(const std::string &fileName)
{
    //Output for 3D cantilver benchmark for fracture

    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    int materialType = 1;
    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ||  this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) {
            if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) {
                this->giveRegion(1)->modifyVoronoiCrossSection(i + 1);
            }
            numberOfLines++;
        }
    }

    //Find boundary node closes to loadnode x = 0.5, y = 0.25, z = 0.25 (boundaryFlag = 2) and support node x = 0, y = 0.25 z = 0.25  (boundaryFlag = 1)
    oofem::FloatArray supportCoords(3);
    supportCoords.at(1) = 0.;
    supportCoords.at(2) = 0.125;
    supportCoords.at(3) = 0.25;
    double distanceSupport = 1.e10;
    int supportNode = 0;

    oofem::FloatArray loadCoords(3);
    loadCoords.at(1) = 0.25;
    loadCoords.at(2) = 0.125;
    loadCoords.at(3) = 0.25;
    double distanceLoad = 1.e10;
    int loadNode = 0;

    //Fine control nodes
    double distance = 0.;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            distance = sqrt(pow(supportCoords.at(1) - coords.at(1), 2.) + pow(supportCoords.at(2) - coords.at(2), 2.) + pow(supportCoords.at(3) - coords.at(3), 2.) );
            if ( distance < distanceSupport ) {
                distanceSupport = distance;
                supportNode = this->giveDelaunayVertex(i + 1)->giveNumber();
            }
            distance = sqrt(pow(loadCoords.at(1) - coords.at(1), 2.) + pow(loadCoords.at(2) - coords.at(2), 2.) + pow(loadCoords.at(3) - coords.at(3), 2.) );
            if ( distance < distanceLoad ) {
                distanceLoad = distance;
                loadNode = this->giveDelaunayVertex(i + 1)->giveNumber();
            }
        }
    }
    printf("Start writing the output\n");

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Mechanical 3D model\n");
    fprintf(outputStream, "NonLinearStatic nmsteps 2 nsteps 1 contextOutputStep 1000 nmodules 2 profileopt 1 lstype 3 smtype 7\n");
    fprintf(outputStream, "nsteps 10 rtolv 1.e-3 reqIterations 100 stiffMode 1 manrmsteps 10 maxiter 200 controllmode 0 stepLength 8.e-6 minsteplength 8.e-6 maxrestarts 0 hpcmode 2 hpc 2 %d 1 hpcw 1 1. lstype 3 smtype 7 donotfixload\n", loadNode);
    fprintf(outputStream, "nsteps 90 rtolv 1.e-3 reqIterations 100 stiffMode 1 manrmsteps 10 maxiter 200 controllmode 0 stepLength 1.2e-5 minsteplength 1.2e-5 maxrestarts 0 hpcmode 2 hpc 2 %d 1 hpcw 1 1. lstype 3 smtype 7\n", loadNode);
    fprintf(outputStream, "vtkxml primvars 1 1 tstep_all domain_all\n");
    fprintf(outputStream, "gpexportmodule vars 5 59 90 84 85 78 tstep_all domain_all\n");
    fprintf(outputStream, "domain 3dLattice\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_output {%d}\n", loadNode);
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 3 nbc 2 nic 0 nltf 1\n", numberOfNodes, numberOfLines);

    //Print Nodes

    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            if ( this->giveDelaunayVertex(i + 1)->giveNumber() == supportNode ) {
                fprintf(outputStream, "node %d coords 3 %e %e %e bc 6 1 1 1 0 0 0\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else if ( this->giveDelaunayVertex(i + 1)->giveNumber() == loadNode ) {
                fprintf(outputStream, "node %d coords 3 %e %e %e load 1 2 bc 6 0 1 1 0 0 0\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else if ( fabs(coords.at(1) - boundaries.at(1) ) < TOL ) {
                fprintf(outputStream, "rigidarmnode %d coords 3 %e %e %e master %d mastermask 6 1 0 0 1 1 1 doftype 6 2 0 0 2 2 2\n", i + 1, coords.at(1), coords.at(2), coords.at(3), supportNode);
            } else if ( fabs(coords.at(1) - boundaries.at(2) ) < TOL ) {
                fprintf(outputStream, "rigidarmnode %d coords 3 %e %e %e master %d mastermask 6 1 0 0 1 1 1 doftype 6 2 0 0 2 2 2\n", i + 1, coords.at(1), coords.at(2), coords.at(3), loadNode);
            } else {
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    //Print Elements

    oofem::FloatArray coordsOne(3);
    oofem::FloatArray coordsTwo(3);
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) { //Elements are inside
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;

            //Check  that notch is
            this->giveDelaunayVertex(nodes.at(1) )->giveCoordinates(coordsOne);
            this->giveDelaunayVertex(nodes.at(2) )->giveCoordinates(coordsTwo);

            //Check if element is attached to support or load node. If yes, set it to be elastic.
            if ( nodes.at(1) == loadNode || nodes.at(2) == loadNode || nodes.at(1) == supportNode || nodes.at(2) == supportNode ) {
                materialType = 3;
            }

            //Introduce Notch
            if ( ( coordsOne.at(1) > 0.125 && coordsTwo.at(1) < 0.125 ) || ( coordsOne.at(1) < 0.125 && coordsTwo.at(1) > 0.125 ) ) {
                if ( ( coordsOne.at(3) < 0.25 ) || ( coordsTwo.at(3) < 0.25 ) ) {
                    materialType = 2;
                }
            }

            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "lattice3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, "\n");
        }
    }
    fprintf(outputStream, "simplecs 1\n");

    fprintf(outputStream, "latticedamage 1 talpha 0. d 0. e 30.e9 e0 100.e-6 stype 3 wf 40.e-6\n");
    fprintf(outputStream, "latticedamage 2 talpha 0. d 0. e 30.e9 e0 1.e-6 stype 3 wf 40.e-6\n");
    fprintf(outputStream, "latticelinearelastic 3 talpha 0. d 0. e 30.e9\n");
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStream, "NodalLoad 2 loadTimeFunction 1 Components 6 1. 0. 0. 0. 0. 0.\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#NODE number %d dof 1 unknown d\n", loadNode);
    fprintf(outputStream, "#LOADLEVEL\n");
    fprintf(outputStream, "##TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");
    return;
}




void
Grid::give3DCantileverSMTMOutput(const std::string &fileName)
{
    //Output for 3D cantilver benchmark for fracture

    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    const std::string fileName1 = fileName + ".sm";
    FILE *outputStreamSM = converter::fopen_or_die(fileName1, "w");

    const std::string fileName2 = fileName + ".tm";
    FILE *outputStreamTM = converter::fopen_or_die(fileName2, "w");




    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfDelaunayNodes, numberOfDelaunayLines;
    int numberOfVoronoiNodes, numberOfVoronoiLines;

    oofem::FloatArray coords;
    int materialType = 1;
    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;


    //*****************************
    //Write the control file
    //*****************************
    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Control file for coupled application\n");
    fprintf(outputStream, "StaggeredProblem nsteps 12000 prob1  \"oofem.in.tm\" prob2 \"oofem.in.sm\" coupling 3 2 1 0 dtf 1\n");
    fprintf(outputStream, "domain heattransfer\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_all element_all\n");
    fprintf(outputStream, "ndofman 0 nelem 0 nmat 0 ncrosssect 0 nbc 0 nic 0 nltf 1\n");
    fprintf(outputStream, "PiecewiseLinFunction 1 npoints 5 t 5 0. 10. 20. 200. 12000. f(t) 5 1. 1. 1. 20. 20.\n");
    fprintf(outputStream, "##%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "##INCLUDE oofem.in.tm\n");
    fprintf(outputStream, "##%%END_CHECK%%\n");


    //Determine the number of Delaunay nodes in the domain
    numberOfDelaunayNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ||  this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfDelaunayNodes++;
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfDelaunayLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) {
            if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) {
                this->giveRegion(1)->modifyVoronoiCrossSection(i + 1);
            }
            numberOfDelaunayLines++;
        }
    }

    //Find boundary node closes to loadnode x = 0.5, y = 0.25, z = 0.25 (boundaryFlag = 2) and support node x = 0, y = 0.25 z = 0.25  (boundaryFlag = 1)
    oofem::FloatArray supportCoords(3);
    supportCoords.at(1) = 0.;
    supportCoords.at(2) = 0.125;
    supportCoords.at(3) = 0.25;
    double distanceSupport = 1.e10;
    int supportNode = 0;

    oofem::FloatArray loadCoords(3);
    loadCoords.at(1) = 0.25;
    loadCoords.at(2) = 0.125;
    loadCoords.at(3) = 0.25;
    double distanceLoad = 1.e10;
    int loadNode = 0;

    double distance = 0.;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            distance = sqrt(pow(supportCoords.at(1) - coords.at(1), 2.) + pow(supportCoords.at(2) - coords.at(2), 2.) + pow(supportCoords.at(3) - coords.at(3), 2.) );
            if ( distance < distanceSupport ) {
                distanceSupport = distance;
                supportNode = this->giveDelaunayVertex(i + 1)->giveNumber();
            }

            distance = sqrt(pow(loadCoords.at(1) - coords.at(1), 2.) + pow(loadCoords.at(2) - coords.at(2), 2.) + pow(loadCoords.at(3) - coords.at(3), 2.) );
            if ( distance < distanceLoad ) {
                distanceLoad = distance;
                loadNode = this->giveDelaunayVertex(i + 1)->giveNumber();
            }
        }
    }


    //****************************************
    //Write the output file for the SM part
    //******************************************
    fprintf(outputStreamSM, "oofem.out.sm\n");
    fprintf(outputStreamSM, "Mechanical 3D model\n");
    fprintf(outputStreamSM, "NonLinearStatic nmsteps 1 nsteps 1 contextOutputStep 10000 nmodules 2 profileopt 1\n");
    fprintf(outputStreamSM, "nsteps 12000 rtolf 1.e-3 stiffMode 1 manrmsteps 10 maxiter 200 controllmode 1 refloadmode 0 minsteplength 1.e-5 ddm 2 %d 1 ddv 1 1 ddltf 2\n", loadNode);
    fprintf(outputStreamSM, "vtkxml primvars 1 1 tstep_step 100 domain_all\n");
    fprintf(outputStreamSM, "gpexportmodule vars 5 59 90 84 85 78 tstep_step 100 domain_all\n");
    fprintf(outputStreamSM, "domain 3dLattice\n");
    fprintf(outputStreamSM, "OutputManager tstep_all dofman_output {%d}\n", loadNode);
    fprintf(outputStreamSM, "ndofman %d nelem %d ncrosssect 1 nmat 3 nbc 2 nic 0 nltf 2\n", numberOfDelaunayNodes, numberOfDelaunayLines);

    //Print Nodes

    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            if ( this->giveDelaunayVertex(i + 1)->giveNumber() == supportNode ) {
                fprintf(outputStreamSM, "node %d coords 3 %e %e %e bc 6 1 1 1 0 0 0\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else if ( this->giveDelaunayVertex(i + 1)->giveNumber() == loadNode ) {
                fprintf(outputStreamSM, "node %d coords 3 %e %e %e load 1 2 bc 6 0 1 1 0 0 0\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else if ( fabs(coords.at(1) - boundaries.at(1) ) < TOL ) {
                fprintf(outputStreamSM, "rigidarmnode %d coords 3 %e %e %e master %d mastermask 6 1 0 0 1 1 1 doftype 6 2 0 0 2 2 2\n", i + 1, coords.at(1), coords.at(2), coords.at(3), supportNode);
            } else if ( fabs(coords.at(1) - boundaries.at(2) ) < TOL ) {
                fprintf(outputStreamSM, "rigidarmnode %d coords 3 %e %e %e master %d mastermask 6 1 0 0 1 1 1 doftype 6 2 0 0 2 2 2\n", i + 1, coords.at(1), coords.at(2), coords.at(3), loadNode);
            } else {
                fprintf(outputStreamSM, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    //Print Elements

    oofem::FloatArray coordsOne(3);
    oofem::FloatArray coordsTwo(3);
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) { //Elements are inside
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;

            //Check  that notch is
            this->giveDelaunayVertex(nodes.at(1) )->giveCoordinates(coordsOne);
            this->giveDelaunayVertex(nodes.at(2) )->giveCoordinates(coordsTwo);

            //Check if element is attached to support or load node. If yes, set it to be elastic.
            if ( nodes.at(1) == loadNode || nodes.at(2) == loadNode || nodes.at(1) == supportNode || nodes.at(2) == supportNode ) {
                materialType = 3;
            }

            //Introduce Notch
            if ( ( coordsOne.at(1) > 0.125 && coordsTwo.at(1) < 0.125 ) || ( coordsOne.at(1) < 0.125 && coordsTwo.at(1) > 0.125 ) ) {
                if ( ( coordsOne.at(3) < 0.25 ) || ( coordsTwo.at(3) < 0.25 ) ) {
                    materialType = 2;
                }
            }

            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamSM, "lattice3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStreamSM, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStreamSM, "\n");
        }
    }
    fprintf(outputStreamSM, "simplecs 1\n");

    fprintf(outputStreamSM, "latticedamage 1 talpha 0. d 0. e 30.e9 e0 100.e-6 stype 3 wf 40.e-6\n");
    fprintf(outputStreamSM, "latticedamage 2 talpha 0. d 0. e 30.e9 e0 1.e-6 stype 3 wf 40.e-6\n");
    fprintf(outputStreamSM, "latticelinearelastic 3 talpha 0. d 0. e 30.e9\n");
    fprintf(outputStreamSM, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStreamSM, "NodalLoad 2 loadTimeFunction 1 Components 6 1. 0. 0. 0. 0. 0.\n");
    fprintf(outputStreamSM, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStreamSM, "PiecewiseLinFunction 2 npoints 3 t 3 0. 9. 12000000 f(t) 3 0. 1.e-4 1.e-4\n");
    fprintf(outputStreamSM, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStreamSM, "#NODE number %d dof 1 unknown d\n", loadNode);
    fprintf(outputStreamSM, "#LOADLEVEL\n");
    fprintf(outputStreamSM, "##TIME\n");
    fprintf(outputStreamSM, "#%%END_CHECK%%\n");




    //****************************************
    //Write the output file for the TM part
    //******************************************

    //Reset material type
    materialType = 1;

    oofem::IntArray crossSectionElements;

    //Determine the number of Delaunay lines in the domain
    numberOfDelaunayLines = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0  ) {
            for ( int n = 0; n < 2; n++ ) {
                this->giveVoronoiVertex(this->giveVoronoiLine(i + 1)->giveLocalVertex(n + 1) )->setPrintFlag(1);
            }
            numberOfDelaunayLines++;
            orderDelaunayCrossSectionVertices(i + 1);
        }
    }


    //Determine the number of Delaunay nodes in the domain
    numberOfVoronoiNodes = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->givePrintFlag() == 1 ) {
            numberOfVoronoiNodes++;
        }
    }


    fprintf(outputStreamTM, "oofem.out.tm\n");
    fprintf(outputStreamTM, "3D coupled benchmark for 2015 CMAME paper\n");
    fprintf(outputStreamTM, "nltransienttransportproblem nsteps 12000 rtol 1.e-3 alpha 1. nsmax 2000 contextOutputStep 3000 nmodules 2 profileopt 1 lumpedcapa deltatfunction 2\n");
    fprintf(outputStreamTM, "vtkxml primvars 1 5 tstep_step 100 domain_all\n");
    fprintf(outputStreamTM, "dofman primvars 1 5 tstep_step 100 domain_all\n");
    fprintf(outputStreamTM, "domain 3dMassLatticeTransport\n");
    fprintf(outputStreamTM, "OutputManager tstep_step 10 dofman_output {");
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->givePrintFlag() == 1 ) {
            this->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
            if ( fabs(coords.at(3) - boundaries.at(5) ) < TOL ) {
                fprintf(outputStreamTM, " %d", i + 1);
            }
        }
    }
    fprintf(outputStreamTM, " }\n");
    fprintf(outputStreamTM, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 1 nic 1 nltf 2\n", numberOfVoronoiNodes, numberOfDelaunayLines);
    int nodeCounter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->givePrintFlag() == 1 ) {
            this->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
            nodeCounter++;
            if ( fabs(coords.at(3) - boundaries.at(5) ) < TOL ) {
                fprintf(outputStreamTM, "node %d coords 3 %e %e %e ic 1 1 bc 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStreamTM, "node %d coords 3 %e %e %e ic 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamTM, "latticemt3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveDelaunayVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStreamTM, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }

            this->giveVoronoiLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            //	    crossSectionElements.printYourself();
            fprintf(outputStreamTM, " couplingflag 1 couplingnumber %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( !( this->giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 0 || this->giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 3 ) ) {
                    fprintf(outputStreamTM, "-1 ");
                    printf("crosssection element not found\n");
                } else {
                    fprintf(outputStreamTM, "%d ", crossSectionElements.at(m + 1) );
                }
            }
            fprintf(outputStreamTM, "mlength 1.e-8");

            fprintf(outputStreamTM, "\n");
        }
    }

    fprintf(outputStreamTM, "simplecs 1\n");
    fprintf(outputStreamTM, "latticetransmat 1 d 1.e3 k 1.e-17 vis 1.e-3 thetas 0.1 thetar 0.0 a 1.e6 m 0.5 contype 1 ctor 0.001\n");
    fprintf(outputStreamTM, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 1.732e6\n");
    fprintf(outputStreamTM, "InitialCondition 1 conditions 1 u 1.732e6\n");
    fprintf(outputStreamTM, "PiecewiseLinFunction 1 npoints 4 t 4 0. 10. 20. 1200000. f(t) 4 1. 1. 0. 0.\n");
    fprintf(outputStreamTM, "PiecewiseLinFunction 2 npoints 5 t 5 0. 10. 20. 200. 12000. f(t) 5 1. 1. 1. 20. 20.\n");
    fprintf(outputStreamTM, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStreamTM, "#LOADLEVEL\n");
    fprintf(outputStreamTM, "##TIME\n");
    fprintf(outputStreamTM, "#%%END_CHECK%%\n");

    return;
}


void
Grid::give3DFPZOutput(const std::string &fileName)
{
    //Output for 3D fracture process zone modelling

    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    int materialType = 1;
    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) {
            numberOfLines++;
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Mechanical 3D model\n");
    fprintf(outputStream, "NonLinearStatic nmsteps 1 nsteps 1 contextOutputStep 1000 nmodules 2 profileopt 1 lstype 3 smtype 7\n");
    fprintf(outputStream, "nsteps 200 rtolv 1.e-3 reqIterations 100 stiffMode 1 manrmsteps 10 maxiter 200 controllmode 0 stepLength 5.e-7 minsteplength 5.e-7 maxrestarts 0 hpcmode 2 hpc 2 %d 31 hpcw 1 1. lstype 3 smtype 7\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "vtkxmllattice primvars 1 1 cellvars 2 60 90 tstep_all domain_all cross 1\n");
    fprintf(outputStream, "gpexportmodule vars 5 59 90 84 85 78 tstep_all domain_all\n");
    fprintf(outputStream, "domain 3dLattice\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_output {%d}\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 1\n", numberOfNodes + 1, numberOfLines);

    int firstFlag = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                fprintf(outputStream, "node %d coords 3 %e %e %e bc 6 1 1 1 1 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    //Periodic control node
    fprintf(outputStream, "node %d coords 3 %e %e %e load 1 2\n", this->giveNumberOfDelaunayVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );

    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ||  this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) { //Elements are inside
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;
            this->giveDelaunayLine(i + 1)->updateMaterial(materialType);
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "lattice3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, "\n");
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveDelaunayVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
            }
            materialType = 1;
            this->giveDelaunayLine(i + 1)->updateMaterial(materialType);
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "lattice3Dboundarytruss %d nodes 3 %d %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), this->giveNumberOfDelaunayVertices() + 1, materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, " location 2 %d %d", location.at(1), location.at(2) );

            fprintf(outputStream, "\n");
        }
    }
    fprintf(outputStream, "latticecs 1 material 1\n");
    fprintf(outputStream, "latticeplastdam 1 d 0 talpha 0. calpha 0. e 50.46e9 a1 0.215 ft 3.e6 fc 30.e6 wf 50.e-6 ahard 1.e-6 angle1 0.5 flow 0.25\n");
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStream, "NodalLoad 2 loadTimeFunction 1 dofs 1 31 Components 1 1.\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#NODE number %d dof 31 unknown d\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "#LOADLEVEL\n");
    fprintf(outputStream, "##TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");
    return;
}


void
Grid::give3DFPZFibreOutput(const std::string &fileName)
{
    //Output for 3D fracture process zone modelling

    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    printf("Starting writing data in output file \n");

    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    int materialType = 1;
    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;
    oofem::FloatArray coordsOne, coordsTwo, centre, radii;
    double radius, distanceOne, distanceTwo;
    int isIn1, isIn2; // bool var to know whether endpoints of element qre inside or outide an inclusion
    double coord_x1, coord_y1, coord_z1, coord_x2, coord_y2, coord_z2;

    //Determine the number of nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    for ( int i = 0; i < this->giveNumberOfReinforcementNode(); i++ ) {
        if ( this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfLines++;
        }
    }
    for ( int i = 0; i < this->giveNumberOfLatticeBeams(); i++ ) {
        if ( this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 0 || this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfLines++;
        }
    }
    for ( int i = 0; i < this->giveNumberOfLatticeLinks(); i++ ) {
        if ( this->giveLatticeLink(i + 1)->giveOutsideFlag() == 0 || this->giveLatticeLink(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfLines++;
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Mechanical 3D periodic cell for fracture\n");
    fprintf(outputStream, "StaticStructural nmsteps 2 nsteps 1 contextOutputStep 1000 nmodules 2 profileopt 1 lstype 4 smtype 8\n");
    fprintf(outputStream, "nsteps 1 rtolv 1.e-3 maxiter 200 stiffmode 1 lstype 4 smtype 8 initialguess 1\n");
    fprintf(outputStream, "nsteps 200 rtolv 1.e-3 maxiter 200 ddm 2 %d 31 ddv 1 4.0e-4 ddltf 3 stiffmode 1 initialguess 1 lstype 4 smtype 8\n", this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode() + 1);
    fprintf(outputStream, "vtkxmllattice primvars 1 1 cellvars 3 60 90 111 tstep_all domain_all cross 1\n");
    fprintf(outputStream, "gpexportmodule vars 5 59 90 84 85 78 tstep_all domain_all\n");
    fprintf(outputStream, "domain 3dLattice\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_output {%d}\n", this->giveNumberOfDelaunayVertices()
            + this->giveNumberOfReinforcementNode() + 1); // a verifier
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 4\n", numberOfNodes + 1, numberOfLines);

    //Print nodes

    int firstFlag = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                fprintf(outputStream, "node %d coords 3 %e %e %e bc 6 1 1 1 1 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    // Print reinforcement nodes

    for ( int i = 0; i < this->giveNumberOfReinforcementNode(); i++ ) {
        if ( this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveReinforcementNode(i + 1)->giveCoordinates(coords);

            fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + this->giveNumberOfDelaunayVertices() + 1, coords.at(1), coords.at(2), coords.at(3) );
        }
    }

    //Periodic control node
    fprintf(outputStream, "node %d coords 3 %e %e %e load 1 2\n", this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );

    printf("Finished writing data for nodes \n");

    //Print Delaunay elements

    int global_index(0);//usefull to have a global numerotation of elements

    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        // take into account of the mesostructure to choose the materials of the element
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);

            //Deal with inclusions
            materialType = 1;
            this->giveDelaunayLine(i + 1)->updateMaterial(1);
            this->giveDelaunayVertex(nodes.at(1) )->giveCoordinates(coordsOne);
            this->giveDelaunayVertex(nodes.at(2) )->giveCoordinates(coordsTwo);


            for ( int m = 0; m < this->giveNumberOfInclusions(); m++ ) {
                //Distinguish between inclusions
                if ( !strcmp(this->giveInclusion(m + 1)->giveClassName(), "InterfaceSphere") ) {
                    ( ( InterfaceSphere * ) this->giveInclusion(m + 1) )->giveCentre(centre);
                    radius = ( ( InterfaceSphere * ) this->giveInclusion(m + 1) )->giveRadius() + ( ( InterfaceSphere * ) this->giveInclusion(m + 1) )->giveITZThickness() / 2.;
                    distanceOne = sqrt(pow(coordsOne.at(1) - centre.at(1), 2.) +
                                       pow(coordsOne.at(2) - centre.at(2), 2.) +
                                       pow(coordsOne.at(3) - centre.at(3), 2.) );
                    distanceTwo = sqrt(pow(coordsTwo.at(1) - centre.at(1), 2.) +
                                       pow(coordsTwo.at(2) - centre.at(2), 2.) +
                                       pow(coordsTwo.at(3) - centre.at(3), 2.) );

                    if ( distanceOne > radius && distanceTwo > radius ) {
                        materialType = 1;
                        this->giveDelaunayLine(i + 1)->updateMaterial(1);
                    } else if ( ( distanceOne > radius && distanceTwo < radius ) || ( distanceOne < radius && distanceTwo > radius ) ) {
                        materialType = 3;
                        this->giveDelaunayLine(i + 1)->updateMaterial(3);
                        break;
                    } else {
                        materialType = 2;
                        this->giveDelaunayLine(i + 1)->updateMaterial(2);
                        break;
                    }
                } else if ( !strcmp(this->giveInclusion(m + 1)->giveClassName(), "Ellipsoid") )   {
                    ( ( Ellipsoid * ) this->giveInclusion(m + 1) )->giveCentre(centre);
                    ( ( Ellipsoid * ) this->giveInclusion(m + 1) )->giveRadii(radii);

                    coord_x1 = coordsOne.at(1);
                    coord_y1 = coordsOne.at(2);
                    coord_z1 = coordsOne.at(3);
                    coord_x2 = coordsTwo.at(1);
                    coord_y2 = coordsTwo.at(2);
                    coord_z2 = coordsTwo.at(3);
                    isIn1 = ( ( Ellipsoid * ) this->giveInclusion(m + 1) )->isInside(coord_x1, coord_y1, coord_z1);
                    isIn2 = ( ( Ellipsoid * ) this->giveInclusion(m + 1) )->isInside(coord_x2, coord_y2, coord_z2);
                    if ( isIn1 == 0 && isIn2 == 0 ) {
                        materialType = 1;     // matrix : we suppose elements to be small enough
                        this->giveDelaunayLine(i + 1)->updateMaterial(1);
                    } else if ( ( isIn1 == 0 && isIn2 == 1 ) || ( isIn2 == 0 && isIn1 == 1 ) || ( isIn1 == 2 || isIn2 == 2 ) ) {
                        materialType = 3;     //ITZ
                        this->giveDelaunayLine(i + 1)->updateMaterial(3);

                        break;
                    } else {
                        materialType = 2;     //aggregate
                        this->giveDelaunayLine(i + 1)->updateMaterial(2);
                        break;
                    }
                }
            }
        }

        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            //    materialType = 1;
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            global_index++;
            fprintf(outputStream, "lattice3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", global_index, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );


            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, "\n");
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveDelaunayVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
            }
            //  materialType = 1;
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);
            global_index++;
            fprintf(outputStream, "lattice3Dboundary %d nodes 3 %d %d %d crossSect 1 mat %d polycoords %d", global_index, nodes.at(1), nodes.at(2), this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode() + 1, materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, " location 2 %d %d", location.at(1), location.at(2) );

            fprintf(outputStream, "\n");
        }
    }

    printf("Finished writing data for delaunay elements \n");

    // print latticebeams


    for ( int i = 0; i < this->giveNumberOfLatticeBeams(); i++ ) {
        if ( this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveLatticeBeam(i + 1)->giveLocalVertices(nodes);
            global_index++;
            fprintf(outputStream, "latticeBeam3D %d nodes 2 %d %d crossSect 1 mat 4 diameter %e", global_index,
                    nodes.at(1) + this->giveNumberOfDelaunayVertices(),
                    nodes.at(2) + this->giveNumberOfDelaunayVertices(),
                    this->giveLatticeBeam(i + 1)->giveDiameter() );
            fprintf(outputStream, "\n");
        } else if ( this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveLatticeBeam(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveReinforcementNode(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->givePeriodicNode();
                }
            }

            global_index++;
            fprintf(outputStream, "latticeBeam3Dboundary %d nodes 3 %d %d %d crossSect 1 mat 4 diameter %e ", global_index,
                    nodes.at(1) + this->giveNumberOfDelaunayVertices(),
                    nodes.at(2) + this->giveNumberOfDelaunayVertices(),
                    this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode() + 1,
                    this->giveLatticeBeam(i + 1)->giveDiameter() );
            fprintf(outputStream, " location 2 %d %d", location.at(1), location.at(2) );
            fprintf(outputStream, "\n");
        }
    }

    printf("Finished writing data for lattice truss (fibre) \n");

    // print latticelinks

    for ( int i = 0; i < this->giveNumberOfLatticeLinks(); i++ ) {
        if ( this->giveLatticeLink(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveLatticeLink(i + 1)->giveLocalVertices(nodes);
            global_index++;
            fprintf(outputStream, "latticeLink3D %d nodes 2 %d %d crossSect 1 mat 5 length %e diameter %e dirvector 3 %e %e %e L_end %e ",
                    global_index,
                    nodes.at(1) + this->giveNumberOfDelaunayVertices(),
                    nodes.at(2),
                    this->giveLatticeLink(i + 1)->giveAssociatedLength(),
                    this->giveLatticeLink(i + 1)->giveDiameter(),
                    ( this->giveLatticeLink(i + 1)->giveDirectionVector() ).at(1),
                    ( this->giveLatticeLink(i + 1)->giveDirectionVector() ).at(2),
                    ( this->giveLatticeLink(i + 1)->giveDirectionVector() ).at(3),
                    this->giveLatticeLink(i + 1)->giveL_end() );
            fprintf(outputStream, "\n");
        } else if ( this->giveLatticeLink(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveLatticeLink(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes

            if ( this->giveReinforcementNode(nodes.at(1) )->giveOutsideFlag() == 1 ) {
                location.at(1) = this->giveReinforcementNode(nodes.at(1) )->giveLocation();
                nodes.at(1) = this->giveReinforcementNode(nodes.at(1) )->givePeriodicNode();
            }
            if ( this->giveDelaunayVertex(nodes.at(2) )->giveOutsideFlag() == 1 ) {
                location.at(2) = this->giveDelaunayVertex(nodes.at(2) )->giveLocation();
                nodes.at(2) = this->giveDelaunayVertex(nodes.at(2) )->givePeriodicNode();
            }


            global_index++;
            fprintf(outputStream, "latticeLink3Dboundary %d nodes 3 %d %d %d crossSect 1 mat 5 length %e  diameter %e dirvector 3 %e %e %e L_end %e", global_index,
                    nodes.at(1) + this->giveNumberOfDelaunayVertices(),
                    nodes.at(2),
                    this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode() + 1,
                    this->giveLatticeLink(i + 1)->giveAssociatedLength(),
                    this->giveLatticeLink(i + 1)->giveDiameter(),
                    this->giveLatticeLink(i + 1)->giveDirectionVector().at(1),
                    this->giveLatticeLink(i + 1)->giveDirectionVector().at(2),
                    this->giveLatticeLink(i + 1)->giveDirectionVector().at(3),
                    this->giveLatticeLink(i + 1)->giveL_end() );

            fprintf(outputStream, " location 2 %d %d", location.at(1), location.at(2) );
            fprintf(outputStream, "\n");
        }
    }

    printf("Finished writing data for lattice links (fibre) \n");


    fprintf(outputStream, "latticecs 1 material 1\n");

    fprintf(outputStream, "latticeplastdam 1 d 0 talpha 0. calpha 0. e 50.46e9 a1 0.215 ft 1.92e6 fc 30.e6 wf 25.e-6 ahard 1.e-3 angle1 0.5 flow 0.25 randvars 2 806 807 randgen 2 4 4\n");
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStream, "NodalLoad 2 loadTimeFunction 1 dofs 6 31 32 33 40 41 42 Components 6 0. -1.5e4 0. 0. 0. 0.\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "PiecewiseLinFunction 2 nPoints 4 t 4  0. 1.1 1.2 201. f(t) 4 1. 1. 1. 1.\n");
    fprintf(outputStream, "PiecewiseLinFunction 3 nPoints 4 t 4  0. 1.1 1.2 201. f(t) 4 0. 0. 0. 1.\n");
    fprintf(outputStream, "InterpolatingFunction 4 name random.dat dim 3\n");
    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#NODE number %d dof 31 unknown d\n", this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode() + 1); // a verifier
    fprintf(outputStream, "##TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");
    return;
}



void
Grid::give3DFibreBenchmarkOutput(const std::string &fileName)
{
    //Output for 3D fracture process zone modelling

    FILE *outputStream = converter::fopen_or_die(fileName, "w");


    printf("Starting writing data in output file \n");

    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    int materialType = 1;
    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;
    oofem::FloatArray coordsOne, coordsTwo, centre, radii;
    double radius, distanceOne, distanceTwo;
    int isIn1, isIn2; // bool var to know whether endpoints of element qre inside or outide an inclusion
    double coord_x1, coord_y1, coord_z1, coord_x2, coord_y2, coord_z2;

    //Determine the number of nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    for ( int i = 0; i < this->giveNumberOfReinforcementNode(); i++ ) {
        if ( this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }


    //Determine the number of lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfLines++;
        }
    }
    for ( int i = 0; i < this->giveNumberOfLatticeBeams(); i++ ) {
        if ( this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 0 || this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfLines++;
        }
    }
    for ( int i = 0; i < this->giveNumberOfLatticeLinks(); i++ ) {
        if ( this->giveLatticeLink(i + 1)->giveOutsideFlag() == 0 || this->giveLatticeLink(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfLines++;
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Mechanical 3D model for fibres and ellipsoide\n");
    fprintf(outputStream, "NonLinearStatic nmsteps 2 nsteps 1 contextOutputStep 1 nmodules 2 profileopt 1\n");
    fprintf(outputStream, "nsteps 2 rtolv 1.e-3 reqIterations 100 stiffMode 1 manrmsteps 10 maxiter 200 controllmode 0 stepLength 1.e-6 minsteplength 1.e-6 maxrestarts 0 hpcmode 2 hpc 2 %d 2 hpcw 1 1. donotfixload\n", this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode() + 1); // a verifier
    fprintf(outputStream, "nsteps 100 rtolv 1.e-3 reqIterations 100 stiffMode 1 manrmsteps 10 maxiter 200 controllmode 0 stepLength 1.e-6 minsteplength 1.e-6 maxrestarts 0 hpcmode 2 hpc 2 %d 2 hpcw 1 1. donotfixload\n", this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode() + 1);   // a verifier
    fprintf(outputStream, "vtkxmldiscrete tstep_all domain_all\n");
    fprintf(outputStream, "gpexportmodule vars 5 59 90 84 85 78 tstep_all domain_all\n");
    fprintf(outputStream, "domain 3dLattice\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_output {%d}\n", this->giveNumberOfDelaunayVertices()
            + this->giveNumberOfReinforcementNode() + 1); // a verifier
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 4 nbc 2 nic 0 nltf 1\n", numberOfNodes + 1, numberOfLines);

    //Print nodes

    int firstFlag = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                fprintf(outputStream, "node %d coords 3 %e %e %e bc 6 1 1 1 1 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    // Print reinforcement nodes

    for ( int i = 0; i < this->giveNumberOfReinforcementNode(); i++ ) {
        if ( this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveReinforcementNode(i + 1)->giveCoordinates(coords);

            fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + this->giveNumberOfDelaunayVertices() + 1, coords.at(1), coords.at(2), coords.at(3) );
        }
    }


    //Periodic control node
    fprintf(outputStream, "node %d coords 3 %e %e %e load 1 2\n", this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );

    printf("Finished writing data for nodes \n");

    //Print Delaunay elements

    int global_index(0);//usefull to have a global numerotation of elements

    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        // take into account weak midplane to choose the materials of the element
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);

            //Deal with inclusions
            materialType = 1;
            this->giveDelaunayLine(i + 1)->updateMaterial(1);
            this->giveDelaunayVertex(nodes.at(1) )->giveCoordinates(coordsOne);
            this->giveDelaunayVertex(nodes.at(2) )->giveCoordinates(coordsTwo);

            //For fibre benchmark elements with cross-sections crossing the midplane should be given special material property
            if ( ( coordsOne.at(2) < 0.5 * specimenDimension.at(2) && coordsTwo.at(2) > 0.5 * specimenDimension.at(2) ) || ( coordsTwo.at(2) < 0.5 * specimenDimension.at(2) && coordsOne.at(2) > 0.5 * specimenDimension.at(2) ) ) {
                materialType = 2;//Weak interface
                this->giveDelaunayLine(i + 1)->updateMaterial(2);
            }
        }



        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            global_index++;
            fprintf(outputStream, "lattice3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", global_index, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );


            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, "\n");
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveDelaunayVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
            }
            //  materialType = 1;
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);
            global_index++;
            fprintf(outputStream, "lattice3Dboundary %d nodes 3 %d %d %d crossSect 1 mat %d polycoords %d", global_index, nodes.at(1), nodes.at(2), this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode() + 1, materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, " location 2 %d %d", location.at(1), location.at(2) );

            fprintf(outputStream, "\n");
        }
    }

    printf("Finished writing data for delaunay elements \n");

    // print latticebeams
    double myPi = 3.14159265;
    for ( int i = 0; i < this->giveNumberOfLatticeBeams(); i++ ) {
        if ( this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveLatticeBeam(i + 1)->giveLocalVertices(nodes);
            global_index++;
            fprintf(outputStream, "latticeBeam3D %d nodes 2 %d %d crossSect 1 mat 3 diameter %e", global_index,
                    nodes.at(1) + this->giveNumberOfDelaunayVertices(),
                    nodes.at(2) + this->giveNumberOfDelaunayVertices(),
                    this->giveLatticeBeam(i + 1)->giveDiameter() );
            fprintf(outputStream, "\n");
        } else if ( this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveLatticeBeam(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveReinforcementNode(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->givePeriodicNode();
                }
            }

            global_index++;
            fprintf(outputStream, "latticeBeam3Dboundary %d nodes 3 %d %d %d crossSect 1 mat 3 diameter %e ", global_index,
                    nodes.at(1) + this->giveNumberOfDelaunayVertices(),
                    nodes.at(2) + this->giveNumberOfDelaunayVertices(),
                    this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode() + 1,
                    this->giveLatticeBeam(i + 1)->giveDiameter() );
            fprintf(outputStream, " location 2 %d %d", location.at(1), location.at(2) );
            fprintf(outputStream, "\n");
        }
    }

    printf("Finished writing data for lattice truss (fibre) \n");

    // print latticelinks

    for ( int i = 0; i < this->giveNumberOfLatticeLinks(); i++ ) {
        if ( this->giveLatticeLink(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveLatticeLink(i + 1)->giveLocalVertices(nodes);
            global_index++;
            fprintf(outputStream, "latticeLink3D %d nodes 2 %d %d crossSect 1 mat 4 length %e diameter %e dirvector 3 %e %e %e L_end %e ",
                    global_index,
                    nodes.at(1) + this->giveNumberOfDelaunayVertices(),
                    nodes.at(2),
                    this->giveLatticeLink(i + 1)->giveAssociatedLength(),
                    this->giveLatticeLink(i + 1)->giveDiameter(),
                    ( this->giveLatticeLink(i + 1)->giveDirectionVector() ).at(1),
                    ( this->giveLatticeLink(i + 1)->giveDirectionVector() ).at(2),
                    ( this->giveLatticeLink(i + 1)->giveDirectionVector() ).at(3),
                    this->giveLatticeLink(i + 1)->giveL_end() );
            fprintf(outputStream, "\n");
        } else if ( this->giveLatticeLink(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveLatticeLink(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes

            if ( this->giveReinforcementNode(nodes.at(1) )->giveOutsideFlag() == 1 ) {
                location.at(1) = this->giveReinforcementNode(nodes.at(1) )->giveLocation();
                nodes.at(1) = this->giveReinforcementNode(nodes.at(1) )->givePeriodicNode();
            }
            if ( this->giveDelaunayVertex(nodes.at(2) )->giveOutsideFlag() == 1 ) {
                location.at(2) = this->giveDelaunayVertex(nodes.at(2) )->giveLocation();
                nodes.at(2) = this->giveDelaunayVertex(nodes.at(2) )->givePeriodicNode();
            }


            global_index++;
            fprintf(outputStream, "latticeLink3Dboundary %d nodes 3 %d %d %d crossSect 1 mat 4 length %e  diameter %e dirvector 3 %e %e %e L_end %e", global_index,
                    nodes.at(1) + this->giveNumberOfDelaunayVertices(),
                    nodes.at(2),
                    this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode() + 1,
                    this->giveLatticeLink(i + 1)->giveAssociatedLength(),
                    this->giveLatticeLink(i + 1)->giveDiameter(),
                    this->giveLatticeLink(i + 1)->giveDirectionVector().at(1),
                    this->giveLatticeLink(i + 1)->giveDirectionVector().at(2),
                    this->giveLatticeLink(i + 1)->giveDirectionVector().at(3),
                    this->giveLatticeLink(i + 1)->giveL_end() );

            fprintf(outputStream, " location 2 %d %d", location.at(1), location.at(2) );
            fprintf(outputStream, "\n");
        }
    }

    printf("Finished writing data for lattice links (fibre) \n");

    fprintf(outputStream, "simplecs 1\n");

    fprintf(outputStream, "latticedamage 1 talpha 0. d 0. e 30.e9 e0 100.e6 stype 3 wf 40.e6\n");
    fprintf(outputStream, "latticedamage 2 talpha 0. d 0. e 30.e9 e0 1.e-12 stype 3 wf 1.e-9\n");
    fprintf(outputStream, "latticelinearelastic 3 d 0. e 200.e9 n 0.\n");
    fprintf(outputStream, "latticeslip 4 talpha 0. d 0. e 30.e11 a1 1 a2 1 t0 4.e6\n");
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStream, "NodalLoad 2 loadTimeFunction 1 Components 6 0. 1. 0. 0. 0. 0.\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode() + 1); // a verifier
    fprintf(outputStream, "#LOADLEVEL\n");
    fprintf(outputStream, "##TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");
    return;
}



void
Grid::give3DWongOutput(const std::string &fileName)
{
    //Template for irregular nonperiodic coupled mechanical transport models. Do not change for applications

    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    const std::string fileName1 = fileName + ".sm";
    FILE *outputStreamSM = converter::fopen_or_die(fileName1, "w");

    const std::string fileName2 = fileName + ".tm";
    FILE *outputStreamTM = converter::fopen_or_die(fileName2, "w");


    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfDelaunayNodes, numberOfDelaunayLines;
    int numberOfVoronoiNodes, numberOfVoronoiLines;

    oofem::FloatArray coords;
    int materialType = 1;

    oofem::IntArray nodes;
    oofem::IntArray crossSectionNodes;
    oofem::IntArray crossSectionElements;
    oofem::IntArray location(2);
    oofem::FloatArray coordsOne, coordsTwo, centre;
    double radius, distanceOne, distanceTwo;

    //*****************************
    //Write the control file
    //*****************************
    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Control file for percolation application\n");
    fprintf(outputStream, "StaggeredProblem nsteps 100 deltat 1. prob1  \"oofem.in.tm\" prob2 \"oofem.in.sm\" contextoutputstep 100 coupling 3 2 1 0\n");
    fprintf(outputStream, "##%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "##INCLUDE oofem.in.tm\n");
    fprintf(outputStream, "##%%END_CHECK%%\n");


    //Determine the number of Delaunay nodes in the domain
    numberOfDelaunayNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfDelaunayNodes++;
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfDelaunayLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfDelaunayLines++;
        }
    }

    //****************************************
    //Write the output file for the SM part
    //******************************************
    fprintf(outputStreamSM, "oofem.out.sm\n");
    fprintf(outputStreamSM, "Mechanical part of percolation model\n");
    fprintf(outputStreamSM, "NonLinearStatic nmsteps 1 nsteps 1 contextOutputStep 1000 nmodules 2 profileopt 1 lstype 3 smtype 7\n");
    fprintf(outputStreamSM, "nsteps 100 rtolv 1.e-3 stiffMode 1 manrmsteps 10 maxiter 200 controllmode 1 minsteplength 1.e-10 maxrestarts 0 lstype 3 smtype 7\n");
    fprintf(outputStreamSM, "vtkxml primvars 1 1 tstep_all domain_all\n");
    fprintf(outputStreamSM, "gpexportmodule vars 5 59 90 84 85 78 tstep_all domain_all\n");
    fprintf(outputStreamSM, "domain 3dLattice\n");
    fprintf(outputStreamSM, "OutputManager tstep_all dofman_output {%d}\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStreamSM, "ndofman %d nelem %d ncrosssect 1 nmat 3 nbc 2 nic 0 nltf 2\n", numberOfDelaunayNodes + 1, numberOfDelaunayLines);

    int firstFlag = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                fprintf(outputStreamSM, "node %d coords 3 %e %e %e bc 6 1 1 1 1 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStreamSM, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }
    //Periodic control node
    fprintf(outputStreamSM, "node %d coords 3 %e %e %e\n", this->giveNumberOfDelaunayVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );

    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);

            //Deal with inclusions
            materialType = 1;
            this->giveDelaunayVertex(nodes.at(1) )->giveCoordinates(coordsOne);
            this->giveDelaunayVertex(nodes.at(2) )->giveCoordinates(coordsTwo);

            for ( int m = 0; m < this->giveNumberOfInclusions(); m++ ) {
                //Distinguish between inclusions
                if ( !strcmp(this->giveInclusion(m + 1)->giveClassName(), "Sphere") ) {
                    ( ( InterfaceSphere * ) this->giveInclusion(m + 1) )->giveCentre(centre);
                    radius = ( ( InterfaceSphere * ) this->giveInclusion(m + 1) )->giveRadius() + ( ( InterfaceSphere * ) this->giveInclusion(m + 1) )->giveITZThickness() / 2.;
                    distanceOne = sqrt(pow(coordsOne.at(1) - centre.at(1), 2.) +
                                       pow(coordsOne.at(2) - centre.at(2), 2.) +
                                       pow(coordsOne.at(3) - centre.at(3), 2.) );
                    distanceTwo = sqrt(pow(coordsTwo.at(1) - centre.at(1), 2.) +
                                       pow(coordsTwo.at(2) - centre.at(2), 2.) +
                                       pow(coordsTwo.at(3) - centre.at(3), 2.) );

                    if ( distanceOne > radius && distanceTwo > radius ) {
                        materialType = 1;
                    } else if ( ( distanceOne > radius && distanceTwo < radius ) || ( distanceOne < radius && distanceTwo > radius ) ) {
                        materialType = 3;
                        break;
                    } else {
                        materialType = 2;
                        break;
                    }
                }
            }
        }

        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);
            //Need to introduce both normal and boundary element. See FPZ output file.
            fprintf(outputStreamSM, "lattice3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStreamSM, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            this->giveDelaunayLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            fprintf(outputStreamSM, " couplingflag 1 couplingnumber %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( giveVoronoiLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    fprintf(outputStreamSM, "%d ", giveVoronoiLine(crossSectionElements.at(m + 1) )->givePeriodicElement() );
                } else {
                    fprintf(outputStreamSM, "%d ", crossSectionElements.at(m + 1) );
                }
            }
            fprintf(outputStreamSM, "bodyloads 1 2\n");
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveDelaunayVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
            }
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStreamSM, "lattice3Dboundary %d nodes 3 %d %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), this->giveNumberOfDelaunayVertices() + 1, materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStreamSM, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            this->giveDelaunayLine(i + 1)->giveCrossSectionElements(crossSectionElements);
            fprintf(outputStreamSM, " couplingflag 1 couplingnumber %d ", crossSectionElements.giveSize() );
            for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                if ( giveVoronoiLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    fprintf(outputStreamSM, "%d ", giveVoronoiLine(crossSectionElements.at(m + 1) )->givePeriodicElement() );
                } else {
                    fprintf(outputStreamSM, "%d ", crossSectionElements.at(m + 1) );
                }
            }
            fprintf(outputStreamSM, " location 2 %d %d ", location.at(1), location.at(2) );
            fprintf(outputStreamSM, "bodyloads 1 2\n");
        }
    }
    fprintf(outputStreamSM, "simplecs 1\n");

    fprintf(outputStreamSM, "latticeplastdam 1 d 0 talpha -5.e-3 e 40.e9 ft 6.5e6 fc 65.e6 wf 15.4e-6 ahard 1.e-3\n");
    fprintf(outputStreamSM, "latticeplastdam 2 d 0 talpha 0. e 100.e9 ft 1.e16 fc 30.e16 wf 15.4e+6 ahard 1.e-3\n");
    fprintf(outputStreamSM, "latticeplastdam 3 d 0 talpha -5.e-3 e 57.1e9 a1 0.01 ft 3.25e6 fc 32.5e6 wf 15.4e-6 ahard 1.e-3\n");

    fprintf(outputStreamSM, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStreamSM, "StructTemperatureLoad 2 loadTimeFunction 2 Components 1 1.0\n");
    fprintf(outputStreamSM, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStreamSM, "PiecewiseLinFunction 2 nPoints 2 t 2 0. 101. f(t) 2 0. 1.\n");

    //Determine the nodes and elements inside the specimen

    //Determine the number of Voronoi nodes in the domain
    numberOfVoronoiNodes = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfVoronoiNodes++;
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfVoronoiLines = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfVoronoiLines++;
        }
    }


    //*******************************************
    //Write the input file of the transport model
    //*******************************************
    fprintf(outputStreamTM, "oofem.out.tm\n");
    fprintf(outputStreamTM, "Transport part of the percolation model\n");
    fprintf(outputStreamTM, "nltransienttransportproblem nsteps 100 deltat 1.0 rtol 0.001 alpha 1. nsmax 200 profileopt 1 contextOutputStep 100 nmodules 0\n");
    fprintf(outputStreamTM, "domain 3dMassLatticeTransport\n");
    fprintf(outputStreamTM, "OutputManager tstep_all dofman_output {%d}\n", giveNumberOfVoronoiVertices() + 1);
    fprintf(outputStreamTM, "ndofman %d nelem %d ncrosssect 1 nmat 3 nbc 2 nic 0 nltf 1\n", numberOfVoronoiNodes + 1, numberOfVoronoiLines);

    firstFlag = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                fprintf(outputStreamTM, "node %d coords 3 %e %e %e bc 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStreamTM, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    //Periodic control node
    fprintf(outputStreamTM, "node %d coords 3 %e %e %e ndofs 3 dofIDmask 3 1 2 3 bc 3 1 1 2\n", this->giveNumberOfVoronoiVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );


    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);

            //Deal with inclusions
            materialType = 1;
            this->giveVoronoiVertex(nodes.at(1) )->giveCoordinates(coordsOne);
            this->giveVoronoiVertex(nodes.at(2) )->giveCoordinates(coordsTwo);

            for ( int m = 0; m < this->giveNumberOfInclusions(); m++ ) {
                //Distinguish between inclusions
                if ( !strcmp(this->giveInclusion(m + 1)->giveClassName(), "Sphere") ) {
                    ( ( InterfaceSphere * ) this->giveInclusion(m + 1) )->giveCentre(centre);
                    radius = ( ( InterfaceSphere * ) this->giveInclusion(m + 1) )->giveRadius() + ( ( InterfaceSphere * ) this->giveInclusion(m + 1) )->giveITZThickness() / 2.;
                    distanceOne = sqrt(pow(coordsOne.at(1) - centre.at(1), 2.) +
                                       pow(coordsOne.at(2) - centre.at(2), 2.) +
                                       pow(coordsOne.at(3) - centre.at(3), 2.) );
                    distanceTwo = sqrt(pow(coordsTwo.at(1) - centre.at(1), 2.) +
                                       pow(coordsTwo.at(2) - centre.at(2), 2.) +
                                       pow(coordsTwo.at(3) - centre.at(3), 2.) );

                    if ( distanceOne > radius && distanceTwo > radius ) {
                        materialType = 1;
                    } else if ( ( distanceOne > radius && distanceTwo < radius ) || ( distanceOne < radius && distanceTwo > radius ) ) {
                        materialType = 3;
                        break;
                    } else {
                        materialType = 2;
                        break;
                    }
                }
            }


            if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) {
                this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

                fprintf(outputStreamTM, "latticemt3D %d nodes 2 %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, 3 * crossSectionNodes.giveSize() );

                for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                    this->giveDelaunayVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                    fprintf(outputStreamTM, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
                }

                this->giveVoronoiLine(i + 1)->giveCrossSectionElements(crossSectionElements);
                fprintf(outputStreamTM, " couplingflag 1 couplingnumber %d ", crossSectionElements.giveSize() );
                for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                    if ( giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 1 ) {
                        fprintf(outputStreamTM, "%d ", giveDelaunayLine(crossSectionElements.at(m + 1) )->givePeriodicElement() );
                    } else {
                        fprintf(outputStreamTM, "%d ", crossSectionElements.at(m + 1) );
                    }
                }
                fprintf(outputStreamTM, "mlength 1.e-8");
                fprintf(outputStreamTM, "\n");
            } //Element is inside
            else if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {  //Element crosses the boundary
                location.zero();
                this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
                //Go through nodes and replace the ones outside with periodic nodes
                for ( int m = 0; m < 2; m++ ) {
                    if ( this->giveVoronoiVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                        location.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->giveLocation();
                        nodes.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                }
                this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

                fprintf(outputStreamTM, "latticemt3Dboundary %d nodes 3 %d %d %d crossSect 1 mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), this->giveNumberOfVoronoiVertices() + 1, materialType, 3 * crossSectionNodes.giveSize() );

                for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                    this->giveDelaunayVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                    fprintf(outputStreamTM, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
                }
                this->giveVoronoiLine(i + 1)->giveCrossSectionElements(crossSectionElements);
                fprintf(outputStreamTM, " couplingflag 1 couplingnumber %d ", crossSectionElements.giveSize() );
                for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                    if ( giveDelaunayLine(crossSectionElements.at(m + 1) )->giveOutsideFlag() == 1 ) {
                        fprintf(outputStreamTM, "%d ", giveDelaunayLine(crossSectionElements.at(m + 1) )->givePeriodicElement() );
                    } else {
                        fprintf(outputStreamTM, "%d ", crossSectionElements.at(m + 1) );
                    }
                }
                fprintf(outputStreamTM, " location 2 %d %d ", location.at(1), location.at(2) );
                fprintf(outputStreamTM, "mlength 1.e-8");
                fprintf(outputStreamTM, "\n");
            }
        }
    }

    fprintf(outputStreamTM, "simplecs 1\n");
    fprintf(outputStreamTM, "latticetransmat 1 d 1.e3 k 1.e-19 vis 1.e-3 thetas 1. thetar 0. contype 0 c 0. ctor 0.001\n");
    fprintf(outputStreamTM, "latticetransmat 2 d 1.e3 k 1.e-22 vis 1.e-3 thetas 1. thetar 0. contype 0 c 0. ctor 0.001\n");
    fprintf(outputStreamTM, "latticetransmat 3 d 1.e3 k 1.e-19 vis 1.e-3 thetas 1. thetar 0. contype 0 c 0. ctor 0.001\n");
    fprintf(outputStreamTM, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.\n");
    fprintf(outputStreamTM, "BoundaryCondition 2 loadTimeFunction 1 prescribedvalue 1.\n");
    fprintf(outputStreamTM, "ConstantFunction 1 f(t) 1.\n");

    fprintf(outputStreamTM, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStreamTM, "#TIME\n");
    fprintf(outputStreamTM, "#REACTION number %d dof 3\n", giveNumberOfVoronoiVertices() + 1);
    fprintf(outputStreamTM, "#%%END_CHECK%%\n");

    return;
}


void
Grid::give3DPeriodicPoreTMOutput(const std::string &fileName)
{
    //Template for irregular periodic mechanical models. Do not change for applications

    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    const std::string fileName1 = fileName + "_lengths.dat";
    FILE *outputStreamLengths = converter::fopen_or_die(fileName1, "w");


    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    int materialType = 1;
    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;
    int numberOfPeriodicLines;
    int numberOfPipeDiameters;

    //Determine the number of Voronoi nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    //Determine the total number and the number of the periodic Voronoi lines in the domain
    numberOfLines = 0;
    numberOfPeriodicLines = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        oofem::IntArray helpNodes(2), locationArray(2), switches(2);
        oofem::FloatArray helpCoords(3);

        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 3 ) {
            this->giveVoronoiLine(i + 1)->giveLocalVertices(helpNodes);
            for ( int m = 0; m < 2; m++ ) {
                giveVoronoiVertex(helpNodes.at(m + 1) )->giveCoordinates(helpCoords);
                locationArray.at(m + 1) = this->giveRegion(1)->giveSwitches(switches, helpCoords);
            }
        }
        //If they do, one of the switches array from the nodes is used to shift the line

        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 3 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 4 ) {
            if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 || ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 3 && locationArray.at(1) == locationArray.at(2) ) ) {
                numberOfPeriodicLines++;
            }
            numberOfLines++;
        }
    }

    //Check if the number of the periodic Voronoi Lines is odd.
    //The generation continues if the number is even, otherwise it exits.
    if ( numberOfPeriodicLines % 2 != 0 ) {
        converter::error("In give3DPeriodicPoreTMOutput: The number of perdiodic lines is not even");
    }

    //Determine the number of the diameter values for the pipes that will have to be generated.
    numberOfPipeDiameters = numberOfLines - numberOfPeriodicLines / 2;

    //Create first random numbers for vertices and lines
    //Create random numbers for vertices
    oofem::FloatArray voronoiVertexRadius(numberOfNodes);
    long randomIntegerOne = this->randomInteger;
    double randomRadius = 0.;

    double gaussianPoreMean, gaussianPoreSTD, gaussianPoreCOV;

    gaussianPoreMean = log(this->poreMean / sqrt(1 + pow(this->poreCOV, 2) ) );
    gaussianPoreSTD = sqrt(log(1 + pow(this->poreCOV, 2) ) );
    gaussianPoreCOV = gaussianPoreSTD / gaussianPoreMean;


    if ( gaussianPoreMean < 0. ) {
        gaussianPoreSTD = -gaussianPoreSTD;
    }

    for ( int i = 0; i < numberOfNodes; i++ ) {
        randomRadius = normalCdfInverse(ran1(& randomIntegerOne), gaussianPoreMean, gaussianPoreSTD);

        //Apply cut-offs if necessary
        if ( randomRadius > poreMax ) {
            randomRadius = poreMax;
            //maxCounter++;
        } else if ( randomRadius < poreMin ) {
            randomRadius = poreMin;
            //minCounter++;
        }
        voronoiVertexRadius.at(i + 1) = randomRadius;
    }
    //Create random numbers for lines
    oofem::FloatArray voronoiLineRadius(numberOfPipeDiameters);
    long randomIntegerTwo = this->randomInteger - 1;


    double gaussianThroatMean, gaussianThroatSTD, gaussianThroatCOV;

    gaussianThroatMean = log(this->throatMean / sqrt(1 + pow(this->throatCOV, 2) ) );
    gaussianThroatSTD = sqrt(log(1 + pow(this->throatCOV, 2) ) );
    gaussianThroatCOV = gaussianThroatSTD / gaussianThroatMean;


    if ( gaussianThroatMean < 0. ) {
        gaussianThroatSTD = -gaussianThroatSTD;
    }

    //Create random numbers for lines
    //oofem::FloatArray voronoiLineRadius(numberOfPipeDiameters);
    //long randomIntegerTwo = this->randomInteger - 1;
    for ( int i = 0; i < numberOfPipeDiameters; i++ ) {
        randomRadius =  normalCdfInverse(ran1(& randomIntegerTwo),  gaussianThroatMean, gaussianThroatSTD);

        //Apply cut-offs if necessary
        if ( randomRadius > throatMax ) {
            randomRadius = throatMax;
        } else if ( randomRadius < throatMin ) {
            randomRadius = throatMin;
        }
        voronoiLineRadius.at(i + 1) = randomRadius;
    }


    /*Strategy to apply random numbers to vertices and elements
     * 1. Sort the line radii (smallest first)
     * 2. Apply the pores radii to pores randomly. This is different to our original idea, but easier to implement.
     * 3. Go trough all elements and find the maximum pore radius for each element. One element has two pores attached to it.
     * 4. Rank the elements according to their maximum pore radius.
     * 5. Apply the smallest element radius to the element with the smallest maximum pore radius using the rank table.
     * Any geometrical violations are not correct. Hopefully there will not be many and will not influence the results too much.
     */

    /*@todo:
     * 1. Check how many violations of the pipe to sphere radii occur.
     */

    //Sort the line radii according to size
    oofem::FloatArray sortedVoronoiLineRadius(numberOfPipeDiameters);
    sortRandomNumbers(sortedVoronoiLineRadius, voronoiLineRadius);

    //Apply the vertex radii to all vertices.
    //This is random as vertices were not sorted and randomly generated.

    int counter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            counter++;
            this->giveVoronoiVertex(i + 1)->setRadius(voronoiVertexRadius.at(counter) );
        }
    }

    //Got through all elements and find largest pore
    oofem::IntArray localVertices(2);
    double radius = 0.;
    oofem::FloatArray minRadius(numberOfPipeDiameters);
    counter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        oofem::IntArray helpNodes(2), locationArray(2), switches(2);
        oofem::FloatArray helpCoords(3);

        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 3 ) {
            this->giveVoronoiLine(i + 1)->giveLocalVertices(helpNodes);
            for ( int m = 0; m < 2; m++ ) {
                giveVoronoiVertex(helpNodes.at(m + 1) )->giveCoordinates(helpCoords);
                locationArray.at(m + 1) = this->giveRegion(1)->giveSwitches(switches, helpCoords);
            }
        }

        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 4 ) {
            counter++;
            minRadius.at(counter) = 1000000.;
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            for ( int m = 0; m < 2; m++ ) {
                //new
                if ( this->giveVoronoiVertex(nodes.at(m + 1) )->giveOutsideFlag() == 2 ) {
                    nodes.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
                radius = this->giveVoronoiVertex(nodes.at(m + 1) )->giveRadius();
                if ( radius < minRadius.at(counter) ) {
                    minRadius.at(counter) = radius;
                }
            }
        } else if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 || ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 3 && locationArray.at(1) == locationArray.at(2) ) ) {
            //A flag was implemented in line.h.
            //This flag shows whether the periodic friend of the line under consideration has been taken into account for the rank vector.
            //flag = 0 if not passed, 1 if passed.
            //Only one of the two periodic friends have to be taken into account.
            if ( this->giveVoronoiLine(i + 1)->givePeriodicElement() > i + 1 ) {
                counter++;
                minRadius.at(counter) = 1000000.;
                this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
                //Go through nodes and replace the ones outside with periodic nodes
                for ( int m = 0; m < 2; m++ ) {
                    if ( this->giveVoronoiVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 || this->giveVoronoiVertex(nodes.at(m + 1) )->giveOutsideFlag() == 2 ) {
                        location.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->giveLocation();
                        nodes.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                    radius = this->giveVoronoiVertex(nodes.at(m + 1) )->giveRadius();
                    if ( radius < minRadius.at(counter) ) {
                        minRadius.at(counter) = radius;
                    }
                }
            }
        }
    }

    int help;
    help = 1;
    //Need to sort it now.
    //Firstly create a rank table
    //Create indexx talble for radii
    oofem::IntArray rankVector(numberOfPipeDiameters);
    createRankTable(rankVector, minRadius);

    //Apply sorted line radii to the elements

    double helpRadius = 0;
    counter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) {
            counter++;
            helpRadius = sortedVoronoiLineRadius.at(rankVector.at(counter) );
            this->giveVoronoiLine(i + 1)->setRadius(helpRadius);
        } else if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 3 ) {
            //Here only the periodic lines that have been used for generating the rank vector are allowed to pass.
            if ( this->giveVoronoiLine(i + 1)->givePeriodicElement() > i + 1 ) {
                counter++;
                helpRadius = sortedVoronoiLineRadius.at(rankVector.at(counter) );
                this->giveVoronoiLine(i + 1)->setRadius(helpRadius);
                //Assign an equal diameter value to the periodic friend.
                this->giveVoronoiLine(this->giveVoronoiLine(i + 1)->givePeriodicElement() )->setRadius(helpRadius);
            }
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Transport discrete 3D model\n");
    fprintf(outputStream, "discretetransportproblem Nsteps %d deltat 1. homogenFlag 1 configurationflag 0  rtol 1.e-3 alpha 1. nsmax 2000 contextOutputStep 100000 profileopt 1 nmodules 0\n", 2 * numberOfNodes);
    fprintf(outputStream, "domain 3dMassLatticeTransport\n");
    fprintf(outputStream, "OutputManager tstep_all element_all dofman_output {%d}\n", this->giveNumberOfVoronoiVertices() + 1);
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 1\n", numberOfNodes + 1, numberOfLines);

    int firstFlag = 0;
    counter = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                counter++;
                fprintf(outputStream, "pore %d coords 3 %e %e %e bc 1 1 rad %.16e\n", i + 1, coords.at(1), coords.at(2), coords.at(3), voronoiVertexRadius.at(counter) );
            } else {
                counter++;
                fprintf(outputStream, "pore %d coords 3 %e %e %e rad %.16e\n", i + 1, coords.at(1), coords.at(2), coords.at(3), voronoiVertexRadius.at(counter) );
            }
        }
    }

    //Periodic control node
    fprintf(outputStream, "node %d coords 3 %e %e %e ndofs 3 dofIDmask 3 1 2 3 bc 3 2 1 1\n", this->giveNumberOfVoronoiVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );


    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            materialType = 1;
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);
            fprintf(outputStream, "latticemt3D_Discrete %d nodes 2 %d %d crossSect 1 mat %d rad %.16e ", i + 1, nodes.at(1), nodes.at(2), materialType, this->giveVoronoiLine(i + 1)->giveRadius() );
            fprintf(outputStream, "\n");
        } else if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {      //Element crosses the boundary
            location.zero();
            this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < 2; m++ ) {
                if ( this->giveVoronoiVertex(nodes.at(m + 1) )->giveOutsideFlag() == 1 ) {
                    location.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->giveLocation();
                    nodes.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1) )->givePeriodicNode();
                }
            }
            materialType = 1;
            this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "latticemt3Dboundary_Discrete %d nodes 3 %d %d %d crossSect 1 mat %d rad %.16e ", i + 1, nodes.at(1), nodes.at(2), this->giveNumberOfVoronoiVertices() + 1, materialType, this->giveVoronoiLine(i + 1)->giveRadius() );
            fprintf(outputStream, " location 2 %d %d", location.at(1), location.at(2) );
            fprintf(outputStream, "\n");
        }
    }

    double maxSuction;
    maxSuction = 2. * 1. / sortedVoronoiLineRadius.at(1);

    fprintf(outputStream, "simplecs 1\n");
    fprintf(outputStream, "latticelammat 1 d 1000. vis 0.001002 rate %e poremean %e pipemin 0.01e-9 poremin 0.1e-9\n", this->deltarad, this->poreMean);
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.\n");
    fprintf(outputStream, "BoundaryCondition 2 loadTimeFunction 1 prescribedvalue 1.\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "PiecewiseLinFunction 2 nPoints 6 t 6 -1. 0.  %d %d %d %d f(t) 6 0. 0. %e %e %e 0.\n", numberOfNodes / 2, numberOfNodes / 2 + 1, numberOfNodes, 2 * numberOfNodes - 1,  2. * ( 1 - .0001 ) / this->throatMean, 2. * ( 1 + .0001 ) / this->throatMean,  2 * maxSuction);
    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", this->giveNumberOfVoronoiVertices() + 1);
    fprintf(outputStream, "#LOADLEVEL\n");
    fprintf(outputStream, "##TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");
    return;
}


void
Grid::give3DSphereOutput(const std::string &fileName)
{
    //Output for hydraulic fracture of sphere.
    //Start with mechanical model

    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords(3);
    oofem::FloatArray coordsOne, coordsTwo, centre;
    int materialType = 1;
    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;
    double radius, distanceOne, distanceTwo;

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ||  this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) {
            if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) {
                this->giveRegion(1)->modifyVoronoiCrossSection(i + 1);
            }
            numberOfLines++;
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Mechanical 3D model\n");
    fprintf(outputStream, "NonLinearStatic nmsteps 1 nsteps 1 contextOutputStep 2000 nmodules 2 profileopt 1 lstype 3 smtype 7\n");
    fprintf(outputStream, "nsteps 1500 rtolv 1.e-3 stiffMode 2 manrmsteps 10 maxiter 200 controllmode 1 lstype 3 smtype 7\n");
    fprintf(outputStream, "vtkxml primvars 1 1 tstep_all domain_all\n");
    fprintf(outputStream, "gpexportmodule vars 5 59 90 84 85 78 tstep_all domain_all\n");
    fprintf(outputStream, "domain 3dLattice\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_output {%d}\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 3 nmat 3 nbc 2 nic 0 nltf 2\n", numberOfNodes, numberOfLines);

    printf("start nodes\n");

    int firstFlag = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                fprintf(outputStream, "node %d coords 3 %e %e %e bc 6 1 1 1 1 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    printf("finished nodes\n");

    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) { //Elements are inside
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);

            //=======================================================
            //Deal with inclusions
            materialType = 1;
            this->giveDelaunayLine(i + 1)->updateMaterial(1);
            this->giveDelaunayVertex(nodes.at(1) )->giveCoordinates(coordsOne);
            this->giveDelaunayVertex(nodes.at(2) )->giveCoordinates(coordsTwo);


            for ( int m = 0; m < this->giveNumberOfInclusions(); m++ ) {
                //Distinguish between inclusions
                if ( !strcmp(this->giveInclusion(m + 1)->giveClassName(), "Sphere") ) {
                    ( ( InterfaceSphere * ) this->giveInclusion(m + 1) )->giveCentre(centre);
                    radius = ( ( InterfaceSphere * ) this->giveInclusion(m + 1) )->giveRadius() + ( ( InterfaceSphere * ) this->giveInclusion(m + 1) )->giveITZThickness() / 2.;
                    distanceOne = sqrt(pow(coordsOne.at(1) - centre.at(1), 2.) +
                                       pow(coordsOne.at(2) - centre.at(2), 2.) +
                                       pow(coordsOne.at(3) - centre.at(3), 2.) );
                    distanceTwo = sqrt(pow(coordsTwo.at(1) - centre.at(1), 2.) +
                                       pow(coordsTwo.at(2) - centre.at(2), 2.) +
                                       pow(coordsTwo.at(3) - centre.at(3), 2.) );

                    if ( distanceOne > radius && distanceTwo > radius ) {
                        materialType = 1;
                        this->giveDelaunayLine(i + 1)->updateMaterial(1);
                    } else if ( ( distanceOne > radius && distanceTwo < radius ) || ( distanceOne < radius && distanceTwo > radius ) ) {
                        materialType = 3;
                        this->giveDelaunayLine(i + 1)->updateMaterial(3);
                        break;
                    } else {
                        materialType = 2;
                        this->giveDelaunayLine(i + 1)->updateMaterial(2);
                        break;
                    }
                }
            }
            //=================================================


            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            if ( materialType == 3 ) {
                fprintf(outputStream, "lattice3D %d nodes 2 %d %d crossSect %d mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, materialType, 3 * crossSectionNodes.giveSize() );

                for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                    this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                    fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
                }
                fprintf(outputStream, " bodyloads 1 2\n");
            } else   {
                fprintf(outputStream, "lattice3D %d nodes 2 %d %d crossSect %d mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, materialType, 3 * crossSectionNodes.giveSize() );

                for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                    this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                    fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
                }
                fprintf(outputStream, "\n");
            }
        }
    }
    fprintf(outputStream, "simplecs 1\n");
    fprintf(outputStream, "simplecs 2\n");
    fprintf(outputStream, "simplecs 3\n");

    fprintf(outputStream, "latticeplastdam 1 d 0 talpha 0. calpha 0. e 30.e9 ft 3.e6 fc 30.e6 wf 50.e-6 ahard 1.e-3\n");
    fprintf(outputStream, "latticeplastdam 2 d 0 talpha 0. calpha 0. e 300.e9 ft 1.e16 fc 30.e16 wf 15.4e+6 ahard 1.e-3\n");
    fprintf(outputStream, "talpha 0. calpha 0.15e-3 e 30.e9 a1 0.001 ft 3.e6 fc 30.e6 wf 50.e-6 ahard 1.e-3\n");
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStream, "StructTemperatureLoad 2 loadTimeFunction 2 Components 1 1.0\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "PiecewiseLinFunction 2 nPoints 2 t 2 0. 1499. f(t) 2 0. 1.\n");
    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "#LOADLEVEL\n");
    fprintf(outputStream, "##TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");
    return;
}


void
Grid::give3DCylinderOutput(const std::string &fileName)
{
    //Output for hydraulic fracture of cylinder.
    //Start with mechanical model

    FILE *outputStream = converter::fopen_or_die(fileName, "w");


    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords(3);
    oofem::FloatArray coordsOne, coordsTwo, line;
    int materialType = 1;
    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;
    double radius, distanceOne, distanceTwo;



    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ||  this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 && ( this->giveDelaunayLine(i + 1) )->delaunayAreaCheck() == 1 ) {
            numberOfLines++;
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 && this->giveRegion(1)->modifyVoronoiCrossSection(i + 1) == 1 )      {
            numberOfLines++;
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Mechanical 3D model of cylinder\n");
    fprintf(outputStream, "NonLinearStatic nsteps 1500 contextOutputStep 2000 nmodules 2 updateelasticstiffnessflag deltatfunction 3 rtolv 1.e-3 stiffMode 1 manrmsteps 10 maxiter 200 controllmode 1 lstype 4 smtype 8\n");
    fprintf(outputStream, "vtkxmllattice primvars 1 1 cellvars 4 46 60 90 111 tstep_all domain_all cross 1\n");
    fprintf(outputStream, "gpexportmodule vars 2 46 139 tstep_all\n");
    fprintf(outputStream, "domain 3dLattice\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_output {%d}\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 3 nmat 4 nbc 2 nic 0 nltf 4 nset 3\n", numberOfNodes, numberOfLines);

    printf("start nodes\n");

    int firstFlag = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                fprintf(outputStream, "node %d coords 3 %e %e %e bc 6 1 1 1 1 0 0\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
                i++;
                this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
                fprintf(outputStream, "node %d coords 3 %e %e %e bc 6 0 1 1 1 0 0\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    printf("finished nodes\n");

    int set1Counter = 0, set2Counter = 0, set3Counter = 0;
    oofem::IntArray set1Temp(numberOfLines);
    set1Temp.zero();

    oofem::IntArray set2Temp(numberOfLines);
    set2Temp.zero();

    oofem::IntArray set3Temp(numberOfLines);
    set3Temp.zero();

    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) && ( this->giveDelaunayLine(i + 1) )->delaunayAreaCheck() == 1 ) { //Elements are inside
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);

            //=======================================================
            //Deal with inclusions
            materialType = 2;
            this->giveDelaunayLine(i + 1)->updateMaterial(2);
            this->giveDelaunayVertex(nodes.at(1) )->giveCoordinates(coordsOne);
            this->giveDelaunayVertex(nodes.at(2) )->giveCoordinates(coordsTwo);


            for ( int m = 0; m < this->giveNumberOfInclusions(); m++ ) {
                //Distinguish between inclusions
                if ( !strcmp(this->giveInclusion(m + 1)->giveClassName(), "InterfaceCylinder") ) {
                    ( ( InterfaceCylinder * ) this->giveInclusion(m + 1) )->giveLine(line);
                    radius = ( ( InterfaceCylinder * ) this->giveInclusion(m + 1) )->giveRadius() + ( ( InterfaceCylinder * ) this->giveInclusion(m + 1) )->giveITZThickness() / 2.;


                    distanceOne = sqrt(pow(coordsOne.at(2) - line.at(2), 2.) +
                                       pow(coordsOne.at(3) - line.at(3), 2.) );
                    distanceTwo = sqrt(pow(coordsTwo.at(2) - line.at(2), 2.) +
                                       pow(coordsTwo.at(3) - line.at(3), 2.) );

                    if ( distanceOne > radius && distanceTwo > radius ) {
                        set1Counter++;
                        set1Temp.at(set1Counter) = i + 1;
                        materialType = 2;
                        this->giveDelaunayLine(i + 1)->updateMaterial(2);
                    } else if ( ( distanceOne > radius && distanceTwo < radius ) || ( distanceOne < radius && distanceTwo > radius ) ) {
                        set3Counter++;
                        set3Temp.at(set3Counter) = i + 1;
                        materialType = 4;
                        this->giveDelaunayLine(i + 1)->updateMaterial(4);
                        break;
                    } else {
                        set2Counter++;
                        set2Temp.at(set2Counter) = i + 1;
                        materialType = 3;
                        this->giveDelaunayLine(i + 1)->updateMaterial(3);
                        break;
                    }
                }
            }
            //=================================================

            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            if ( materialType == 4 ) {
                fprintf(outputStream, "lattice3D %d nodes 2 %d %d crossSect %d mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType - 1, materialType, 3 * crossSectionNodes.giveSize() );

                for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                    this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                    fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
                }
                fprintf(outputStream, " bodyloads 1 2\n");
            } else   {
                fprintf(outputStream, "lattice3D %d nodes 2 %d %d crossSect %d mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType - 1, materialType, 3 * crossSectionNodes.giveSize() );

                for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                    this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                    fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
                }
                fprintf(outputStream, "\n");
            }
        }
    }

    oofem::IntArray set1(set1Counter);
    set1.zero();

    oofem::IntArray set2(set2Counter);
    set2.zero();

    oofem::IntArray set3(set3Counter);
    set3.zero();

    for (int i = 0; i < set1Counter; i++) {
        set1.at(i + 1) = set1Temp.at(i + 1);
    }

    for (int i = 0; i < set2Counter; i++) {
        set2.at(i + 1) = set2Temp.at(i + 1);
    }

    for (int i = 0; i < set3Counter; i++) {
        set3.at(i + 1) = set3Temp.at(i + 1);
    }

    fprintf(outputStream, "latticecs 1 material 2\n");
    fprintf(outputStream, "latticecs 2 material 3\n");
    fprintf(outputStream, "latticecs 3 material 4\n");

    fprintf(outputStream, "mps 1 d 0. lattice a1 0.297 a2 1.e-12 talpha 0. referencetemperature 296. mode 1 q1 1.26403980927102E-11 q2 8.99822274579516E-11 q3 1.63092787267537E-12 q4 3.51199279469883E-12  stiffnessfactor 1. timefactor 1. lambda0 1. begoftimeofinterest 1.e-8 endoftimeofinterest 100000. relMatAge 28. CoupledAnalysisType 0\n");
    fprintf(outputStream, "latticeplasticitydamageviscoelastic 2 d 0 talpha 0. viscomat 1 calpha 0. e 45.91e9 a1 0.297 a2 1.e-12 ft 2.35e6 fc 30.e6 wf 20.e-6 ahard 1.e-3 angle1 0.5 flow 0.25 timefactor 1. iter 100 tol 1.e-6 timedepfracturing fcm28 30 fib_s 0.25 stiffnessfactor 1. randvars 2 806 807 randgen 2 4 4\n");
    fprintf(outputStream, "latticelinearelastic 3 d 0 talpha 0. calpha 0. e 300e10\n");
    fprintf(outputStream, "latticelinearelastic 4 d 0 talpha 0. e 45.91e9 a1 0.001 calpha 0.15e-3\n");

    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStream, "StructTemperatureLoad 2 loadTimeFunction 2 Components 1 1.0\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "PiecewiseLinFunction 2 nPoints 2 t 2 0. 47600. f(t) 2 0. 1.\n");
    fprintf(outputStream, "PiecewiseLinFunction 3 nPoints 2 t 2 0. 1499. f(t) 2 31.75 31.75\n");
    fprintf(outputStream, "InterpolatingFunction 4 name random.dat dim 3\n");

    //Print set1
    fprintf(outputStream, "Set 1 elements %d", set1Counter);
    for (int i = 0; i < set1Counter; i++) {
        fprintf(outputStream, " %d", set1.at(i + 1) );
    }
    fprintf(outputStream, "\n");

    //Print set2
    fprintf(outputStream, "Set 2 elements %d", set2Counter);
    for (int i = 0; i < set2Counter; i++) {
        fprintf(outputStream, " %d", set2.at(i + 1) );
    }
    fprintf(outputStream, "\n");

    //Print set2
    fprintf(outputStream, "Set 3 elements %d", set3Counter);
    for (int i = 0; i < set3Counter; i++) {
        fprintf(outputStream, " %d", set3.at(i + 1) );
    }
    fprintf(outputStream, "\n");



    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "#LOADLEVEL\n");
    fprintf(outputStream, "##TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");
    return;
}


void
Grid::give3DTensionOutput(const std::string &fileName)
{
    //Output for direct tension.
    //Start with mechanical model

    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords(3);
    oofem::FloatArray coordsOne, coordsTwo, line;
    int materialType = 1;
    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;
    double radius, distanceOne, distanceTwo;


    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ||  this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 && ( this->giveDelaunayLine(i + 1) )->delaunayAreaCheck() == 1 ) {
            numberOfLines++;
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 && this->giveRegion(1)->modifyVoronoiCrossSection(i + 1) == 1 )      {
            numberOfLines++;
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Mechanical 3D model of cylinder\n");
    fprintf(outputStream, "NonLinearStatic nmsteps 1 nsteps 1 contextOutputStep 2000 nmodules 2 lstype 3 smtype 7\n");
    fprintf(outputStream, "nsteps 200 rtolv 1.e-3 stiffMode 1 manrmsteps 10 maxiter 200 controllmode 1 lstype 3 smtype 7\n");
    fprintf(outputStream, "vtkxmllattice primvars 1 1 tstep_all domain_all cross 1 cellvars 2 60 90\n");
    fprintf(outputStream, "gpexportmodule vars 2 46 139 tstep_all\n");
    fprintf(outputStream, "domain 3dLattice\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_output { 3 4 }\n");
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 2 nmat 2 nbc 2 nic 0 nltf 2 nset 5\n", numberOfNodes, numberOfLines);

    printf("start nodes\n");

    int supportNode;
    int loadNode;
    int firstFlag = 0;
    int bottomSetCounter = 0;
    oofem::IntArray bottomSetTemp(this->giveNumberOfDelaunayVertices() );
    int topSetCounter = 0;
    oofem::IntArray topSetTemp(this->giveNumberOfDelaunayVertices() );

    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);


            //Create support and load nodes
            if ( fabs(coords.at(1) - 0.025) < giveTol() && fabs(coords.at(2) - 0.0) < giveTol() && fabs(coords.at(3) - 0.005) < giveTol() ) {
                supportNode = i + 1;
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else if ( fabs(coords.at(1) - 0.025) < giveTol() && fabs(coords.at(2) - 0.075) < giveTol() && fabs(coords.at(3) - 0.005) < giveTol() )                  {
                loadNode = i + 1;
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else if ( fabs(coords.at(2) - 0.0) < giveTol() )          {
                bottomSetCounter++;
                bottomSetTemp.at(bottomSetCounter) = i + 1;
                fprintf(outputStream, "rigidarmnode %d coords 3 %e %e %e master %d mastermask 6 0 1 0 0 0 0 doftype 6 0 2 0 0 0 0\n", i + 1, coords.at(1), coords.at(2), coords.at(3), supportNode);
            } else if ( fabs(coords.at(2) - 0.075) < giveTol() )          {
                topSetCounter++;
                topSetTemp.at(topSetCounter) = i + 1;
                fprintf(outputStream, "rigidarmnode %d coords 3 %e %e %e master %d mastermask 6 0 1 0 0 0 0 doftype 6 0 2 0 0 0 0\n", i + 1, coords.at(1), coords.at(2), coords.at(3), loadNode);
            } else   {
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    oofem::IntArray bottomSet(bottomSetCounter);
    bottomSet.zero();
    for (int i = 0; i < bottomSetCounter; i++) {
        bottomSet.at(i + 1) = bottomSetTemp.at(i + 1);
    }

    oofem::IntArray topSet(topSetCounter);
    topSet.zero();
    for (int i = 0; i < topSetCounter; i++) {
        topSet.at(i + 1) = topSetTemp.at(i + 1);
    }

    printf("finished nodes\n");

    int set1Counter = 0, set2Counter = 0, set3Counter = 0;
    oofem::IntArray set1Temp(numberOfLines);
    set1Temp.zero();

    oofem::IntArray set2Temp(numberOfLines);
    set2Temp.zero();

    oofem::IntArray set3Temp(numberOfLines);
    set3Temp.zero();

    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) && ( this->giveDelaunayLine(i + 1) )->delaunayAreaCheck() == 1 ) { //Elements are inside
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);

            //=======================================================
            //Deal with inclusions
            materialType = 1;
            //	this->giveDelaunayLine(i + 1)->updateMaterial(1);
            this->giveDelaunayVertex(nodes.at(1) )->giveCoordinates(coordsOne);
            this->giveDelaunayVertex(nodes.at(2) )->giveCoordinates(coordsTwo);


            if ( coordsOne.at(2) < 0.0125 && coordsTwo.at(2) < 0.0125 ) {
                set2Counter++;
                set2Temp.at(set2Counter) = i + 1;
                materialType = 2;
                this->giveDelaunayLine(i + 1)->updateMaterial(2);
            } else if ( (coordsOne.at(2) >= 0.0125 && coordsOne.at(2) <= 0.0625) || (coordsTwo.at(2) >= 0.0125 && coordsTwo.at(2) <= 0.0625) )    {
                set1Counter++;
                set1Temp.at(set1Counter) = i + 1;
                materialType = 1;
                this->giveDelaunayLine(i + 1)->updateMaterial(1);
            } else if ( coordsOne.at(2) > 0.0625 && coordsTwo.at(2) > 0.0625 ) {
                set3Counter++;
                set3Temp.at(set3Counter) = i + 1;
                materialType = 2;
                this->giveDelaunayLine(i + 1)->updateMaterial(2);
            }

            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "lattice3D %d nodes 2 %d %d crossSect %d mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, "\n");
        }
    }

    oofem::IntArray set1(set1Counter);
    set1.zero();

    oofem::IntArray set2(set2Counter);
    set2.zero();

    oofem::IntArray set3(set3Counter);
    set3.zero();

    for (int i = 0; i < set1Counter; i++) {
        set1.at(i + 1) = set1Temp.at(i + 1);
    }

    for (int i = 0; i < set2Counter; i++) {
        set2.at(i + 1) = set2Temp.at(i + 1);
    }

    for (int i = 0; i < set3Counter; i++) {
        set3.at(i + 1) = set3Temp.at(i + 1);
    }

    fprintf(outputStream, "latticecs 1 material 1\n");
    fprintf(outputStream, "latticecs 2 material 2\n");

    fprintf(outputStream, "latticeplastdam 1 d 0 talpha 0. calpha 0. e 50.46e9 a1 0.215 ft 2.44e6 fc 30.e6 wf 50.e-6 ahard 1.e-3 angle1 0.5 flow 0.25\n");
    fprintf(outputStream, "latticelinearelastic 2 d 0 talpha 0. calpha 0. e 50.46e9 a1 0.215\n");

    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 dofs 6 1 2 3 4 5 6 values 6 0 0 0 0 0 0 set 4\n");
    fprintf(outputStream, "BoundaryCondition 2 loadTimeFunction 2 dofs 6 1 2 3 4 5 6 values 6 0 4.e-4 0 0 0 0 set 5\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "PiecewiseLinFunction 2 nPoints 2 t 2 0. 199. f(t) 2 0. 1.\n");

    //Print set1
    fprintf(outputStream, "Set 1 elements %d", set1Counter);
    for (int i = 0; i < set1Counter; i++) {
        fprintf(outputStream, " %d", set1.at(i + 1) );
    }
    fprintf(outputStream, "\n");

    //Print set2
    fprintf(outputStream, "Set 2 elements %d", set2Counter);
    for (int i = 0; i < set2Counter; i++) {
        fprintf(outputStream, " %d", set2.at(i + 1) );
    }
    fprintf(outputStream, "\n");

    //Print set2
    fprintf(outputStream, "Set 3 elements %d", set3Counter);
    for (int i = 0; i < set3Counter; i++) {
        fprintf(outputStream, " %d", set3.at(i + 1) );
    }
    fprintf(outputStream, "\n");

    //Print set4
    fprintf(outputStream, "Set 4 nodes 1 %d\n", supportNode);
    fprintf(outputStream, "Set 5 nodes 1 %d\n", loadNode);
    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#NODE number 4 dof 2 unknown d\n");
    fprintf(outputStream, "#REACTION number 3 dof 2 unknown d\n");
    fprintf(outputStream, "#TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");
    return;
}

void
Grid::give3DPeriodicTetraSMOutput(const std::string &fileName)
{
    FILE *outputStream = converter::fopen_or_die(fileName, "w");


    oofem::IntArray periodicFlag(3);
    this->givePeriodicityFlag(periodicFlag);

    oofem::FloatArray boundaries(6);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfNodes, numberOfTetras;
    oofem::FloatArray coords, coordTest;
    int materialType = 1;
    oofem::IntArray nodes, location(4);
    oofem::IntArray newNodes, imageNodeList, mirrorNodeList;

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        //Determine which nodes are mirror images
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coordTest);

            if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0 ) {
                //Periodicity in x
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                //Periodicity in y
                if ( fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in z
                if ( fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                //Periodicity in x y
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol()  ||
                     fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in x z
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol()  ||
                     fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in y z
                if ( fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol()  ||
                     fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in x y z
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol()  ||
                     fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol()  ||
                     fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            }
        }
    }

    numberOfNodes = mirrorNodeList.giveSize() + 1; //add control node

    //Determine the number of Delaunay tetrahedra in the domain
    numberOfTetras = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayTetras(); i++ ) {
        if (  ( this->giveDelaunayTetra(i + 1) )->giveOutsideFlag() == 0 ) {
            numberOfTetras++;
        }
    }

    //OOFEM INPUT BEGINS
    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "3D RVE - periodic boundary conditions with LTRSpaceBoundary elements\n");
    fprintf(outputStream, "StaticStructural nsteps 1 nmodules 0 lstype 3 smtype 7\n");
    fprintf(outputStream, "#vtkxmlperiodic primvars 1 1 vars 3 4 1 12 cellvars 1 46 tstep_all domain_all stype 0\n");
    fprintf(outputStream, "#vtkxmlperiodic ipvars 3 4 1 12 tstep_all domain_all\n");
    fprintf(outputStream, "#gpexportmodule vars 3 4 1 12 tstep_all domain_all\n");
    fprintf(outputStream, "#matlab tstep_all mesh data specials reactionforces integrationpoints internalvars 3 4 1 12\n");
    fprintf(outputStream, "domain 3d\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_all element_all\n");
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 1 nset 1\n", numberOfNodes, numberOfTetras);

    int firstFlag = 0;

    //node output
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);

            if ( mirrorNodeList.isEmpty() == 0 ) {
                if ( mirrorNodeList.contains(i + 1) ) {
                    if ( firstFlag == 0 ) {
                        firstFlag = 1;
                        fprintf(outputStream, "node %d coords 3 %e %e %e bc 3 1 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
                    } else {
                        fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
                    }
                }
            }
        }
    }

    //control node
    if ( this->macroType == _Truss ) {
        fprintf(outputStream, "node %d coords 3 %e %e %e dofidmask 1 1 bc 1 2\n", this->giveNumberOfDelaunayVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );
    } else if ( this->macroType == _Membrane ) {
        fprintf(outputStream, "node %d coords 3 %e %e %e dofidmask 4 1 2 4 5 bc 4 2 2 2 2\n", this->giveNumberOfDelaunayVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );
    } else if ( this->macroType == _Beam ) {
        fprintf(outputStream, "node %d coords 3 %e %e %e dofidmask 3 1 7 10 bc 3 2 2 2\n", this->giveNumberOfDelaunayVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );
    } else if ( this->macroType == _Plate ) {
        fprintf(outputStream, "node %d coords 3 %e %e %e dofidmask 10 1 2 4 5 7 8 10 11 12 13 bc 10 2 2 2 2 2 2 2 2 2 2\n", this->giveNumberOfDelaunayVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );
    } else if ( this->macroType == _3dVoigt ) {
        fprintf(outputStream, "node %d coords 3 %e %e %e dofidmask 6 1 2 3 4 5 6 bc 6 2 2 2 2 2 2\n", this->giveNumberOfDelaunayVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );
    } else if ( this->macroType == _3d ) {
        fprintf(outputStream, "node %d coords 3 %e %e %e dofidmask 9 1 2 3 4 5 6 7 8 9 bc 9 2 2 2 2 2 2 2 2 2\n", this->giveNumberOfDelaunayVertices() + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );
    }

    //element output
    int boundaryFlag = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayTetras(); i++ ) {
        //First plot all of them
        if ( this->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            location.zero();
            boundaryFlag = 0;
            materialType = 1;
            this->giveDelaunayTetra(i + 1)->giveLocalVertices(nodes);
            newNodes = nodes;

            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < nodes.giveSize(); m++ ) {
                giveDelaunayVertex(nodes.at(m + 1) )->giveCoordinates(coords);
                if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0 ) {
                    //Periodicity in x
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol() ) { //x=xmax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                    //Periodicity in y
                    if ( fabs(coords.at(2) - boundaries.at(4) ) < giveTol() ) { //y=ymax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in z
                    if ( fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //z=zmax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                    //Periodicity in x y
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol()  ||
                         fabs(coords.at(2) - boundaries.at(4) ) < giveTol() ) { //x=xmax or y=ymax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in x z
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol()  ||
                         fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //x=xmax or z=zmax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in y z
                    if ( fabs(coords.at(2) - boundaries.at(4) ) < giveTol()  ||
                         fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //y=ymax or z=zmax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in x y z
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol()  ||
                         fabs(coords.at(2) - boundaries.at(4) ) < giveTol()  ||
                         fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //x=xmax or y=ymax or z=zmax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                }
            }

            if ( boundaryFlag == 0 ) {
                fprintf(outputStream, "ltrspace %d nodes 4 %d %d %d %d crossSect 1 mat %d", i + 1, nodes.at(1), nodes.at(2), nodes.at(3), nodes.at(4), materialType);
                fprintf(outputStream, "\n");
            } else {
                fprintf(outputStream, "%s %d nodes 5 %d %d %d %d %d crossSect 1 mat %d", boundElemName, i + 1, newNodes.at(1), newNodes.at(2), newNodes.at(3), newNodes.at(4), this->giveNumberOfDelaunayVertices() + 1, materialType);
                fprintf(outputStream, " location 4 %d %d %d %d", location.at(1), location.at(2), location.at(3), location.at(4) );
                fprintf(outputStream, "\n");
            }
        }
    }

    fprintf(outputStream, "simplecs 1\n");
    fprintf(outputStream, "isole 1 talpha 0. d 0. e 30.e9 n 0.\n");
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    if ( this->macroType == _Truss ) {
        fprintf(outputStream, "BoundaryCondition 2 loadTimeFunction 1 dofs 1 1 values 1 0 set 1\n");
    } else if ( this->macroType == _Membrane ) {
        fprintf(outputStream, "BoundaryCondition 2 loadTimeFunction 1 dofs 4 1 2 4 5  values 4 0 0 0 0 set 1\n");
    } else if ( this->macroType == _Beam ) {
        fprintf(outputStream, "BoundaryCondition 2 loadTimeFunction 1 dofs 3 1 7 10 values 3 0 0 0 set 1\n");
    } else if ( this->macroType == _Plate ) {
        fprintf(outputStream, "BoundaryCondition 2 loadTimeFunction 1 dofs 10 1 2 4 5 7 8 10 11 12 13 values 10 0 0 0 0 0 0 0 0 0 0 set 1\n");
    } else if ( this->macroType == _3dVoigt ) {
        fprintf(outputStream, "BoundaryCondition 2 loadTimeFunction 1 dofs 6 1 2 3 4 5 6 values 6 0 0 0 0 0 0 set 1\n");
    } else if ( this->macroType == _3d ) {
        fprintf(outputStream, "BoundaryCondition 2 loadTimeFunction 1 dofs 9 1 2 3 4 5 6 7 8 9 values 9 0 0 0 0 0 0 0 0 0 set 1\n");
    }
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "set 1 nodes 1 %d\n", this->giveNumberOfDelaunayVertices() + 1);
    return;
}



void
Grid::give3DTetraSMOutput(const std::string &fileName)
{
    //Output for 3D fracture process zone modelling

    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    oofem::FloatArray boundaries(3);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfNodes, numberOfTetras;
    oofem::FloatArray coords;
    int materialType = 1;
    oofem::IntArray nodes, location(2);

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ||  this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfNodes++;
        }
    }

    //Determine the number of Delaunay tetrahedra in the domain
    numberOfTetras = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayTetras(); i++ ) {
        if (  ( this->giveDelaunayTetra(i + 1) )->giveOutsideFlag() == 0 ) {
            numberOfTetras++;
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Mechanical 3D model\n");
    fprintf(outputStream, "NonLinearStatic nmsteps 1 nsteps 1 contextOutputStep 1000 nmodules 2 profileopt 1 lstype 3 smtype 7\n");
    fprintf(outputStream, "nsteps 200 rtolv 1.e-3 reqIterations 100 stiffMode 1 manrmsteps 10 maxiter 200 controllmode 0 stepLength 5.e-7 minsteplength 5.e-7 maxrestarts 0 hpcmode 2 hpc 2 X 2 hpcw 1 1. lstype 3 smtype 7\n");
    fprintf(outputStream, "vtkxml primvars 1 1 tstep_all domain_all\n");
    fprintf(outputStream, "gpexportmodule vars 5 59 90 84 85 78 tstep_all domain_all\n");
    fprintf(outputStream, "domain 3d\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_output {%d}\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 1\n", numberOfNodes, numberOfTetras);

    int firstFlag = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ||  this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            if ( firstFlag == 0 ) {
                firstFlag = 1;
                fprintf(outputStream, "node %d coords 3 %e %e %e bc 3 1 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    for ( int i = 0; i < this->giveNumberOfDelaunayTetras(); i++ ) {
        //First plot all of them
        if ( this->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveDelaunayTetra(i + 1)->giveLocalVertices(nodes);
            materialType = 1;
            fprintf(outputStream, "ltrspace %d nodes 4 %d %d %d %d crossSect 1 mat %d", i + 1, nodes.at(1), nodes.at(2), nodes.at(3), nodes.at(4), materialType);
            fprintf(outputStream, "\n");
        }
    }
    fprintf(outputStream, "simplecs 1\n");
    fprintf(outputStream, "isole 1 talpha 0. d 0. e 30.e9 n 0.\n");
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0.0\n");
    fprintf(outputStream, "NodalLoad 2 loadTimeFunction 1 Components 6 0. 1. 0. 0. 0. 0.\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", this->giveNumberOfDelaunayVertices() + 1);
    fprintf(outputStream, "#LOADLEVEL\n");
    fprintf(outputStream, "##TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");
    return;
}

void
Grid::give3DRCSMOutput(const std::string &fileName)
/**
 * Output for random mesh of a 3D reinforced concrete RVE.  Concrete is modelled with tetrahedras (LTRSpace),
 * reinforcement is modelled with beam elements (LIBeam3d), and the interface is modelled with pointwise
 * interface elements (IntElPoint).
 *
 * @author: Adam Sciegaj
 */
{
    FILE *outputStream = converter::fopen_or_die(fileName, "w");



    printf("Processing nodes\n");

    oofem::FloatArray boundaries(6);
    this->giveRegion(1)->defineBoundaries(boundaries);

    int numberOfNodes, numberOfTetras, numberOfBeams, numberOfLinks;
    oofem::FloatArray coords;
    oofem::IntArray nodes;

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    for ( int i = 0; i < this->giveNumberOfReinforcementNode(); i++ ) {
        if ( this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 0 || this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfNodes++;
        }
    }

    //add reference nodes (1 per fibre)
    numberOfNodes +=  converter::size1(fibreList);

    //Determine the number of Delaunay tetrahedra in the domain
    numberOfTetras = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayTetras(); i++ ) {
        if (  this->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfTetras++;
        }
    }

    numberOfBeams = 0;
    for ( int i = 0; i < this->giveNumberOfLatticeBeams(); i++ ) {
        if ( this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 0 || this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfBeams++;
        }
    }

    numberOfLinks = 0;
    for ( int i = 0; i < this->giveNumberOfLatticeLinks(); i++ ) {
        if ( this->giveLatticeLink(i + 1)->giveOutsideFlag() == 0 || this->giveLatticeLink(i + 1)->giveOutsideFlag() == 3 ) {
            numberOfLinks++;
        }
    }

    int numberOfElem = numberOfTetras + numberOfBeams + numberOfLinks;

    //OOFEM INPUT BEGINS
    //potentially different cross section, material and bond properties for each fiber
    fprintf(outputStream, "rve.out\n");
    fprintf(outputStream, "Reinforced concrete RVE in 3D\n");
    fprintf(outputStream, "StaticStructural nsteps 1 initialguess 1 nmodules 4 lstype 3 smtype 7\n");
    fprintf(outputStream, "vtkxml tstep_all primvars 1 1 vars 2 4 1 stype 0 regionsets 1 1\n");
    fprintf(outputStream, "vtkxml tstep_all primvars 1 1 vars 2 7 8 stype 0 regionsets 1 2\n");
    fprintf(outputStream, "vtkxml tstep_all ipvars 2 98 99 regionsets 1 3\n");
    fprintf(outputStream, "matlab tstep_all mesh data reactionforces integrationpoints internalvars 7 4 1 12 7 8 98 99\n");
    fprintf(outputStream, "domain 3d\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_all element_all\n");
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect %d nmat %d nbc 2 nic 0 nltf 1 nset %d\n", numberOfNodes, numberOfElem,
            2 * converter::size1(fibreList) + 1,  2 * converter::size1(fibreList) + 1, 3 + 2 * converter::size1(fibreList) + 1);

    int firstFlag = 0;

    //node output - tetrahedra
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);

            if ( firstFlag == 0 ) {
                firstFlag = 1;
                fprintf(outputStream, "node %d coords 3 %e %e %e bc 3 1 1 1\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            } else {
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    //node output - beams
    oofem::IntArray L2Gmap;
    L2Gmap.resize(this->giveNumberOfReinforcementNode() );
    L2Gmap.zero();
    for ( int i = 0; i < this->giveNumberOfReinforcementNode(); i++ ) {
        if ( this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 0 || this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveReinforcementNode(i + 1)->giveCoordinates(coords);

            L2Gmap.at(i + 1) = this->giveNumberOfDelaunayVertices() + i + 1;
            fprintf(outputStream, "node %d coords 3 %e %e %e \n", L2Gmap.at(i + 1), coords.at(1), coords.at(2), coords.at(3) );
        }
    }

    //reference nodes for 3d beams. Assuming only horizontal reinforcement
    oofem::FloatMatrix R(3, 3);
    R.zero();
    R.at(1, 2) = -1;
    R.at(2, 1) = 1;
    R.at(3, 3) = 1;
    oofem::IntArray refNodeNumbers(this->giveNumberOfFibres() );
    oofem::FloatMatrix referenceNodes;
    referenceNodes.resize(this->giveNumberOfFibres(), 3);
    for ( int i = 0; i < this->giveNumberOfFibres(); i++) {
        oofem::FloatArray dirVec = this->giveFibre(i + 1)->giveDirVector();
        dirVec.rotatedWith(R, 'n');
        int numReinfNode = this->giveFibre(i + 1)->giveNumberReinforcementNode(1);
        oofem::FloatArray refNodeCoords = * this->giveReinforcementNode(numReinfNode)->giveCoordinates();
        refNodeCoords.add(dirVec);
        refNodeNumbers.at(i + 1) = i + 1 + this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode();
        referenceNodes.at(i + 1, 1) = refNodeCoords.at(1);
        referenceNodes.at(i + 1, 2) = refNodeCoords.at(2);
        referenceNodes.at(i + 1, 3) = refNodeCoords.at(3);

        fprintf(outputStream, "node %d coords 3 %e %e %e \n", refNodeNumbers.at(i + 1), refNodeCoords.at(1), refNodeCoords.at(2), refNodeCoords.at(3) );
    }

    printf("Finished writing node data\n");

    //element output - tetrahedra
    oofem::IntArray set1;
    for ( int i = 0; i < this->giveNumberOfDelaunayTetras(); i++ ) {
        //First plot all of them
        if ( this->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            this->giveDelaunayTetra(i + 1)->giveLocalVertices(nodes);

            fprintf(outputStream, "ltrspace %d nodes 4 %d %d %d %d ", i + 1, nodes.at(1), nodes.at(2), nodes.at(3), nodes.at(4) );
            fprintf(outputStream, "\n");
            set1.followedBy(i + 1);
        }
    }
    printf("Finished writing Delaunay element (tetrahedra) data\n");

    //element output - beams
    oofem::IntArray set2, beamElL2Gmap;
    beamElL2Gmap.resize(this->giveNumberOfLatticeBeams() );
    beamElL2Gmap.zero();
    for ( int i = 0; i < this->giveNumberOfLatticeBeams(); i++ ) {
        //First plot all of them
        if ( this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 0 || this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 2 ) { //Elements are inside
            this->giveLatticeBeam(i + 1)->giveLocalVertices(nodes);

            //find the reference node
            int refNode = 0;
            oofem::FloatArray dirVec = this->giveLatticeBeam(i + 1)->giveDirectionVector();
            for ( int j = 0; j < this->giveNumberOfFibres(); j++ ) {
                oofem::FloatArray perpVec(3), firstNodeCoords(3);
                int numReinfNode = this->giveFibre(j + 1)->giveNumberReinforcementNode(1);
                firstNodeCoords = * this->giveReinforcementNode(numReinfNode)->giveCoordinates();
                perpVec.at(1) = referenceNodes.at(j + 1, 1) - firstNodeCoords.at(1);
                perpVec.at(2) = referenceNodes.at(j + 1, 2) - firstNodeCoords.at(2);
                perpVec.at(3) = referenceNodes.at(j + 1, 3) - firstNodeCoords.at(3);
                double sp = dirVec.dotProduct(perpVec, 3);

                if ( fabs(sp) < this->giveTol() ) {
                    refNode = refNodeNumbers.at(j + 1);
                    break;
                }
            }

            fprintf(outputStream, "libeam3d %d nodes 2 %d %d refnode %d", this->giveNumberOfDelaunayTetras() + i + 1, L2Gmap.at(nodes.at(1) ), L2Gmap.at(nodes.at(2) ), refNode);
            fprintf(outputStream, "\n");
            set2.followedBy(this->giveNumberOfDelaunayTetras() + i + 1);
            beamElL2Gmap.at(i + 1) = this->giveNumberOfDelaunayTetras() + i + 1;
        }
    }
    printf("Finished writing reinforcement data\n");

    //element output - interface
    oofem::IntArray set3, intElL2Gmap;
    intElL2Gmap.resize(this->giveNumberOfLatticeLinks() );
    intElL2Gmap.zero();
    for ( int i = 0; i < this->giveNumberOfLatticeLinks(); i++ ) {
        int linkVer1 = this->giveLatticeLink(i + 1)->giveLocalVertex(1);
        int linkVer2 = this->giveLatticeLink(i + 1)->giveLocalVertex(2);
        //recalculate normal vector
        oofem::FloatArray normalVec(3);
        oofem::FloatArray beamNode = * this->giveReinforcementNode(linkVer1)->giveCoordinates();
        oofem::FloatArray tetraNode = * this->giveDelaunayVertex(linkVer2)->giveCoordinates();
        normalVec.add(beamNode);
        normalVec.subtract(tetraNode);
        normalVec.normalize();
        double segmentLength(1);
        if ( this->giveLatticeLink(i + 1)->giveOutsideFlag() == 0 || this->giveLatticeLink(i + 1)->giveOutsideFlag() == 3 ) {
            if ( this->giveLatticeLink(i + 1)->giveOutsideFlag() == 3 ) {
                //shorten the associated length of the boundary links
                segmentLength = this->giveLatticeLink(i + 1)->giveAssociatedLength() / 2;
            } else {
                segmentLength = this->giveLatticeLink(i + 1)->giveAssociatedLength();
            }

            //             fprintf( outputStream, "intelpoint %d nodes 2 %d %d normal 3 %e %e %e length %e\n", this->giveNumberOfDelaunayTetras() + this->giveNumberOfLatticeBeams() + i + 1, linkVer2, L2Gmap.at(linkVer1), normalVec.at(1), normalVec.at(2), normalVec.at(3), segmentLength );
            fprintf(outputStream, "intelpoint %d nodes 2 %d %d normal 3 0 0 1 length %e\n", this->giveNumberOfDelaunayTetras() + this->giveNumberOfLatticeBeams() + i + 1, linkVer2, L2Gmap.at(linkVer1), this->giveLatticeLink(i + 1)->giveAssociatedLength() );
            set3.followedBy(this->giveNumberOfDelaunayTetras() + this->giveNumberOfLatticeBeams() + i + 1);
            intElL2Gmap.at(i + 1) = this->giveNumberOfDelaunayTetras() + this->giveNumberOfLatticeBeams() + i + 1;
        }
    }
    printf("Finished writing interface data\n");

    //split into sets (along fibres)
    std::vector < std::vector < int >> beamSets(this->giveNumberOfFibres() );
    std::vector < std::vector < int >> interfaceSets(this->giveNumberOfFibres() );
    oofem::IntArray beamsPerFibre(this->giveNumberOfFibres() ), firstNodeinFibres(this->giveNumberOfFibres() );
    for ( int i = 0; i < this->giveNumberOfFibres(); i++ ) {
        beamsPerFibre.at(i + 1) = this->giveFibre(i + 1)->NbOfReinfNodes() - 1;
        int beamNo(0);
        if ( i != 0 ) {
            for ( int r = 1; r <= i; r++ ) {
                beamNo += beamsPerFibre.at(r);
            }
        }
        //first split reinforcement
        for ( int j = 0; j < this->giveFibre(i + 1)->NbOfReinfNodes() - 1; j++ ) {
            if ( this->giveLatticeBeam(beamNo + j + 1)->giveOutsideFlag() == 0 || this->giveLatticeBeam(beamNo + j + 1)->giveOutsideFlag() == 2 ) {
                ( beamSets [ i ] ).push_back(beamNo + j + 1);
                oofem::IntArray beamNodes;
                this->giveLatticeBeam(beamNo + j + 1)->giveLocalVertices(beamNodes);
                //for the beam elements, get the corresponding link elements
                for ( int k = 0; k < this->giveNumberOfLatticeLinks(); k++ ) {
                    oofem::IntArray linkNodes;
                    this->giveLatticeLink(k + 1)->giveLocalVertices(linkNodes);
                    if ( beamNodes.at(1) == linkNodes.at(1)  ) {
                        ( interfaceSets [ i ] ).push_back(k + 1);
                    } else if ( this->giveLatticeLink(k + 1)->giveOutsideFlag() == 3 && beamNodes.at(2) == linkNodes.at(1) ) {
                        //include the last boundary link
                        ( interfaceSets [ i ] ).push_back(k + 1);
                    }
                }
            }
        }
        firstNodeinFibres.at(i + 1) = L2Gmap.at(this->giveLatticeLink(interfaceSets [ i ] [ 0 ])->giveLocalVertex(1) );
    }
    printf("Finished splitting into sets\n");

    fprintf(outputStream, "SimpleCS 1 material 1 set 1\n");
    for ( int j = 0; j < this->giveNumberOfFibres(); j++ ) {
        double diam = this->giveFibre(j + 1)->giveDiameter();
        double sideLength = sqrt(M_PI * diam * diam * 0.25);

        fprintf(outputStream, "FiberedCS %d fibermaterials 1 %d thicks 1 %e widths 1 %e thick %e width %e fiberycentrecoords 1 %e fiberzcentrecoords 1 %e set %d\n", 2 * ( j + 1 ), 2 * ( j + 1 ), sideLength, sideLength, sideLength, sideLength, sideLength / 2, sideLength / 2, 2 * ( j + 2 ) );
        fprintf(outputStream, "InterfaceCS %d thickness %e material %d set %d\n", 2 * ( j + 1 ) + 1, M_PI * diam, 2 * ( j + 1 ) + 1, 2 * ( j + 2 ) + 1);
    }
    fprintf(outputStream, "IsoLE 1 talpha 0 d 1.0 e 30.e9 n 0.2\n");
    for ( int j = 0; j < this->giveNumberOfFibres(); j++ ) {
        double s3 = 5. / 24 * this->giveFibre(j + 1)->giveDiameter() + 7. / 3000; //approximate interpolation
        fprintf(outputStream, "MisesMat %d d 1.0 E 2e11 n 0.3 sig0 500e6 H 8.456659619e8 omega_crit 0 a 0 tAlpha 1.0\n", 2 * ( j + 1 ) );
        fprintf(outputStream, "BondCEB %d kn 6e12 ks 6.135e10 s1 0.001 s2 0.002 s3 %e taumax 1.541e7 tauf 6.164e6\n", 2 * ( j + 1 ) + 1, s3);
    }
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 dofs 3 1 2 3 values 3 0 0 0\n");
    fprintf(outputStream, "BoundaryCondition 2 loadTimeFunction 1 dofs 3 4 5 6 values 2 0 0 set 2\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "set 1 elements %d ", set1.giveSize() );
    for ( int i = 0; i < set1.giveSize(); i++ ) {
        fprintf(outputStream, "%d ", set1.at(i + 1) );
    }
    fprintf(outputStream, "\n");
    fprintf(outputStream, "set 2 elements %d ", set2.giveSize() );
    for ( int i = 0; i < set2.giveSize(); i++ ) {
        fprintf(outputStream, "%d ", set2.at(i + 1) );
    }
    fprintf(outputStream, "\n");
    fprintf(outputStream, "set 3 elements %d ", set3.giveSize() );
    for ( int i = 0; i < set3.giveSize(); i++ ) {
        fprintf(outputStream, "%d ", set3.at(i + 1) );
    }
    fprintf(outputStream, "\n");
    for ( int i = 0; i < this->giveNumberOfFibres(); i++ ) {
        fprintf(outputStream, "set %d elements %zu ", 3 + 2 * ( i + 1 ) - 1, beamSets [ i ].size() );
        for ( int j = 0; j < beamSets [ i ].size(); j++ ) {
            fprintf(outputStream, "%d ", beamElL2Gmap.at(beamSets [ i ] [ j ]) );
        }
        fprintf(outputStream, "\n");
        fprintf(outputStream, "set %d elements %zu ", 3 + 2 * ( i + 1 ), interfaceSets [ i ].size() );
        for ( int j = 0; j < interfaceSets [ i ].size(); j++ ) {
            fprintf(outputStream, "%d ", intElL2Gmap.at(interfaceSets [ i ] [ j ]) );
        }
        fprintf(outputStream, "\n");
    }
    fprintf(outputStream, "set %d nodes %d ", 3 + 2 * converter::size1(fibreList) + 1, converter::size1(fibreList) );
    for ( int i = 0; i < this->giveNumberOfFibres(); i++ ) {
        fprintf(outputStream, "%d ", firstNodeinFibres.at(i + 1) );
    }
    fprintf(outputStream, "\n");


    return;
}






void
Grid::give3DRCPeriodicSMOutput(const std::string &fileName)
/**
 * Output for periodic mesh of a 3D reinforced concrete RVE.  Concrete is modelled with tetrahedras (LTRSpace and LTRSpaceBoundary),
 * reinforcement is modelled with beam elements (LIBeam3d and LIBeam3dBoundary), and the interface is modelled with link elements (BondLink3d).
 *
 * @authors: Adam Sciegaj, Peter Grassl
 */
{
    FILE *outputStream = converter::fopen_or_die(fileName, "w");



    printf("Processing nodes\n");

    oofem::IntArray periodicFlag(3);
    this->givePeriodicityFlag(periodicFlag);

    oofem::FloatArray boundaries(6);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfNodes, numberOfTetras, numberOfBeams, numberOfLinks;
    oofem::FloatArray coords, coordTest;
    int materialType = 1;
    oofem::IntArray nodes, location(4);
    oofem::IntArray newNodes, imageNodeList, mirrorNodeList;
    oofem::IntArray reinfImageNodeList, reinfMirrorNodeList;

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        //Determine which nodes are on mirror/image boundary
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coordTest);

            if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0 ) {
                //Periodicity in x
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                //Periodicity in y
                if ( fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in z
                if ( fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                //Periodicity in x y
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol()  ||
                     fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in x z
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol()  ||
                     fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in y z
                if ( fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol()  ||
                     fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in x y z
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol()  ||
                     fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol()  ||
                     fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            }
        }
    }

    for ( int i = 0; i < this->giveNumberOfReinforcementNode(); i++ ) {
        this->giveReinforcementNode(i + 1)->giveCoordinates(coordTest);
        if ( this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 0 || this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 2 ) {
            //Determine if on mirror or image boundary
            this->giveReinforcementNode(i + 1)->giveCoordinates(coordTest);

            if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0 ) {
                //Periodicity in x
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol() ) {
                    if ( reinfImageNodeList.contains(i + 1) == 0 ) {
                        reinfImageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( reinfMirrorNodeList.contains(i + 1) == 0 ) {
                        reinfMirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                //Periodicity in y
                if ( fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol() ) {
                    if ( reinfImageNodeList.contains(i + 1) == 0 ) {
                        reinfImageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( reinfMirrorNodeList.contains(i + 1) == 0 ) {
                        reinfMirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in z
                if ( fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( reinfImageNodeList.contains(i + 1) == 0 ) {
                        reinfImageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( reinfMirrorNodeList.contains(i + 1) == 0 ) {
                        reinfMirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                //Periodicity in x y
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol()  ||
                     fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol() ) {
                    if ( reinfImageNodeList.contains(i + 1) == 0 ) {
                        reinfImageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( reinfMirrorNodeList.contains(i + 1) == 0 ) {
                        reinfMirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in x z
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol()  ||
                     fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( reinfImageNodeList.contains(i + 1) == 0 ) {
                        reinfImageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( reinfMirrorNodeList.contains(i + 1) == 0 ) {
                        reinfMirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in y z
                if ( fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol()  ||
                     fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( reinfImageNodeList.contains(i + 1) == 0 ) {
                        reinfImageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( reinfMirrorNodeList.contains(i + 1) == 0 ) {
                        reinfMirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in x y z
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol()  ||
                     fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol()  ||
                     fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( reinfImageNodeList.contains(i + 1) == 0 ) {
                        reinfImageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( reinfMirrorNodeList.contains(i + 1) == 0 ) {
                        reinfMirrorNodeList.followedBy(i + 1);
                    }
                }
            }
        }
    }

    //total number of nodes: number of master tetra nodes + number of master reinforcement nodes +
    // reference nodes (assuming straight reinforcement -- 1 per fibre) + control node
    numberOfNodes = mirrorNodeList.giveSize() + reinfMirrorNodeList.giveSize() + converter::size1(fibreList) + 1;

    //Determine the number of Delaunay tetrahedra in the domain
    numberOfTetras = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayTetras(); i++ ) {
        if (  this->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfTetras++;
        }
    }

    numberOfBeams = 0;
    for ( int i = 0; i < this->giveNumberOfLatticeBeams(); i++ ) {
        if ( this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 0 || this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfBeams++;
        }
    }

    numberOfLinks = 0;
    for ( int i = 0; i < this->giveNumberOfLatticeLinks(); i++ ) {
        int linkVer1 = this->giveLatticeLink(i + 1)->giveLocalVertex(1);
        int linkVer2 = this->giveLatticeLink(i + 1)->giveLocalVertex(2);
        if ( reinfMirrorNodeList.contains(linkVer1) && mirrorNodeList.contains(linkVer2) ) {
            numberOfLinks++;
        }
    }

    int numberOfElem = numberOfTetras + numberOfBeams + numberOfLinks;

    int controlNode = this->giveNumberOfFibres() + this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode() + 1;

    //OOFEM INPUT BEGINS
    //potentially different cross section, material and bond properties for each fiber
    fprintf(outputStream, "rve.out\n");
    fprintf(outputStream, "Periodic reinforced concrete RVE in 3D\n");
    fprintf(outputStream, "StaticStructural nsteps 200 nmodules 3 initialguess 1 lstype 3 smtype 7 stiffmode 2 rtolv 1.e-3 maxiter 10000\n");
    //fprintf(outputStream, "NonLinearStatic nmsteps 1 nsteps 1 contextOutputStep 1000 nmodules 3 profileopt 1 lstype 3 smtype 7\nnsteps 200 rtolv 1.e-3 reqIterations 100 stiffMode 1 manrmsteps 10 maxiter 200 controllmode 0 stepLength 1.e-6 minsteplength 1.e-6 maxrestarts 0 hpcmode 2 hpc 2 %d 1 hpcw 1 1. lstype 3 smtype 7\n", controlNode);
    //fprintf(outputStream, "StaticStructural nsteps 10 solverType \"calm\" stepLength 5.e-4 Psi 0 hpcmode 2 hpc 2 %d 1 hpcw 1 1 nmodules 3 initialguess 1 lstype 0 smtype 1 stiffmode 2\n", controlNode);
    fprintf(outputStream, "vtkxmlperiodic tstep_all primvars 1 1 vars 2 4 1 stype 0 regionsets 1 1\n");
    fprintf(outputStream, "vtkxmlperiodic tstep_all primvars 1 1 vars 2 7 8 stype 0 regionsets 1 2\n");
    fprintf(outputStream, "vtkxmlperiodic tstep_all primvars 1 1 ipvars 2 98 99 regionsets 1 3\n");
    //    fprintf(outputStream, "matlab tstep_all mesh data reactionforces integrationpoints internalvars 7 4 1 12 7 8 98 99\n");
    fprintf(outputStream, "domain 3d\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_output {%d}\n", controlNode);
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect %d nmat %d nbc 4 nic 0 nltf 1 nset %d\n", numberOfNodes, numberOfElem,
            2 * converter::size1(fibreList) + 1,  2 * converter::size1(fibreList) + 1, 7 + 2 * converter::size1(fibreList) );

    int firstFlag = 0;

    //Find three nodes which can be used for constraining the specimen.
    //Translations can be fixed by setting DOFs of first node to zero.
    int firstNode = 0, secondNode = 0, thirdNode = 0;
    oofem::FloatArray firstNodeCoords(3);

    //node output - tetrahedra
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);

            if ( mirrorNodeList.isEmpty() == 0 ) {
                if ( mirrorNodeList.contains(i + 1) ) {
                    fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
                }
            }
        }
    }

    //Define nodes for supports to block rigid body rotations
    //These nodes are control nodes defined in the mesh.in file
    firstNode = 9;
    secondNode = 10;
    thirdNode = 11;

    //node output - beams
    oofem::IntArray L2Gmap;
    L2Gmap.resize(this->giveNumberOfReinforcementNode() );
    L2Gmap.zero();
    for ( int i = 0; i < this->giveNumberOfReinforcementNode(); i++ ) {
        if ( this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 0 || this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveReinforcementNode(i + 1)->giveCoordinates(coords);

            if ( reinfMirrorNodeList.isEmpty() == 0 ) {
                if ( reinfMirrorNodeList.contains(i + 1) ) {
                    L2Gmap.at(i + 1) = this->giveNumberOfDelaunayVertices() + i + 1;
                    fprintf(outputStream, "node %d coords 3 %e %e %e\n", L2Gmap.at(i + 1), coords.at(1), coords.at(2), coords.at(3) );
                }
            }
        }
    }

    //reference nodes for 3d beams. Assuming only straight reinforcement bars.
    //Therefore, only one refnode per fibre. For curved reinforcement, each beam
    //element would require a seperate reference node.
    oofem::FloatMatrix R(3, 3);
    R.zero();
    R.at(1, 2) = -1;
    R.at(2, 1) = 1;
    R.at(3, 3) = 1;
    oofem::IntArray refNodeNumbers(this->giveNumberOfFibres() );
    oofem::FloatMatrix referenceNodes;
    referenceNodes.resize(this->giveNumberOfFibres(), 3);
    for ( int i = 0; i < this->giveNumberOfFibres(); i++) {
        oofem::FloatArray dirVec = this->giveFibre(i + 1)->giveDirVector();
        dirVec.rotatedWith(R, 'n');
        int numReinfNode = this->giveFibre(i + 1)->giveNumberReinforcementNode(1);
        oofem::FloatArray refNodeCoords = * this->giveReinforcementNode(numReinfNode)->giveCoordinates();
        refNodeCoords.add(dirVec);
        refNodeNumbers.at(i + 1) = i + 1 + this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode();
        referenceNodes.at(i + 1, 1) = refNodeCoords.at(1);
        referenceNodes.at(i + 1, 2) = refNodeCoords.at(2);
        referenceNodes.at(i + 1, 3) = refNodeCoords.at(3);

        fprintf(outputStream, "node %d coords 3 %e %e %e \n", refNodeNumbers.at(i + 1), refNodeCoords.at(1), refNodeCoords.at(2), refNodeCoords.at(3) );
    }


    //control node
    fprintf(outputStream, "node %d coords 3 %e %e %e\n", controlNode, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );

    // if (this->macroType == _Truss) {
    //     fprintf( outputStream, "node %d coords 3 %e %e %e dofidmask 1 1 \n", refNodeNumbers.at(this->giveNumberOfFibres()) + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );
    // } else if (this->macroType == _Membrane ) {
    //     fprintf( outputStream, "node %d coords 3 %e %e %e dofidmask 4 1 2 4 5 \n", refNodeNumbers.at(this->giveNumberOfFibres()) + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );
    // } else if (this->macroType == _Beam) {
    //     fprintf( outputStream, "node %d coords 3 %e %e %e dofidmask 3 1 7 10 \n", refNodeNumbers.at(this->giveNumberOfFibres()) + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );
    // } else if (this->macroType == _Plate) {
    //     fprintf( outputStream, "node %d coords 3 %e %e %e dofidmask 10 1 2 4 5 7 8 10 11 12 13 \n", refNodeNumbers.at(this->giveNumberOfFibres()) + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );
    // } else if (this->macroType == _3dVoigt) {
    //     fprintf( outputStream, "node %d coords 3 %e %e %e dofidmask 6 1 2 3 4 5 6 \n", refNodeNumbers.at(this->giveNumberOfFibres()) + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );
    // } else if (this->macroType == _3d) {
    //     fprintf( outputStream, "node %d coords 3 %e %e %e dofidmask 9 1 2 3 4 5 6 7 8 9 \n", refNodeNumbers.at(this->giveNumberOfFibres()) + 1, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );
    // }

    printf("Finished writing node data\n");

    //element output - tetrahedra
    int boundaryFlag = 0;
    oofem::IntArray set1;
    for ( int i = 0; i < this->giveNumberOfDelaunayTetras(); i++ ) {
        //First plot all of them
        if ( this->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            location.zero();
            boundaryFlag = 0;
            materialType = 1;
            this->giveDelaunayTetra(i + 1)->giveLocalVertices(nodes);
            newNodes = nodes;

            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < nodes.giveSize(); m++ ) {
                giveDelaunayVertex(nodes.at(m + 1) )->giveCoordinates(coords);
                if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0 ) {
                    //Periodicity in x
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol() ) { //x=xmax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                    //Periodicity in y
                    if ( fabs(coords.at(2) - boundaries.at(4) ) < giveTol() ) { //y=ymax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in z
                    if ( fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //z=zmax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                    //Periodicity in x y
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol()  ||
                         fabs(coords.at(2) - boundaries.at(4) ) < giveTol() ) { //x=xmax or y=ymax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in x z
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol()  ||
                         fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //x=xmax or z=zmax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in y z
                    if ( fabs(coords.at(2) - boundaries.at(4) ) < giveTol()  ||
                         fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //y=ymax or z=zmax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in x y z
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol()  ||
                         fabs(coords.at(2) - boundaries.at(4) ) < giveTol()  ||
                         fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //x=xmax or y=ymax or z=zmax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                }
            }

            if ( boundaryFlag == 0 ) {
                fprintf(outputStream, "ltrspace %d nodes 4 %d %d %d %d ", i + 1, nodes.at(1), nodes.at(2), nodes.at(3), nodes.at(4) );
                fprintf(outputStream, "\n");
                set1.followedBy(i + 1);
            } else {
                fprintf(outputStream, "%s %d nodes 5 %d %d %d %d %d ", boundElemName, i + 1, newNodes.at(1), newNodes.at(2), newNodes.at(3), newNodes.at(4), controlNode);
                fprintf(outputStream, " location 4 %d %d %d %d", location.at(1), location.at(2), location.at(3), location.at(4) );
                fprintf(outputStream, "\n");
                set1.followedBy(i + 1);
            }
        }
    }
    printf("Finished writing Delaunay element (tetrahedra) data\n");

    //element output - beams
    boundaryFlag = 0;
    oofem::IntArray locationBeam(2), set2, beamElL2Gmap;
    beamElL2Gmap.resize(this->giveNumberOfLatticeBeams() );
    beamElL2Gmap.zero();
    for ( int i = 0; i < this->giveNumberOfLatticeBeams(); i++ ) {
        //First plot all of them
        if ( this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 0 || this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 2 ) { //Elements are inside
            locationBeam.zero();
            boundaryFlag = 0;
            this->giveLatticeBeam(i + 1)->giveLocalVertices(nodes);
            newNodes = nodes;

            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < nodes.giveSize(); m++ ) {
                giveReinforcementNode(nodes.at(m + 1) )->giveCoordinates(coords);
                if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0 ) {
                    //Periodicity in x
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol() ) { //x=xmax
                        boundaryFlag = 1;
                        locationBeam.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                    //Periodicity in y
                    if ( fabs(coords.at(2) - boundaries.at(4) ) < giveTol() ) { //y=ymax
                        boundaryFlag = 1;
                        locationBeam.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in z
                    if ( fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //z=zmax
                        boundaryFlag = 1;
                        locationBeam.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                    //Periodicity in x y
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol()  ||
                         fabs(coords.at(2) - boundaries.at(4) ) < giveTol() ) { //x=xmax or y=ymax
                        boundaryFlag = 1;
                        locationBeam.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in x z
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol()  ||
                         fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //x=xmax or z=zmax
                        boundaryFlag = 1;
                        locationBeam.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in y z
                    if ( fabs(coords.at(2) - boundaries.at(4) ) < giveTol()  ||
                         fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //y=ymax or z=zmax
                        boundaryFlag = 1;
                        locationBeam.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in x y z
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol()  ||
                         fabs(coords.at(2) - boundaries.at(4) ) < giveTol()  ||
                         fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //x=xmax or y=ymax or z=zmax
                        boundaryFlag = 1;
                        locationBeam.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->givePeriodicNode();
                    }
                }
            }

            //find the reference node
            int refNode = 0;
            oofem::FloatArray dirVec = this->giveLatticeBeam(i + 1)->giveDirectionVector();
            for ( int j = 0; j < this->giveNumberOfFibres(); j++ ) {
                oofem::FloatArray perpVec(3), firstNodeCoords(3);
                int numReinfNode = this->giveFibre(j + 1)->giveNumberReinforcementNode(1);
                firstNodeCoords = * this->giveReinforcementNode(numReinfNode)->giveCoordinates();
                perpVec.at(1) = referenceNodes.at(j + 1, 1) - firstNodeCoords.at(1);
                perpVec.at(2) = referenceNodes.at(j + 1, 2) - firstNodeCoords.at(2);
                perpVec.at(3) = referenceNodes.at(j + 1, 3) - firstNodeCoords.at(3);
                double sp = dirVec.dotProduct(perpVec, 3);

                if ( fabs(sp) < this->giveTol() ) {
                    refNode = refNodeNumbers.at(j + 1);
                    break;
                }
            }

            if ( boundaryFlag == 0 ) {
                fprintf(outputStream, "libeam3d %d nodes 2 %d %d refnode %d", this->giveNumberOfDelaunayTetras() + i + 1, L2Gmap.at(nodes.at(1) ), L2Gmap.at(nodes.at(2) ), refNode);
                fprintf(outputStream, "\n");
                set2.followedBy(this->giveNumberOfDelaunayTetras() + i + 1);
                beamElL2Gmap.at(i + 1) = this->giveNumberOfDelaunayTetras() + i + 1;
            } else {
                fprintf(outputStream, "%s %d nodes 3 %d %d %d refnode %d", boundBeamElemName, this->giveNumberOfDelaunayTetras() + i + 1, L2Gmap.at(newNodes.at(1) ), L2Gmap.at(newNodes.at(2) ), controlNode, refNode);
                fprintf(outputStream, " location 2 %d %d", locationBeam.at(1), locationBeam.at(2) );
                fprintf(outputStream, "\n");
                set2.followedBy(this->giveNumberOfDelaunayTetras() + i + 1);
                beamElL2Gmap.at(i + 1) = this->giveNumberOfDelaunayTetras() + i + 1;
            }
        }
    }
    printf("Finished writing reinforcement data\n");

    //element output - interface
    oofem::IntArray set3, intElL2Gmap;
    intElL2Gmap.resize(this->giveNumberOfLatticeLinks() );
    intElL2Gmap.zero();
    oofem::FloatArray direction(3);
    for ( int i = 0; i < this->giveNumberOfLatticeLinks(); i++ ) {
        int linkVer1 = this->giveLatticeLink(i + 1)->giveLocalVertex(1);
        int linkVer2 = this->giveLatticeLink(i + 1)->giveLocalVertex(2);
        //recalculate normal vector
        direction = this->giveLatticeLink(i + 1)->giveDirectionVector();
        oofem::FloatArray normalVec(3);
        oofem::FloatArray beamNode = * this->giveReinforcementNode(linkVer1)->giveCoordinates();
        oofem::FloatArray tetraNode = * this->giveDelaunayVertex(linkVer2)->giveCoordinates();
        normalVec.add(beamNode);
        normalVec.subtract(tetraNode);
        normalVec.normalize();
        if ( reinfMirrorNodeList.contains(linkVer1) && mirrorNodeList.contains(linkVer2) ) { //to make sure that the link on master surface is not omitted
            //             fprintf( outputStream, "intelpoint %d nodes 2 %d %d normal 3 %e %e %e length %e\n", this->giveNumberOfDelaunayTetras() + this->giveNumberOfLatticeBeams() + i + 1, linkVer2, L2Gmap.at(linkVer1), normalVec.at(1), normalVec.at(2), normalVec.at(3), this->giveLatticeLink(i+1)->giveAssociatedLength() );

            fprintf(outputStream, "bondlink3d %d nodes 2 %d %d dirvector 3 %e %e %e length %e length_end %e diameter %e\n", this->giveNumberOfDelaunayTetras() + this->giveNumberOfLatticeBeams() + i + 1, L2Gmap.at(linkVer1), linkVer2, direction.at(1), direction.at(2), direction.at(3), this->giveLatticeLink(i + 1)->giveAssociatedLength(), this->giveLatticeLink(i + 1)->giveL_end(), this->giveLatticeLink(i + 1)->giveDiameter() );
            set3.followedBy(this->giveNumberOfDelaunayTetras() + this->giveNumberOfLatticeBeams() + i + 1);
            intElL2Gmap.at(i + 1) = this->giveNumberOfDelaunayTetras() + this->giveNumberOfLatticeBeams() + i + 1;
        }
    }
    printf("Finished writing interface data\n");

    //split into sets (along fibres)
    std::vector < std::vector < int >> beamSets(this->giveNumberOfFibres() );
    std::vector < std::vector < int >> interfaceSets(this->giveNumberOfFibres() );
    oofem::IntArray beamsPerFibre(this->giveNumberOfFibres() ), firstNodeinFibres(this->giveNumberOfFibres() );
    for ( int i = 0; i < this->giveNumberOfFibres(); i++ ) {
        beamsPerFibre.at(i + 1) = this->giveFibre(i + 1)->NbOfReinfNodes() - 1;
        int beamNo(0);
        if ( i != 0 ) {
            for ( int r = 1; r <= i; r++ ) {
                beamNo += beamsPerFibre.at(r);
            }
        }
        //first split reinforcement
        for ( int j = 0; j < this->giveFibre(i + 1)->NbOfReinfNodes() - 1; j++ ) {
            if ( this->giveLatticeBeam(beamNo + j + 1)->giveOutsideFlag() == 0 || this->giveLatticeBeam(beamNo + j + 1)->giveOutsideFlag() == 2 ) {
                ( beamSets [ i ] ).push_back(beamNo + j + 1);
                oofem::IntArray beamNodes;
                this->giveLatticeBeam(beamNo + j + 1)->giveLocalVertices(beamNodes);
                //for the beam elements, get the corresponding link elements
                for ( int k = 0; k < this->giveNumberOfLatticeLinks(); k++ ) {
                    oofem::IntArray linkNodes;
                    this->giveLatticeLink(k + 1)->giveLocalVertices(linkNodes);
                    if ( beamNodes.at(1) == linkNodes.at(1)  ) {
                        ( interfaceSets [ i ] ).push_back(k + 1);
                    }
                }
            }
        }
        firstNodeinFibres.at(i + 1) = L2Gmap.at(this->giveLatticeLink(interfaceSets [ i ] [ 0 ])->giveLocalVertex(1) );
    }
    printf("Finished splitting into sets\n");

    fprintf(outputStream, "SimpleCS 1 material 1 set 1\n");
    for ( int j = 0; j < this->giveNumberOfFibres(); j++ ) {
        double diam = this->giveFibre(j + 1)->giveDiameter();

        fprintf(outputStream, "FiberedCS %d ", 2 * ( j + 1 ) );
        fprintf(outputStream, "fibermaterials 16 %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d ", 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ) );
        fprintf(outputStream, "thicks 16 %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e ", 0.25 * diam, 0.25 * diam, 0.25 * diam, 0.25 * diam, 0.125 * diam, 0.25 * diam, 0.25 * diam, 0.125 * diam, 0.25 * diam, 0.25 * diam, 0.25 * diam, 0.25 * diam, 0.125 * diam, 0.25 * diam, 0.25 * diam, 0.125 * diam);
        fprintf(outputStream, "widths 16 %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e ", 0.25 * diam, 0.25 * diam, 0.25 * diam, 0.25 * diam, 0.125 * diam, 0.25 * diam, 0.25 * diam, 0.125 * diam, 0.25 * diam, 0.25 * diam, 0.25 * diam, 0.25 * diam, 0.125 * diam, 0.25 * diam, 0.25 * diam, 0.125 * diam);
        fprintf(outputStream, "thick %e width %e ", diam, diam);
        fprintf(outputStream, "fiberycentrecoords 16 %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e ", 0.375 * diam, 0.125 * diam, -0.125 * diam, -0.375 * diam, 0.3125 * diam, 0.125 * diam, -0.125 * diam, -0.3125 * diam, 0.375 * diam, 0.125 * diam, -0.125 * diam, -0.375 * diam, 0.3125 * diam, 0.125 * diam, -0.125 * diam, -0.3125 * diam);
        fprintf(outputStream, "fiberzcentrecoords 16 %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e ", -0.125 * diam, -0.125 * diam, -0.125 * diam, -0.125 * diam, -0.3125 * diam, -0.375 * diam, -0.375 * diam, -0.3125 * diam, 0.125 * diam, 0.125 * diam, 0.125 * diam, 0.125 * diam, 0.3125 * diam, 0.375 * diam, 0.375 * diam, 0.3125 * diam);
        fprintf(outputStream, "set %d\n", 2 * ( j + 2 ) );
        fprintf(outputStream, "SimpleCS %d material %d set %d\n", 2 * ( j + 1 ) + 1, 2 * ( j + 1 ) + 1, 2 * ( j + 2 ) + 1);
    }
    fprintf(outputStream, "idm1 1 talpha 0 d 0.0 e 30.e9 e0 100.e-6 wf 33.e-6 n 0.2 equivstraintype 1\n");
    for ( int j = 0; j < this->giveNumberOfFibres(); j++ ) {
        double s3 = 5. / 24 * this->giveFibre(j + 1)->giveDiameter() + 7. / 3000; //approximate interpolation
        fprintf(outputStream, "MisesMat %d d 0.0 E 2.e11 n 0.3 sig0 500e6 H 0. omega_crit 0 a 0 tAlpha 0.0 yieldtol 1.e-10\n", 2 * ( j + 1 ) );
        fprintf(outputStream, "linkslip %d talpha 0. d 0. kn 1.e12 a1 1000 t0 13.7e6 type 1 s1 1.e-3 alpha 0.4\n", 2 * ( j + 1 ) + 1);
    }

    if ( this->macroType == _Beam ) {
        //     fprintf(outputStream, "BoundaryCondition 3 loadTimeFunction 1 dofs 3 1 7 10 values 3 0 0 0 set %d\n", 3 + 2*fibreList.size1()+1);
        fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 dofs 3 1 7 10 values 3 2.6e-3 0 0 set %d\n", 3 + 2 * converter::size1(fibreList) + 1);
    } else if ( this->macroType == _Plate )     {
        fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 dofs 7 1 7 8 10 11 12 13 values 7 2.e-3 0 0 0 0 0 0 set %d\n", 3 + 2 * converter::size1(fibreList) + 1);
    } else   {
        printf("output for this macrotype not yet implemented in grid.C\n");
    }


    fprintf(outputStream, "BoundaryCondition 2 loadTimeFunction 1 dofs 3 1 2 3 values 3 0 0 0 set %d\n", 3 + 2 * converter::size1(fibreList) + 2);
    fprintf(outputStream, "BoundaryCondition 3 loadTimeFunction 1 dofs 2 2 3 values 2 0 0 set %d\n", 3 + 2 * converter::size1(fibreList) + 3);
    fprintf(outputStream, "BoundaryCondition 4 loadTimeFunction 1 dofs 1 3 values 1 0 set %d\n", 3 + 2 * converter::size1(fibreList) + 4);
    fprintf(outputStream, "PiecewiseLinFunction 1 nPoints 2 t 2 0. 200. f(t) 2 0. 1.\n");
    fprintf(outputStream, "set 1 elements %d ", set1.giveSize() );
    for ( int i = 0; i < set1.giveSize(); i++ ) {
        fprintf(outputStream, "%d ", set1.at(i + 1) );
    }
    fprintf(outputStream, "\n");
    fprintf(outputStream, "set 2 elements %d ", set2.giveSize() );
    for ( int i = 0; i < set2.giveSize(); i++ ) {
        fprintf(outputStream, "%d ", set2.at(i + 1) );
    }
    fprintf(outputStream, "\n");
    fprintf(outputStream, "set 3 elements %d ", set3.giveSize() );
    for ( int i = 0; i < set3.giveSize(); i++ ) {
        fprintf(outputStream, "%d ", set3.at(i + 1) );
    }
    fprintf(outputStream, "\n");
    for ( int i = 0; i < this->giveNumberOfFibres(); i++ ) {
        fprintf(outputStream, "set %d elements %d ", 3 + 2 * ( i + 1 ) - 1, beamSets [ i ].size() );
        for ( int j = 0; j < beamSets [ i ].size(); j++ ) {
            fprintf(outputStream, "%d ", beamElL2Gmap.at(beamSets [ i ] [ j ]) );
        }
        fprintf(outputStream, "\n");
        fprintf(outputStream, "set %d elements %d ", 3 + 2 * ( i + 1 ), interfaceSets [ i ].size() );
        for ( int j = 0; j < interfaceSets [ i ].size(); j++ ) {
            fprintf(outputStream, "%d ", intElL2Gmap.at(interfaceSets [ i ] [ j ]) );
        }
        fprintf(outputStream, "\n");
    }
    fprintf(outputStream, "set %d nodes 1 %d \n", 3 + 2 * converter::size1(fibreList) + 1, controlNode);
    fprintf(outputStream, "set %d nodes 1 %d \n", 3 + 2 * converter::size1(fibreList) + 2, firstNode);
    fprintf(outputStream, "set %d nodes 1 %d \n", 3 + 2 * converter::size1(fibreList) + 3, secondNode);
    fprintf(outputStream, "set %d nodes 1 %d \n", 3 + 2 * converter::size1(fibreList) + 4, thirdNode);

    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#REACTION number %d dof 1\n", controlNode);
    fprintf(outputStream, "#NODE number %d dof 1 unknown d\n", controlNode);
    fprintf(outputStream, "#TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");

    return;
}

void
Grid::give3DRCPeriodicSMOutput2(const std::string &fileName)
/**
 * Alternative (2) output for periodic mesh of a 3D reinforced concrete RVE.  Concrete is modelled with tetrahedras (LTRSpace and LTRSpaceBoundary),
 * reinforcement is modelled with beam elements (LIBeam3d and LIBeam3dBoundary), and the interface is modelled with link elements (BondLink3d). Four link elements are used to anchor the rebar node in the tetrahedron.
 *
 * @authors: Adam Sciegaj, Peter Grassl
 */
{
    FILE *outputStream = converter::fopen_or_die(fileName, "w");


    printf("Processing nodes\n");

    oofem::IntArray periodicFlag(3);
    this->givePeriodicityFlag(periodicFlag);

    oofem::FloatArray boundaries(6);
    this->giveRegion(1)->defineBoundaries(boundaries);
    oofem::FloatArray specimenDimension(3);
    specimenDimension.at(1) = boundaries.at(2) - boundaries.at(1);
    specimenDimension.at(2) = boundaries.at(4) - boundaries.at(3);
    specimenDimension.at(3) = boundaries.at(6) - boundaries.at(5);

    int numberOfNodes, numberOfTetras, numberOfBeams, numberOfLinks;
    oofem::FloatArray coords, coordTest;
    int materialType = 1;
    oofem::IntArray nodes, location(4);
    oofem::IntArray newNodes, imageNodeList, mirrorNodeList;
    oofem::IntArray reinfImageNodeList, reinfMirrorNodeList;

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        //Determine which nodes are on mirror/image boundary
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coordTest);

            if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0 ) {
                //Periodicity in x
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                //Periodicity in y
                if ( fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in z
                if ( fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                //Periodicity in x y
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol()  ||
                     fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in x z
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol()  ||
                     fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in y z
                if ( fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol()  ||
                     fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in x y z
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol()  ||
                     fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol()  ||
                     fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( imageNodeList.contains(i + 1) == 0 ) {
                        imageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( mirrorNodeList.contains(i + 1) == 0 ) {
                        mirrorNodeList.followedBy(i + 1);
                    }
                }
            }
        }
    }

    for ( int i = 0; i < this->giveNumberOfReinforcementNode(); i++ ) {
        this->giveReinforcementNode(i + 1)->giveCoordinates(coordTest);
        if ( this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 0 || this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 2 ) {
            //Determine if on mirror or image boundary
            this->giveReinforcementNode(i + 1)->giveCoordinates(coordTest);

            if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0 ) {
                //Periodicity in x
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol() ) {
                    if ( reinfImageNodeList.contains(i + 1) == 0 ) {
                        reinfImageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( reinfMirrorNodeList.contains(i + 1) == 0 ) {
                        reinfMirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                //Periodicity in y
                if ( fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol() ) {
                    if ( reinfImageNodeList.contains(i + 1) == 0 ) {
                        reinfImageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( reinfMirrorNodeList.contains(i + 1) == 0 ) {
                        reinfMirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in z
                if ( fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( reinfImageNodeList.contains(i + 1) == 0 ) {
                        reinfImageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( reinfMirrorNodeList.contains(i + 1) == 0 ) {
                        reinfMirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                //Periodicity in x y
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol()  ||
                     fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol() ) {
                    if ( reinfImageNodeList.contains(i + 1) == 0 ) {
                        reinfImageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( reinfMirrorNodeList.contains(i + 1) == 0 ) {
                        reinfMirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in x z
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol()  ||
                     fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( reinfImageNodeList.contains(i + 1) == 0 ) {
                        reinfImageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( reinfMirrorNodeList.contains(i + 1) == 0 ) {
                        reinfMirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in y z
                if ( fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol()  ||
                     fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( reinfImageNodeList.contains(i + 1) == 0 ) {
                        reinfImageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( reinfMirrorNodeList.contains(i + 1) == 0 ) {
                        reinfMirrorNodeList.followedBy(i + 1);
                    }
                }
            } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                //Periodicity in x y z
                if ( fabs(coordTest.at(1) - boundaries.at(2) ) < giveTol()  ||
                     fabs(coordTest.at(2) - boundaries.at(4) ) < giveTol()  ||
                     fabs(coordTest.at(3) - boundaries.at(6) ) < giveTol() ) {
                    if ( reinfImageNodeList.contains(i + 1) == 0 ) {
                        reinfImageNodeList.followedBy(i + 1);
                    }
                } else {
                    if ( reinfMirrorNodeList.contains(i + 1) == 0 ) {
                        reinfMirrorNodeList.followedBy(i + 1);
                    }
                }
            }
        }
    }

    //total number of nodes: number of master tetra nodes + number of master reinforcement nodes +
    // reference nodes (assuming straight reinforcement -- 1 per fibre) + control node
    numberOfNodes = mirrorNodeList.giveSize() + reinfMirrorNodeList.giveSize() + converter::size1(fibreList) + 1;

    //Determine the number of Delaunay tetrahedra in the domain
    numberOfTetras = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayTetras(); i++ ) {
        if (  this->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfTetras++;
        }
    }

    numberOfBeams = 0;
    for ( int i = 0; i < this->giveNumberOfLatticeBeams(); i++ ) {
        if ( this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 0 || this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfBeams++;
        }
    }

    numberOfLinks = 0;
    for ( int i = 0; i < this->giveNumberOfLatticeLinks(); i++ ) {
        int linkVer1 = this->giveLatticeLink(i + 1)->giveLocalVertex(1);
        int linkVer2 = this->giveLatticeLink(i + 1)->giveLocalVertex(2);
        if ( reinfMirrorNodeList.contains(linkVer1) && mirrorNodeList.contains(linkVer2) ) {
            numberOfLinks++;
        }
    }

    //Sum up all elements. Here take each lattice link times 4
    int numberOfElem = numberOfTetras + numberOfBeams + 4. * numberOfLinks;

    int controlNode = this->giveNumberOfFibres() + this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode() + 1;

    //OOFEM INPUT BEGINS
    //potentially different cross section, material and bond properties for each fiber
    fprintf(outputStream, "rve.out\n");
    fprintf(outputStream, "Periodic reinforced concrete RVE in 3D\n");
    fprintf(outputStream, "StaticStructural nsteps 200 nmodules 3 initialguess 1 lstype 3 smtype 7 stiffmode 2 rtolv 1.e-3 maxiter 10000\n");
   fprintf(outputStream, "vtkxmlperiodic tstep_all primvars 1 1 vars 2 4 1 stype 0 regionsets 1 1\n");
    fprintf(outputStream, "vtkxmlperiodic tstep_all primvars 1 1 vars 2 7 8 stype 0 regionsets 1 2\n");
    fprintf(outputStream, "vtkxmlperiodic tstep_all primvars 1 1 ipvars 2 98 99 regionsets 1 3\n");
    //    fprintf(outputStream, "matlab tstep_all mesh data reactionforces integrationpoints internalvars 7 4 1 12 7 8 98 99\n");
    fprintf(outputStream, "domain 3d\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_output {%d}\n", controlNode);
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect %d nmat %d nbc 4 nic 0 nltf 1 nset %d\n", numberOfNodes, numberOfElem,
            2 * converter::size1(fibreList) + 1,  2 * converter::size1(fibreList) + 1, 7 + 2 * converter::size1(fibreList) );

    int firstFlag = 0;

    //Find three nodes which can be used for constraining the specimen.
    //Translations can be fixed by setting DOFs of first node to zero.
    int firstNode = 0, secondNode = 0, thirdNode = 0;
    oofem::FloatArray firstNodeCoords(3);

    //node output - tetrahedra
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);

            if ( mirrorNodeList.isEmpty() == 0 ) {
                if ( mirrorNodeList.contains(i + 1) ) {
                    fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
                }
            }
        }
    }

    //Define nodes for supports to block rigid body rotations
    //These nodes are control nodes defined in the mesh.in file
    firstNode = 9;
    secondNode = 10;
    thirdNode = 11;

    //node output - beams
    oofem::IntArray L2Gmap;
    L2Gmap.resize(this->giveNumberOfReinforcementNode() );
    L2Gmap.zero();
    for ( int i = 0; i < this->giveNumberOfReinforcementNode(); i++ ) {
        if ( this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 0 || this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveReinforcementNode(i + 1)->giveCoordinates(coords);

            if ( reinfMirrorNodeList.isEmpty() == 0 ) {
                if ( reinfMirrorNodeList.contains(i + 1) ) {
                    L2Gmap.at(i + 1) = this->giveNumberOfDelaunayVertices() + i + 1;
                    fprintf(outputStream, "node %d coords 3 %e %e %e\n", L2Gmap.at(i + 1), coords.at(1), coords.at(2), coords.at(3) );
                }
            }
        }
    }

    //reference nodes for 3d beams. Assuming only straight reinforcement bars.
    //Therefore, only one refnode per fibre. For curved reinforcement, each beam
    //element would require a seperate reference node.
    oofem::FloatMatrix R(3, 3);
    R.zero();
    R.at(1, 2) = -1;
    R.at(2, 1) = 1;
    R.at(3, 3) = 1;
    oofem::IntArray refNodeNumbers(this->giveNumberOfFibres() );
    oofem::FloatMatrix referenceNodes;
    referenceNodes.resize(this->giveNumberOfFibres(), 3);
    for ( int i = 0; i < this->giveNumberOfFibres(); i++) {
        oofem::FloatArray dirVec = this->giveFibre(i + 1)->giveDirVector();
        dirVec.rotatedWith(R, 'n');
        int numReinfNode = this->giveFibre(i + 1)->giveNumberReinforcementNode(1);
        oofem::FloatArray refNodeCoords = * this->giveReinforcementNode(numReinfNode)->giveCoordinates();
        refNodeCoords.add(dirVec);
        refNodeNumbers.at(i + 1) = i + 1 + this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode();
        referenceNodes.at(i + 1, 1) = refNodeCoords.at(1);
        referenceNodes.at(i + 1, 2) = refNodeCoords.at(2);
        referenceNodes.at(i + 1, 3) = refNodeCoords.at(3);

        fprintf(outputStream, "node %d coords 3 %e %e %e \n", refNodeNumbers.at(i + 1), refNodeCoords.at(1), refNodeCoords.at(2), refNodeCoords.at(3) );
    }


    //control node
    fprintf(outputStream, "node %d coords 3 %e %e %e\n", controlNode, specimenDimension.at(1), specimenDimension.at(2), specimenDimension.at(3) );

    printf("Finished writing node data\n");

    //element output - tetrahedra
    int boundaryFlag = 0;
    oofem::IntArray set1;
    for ( int i = 0; i < this->giveNumberOfDelaunayTetras(); i++ ) {
        //First plot all of them
        if ( this->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 0 ) { //Elements are inside
            location.zero();
            boundaryFlag = 0;
            materialType = 1;
            this->giveDelaunayTetra(i + 1)->giveLocalVertices(nodes);
            newNodes = nodes;

            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < nodes.giveSize(); m++ ) {
                giveDelaunayVertex(nodes.at(m + 1) )->giveCoordinates(coords);
                if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0 ) {
                    //Periodicity in x
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol() ) { //x=xmax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                    //Periodicity in y
                    if ( fabs(coords.at(2) - boundaries.at(4) ) < giveTol() ) { //y=ymax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in z
                    if ( fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //z=zmax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                    //Periodicity in x y
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol()  ||
                         fabs(coords.at(2) - boundaries.at(4) ) < giveTol() ) { //x=xmax or y=ymax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in x z
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol()  ||
                         fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //x=xmax or z=zmax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in y z
                    if ( fabs(coords.at(2) - boundaries.at(4) ) < giveTol()  ||
                         fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //y=ymax or z=zmax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in x y z
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol()  ||
                         fabs(coords.at(2) - boundaries.at(4) ) < giveTol()  ||
                         fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //x=xmax or y=ymax or z=zmax
                        boundaryFlag = 1;
                        location.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveDelaunayVertex(nodes.at(m + 1) )->givePeriodicNode();
                    }
                }
            }

            if ( boundaryFlag == 0 ) {
                fprintf(outputStream, "ltrspace %d nodes 4 %d %d %d %d ", i + 1, nodes.at(1), nodes.at(2), nodes.at(3), nodes.at(4) );
                fprintf(outputStream, "\n");
                set1.followedBy(i + 1);
            } else {
                fprintf(outputStream, "%s %d nodes 5 %d %d %d %d %d ", boundElemName, i + 1, newNodes.at(1), newNodes.at(2), newNodes.at(3), newNodes.at(4), controlNode);
                fprintf(outputStream, " location 4 %d %d %d %d", location.at(1), location.at(2), location.at(3), location.at(4) );
                fprintf(outputStream, "\n");
                set1.followedBy(i + 1);
            }
        }
    }
    printf("Finished writing Delaunay element (tetrahedra) data\n");

    //element output - beams
    boundaryFlag = 0;
    oofem::IntArray locationBeam(2), set2, beamElL2Gmap;
    beamElL2Gmap.resize(this->giveNumberOfLatticeBeams() );
    beamElL2Gmap.zero();
    for ( int i = 0; i < this->giveNumberOfLatticeBeams(); i++ ) {
        //First plot all of them
        if ( this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 0 || this->giveLatticeBeam(i + 1)->giveOutsideFlag() == 2 ) { //Elements are inside
            locationBeam.zero();
            boundaryFlag = 0;
            this->giveLatticeBeam(i + 1)->giveLocalVertices(nodes);
            newNodes = nodes;

            //Go through nodes and replace the ones outside with periodic nodes
            for ( int m = 0; m < nodes.giveSize(); m++ ) {
                giveReinforcementNode(nodes.at(m + 1) )->giveCoordinates(coords);
                if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 0 ) {
                    //Periodicity in x
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol() ) { //x=xmax
                        boundaryFlag = 1;
                        locationBeam.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                    //Periodicity in y
                    if ( fabs(coords.at(2) - boundaries.at(4) ) < giveTol() ) { //y=ymax
                        boundaryFlag = 1;
                        locationBeam.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in z
                    if ( fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //z=zmax
                        boundaryFlag = 1;
                        locationBeam.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 0 ) {
                    //Periodicity in x y
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol()  ||
                         fabs(coords.at(2) - boundaries.at(4) ) < giveTol() ) { //x=xmax or y=ymax
                        boundaryFlag = 1;
                        locationBeam.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 0 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in x z
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol()  ||
                         fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //x=xmax or z=zmax
                        boundaryFlag = 1;
                        locationBeam.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 0 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in y z
                    if ( fabs(coords.at(2) - boundaries.at(4) ) < giveTol()  ||
                         fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //y=ymax or z=zmax
                        boundaryFlag = 1;
                        locationBeam.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->givePeriodicNode();
                    }
                } else if ( periodicityFlag.at(1) == 1 && periodicityFlag.at(2) == 1 && periodicityFlag.at(3) == 1 ) {
                    //Periodicity in x y z
                    if ( fabs(coords.at(1) - boundaries.at(2) ) < giveTol()  ||
                         fabs(coords.at(2) - boundaries.at(4) ) < giveTol()  ||
                         fabs(coords.at(3) - boundaries.at(6) ) < giveTol() ) { //x=xmax or y=ymax or z=zmax
                        boundaryFlag = 1;
                        locationBeam.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->giveLocation();
                        newNodes.at(m + 1) = this->giveReinforcementNode(nodes.at(m + 1) )->givePeriodicNode();
                    }
                }
            }

            //find the reference node
            int refNode = 0;
            oofem::FloatArray dirVec = this->giveLatticeBeam(i + 1)->giveDirectionVector();
            for ( int j = 0; j < this->giveNumberOfFibres(); j++ ) {
                oofem::FloatArray perpVec(3), firstNodeCoords(3);
                int numReinfNode = this->giveFibre(j + 1)->giveNumberReinforcementNode(1);
                firstNodeCoords = * this->giveReinforcementNode(numReinfNode)->giveCoordinates();
                perpVec.at(1) = referenceNodes.at(j + 1, 1) - firstNodeCoords.at(1);
                perpVec.at(2) = referenceNodes.at(j + 1, 2) - firstNodeCoords.at(2);
                perpVec.at(3) = referenceNodes.at(j + 1, 3) - firstNodeCoords.at(3);
                double sp = dirVec.dotProduct(perpVec, 3);

                if ( fabs(sp) < this->giveTol() ) {
                    refNode = refNodeNumbers.at(j + 1);
                    break;
                }
            }

            if ( boundaryFlag == 0 ) {
                fprintf(outputStream, "libeam3d %d nodes 2 %d %d refnode %d", this->giveNumberOfDelaunayTetras() + i + 1, L2Gmap.at(nodes.at(1) ), L2Gmap.at(nodes.at(2) ), refNode);
                fprintf(outputStream, "\n");
                set2.followedBy(this->giveNumberOfDelaunayTetras() + i + 1);
                beamElL2Gmap.at(i + 1) = this->giveNumberOfDelaunayTetras() + i + 1;
            } else {
                fprintf(outputStream, "%s %d nodes 3 %d %d %d refnode %d", boundBeamElemName, this->giveNumberOfDelaunayTetras() + i + 1, L2Gmap.at(newNodes.at(1) ), L2Gmap.at(newNodes.at(2) ), controlNode, refNode);
                fprintf(outputStream, " location 2 %d %d", locationBeam.at(1), locationBeam.at(2) );
                fprintf(outputStream, "\n");
                set2.followedBy(this->giveNumberOfDelaunayTetras() + i + 1);
                beamElL2Gmap.at(i + 1) = this->giveNumberOfDelaunayTetras() + i + 1;
            }
        }
    }
    printf("Finished writing reinforcement data\n");

    //element output - interface
    oofem::IntArray set3, intElL2Gmap;
    intElL2Gmap.resize(this->giveNumberOfLatticeLinks() );
    intElL2Gmap.zero();
    oofem::FloatArray direction(3);
    for ( int i = 0; i < this->giveNumberOfLatticeLinks(); i++ ) {
        int linkVer1 = this->giveLatticeLink(i + 1)->giveLocalVertex(1);
        int linkVer2 = this->giveLatticeLink(i + 1)->giveLocalVertex(2);
        //recalculate normal vector
        direction = this->giveLatticeLink(i + 1)->giveDirectionVector();
        oofem::FloatArray normalVec(3);
        oofem::FloatArray beamNodeCoords = * this->giveReinforcementNode(linkVer1)->giveCoordinates();
        //        oofem::FloatArray tetraNodeCoords = *this->giveDelaunayVertex(linkVer2)->giveCoordinates();
        //        normalVec.add(beamNodeCoords); normalVec.subtract(tetraNodeCoords); normalVec.normalize();

        //Find tetra in which rebar node is and create for links to nodes on vertices
        //Get the tetras connected to the tetra node.
        //Then check for each tetra if rebar node is inside.

        oofem::FloatArray tetraCoords(12);
        oofem::FloatArray barycentres(3);
        oofem::IntArray localTetras;
        int targetTetra;
        this->giveDelaunayVertex(linkVer2)->giveLocalTetras(localTetras);
        for (int k = 1; k <= localTetras.giveSize(); k++) {
            this->giveDelaunayTetra(localTetras.at(k) )->giveCoordinates(tetraCoords);
            giveTetrahedronBarycentres(barycentres, tetraCoords, beamNodeCoords);
            if ( barycentres.at(1) > 0 && barycentres.at(2) > 0 & barycentres.at(3) > 0 ) {
                targetTetra = localTetras.at(k);
                break;//There can only be one target tetra
            }
        }
        oofem::IntArray targetTetraNodes(4);
        this->giveDelaunayTetra(targetTetra)->giveLocalVertices(targetTetraNodes);

        if ( reinfMirrorNodeList.contains(linkVer1) && mirrorNodeList.contains(linkVer2) ) { //to make sure that the link on master surface is not omitted
            for (int k = 1; k <= 4; k++) {
                fprintf(outputStream, "bondlink3d %d nodes 2 %d %d dirvector 3 %e %e %e length %e length_end %e diameter %e\n", this->giveNumberOfDelaunayTetras() + this->giveNumberOfLatticeBeams() + 4 * i + k, L2Gmap.at(linkVer1), targetTetraNodes.at(k), direction.at(1), direction.at(2), direction.at(3), this->giveLatticeLink(i + 1)->giveAssociatedLength(), this->giveLatticeLink(i + 1)->giveL_end(), this->giveLatticeLink(i + 1)->giveDiameter() );
                set3.followedBy(this->giveNumberOfDelaunayTetras() + this->giveNumberOfLatticeBeams() + i + 1);
                intElL2Gmap.at(i + 1) = this->giveNumberOfDelaunayTetras() + this->giveNumberOfLatticeBeams() + i + 1;
            }
        }
    }
    printf("Finished writing interface data\n");

    //split into sets (along fibres)
    std::vector < std::vector < int >> beamSets(this->giveNumberOfFibres() );
    std::vector < std::vector < int >> interfaceSets(this->giveNumberOfFibres() );
    oofem::IntArray beamsPerFibre(this->giveNumberOfFibres() ), firstNodeinFibres(this->giveNumberOfFibres() );
    for ( int i = 0; i < this->giveNumberOfFibres(); i++ ) {
        beamsPerFibre.at(i + 1) = this->giveFibre(i + 1)->NbOfReinfNodes() - 1;
        int beamNo(0);
        if ( i != 0 ) {
            for ( int r = 1; r <= i; r++ ) {
                beamNo += beamsPerFibre.at(r);
            }
        }
        //first split reinforcement
        for ( int j = 0; j < this->giveFibre(i + 1)->NbOfReinfNodes() - 1; j++ ) {
            if ( this->giveLatticeBeam(beamNo + j + 1)->giveOutsideFlag() == 0 || this->giveLatticeBeam(beamNo + j + 1)->giveOutsideFlag() == 2 ) {
                ( beamSets [ i ] ).push_back(beamNo + j + 1);
                oofem::IntArray beamNodes;
                this->giveLatticeBeam(beamNo + j + 1)->giveLocalVertices(beamNodes);
                //for the beam elements, get the corresponding link elements
                for ( int k = 0; k < this->giveNumberOfLatticeLinks(); k++ ) {
                    oofem::IntArray linkNodes;
                    this->giveLatticeLink(k + 1)->giveLocalVertices(linkNodes);
                    if ( beamNodes.at(1) == linkNodes.at(1)  ) {
                        ( interfaceSets [ i ] ).push_back(k + 1);
                    }
                }
            }
        }
        firstNodeinFibres.at(i + 1) = L2Gmap.at(this->giveLatticeLink(interfaceSets [ i ] [ 0 ])->giveLocalVertex(1) );
    }
    printf("Finished splitting into sets\n");

    fprintf(outputStream, "SimpleCS 1 material 1 set 1\n");
    for ( int j = 0; j < this->giveNumberOfFibres(); j++ ) {
        double diam = this->giveFibre(j + 1)->giveDiameter();

        fprintf(outputStream, "FiberedCS %d ", 2 * ( j + 1 ) );
        fprintf(outputStream, "fibermaterials 16 %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d ", 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ), 2 * ( j + 1 ) );
        fprintf(outputStream, "thicks 16 %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e ", 0.25 * diam, 0.25 * diam, 0.25 * diam, 0.25 * diam, 0.125 * diam, 0.25 * diam, 0.25 * diam, 0.125 * diam, 0.25 * diam, 0.25 * diam, 0.25 * diam, 0.25 * diam, 0.125 * diam, 0.25 * diam, 0.25 * diam, 0.125 * diam);
        fprintf(outputStream, "widths 16 %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e ", 0.25 * diam, 0.25 * diam, 0.25 * diam, 0.25 * diam, 0.125 * diam, 0.25 * diam, 0.25 * diam, 0.125 * diam, 0.25 * diam, 0.25 * diam, 0.25 * diam, 0.25 * diam, 0.125 * diam, 0.25 * diam, 0.25 * diam, 0.125 * diam);
        fprintf(outputStream, "thick %e width %e ", diam, diam);
        fprintf(outputStream, "fiberycentrecoords 16 %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e ", 0.375 * diam, 0.125 * diam, -0.125 * diam, -0.375 * diam, 0.3125 * diam, 0.125 * diam, -0.125 * diam, -0.3125 * diam, 0.375 * diam, 0.125 * diam, -0.125 * diam, -0.375 * diam, 0.3125 * diam, 0.125 * diam, -0.125 * diam, -0.3125 * diam);
        fprintf(outputStream, "fiberzcentrecoords 16 %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e ", -0.125 * diam, -0.125 * diam, -0.125 * diam, -0.125 * diam, -0.3125 * diam, -0.375 * diam, -0.375 * diam, -0.3125 * diam, 0.125 * diam, 0.125 * diam, 0.125 * diam, 0.125 * diam, 0.3125 * diam, 0.375 * diam, 0.375 * diam, 0.3125 * diam);
        fprintf(outputStream, "set %d\n", 2 * ( j + 2 ) );
        fprintf(outputStream, "SimpleCS %d material %d set %d\n", 2 * ( j + 1 ) + 1, 2 * ( j + 1 ) + 1, 2 * ( j + 2 ) + 1);
    }
    fprintf(outputStream, "idm1 1 talpha 0 d 0.0 e 30.e9 e0 100.e-6 wf 33.e-6 n 0.2 equivstraintype 1\n");
    for ( int j = 0; j < this->giveNumberOfFibres(); j++ ) {
        double s3 = 5. / 24 * this->giveFibre(j + 1)->giveDiameter() + 7. / 3000; //approximate interpolation
        fprintf(outputStream, "MisesMat %d d 0.0 E 2.e11 n 0.3 sig0 500e6 H 0. omega_crit 0 a 0 tAlpha 0.0 yieldtol 1.e-10\n", 2 * ( j + 1 ) );
        fprintf(outputStream, "linkslip %d talpha 0. d 0. kn 1.e12 a1 1000 t0 13.7e6 type 1 s1 1.e-3 alpha 0.4\n", 2 * ( j + 1 ) + 1);
    }

    if ( this->macroType == _Beam ) {
        fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 dofs 3 1 7 10 values 3 2.6e-3 0 0 set %d\n", 3 + 2 * converter::size1(fibreList) + 1);
    } else if ( this->macroType == _Plate )     {
        fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 dofs 7 1 7 8 10 11 12 13 values 7 2.e-3 0 0 0 0 0 0 set %d\n", 3 + 2 * converter::size1(fibreList) + 1);
    } else   {
        printf("output for this macrotype %s not yet implemented in grid.C\n", this - macroType);
    }

    fprintf(outputStream, "BoundaryCondition 2 loadTimeFunction 1 dofs 3 1 2 3 values 3 0 0 0 set %d\n", 3 + 2 * converter::size1(fibreList) + 2);
    fprintf(outputStream, "BoundaryCondition 3 loadTimeFunction 1 dofs 2 2 3 values 2 0 0 set %d\n", 3 + 2 * converter::size1(fibreList) + 3);
    fprintf(outputStream, "BoundaryCondition 4 loadTimeFunction 1 dofs 1 3 values 1 0 set %d\n", 3 + 2 * converter::size1(fibreList) + 4);

    fprintf(outputStream, "PiecewiseLinFunction 1 nPoints 2 t 2 0. 200. f(t) 2 0. 1.\n");
    fprintf(outputStream, "set 1 elements %d ", set1.giveSize() );
    for ( int i = 0; i < set1.giveSize(); i++ ) {
        fprintf(outputStream, "%d ", set1.at(i + 1) );
    }
    fprintf(outputStream, "\n");
    fprintf(outputStream, "set 2 elements %d ", set2.giveSize() );
    for ( int i = 0; i < set2.giveSize(); i++ ) {
        fprintf(outputStream, "%d ", set2.at(i + 1) );
    }
    fprintf(outputStream, "\n");
    fprintf(outputStream, "set 3 elements %d ", set3.giveSize() );
    for ( int i = 0; i < set3.giveSize(); i++ ) {
        fprintf(outputStream, "%d ", set3.at(i + 1) );
    }
    fprintf(outputStream, "\n");
    for ( int i = 0; i < this->giveNumberOfFibres(); i++ ) {
        fprintf(outputStream, "set %d elements %d ", 3 + 2 * ( i + 1 ) - 1, beamSets [ i ].size() );
        for ( int j = 0; j < beamSets [ i ].size(); j++ ) {
            fprintf(outputStream, "%d ", beamElL2Gmap.at(beamSets [ i ] [ j ]) );
        }
        fprintf(outputStream, "\n");
        fprintf(outputStream, "set %d elements %d ", 3 + 2 * ( i + 1 ), interfaceSets [ i ].size() );
        for ( int j = 0; j < interfaceSets [ i ].size(); j++ ) {
            fprintf(outputStream, "%d ", intElL2Gmap.at(interfaceSets [ i ] [ j ]) );
        }
        fprintf(outputStream, "\n");
    }
    fprintf(outputStream, "set %d nodes 1 %d \n", 3 + 2 * converter::size1(fibreList) + 1, controlNode);
    fprintf(outputStream, "set %d nodes 1 %d \n", 3 + 2 * converter::size1(fibreList) + 2, firstNode);
    fprintf(outputStream, "set %d nodes 1 %d \n", 3 + 2 * converter::size1(fibreList) + 3, secondNode);
    fprintf(outputStream, "set %d nodes 1 %d \n", 3 + 2 * converter::size1(fibreList) + 4, thirdNode);

    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#REACTION number %d dof 1\n", controlNode);
    fprintf(outputStream, "#NODE number %d dof 1 unknown d\n", controlNode);
    fprintf(outputStream, "#TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");

    printf("Finished the output\n");

    return;
}


void
Grid::give3DGopShaOutput(const std::string &fileName)
{
    //Output for GopSha experiment (Direct tension with two small notches.
    //Only mechanical model


    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords(3);
    oofem::FloatArray coordsOne, coordsTwo, line;
    int materialType = 1;
    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;
    double radius, distanceOne, distanceTwo;

    int supportNode = -1;
    int loadNode = -1;
    int notchOneBottomNode = -1;
    int notchOneTopNode = -1;
    int notchTwoBottomNode = -1;
    int notchTwoTopNode = -1;
    int controlMidBottomNode = -1;
    int controlMidTopNode = -1;

    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ||  this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfNodes++;

            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);

            //Create support and load nodes
            if ( fabs(coords.at(1) - 0.038) < TOL && fabs(coords.at(2) - 0.0) < TOL && fabs(coords.at(3) - 0.019) < TOL ) {
                supportNode = i + 1;
            } else if ( fabs(coords.at(1) - 0.038) < TOL && fabs(coords.at(2) - 0.305) < TOL && fabs(coords.at(3) - 0.019) < TOL )                  {
                loadNode = i + 1;
            } else if ( fabs(coords.at(1) - 0.0) < TOL && fabs(coords.at(2) - 0.146) < TOL && fabs(coords.at(3) - 0.019) < TOL )                  {
                notchOneBottomNode = i + 1;
            } else if ( fabs(coords.at(1) - 0.0) < TOL && fabs(coords.at(2) - 0.159) < TOL && fabs(coords.at(3) - 0.019) < TOL )                  {
                notchOneTopNode = i + 1;
            } else if ( fabs(coords.at(1) - 0.076) < TOL && fabs(coords.at(2) - 0.146) < TOL && fabs(coords.at(3) - 0.019) < TOL )                  {
                notchTwoBottomNode = i + 1;
            } else if ( fabs(coords.at(1) - 0.076) < TOL && fabs(coords.at(2) - 0.159) < TOL && fabs(coords.at(3) - 0.019) < TOL )                  {
                notchTwoTopNode = i + 1;
            } else if ( fabs(coords.at(1) - 0.038) < TOL && fabs(coords.at(2) - 0.111) < TOL && fabs(coords.at(3) - 0.019) < TOL )                  {
                controlMidBottomNode = i + 1;
            } else if ( fabs(coords.at(1) - 0.038) < TOL && fabs(coords.at(2) - 0.194) < TOL && fabs(coords.at(3) - 0.019) < TOL )                  {
                controlMidTopNode = i + 1;
            }
        }
    }

    //Determine the number of Delaunay lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
        this->giveDelaunayVertex(nodes.at(1) )->giveCoordinates(coordsOne);
        this->giveDelaunayVertex(nodes.at(2) )->giveCoordinates(coordsTwo);

        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 && ( this->giveDelaunayLine(i + 1) )->delaunayAreaCheck() == 1 ) {
            numberOfLines++;
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 && this->giveRegion(1)->modifyVoronoiCrossSection(i + 1) == 1 )      {
            numberOfLines++;
        }
    }

    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Mechanical 3D model of GopSha experiment\n");
    fprintf(outputStream, "Nonlinearstatic nmsteps 3 nsteps 1 contextOutputStep 2000 nmodules 2 lstype 3 smtype 7\n");
    fprintf(outputStream, "nsteps 100 rtolv 1.e-3 stiffMode 1 manrmsteps 10 maxiter 400 controllmode 0 steplength 2.e-6 minsteplength 2.e-6 hpcmode 2 hpc 12 %d 2 %d 2 %d 2 %d 2 %d 2 %d 2 hpcw 6 -1 1 -1 1 -1 1 donotfixload lstype 3 smtype 7\n", notchOneBottomNode, notchOneTopNode, notchTwoBottomNode, notchTwoTopNode, controlMidBottomNode, controlMidTopNode);
    fprintf(outputStream, "nsteps 50 rtolv 1.e-3 stiffMode 1 manrmsteps 10 maxiter 400 controllmode 0 steplength 2.e-6 minsteplength 2.e-6 hpcmode 2 hpc 12 %d 2 %d 2 %d 2 %d 2 %d 2 %d 2 hpcw 6 -1 1 -1 1 -1 1 donotfixload lstype 3 smtype 7\n", notchOneBottomNode, notchOneTopNode, notchTwoBottomNode, notchTwoTopNode, controlMidBottomNode, controlMidTopNode);
    fprintf(outputStream, "nsteps 50 rtolv 1.e-3 stiffMode 1 manrmsteps 10 maxiter 400 controllmode 0 steplength 2.e-6 minsteplength 2.e-6 hpcmode 2 hpc 12 %d 2 %d 2 %d 2 %d 2 %d 2 %d 2 hpcw 6 -1 1 -1 1 -1 1 donotfixload lstype 3 smtype 7\n", notchOneBottomNode, notchOneTopNode, notchTwoBottomNode, notchTwoTopNode, controlMidBottomNode, controlMidTopNode);
    fprintf(outputStream, "vtkxmllattice primvars 1 1 tstep_all domain_all cross 1 cellvars 4 46 60 90 111\n");
    fprintf(outputStream, "gpexportmodule vars 2 46 139 tstep_all\n");
    fprintf(outputStream, "domain 3dLattice\n");
    fprintf(outputStream, "OutputManager tstep_all dofman_output { %d %d %d %d %d %d %d %d }\n", loadNode, supportNode, notchOneBottomNode, notchOneTopNode, notchTwoBottomNode, notchTwoTopNode, controlMidBottomNode, controlMidTopNode);
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 3 nmat 3 nbc 3 nic 0 nltf 2 nset 5\n", numberOfNodes, numberOfLines);

    printf("start nodes\n");

    int firstFlag = 0;
    int bottomSetCounter = 0;
    oofem::IntArray bottomSetTemp(this->giveNumberOfDelaunayVertices() );
    int topSetCounter = 0;
    oofem::IntArray topSetTemp(this->giveNumberOfDelaunayVertices() );

    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);

            if ( fabs(coords.at(2) - 0.0) < giveTol() && i + 1 != supportNode ) {
                bottomSetCounter++;
                bottomSetTemp.at(bottomSetCounter) = i + 1;
                fprintf(outputStream, "rigidarmnode %d coords 3 %e %e %e master %d mastermask 6 0 1 0 0 0 1 doftype 6 0 2 0 0 0 2\n", i + 1, coords.at(1), coords.at(2), coords.at(3), supportNode);
            } else if ( fabs(coords.at(2) - 0.305) < giveTol() && i + 1 != loadNode )            {
                topSetCounter++;
                topSetTemp.at(topSetCounter) = i + 1;
                fprintf(outputStream, "rigidarmnode %d coords 3 %e %e %e master %d mastermask 6 0 1 0 0 0 1 doftype 6 0 2 0 0 0 2\n", i + 1, coords.at(1), coords.at(2), coords.at(3), loadNode);
            } else   {
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }

    oofem::IntArray bottomSet(bottomSetCounter);
    bottomSet.zero();
    for (int i = 0; i < bottomSetCounter; i++) {
        bottomSet.at(i + 1) = bottomSetTemp.at(i + 1);
    }

    oofem::IntArray topSet(topSetCounter);
    topSet.zero();
    for (int i = 0; i < topSetCounter; i++) {
        topSet.at(i + 1) = topSetTemp.at(i + 1);
    }

    printf("finished nodes\n");

    int set1Counter = 0, set2Counter = 0, set3Counter = 0;
    oofem::IntArray set1Temp(numberOfLines);
    set1Temp.zero();

    oofem::IntArray set2Temp(numberOfLines);
    set2Temp.zero();

    oofem::IntArray set3Temp(numberOfLines);
    set3Temp.zero();

    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        //Check if element might be in the notch area. For this, nodes 1 or nodes 2 must be in the 13 mm notch length.
        //and cross the middle

        if ( ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) && ( this->giveDelaunayLine(i + 1) )->delaunayAreaCheck() == 1 ) { //Elements are inside
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);

            //=======================================================
            //Deal with different materials (elastic, concrete, notch)
            materialType = 1;
            //	this->giveDelaunayLine(i + 1)->updateMaterial(1);
            this->giveDelaunayVertex(nodes.at(1) )->giveCoordinates(coordsOne);
            this->giveDelaunayVertex(nodes.at(2) )->giveCoordinates(coordsTwo);


            if ( coordsOne.at(2) < 0.1025 && coordsTwo.at(2) < 0.1025 ) {
                set2Counter++;
                set2Temp.at(set2Counter) = i + 1;
                materialType = 2;
                this->giveDelaunayLine(i + 1)->updateMaterial(2);
            } else if ( coordsOne.at(2) >= 0.1025 && coordsOne.at(2) <= 0.2025 || coordsTwo.at(2) >= 0.1025 && coordsTwo.at(2) <= 0.2025 )    {
                //Notch
                if ( ( ( coordsOne.at(2) > 0.1525 - 0.003 / 2. && coordsTwo.at(2) < 0.1525 + 0.003 / 2. ) || ( coordsOne.at(2) < 0.1525 + 0.003 / 2. && coordsTwo.at(2) > 0.1525 - 0.003 / 2. ) ) &&
                     ( ( coordsOne.at(1) < 0.013 && coordsTwo.at(1) < 0.013 ) || ( coordsOne.at(1) > 0.063 && coordsTwo.at(1) > 0.063 ) ) ) {
                    set3Counter++;
                    set3Temp.at(set3Counter) = i + 1;
                    materialType = 3;
                    this->giveDelaunayLine(i + 1)->updateMaterial(3);
                } else {
                    set1Counter++;
                    set1Temp.at(set1Counter) = i + 1;
                    materialType = 1;
                    this->giveDelaunayLine(i + 1)->updateMaterial(1);
                }
            } else if ( coordsOne.at(2) > 0.2025 && coordsTwo.at(2) > 0.2025 ) {
                set2Counter++;
                set2Temp.at(set2Counter) = i + 1;
                materialType = 2;
                this->giveDelaunayLine(i + 1)->updateMaterial(2);
            }

            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "lattice3D %d nodes 2 %d %d crossSect %d mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, "\n");
        }
    }

    oofem::IntArray set1(set1Counter);
    set1.zero();

    oofem::IntArray set2(set2Counter);
    set2.zero();

    oofem::IntArray set3(set3Counter);
    set3.zero();

    for (int i = 0; i < set1Counter; i++) {
        set1.at(i + 1) = set1Temp.at(i + 1);
    }

    for (int i = 0; i < set2Counter; i++) {
        set2.at(i + 1) = set2Temp.at(i + 1);
    }

    for (int i = 0; i < set3Counter; i++) {
        set3.at(i + 1) = set3Temp.at(i + 1);
    }

    fprintf(outputStream, "latticecs 1 material 1\n");
    fprintf(outputStream, "latticecs 2 material 2\n");
    fprintf(outputStream, "latticecs 3 material 3\n");

    fprintf(outputStream, "latticeplastdam 1 d 0 talpha 0. calpha 0. e 50.46e9 a1 0.215 ft 2.44e6 fc 30.e6 wf 50.e-6 ahard 1.e-3 angle1 0.5 flow 0.25 randvars 2 806 807 randgen 2 2 2\n");
    fprintf(outputStream, "latticelinearelastic 2 d 0 talpha 0. calpha 0. e 50.46e11 a1 0.215\n");
    fprintf(outputStream, "latticelinearelastic 3 d 0 talpha 0. calpha 0. e 50.46 a1 0.215\n");
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 dofs 6 1 2 3 4 5 6 values 6 0 0 0 0 0 0 set 4\n");
    fprintf(outputStream, "BoundaryCondition 2 loadTimeFunction 1 dofs 4 1 3 4 5 values 4 0 0 0 0 set 5\n");
    fprintf(outputStream, "NodalLoad 3 loadTimeFunction 1 dofs 6 1 2 3 4 5 6 Components 6 0. 1. 0. 0. 0. 0. set 5\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "InterpolatingFunction 2 name random.dat dim 3\n");

    //Print set1
    fprintf(outputStream, "Set 1 elements %d", set1Counter);
    for (int i = 0; i < set1Counter; i++) {
        fprintf(outputStream, " %d", set1.at(i + 1) );
    }
    fprintf(outputStream, "\n");

    //Print set2
    fprintf(outputStream, "Set 2 elements %d", set2Counter);
    for (int i = 0; i < set2Counter; i++) {
        fprintf(outputStream, " %d", set2.at(i + 1) );
    }
    fprintf(outputStream, "\n");

    //Print set2
    fprintf(outputStream, "Set 3 elements %d", set3Counter);
    for (int i = 0; i < set3Counter; i++) {
        fprintf(outputStream, " %d", set3.at(i + 1) );
    }
    fprintf(outputStream, "\n");

    //Print set4
    fprintf(outputStream, "Set 4 nodes 1 %d\n", supportNode);
    fprintf(outputStream, "Set 5 nodes 1 %d\n", loadNode);

    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", loadNode);
    fprintf(outputStream, "#REACTION number %d dof 2 unknown d\n", supportNode);
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", notchOneBottomNode);
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", notchOneTopNode);
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", notchTwoBottomNode);
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", notchTwoTopNode);
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", controlMidBottomNode);
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", controlMidTopNode);
    fprintf(outputStream, "#TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");
    return;
}

void
Grid::give3DKupferOutput(const std::string &fileName)
{
    //Output for Kupfer experiment (Biaxial compression with stress ratios.
    //Only mechanical model

    FILE *outputStream = converter::fopen_or_die(fileName, "w");

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords(3);
    oofem::FloatArray coordsOne, coordsTwo, line;
    int materialType = 1;
    oofem::IntArray nodes, location(2);
    oofem::IntArray crossSectionNodes;
    double radius, distanceOne, distanceTwo;

    int supportNodeOne = -1;
    int loadNodeOne = -1;

    int supportNodeTwo = -1;
    int loadNodeTwo = -1;


    //Determine the number of Delaunay nodes in the domain
    numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 ||  this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfNodes++;

            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);

            //Create support and load nodes
            if ( fabs(coords.at(1) - 0.05) < TOL && fabs(coords.at(2) - 0.0) < TOL && fabs(coords.at(3) - 0.005) < TOL ) {
                supportNodeOne = i + 1;
            } else if ( fabs(coords.at(1) - 0.05) < TOL && fabs(coords.at(2) - 0.1) < TOL && fabs(coords.at(3) - 0.005) < TOL )                  {
                loadNodeOne = i + 1;
            } else if ( fabs(coords.at(1) - 0.0) < TOL && fabs(coords.at(2) - 0.05) < TOL && fabs(coords.at(3) - 0.005) < TOL )                  {
                supportNodeTwo = i + 1;
            } else if ( fabs(coords.at(1) - 0.1) < TOL && fabs(coords.at(2) - 0.05) < TOL && fabs(coords.at(3) - 0.005) < TOL )                  {
                loadNodeTwo = i + 1;
            }
        }
    }
    //Determine the number of Delaunay lines in the domain
    numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
        this->giveDelaunayVertex(nodes.at(1) )->giveCoordinates(coordsOne);
        this->giveDelaunayVertex(nodes.at(2) )->giveCoordinates(coordsTwo);

        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 && ( this->giveDelaunayLine(i + 1) )->delaunayAreaCheck() == 1 ) {
            numberOfLines++;
        } else if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 && this->giveRegion(1)->modifyVoronoiCrossSection(i + 1) == 1 )        {
            numberOfLines++;
        }
    }
    fprintf(outputStream, "oofem.out\n");
    fprintf(outputStream, "Mechanical 3D model of Kupfer experiment\n");
    fprintf(outputStream, "NlDEIDynamic nsteps 1 dumpcoef 0. deltat 0.1 reduct 0.5 nmodules 2\n");
    fprintf(outputStream, "#Nonlinearstatic nmsteps 1 nsteps 1 contextOutputStep 2000 nmodules 2 lstype 4 smtype 8\n");
    fprintf(outputStream, "#nsteps 100 rtolv 1.e-3 stiffMode 2 maxiter 400 controllmode 0 steplength 1.e-5 minsteplength 1.e-5 hpcmode 2 hpc 2 %d 2 hpcw 1 -1 donotfixload lstype 4 smtype 8\n", loadNodeOne);
    fprintf(outputStream, "vtkxmllattice primvars 1 1 tstep_step 1000 domain_all cross 1 cellvars 4 46 60 90 111 cross 1\n");
    fprintf(outputStream, "gpexportmodule vars 2 46 139 tstep_step 1000\n");
    fprintf(outputStream, "domain 3dLattice\n");
    fprintf(outputStream, "OutputManager tstep_step 1000 dofman_output { %d %d %d %d }\n", loadNodeOne, supportNodeOne, loadNodeTwo, supportNodeTwo);
    fprintf(outputStream, "ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 5 nic 0 nltf 3 nset 4\n", numberOfNodes, numberOfLines);

    printf("start nodes\n");

    int firstFlag = 0;

    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
            //First corners then facets
            if ( fabs(coords.at(1) - 0.0) < giveTol() && fabs(coords.at(2) - 0.0) < giveTol() ) {
                fprintf(outputStream, "node %d coords 3 %e %e %e dofidmask 6 1 2 3 4 5 6 mastermask 6 %d %d 0 0 0 0 doftype 6 1 1 0 0 0 0\n", i + 1, coords.at(1), coords.at(2), coords.at(3), supportNodeTwo, supportNodeOne);
            } else if ( fabs(coords.at(1) - 0.0) < giveTol() && fabs(coords.at(2) - 0.1) < giveTol() )              {
                fprintf(outputStream, "node %d coords 3 %e %e %e dofidmask 6 1 2 3 4 5 6 mastermask 6 %d %d 0 0 0 0 doftype 6 1 1 0 0 0 0\n", i + 1, coords.at(1), coords.at(2), coords.at(3), supportNodeTwo, loadNodeOne);
            } else if ( fabs(coords.at(1) - 0.1) < giveTol() && fabs(coords.at(2) - 0.1) < giveTol() )              {
                fprintf(outputStream, "node %d coords 3 %e %e %e dofidmask 6 1 2 3 4 5 6 mastermask 6 %d %d 0 0 0 0 doftype 6 1 1 0 0 0 0\n", i + 1, coords.at(1), coords.at(2), coords.at(3), loadNodeTwo, loadNodeOne);
            } else if ( fabs(coords.at(1) - 0.1) < giveTol() && fabs(coords.at(2) - 0.0) < giveTol() )              {
                fprintf(outputStream, "node %d coords 3 %e %e %e dofidmask 6 1 2 3 4 5 6 mastermask 6 %d %d 0 0 0 0 doftype 6 1 1 0 0 0 0\n", i + 1, coords.at(1), coords.at(2), coords.at(3), loadNodeTwo, supportNodeOne);
            } else if ( fabs(coords.at(2) - 0.0) < giveTol() && i + 1 != supportNodeOne )           {
                fprintf(outputStream, "node %d coords 3 %e %e %e dofidmask 6 1 2 3 4 5 6 mastermask 6 0 %d 0 0 0 0 doftype 6 0 1 0 0 0 0\n", i + 1, coords.at(1), coords.at(2), coords.at(3), supportNodeOne);
            } else if ( fabs(coords.at(2) - 0.1) < giveTol() && i + 1 != loadNodeOne )            {
                fprintf(outputStream, "node %d coords 3 %e %e %e dofidmask 6 1 2 3 4 5 6  mastermask 6 0 %d 0 0 0 0 doftype 6 0 1 0 0 0 0\n", i + 1, coords.at(1), coords.at(2), coords.at(3), loadNodeOne);
            } else if ( fabs(coords.at(1) - 0.0) < giveTol() && i + 1 != supportNodeTwo )            {
                fprintf(outputStream, "node %d coords 3 %e %e %e dofidmask 6 1 2 3 4 5 6  mastermask 6 %d 0 0 0 0 0 doftype 6 1 0 0 0 0 0\n", i + 1, coords.at(1), coords.at(2), coords.at(3), supportNodeTwo);
            } else if ( fabs(coords.at(1) - 0.1) < giveTol() && i + 1 != loadNodeTwo )            {
                fprintf(outputStream, "node %d coords 3 %e %e %e dofidmask 6 1 2 3 4 5 6  mastermask 6 %d 0 0 0 0 0 doftype 6 1 0 0 0 0 0\n", i + 1, coords.at(1), coords.at(2), coords.at(3), loadNodeTwo);
            } else   {
                fprintf(outputStream, "node %d coords 3 %e %e %e\n", i + 1, coords.at(1), coords.at(2), coords.at(3) );
            }
        }
    }


    printf("finished nodes\n");

    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        //and cross the middle

        if ( ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) && ( this->giveDelaunayLine(i + 1) )->delaunayAreaCheck() == 1 ) { //Elements are inside
            this->giveDelaunayLine(i + 1)->giveLocalVertices(nodes);

            //Only one material at the moment
            materialType = 1;
            this->giveDelaunayLine(i + 1)->updateMaterial(1);
            this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

            fprintf(outputStream, "lattice3D %d nodes 2 %d %d crossSect %d mat %d polycoords %d", i + 1, nodes.at(1), nodes.at(2), materialType, materialType, 3 * crossSectionNodes.giveSize() );

            for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                this->giveVoronoiVertex(crossSectionNodes.at(m + 1) )->giveCoordinates(coords);
                fprintf(outputStream, " %e %e %e", coords.at(1), coords.at(2), coords.at(3) );
            }
            fprintf(outputStream, "\n");
        }
    }

    fprintf(outputStream, "latticecs 1 material 1\n");
    fprintf(outputStream, "latticeplastdam 1 d 2500 talpha 0. calpha 0. e 45.91e9 a1 0.297 a2 1.e-12 ft 2.35e6 fc 30.e6 \
wf 20.e-6 ahard 1.e-3 angle1 0.5 flow 0.25 randvars 2 806 807 randgen 2 3 3\n");
    fprintf(outputStream, "#latticelinearelastic 1 d 0 talpha 0. calpha 0. e 50.46e9 a1 0.215\n");
    fprintf(outputStream, "BoundaryCondition 1 loadTimeFunction 1 dofs 2 2 3 values 2 0 0 set 1\n");
    fprintf(outputStream, "BoundaryCondition 2 loadTimeFunction 1 dofs 1 3 values 1 0 set 2\n");
    fprintf(outputStream, "BoundaryCondition 3 loadTimeFunction 1 dofs 2 1 3 values 2 0 0 set 3\n");
    fprintf(outputStream, "BoundaryCondition 4 loadTimeFunction 1 dofs 1 3 values 1 0 set 4\n");
    fprintf(outputStream, "BoundaryCondition 5 loadTimeFunction 2 dofs 1 2 values 1 -1.e-3 set 2\n");
    fprintf(outputStream, "#NodalLoad 5 loadTimeFunction 2 dofs 6 1 2 3 4 5 6 Components 6 0. 1. 0. 0. 0. 0. set 2\n");
    fprintf(outputStream, "#NodalLoad 6 loadTimeFunction 2 dofs 6 1 2 3 4 5 6 Components 6 0. 0. 0. 0. 0. 0. set 4\n");
    fprintf(outputStream, "ConstantFunction 1 f(t) 1.\n");
    fprintf(outputStream, "PiecewiseLinFunction 2 t 2 0 0.1 f(t) 2  0 1.\n");
    fprintf(outputStream, "InterpolatingFunction 3 name random.dat dim 3\n");

    //Print set4
    fprintf(outputStream, "Set 1 nodes 1 %d\n", supportNodeOne);
    fprintf(outputStream, "Set 2 nodes 1 %d\n", loadNodeOne);
    fprintf(outputStream, "Set 3 nodes 1 %d\n", supportNodeTwo);
    fprintf(outputStream, "Set 4 nodes 1 %d\n", loadNodeTwo);

    fprintf(outputStream, "#%%BEGIN_CHECK%%\n");
    fprintf(outputStream, "#NODE number %d dof 2 unknown d\n", loadNodeOne);
    fprintf(outputStream, "#REACTION number %d dof 2 unknown d\n", supportNodeOne);
    fprintf(outputStream, "#NODE number %d dof 1 unknown d\n", loadNodeTwo);
    fprintf(outputStream, "#REACTION number %d dof 1 unknown d\n", supportNodeTwo);
    fprintf(outputStream, "#TIME\n");
    fprintf(outputStream, "#%%END_CHECK%%\n");
    return;
}

void
Grid::give3DImranOutput(const std::string &fileName)
{
    return;
}


void
Grid::give3DNotchOutput(const std::string &fileName)
{
    return;
}



void
Grid::givePOVOutput(const std::string &fileName)
{
    //    char fileName1 [ MAX_FILENAME_LENGTH + 10 ];

    this->giveVoronoiPOVOutput(fileName);

    this->giveDelaunayPOVOutput(fileName);
}


void
Grid::giveVoronoiPOVOutput(const std::string &fileName)
{
    //    FILE* outputStream = converter::fopen_or_die(fileName, "w");

    const std::string fileName1 = fileName + ".vor.line.pov";
    FILE *outputStream1 = converter::fopen_or_die(fileName1, "w");

    const std::string fileName2 = fileName + ".vor.cross.pov";
    FILE *outputStream2 = converter::fopen_or_die(fileName2, "w");


    oofem::FloatArray coord1, coord2;
    oofem::IntArray nodes, lines;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 || giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 ) {
            giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
            giveVoronoiVertex(nodes.at(1) )->giveCoordinates(coord1);
            giveVoronoiVertex(nodes.at(2) )->giveCoordinates(coord2);
            fprintf(outputStream1, "cylinder{<%e,%e,%e>,<%e,%e,%e>,r}\n", coord1.at(1), coord1.at(2), coord1.at(3), coord2.at(1), coord2.at(2), coord2.at(3) );
            giveVoronoiLine(i + 1)->giveCrossSectionElements(lines);
            for ( int m = 0; m < lines.giveSize(); m++ ) {
                giveDelaunayLine(lines.at(m + 1) )->giveLocalVertices(nodes);
                giveDelaunayVertex(nodes.at(1) )->giveCoordinates(coord1);
                giveDelaunayVertex(nodes.at(2) )->giveCoordinates(coord2);
                fprintf(outputStream2, "cylinder{<%e,%e,%e>,<%e,%e,%e>,r}\n", coord1.at(1), coord1.at(2), coord1.at(3), coord2.at(1), coord2.at(2), coord2.at(3) );
            }
        }
    }
}


void
Grid::giveDelaunayPOVOutput(const std::string &fileName)
{
    const std::string fileName1 = fileName + ".del.line.pov";
    FILE *outputStream1 = converter::fopen_or_die(fileName1, "w");

    const std::string fileName2 = fileName + ".del.cross.pov";
    FILE *outputStream2 = converter::fopen_or_die(fileName2, "w");

    const std::string fileName3 = fileName + ".del.crack.pov";
    FILE *outputStream3 = converter::fopen_or_die(fileName2, "w");


    oofem::FloatArray coord1, coord2;
    oofem::IntArray nodes, lines;
    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        if ( giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || giveDelaunayLine(i + 1)->giveOutsideFlag() == 2 ) {
            giveDelaunayLine(i + 1)->giveLocalVertices(nodes);
            giveDelaunayVertex(nodes.at(1) )->giveCoordinates(coord1);
            giveDelaunayVertex(nodes.at(2) )->giveCoordinates(coord2);
            fprintf(outputStream1, "cylinder{<%e,%e,%e>,<%e,%e,%e>,r}\n", coord1.at(1), coord1.at(2), coord1.at(3), coord2.at(1), coord2.at(2), coord2.at(3) );
            giveDelaunayLine(i + 1)->giveCrossSectionElements(lines);
            for ( int m = 0; m < lines.giveSize(); m++ ) {
                giveVoronoiLine(lines.at(m + 1) )->giveLocalVertices(nodes);
                giveVoronoiVertex(nodes.at(1) )->giveCoordinates(coord1);
                giveVoronoiVertex(nodes.at(2) )->giveCoordinates(coord2);
                fprintf(outputStream2, "cylinder{<%e,%e,%e>,<%e,%e,%e>,r}\n", coord1.at(1), coord1.at(2), coord1.at(3), coord2.at(1), coord2.at(2), coord2.at(3) );
            }
            giveDelaunayLine(i + 1)->giveCrossSectionVertices(nodes);
            //Not working
            fprintf(outputStream3, "polygon{%d", nodes.giveSize() + 1);
            for ( int k = 0; k < nodes.giveSize(); k++ ) {
                giveVoronoiVertex(nodes.at(k + 1) )->giveCoordinates(coord1);
                fprintf(outputStream3, ",<%e,%e,%e>", coord1.at(1), coord1.at(2), coord1.at(3) );
            }
            //Close it
            giveVoronoiVertex(nodes.at(1) )->giveCoordinates(coord1);
            fprintf(outputStream3, ",<%e,%e,%e>", coord1.at(1), coord1.at(2), coord1.at(3) );
            fprintf(outputStream3, "}\n");
        }
    }
}


void
Grid::giveVoronoiCellVTKOutput(FILE *outputStream)
{
    //Here should go the code for the Voronoi nodes

    oofem::FloatArray coords;
    oofem::IntArray nodes;
    oofem::IntArray crossSectionNodes;
    Vertex *delaunayVertex;
    Line *delaunayLine;

    FILE *lengthFile;

    //Cells for Delaunay nodes in the specimen
    oofem::IntArray localDelaunayLines;
    int numberOfCells = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveLocalLines(localDelaunayLines);
            numberOfCells += localDelaunayLines.giveSize();
        }
    }

    //Nodes for Voronoi nodes in the specimen
    int numberOfNodes = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 || this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfNodes++;
        }
    }

    //Start with VTK output
    fprintf(outputStream, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(outputStream, "<UnstructuredGrid>\n");
    fprintf(outputStream, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", numberOfNodes, numberOfCells);
    // export nodes in region as vtk vertices

    fprintf(outputStream, "<Points>\n");
    fprintf(outputStream, "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");

    int nodeCounter = 0;
    oofem::IntArray nodeConverter(this->giveNumberOfVoronoiVertices() );
    for ( int i = 0; i < this->giveNumberOfVoronoiVertices(); i++ ) {
        if ( this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 0 || this->giveVoronoiVertex(i + 1)->giveOutsideFlag() == 2 ) {
            nodeCounter++;
            nodeConverter.at(i + 1) = nodeCounter;
            this->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
            for ( int i = 0; i < 3; i++ ) {
                fprintf(outputStream, "%e ", coords.at(i + 1) );
            }
            fprintf(outputStream, "\n");
        }
    }

    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "</Points>\n");

    fprintf(outputStream, "<Cells>\n");
    fprintf(outputStream, " <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");

    oofem::IntArray crossSectionVertices;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveLocalLines(localDelaunayLines);
            for ( int m = 0; m < localDelaunayLines.giveSize(); m++ ) {
                this->giveDelaunayLine(localDelaunayLines.at(m + 1) )->giveCrossSectionVertices(crossSectionVertices);
                for ( int k = 0; k < crossSectionVertices.giveSize(); k++ ) {
                    fprintf(outputStream, "%d ", nodeConverter.at(crossSectionVertices.at(k + 1) ) - 1);
                }
                fprintf(outputStream, "\n");
            }
        }
    }
    fprintf(outputStream, "</DataArray>\n");

    int offsetCounter = 0;
    fprintf(outputStream, " <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        if ( this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayVertex(i + 1)->giveLocalLines(localDelaunayLines);
            for ( int m = 0; m < localDelaunayLines.giveSize(); m++ ) {
                this->giveDelaunayLine(localDelaunayLines.at(m + 1) )->giveCrossSectionVertices(crossSectionVertices);
                offsetCounter += crossSectionVertices.giveSize();
                fprintf(outputStream, "%d ", offsetCounter);
            }
        }
    }
    fprintf(outputStream, "\n");
    fprintf(outputStream, "</DataArray>\n");

    fprintf(outputStream, " <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for ( int i = 0; i < numberOfCells; i++ ) {
        fprintf(outputStream, "7 ");
    }
    fprintf(outputStream, "\n");
    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "</Cells>\n");
    fprintf(outputStream, "</Piece>\n");
    fprintf(outputStream, "</UnstructuredGrid>\n</VTKFile>");

    return;
}


void
Grid::giveTetraElementVTKOutput(FILE *outputStream)
{
    printf("In Tetra output\n");

    oofem::FloatArray coords;
    oofem::IntArray nodes;
    oofem::IntArray crossSectionNodes;
    Vertex *delaunayVertex;
    Line *delaunayLine;

    FILE *lengthFile;

    //Each Voronoi node corresponds to a Delaunay Cell.
    //Check if one of the Delaunay nodes corresponding to the cell is inside.
    //If yes, this is a Delaunay cell that should be printed.
    int numberOfTetras = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayTetras(); i++ ) {
        if ( this->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 2 ) {
            numberOfTetras++;
        }
    }


    oofem::IntArray localTetras;
    int numberOfVertices = 0;
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        this->giveDelaunayVertex(i + 1)->giveLocalTetras(localTetras);
        for ( int k = 0; k < localTetras.giveSize(); k++) {
            if ( this->giveDelaunayTetra(localTetras.at(k + 1) )->giveOutsideFlag() == 0 || this->giveDelaunayTetra(localTetras.at(k + 1) )->giveOutsideFlag() == 2 ) {
                numberOfVertices++;
                break;
            }
        }
    }

    //Start with VTK output
    fprintf(outputStream, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(outputStream, "<UnstructuredGrid>\n");
    fprintf(outputStream, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", numberOfVertices, numberOfTetras);
    // export nodes in region as vtk vertices

    fprintf(outputStream, "<Points>\n");
    fprintf(outputStream, "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");

    int vertexCounter = 0;
    oofem::IntArray vertexConverter(this->giveNumberOfDelaunayVertices() );
    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        this->giveDelaunayVertex(i + 1)->giveLocalTetras(localTetras);
        for ( int k = 0; k < localTetras.giveSize(); k++) {
            if ( this->giveDelaunayTetra(localTetras.at(k + 1) )->giveOutsideFlag() == 0 || this->giveDelaunayTetra(localTetras.at(k + 1) )->giveOutsideFlag() == 2 ) {
                vertexCounter++;
                vertexConverter.at(i + 1) = vertexCounter;
                this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
                for ( int m = 0; m < 3; m++ ) {
                    fprintf(outputStream, "%e ", coords.at(m + 1) );
                }
                fprintf(outputStream, "\n");
                break;
            }
        }
    }

    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "</Points>\n");

    fprintf(outputStream, "<Cells>\n");
    fprintf(outputStream, " <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");

    oofem::IntArray localDelaunayVertices;
    for ( int i = 0; i < this->giveNumberOfDelaunayTetras(); i++ ) {
        if ( this->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayTetra(i + 1)->giveLocalVertices(localDelaunayVertices);
            for ( int m = 0; m < localDelaunayVertices.giveSize(); m++ ) {
                fprintf(outputStream, "%d ", vertexConverter.at(localDelaunayVertices.at(m + 1) ) - 1);
            }
        }
    }
    fprintf(outputStream, "</DataArray>\n");

    int offsetCounter = 0;
    fprintf(outputStream, " <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    for ( int i = 0; i < this->giveNumberOfDelaunayTetras(); i++ ) {
        if ( this->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 2 ) {
            this->giveDelaunayTetra(i + 1)->giveLocalVertices(localDelaunayVertices);
            offsetCounter += localDelaunayVertices.giveSize();
            fprintf(outputStream, "%d ", offsetCounter);
        }
    }

    fprintf(outputStream, "</DataArray>\n");

    fprintf(outputStream, " <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for ( int i = 0; i < this->giveNumberOfDelaunayTetras(); i++ ) {
        if ( this->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayTetra(i + 1)->giveOutsideFlag() == 2 ) {
            fprintf(outputStream, "10 ");
        }
    }
    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "</Cells>\n");
    fprintf(outputStream, "</Piece>\n");
    fprintf(outputStream, "</UnstructuredGrid>\n</VTKFile>");

    return;
}


void
Grid::giveDelaunayElementVTKOutput(FILE *outputStream)
{
    //Here should go the code for the Voronoi nodes

    //This should be extended to write out also material properties.
    //In this way, it would be possible to show the mesh with projected properties without having to run any analyses.

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    oofem::IntArray nodes;
    oofem::IntArray crossSectionNodes;
    Vertex *delaunayVertex;
    Line *delaunayLine;

    FILE *lengthFile;


    fprintf(outputStream, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(outputStream, "<UnstructuredGrid>\n");

    numberOfNodes = 0;

    //This needs to be changed for periodic cells. If elements cross the boundary then the nodes should be included as well.
    //Thus, it needs to be written differently to allow for the case of crossing nodes
    //Also, the structure should be changed so that info on material or cross-sections could be written which originate from underlying meso-structure.


    for ( int inode = 0; inode < this->giveNumberOfDelaunayVertices(); inode++ ) {
        if ( this->giveDelaunayVertex(inode + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(inode + 1)->giveOutsideFlag() == 2 ) {
            numberOfNodes++;
        }
    }

    numberOfLines = 0;
    for ( int iline = 0; iline < this->giveNumberOfDelaunayLines(); iline++ ) {
        if ( this->giveDelaunayLine(iline + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(iline + 1)->giveOutsideFlag() == 3 ) {
            numberOfLines++;
        }
    }

    printf("numberOfLines = %d\n", numberOfLines);

    fprintf(outputStream, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", numberOfNodes, numberOfLines);
    // export nodes in region as vtk vertices

    fprintf(outputStream, "<Points>\n <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"> \n");

    oofem::IntArray nodeConverter(this->giveNumberOfDelaunayVertices() );
    int newNodeNumber = 0;
    for ( int inode = 0; inode < this->giveNumberOfDelaunayVertices(); inode++ ) {
        if ( this->giveDelaunayVertex(inode + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(inode + 1)->giveOutsideFlag() == 2 ) {
            newNodeNumber++;
            nodeConverter.at(inode + 1) = newNodeNumber;
            delaunayVertex = this->giveDelaunayVertex(inode + 1);
            delaunayVertex->giveCoordinates(coords);
            for ( int i = 1; i <= 3; i++ ) {
                fprintf(outputStream, "%e ", coords.at(i) );
            }
            fprintf(outputStream, "\n");
        }
    }

    fprintf(outputStream, "</DataArray>\n</Points>\n");

    fprintf(outputStream, "<Cells>\n");
    // output the connectivity data
    fprintf(outputStream, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    //

    for ( int i = 0; i <  this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 0 || this->giveDelaunayLine(i + 1)->giveOutsideFlag() == 3 ) {
            delaunayLine = this->giveDelaunayLine(i + 1);
            delaunayLine->giveLocalVertices(nodes);
            delaunayLine->giveCrossSectionVertices(crossSectionNodes);
            fprintf(outputStream, "%d %d ", nodeConverter.at(nodes.at(1) ) - 1, nodeConverter.at(nodes.at(2) ) - 1);
        }
    }
    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");

    for ( int i = 0; i < numberOfLines; i++ ) {
        fprintf(outputStream, "%d ", 2 * ( i + 1 ) );
    }

    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"> ");
    for ( int cell = 0; cell < numberOfLines; cell++ ) {
        fprintf(outputStream, "%d ", 3);
    }
    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "</Cells>\n");
    fprintf(outputStream, "</Piece>\n");
    fprintf(outputStream, "</UnstructuredGrid>\n</VTKFile>");

    return;
}



void
Grid::giveVoronoiCrossSectionVTKOutput(FILE *outputStream)
{
}





void
Grid::giveVoronoiElementVTKOutput(FILE *outputStream)
{
    //Here should go the code for the Voronoi nodes

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    oofem::IntArray nodes;
    oofem::IntArray crossSectionNodes;
    Vertex *voronoiVertex;
    Line *voronoiLine;

    FILE *lengthFile;


    fprintf(outputStream, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(outputStream, "<UnstructuredGrid>\n");

    numberOfNodes = 0;

    for ( int inode = 0; inode < this->giveNumberOfVoronoiVertices(); inode++ ) {
        if ( this->giveVoronoiVertex(inode + 1)->giveOutsideFlag() == 0 || this->giveVoronoiVertex(inode + 1)->giveOutsideFlag() == 1 || this->giveVoronoiVertex(inode + 1)->giveOutsideFlag() == 2 ) {
            numberOfNodes++;
        }
    }

    numberOfLines = 0;
    for ( int iline = 0; iline < this->giveNumberOfVoronoiLines(); iline++ ) {
        if ( this->giveVoronoiLine(iline + 1)->giveOutsideFlag() == 0 || this->giveVoronoiLine(iline + 1)->giveOutsideFlag() == 2 || this->giveVoronoiLine(iline + 1)->giveOutsideFlag() == 3 ) {
            numberOfLines++;
        }
    }


    fprintf(outputStream, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", numberOfNodes, numberOfLines);
    // export nodes in region as vtk vertices

    fprintf(outputStream, "<Points>\n <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");

    int nodeCounter = 0;
    oofem::IntArray nodeConverter(this->giveNumberOfVoronoiVertices() );
    for ( int inode = 0; inode < this->giveNumberOfVoronoiVertices(); inode++ ) {
        if ( this->giveVoronoiVertex(inode + 1)->giveOutsideFlag() == 0 || this->giveVoronoiVertex(inode + 1)->giveOutsideFlag() == 1 || this->giveVoronoiVertex(inode + 1)->giveOutsideFlag() == 2 ) {
            nodeCounter++;
            nodeConverter.at(inode + 1) = nodeCounter;
            voronoiVertex = this->giveVoronoiVertex(inode + 1);
            voronoiVertex->giveCoordinates(coords);
            for ( int i = 1; i <= 3; i++ ) {
                fprintf(outputStream, "%e ", coords.at(i) );
            }
            fprintf(outputStream, "\n");
        }
    }

    fprintf(outputStream, "</DataArray>\n</Points>\n");

    fprintf(outputStream, "<Cells>\n");
    // output the connectivity data
    fprintf(outputStream, " <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    //

    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        if ( this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 0 ||  this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 2 || this->giveVoronoiLine(i + 1)->giveOutsideFlag() == 3 ) {
            voronoiLine = this->giveVoronoiLine(i + 1);
            voronoiLine->giveLocalVertices(nodes);
            voronoiLine->giveCrossSectionVertices(crossSectionNodes);
            fprintf(outputStream, "%d %d ",  nodeConverter.at(nodes.at(1) ) - 1, nodeConverter.at(nodes.at(2) ) - 1);
        }
    }
    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, " <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"> ");

    for ( int i = 0; i <  numberOfLines; i++ ) {
        fprintf(outputStream, "%d ", 2 * ( i + 1 ) );
    }

    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, " <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"> ");
    for ( int cell = 0; cell <  numberOfLines; cell++ ) {
        fprintf(outputStream, "%d ", 3);
    }
    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "</Cells>\n");
    fprintf(outputStream, "</Piece>\n");
    fprintf(outputStream, "</UnstructuredGrid>\n</VTKFile>");

    return;
}



void
Grid::giveDelaunayCrossSectionVTKOutput(FILE *outputStream)
{
}



void
Grid::createRankTable(oofem::IntArray &rankVector, oofem::FloatArray &randomNumbers) {
    int size = randomNumbers.giveSize();
    rankVector.resize(size);

    double *random = new double [ size ];
    for ( int i = 0; i < size; i++ ) {
        random [ i ] = randomNumbers.at(i + 1);
    }

    int *conversion = new int [ size ];
    int *rank = new int [ size ];
    double *sortedRandom = new double [ size ];

    //Test the indexx-function

    //Sort the field using indexx
    conversion--;
    random--;
    rank--;

    indexx(size, random, conversion);

    for ( int i = 1; i <= size; i++ ) {
        rank [ conversion [ i ] ] = i;
    }

    rank++;
    conversion++;
    random++;

    for ( int i = 0; i < size; i++ ) {
        rankVector.at(i + 1) = rank [ i ];
    }

    delete[] random;
    delete[] conversion;
    delete[] rank;
    delete[] sortedRandom;

    return;
}


void
Grid::sortRandomNumbers(oofem::FloatArray &sortedRandomNumbers, oofem::FloatArray &randomNumbers) {
    int size = randomNumbers.giveSize();
    sortedRandomNumbers.resize(size);
    double *random = new double [ size ];
    for ( int i = 0; i < size; i++ ) {
        random [ i ] = randomNumbers.at(i + 1);
    }

    int *conversion = new int [ size ];
    double *sortedRandom = new double [ size ];

    //Test the indexx-function

    //Sort the field using indexx
    conversion--;
    random--;
    sortedRandom--;

    indexx(size, random, conversion);
    for ( int i = 1; i <= size; i++ ) {
        sortedRandom [ i ] = random [ conversion [ i ] ];
    }

    random++;
    sortedRandom++;
    conversion++;

    for ( int i = 0; i < size; i++ ) {
        sortedRandomNumbers.at(i + 1) = sortedRandom [ i ];
    }

    printf("Sorting done\n");

    delete[] sortedRandom;
    delete[] conversion;
    delete[] random;

    return;
}




double Grid::ran1(long *idum)
{
    long k;
    static long iy = 0;
    static long iv[ NTAB ];
    float temp;

    if ( * idum <= 0 || !iy ) {
        if ( -( * idum ) < 1 ) {
            * idum = 1;
        } else {
            * idum = -( * idum );
        }

        for ( int j = NTAB + 7; j >= 0; j-- ) {
            k = ( * idum ) / IQ;
            * idum = IA * ( * idum - k * IQ ) - IR * k;
            if ( * idum < 0 ) {
                * idum += IM;
            }

            if ( j < NTAB ) {
                iv [ j ] = * idum;
            }
        }

        iy = iv [ 0 ];
    }

    k = ( * idum ) / IQ;
    * idum = IA * ( * idum - k * IQ ) - IR * k;
    if ( * idum < 0 ) {
        * idum += IM;
    }

    int j = iy / NDIV;
    iy = iv [ j ];
    iv [ j ] = * idum;
    if ( ( temp = AM * iy ) > RNMX ) {
        return RNMX;
    } else {
        return temp;
    }
}


double Grid::normalCdfInverse(double cdf, double a, double b)
{
    double x;
    double x2;
    if ( cdf < 0.0 || 1.0 < cdf ) {
        printf("Error: in normalCdfInverse. Values outside range 0-1.\n");
        exit;
    }

    x2 = normal01CdfInverse(cdf);
    x = exp(a + b * x2);

    return x;
}

double Grid::normal01CdfInverse(double p)
{
    double a[ 8 ] = {
        3.3871328727963666080,     1.3314166789178437745e+2,
        1.9715909503065514427e+3,  1.3731693765509461125e+4,
        4.5921953931549871457e+4,  6.7265770927008700853e+4,
        3.3430575583588128105e+4,  2.5090809287301226727e+3
    };
    double b[ 8 ] = {
        1.0,                       4.2313330701600911252e+1,
        6.8718700749205790830e+2,  5.3941960214247511077e+3,
        2.1213794301586595867e+4,  3.9307895800092710610e+4,
        2.8729085735721942674e+4,  5.2264952788528545610e+3
    };
    double c[ 8 ] = {
        1.42343711074968357734,     4.63033784615654529590,
        5.76949722146069140550,     3.64784832476320460504,
        1.27045825245236838258,     2.41780725177450611770e-1,
        2.27238449892691845833e-2,  7.74545014278341407640e-4
    };
    double const1 = 0.180625;
    double const2 = 1.6;
    double d[ 8 ] = {
        1.0,                        2.05319162663775882187,
        1.67638483018380384940,     6.89767334985100004550e-1,
        1.48103976427480074590e-1,  1.51986665636164571966e-2,
        5.47593808499534494600e-4,  1.05075007164441684324e-9
    };
    double e[ 8 ] = {
        6.65790464350110377720,     5.46378491116411436990,
        1.78482653991729133580,     2.96560571828504891230e-1,
        2.65321895265761230930e-2,  1.24266094738807843860e-3,
        2.71155556874348757815e-5,  2.01033439929228813265e-7
    };
    double f[ 8 ] = {
        1.0,                        5.99832206555887937690e-1,
        1.36929880922735805310e-1,  1.48753612908506148525e-2,
        7.86869131145613259100e-4,  1.84631831751005468180e-5,
        1.42151175831644588870e-7,  2.04426310338993978564e-15
    };
    double q;
    double r;
    double split1 = 0.425;
    double split2 = 5.0;
    double value;

    if ( p <= 0.0 ) {
        value = -HUGE_VAL;
        return value;
    }

    if ( 1.0 <= p ) {
        value = HUGE_VAL;
        return value;
    }

    q = p - 0.5;

    if ( fabs(q) <= split1 ) {
        r = const1 - q * q;
        value = q * dpolyValue(8, a, r) / dpolyValue(8, b, r);
    } else {
        if ( q < 0.0 ) {
            r = p;
        } else {
            r = 1.0 - p;
        }

        if ( r <= 0.0 ) {
            value = -1.0;
            printf("LocalGaussianRandomGenerator :: normal01CdfInverse - r < 0.0!");
            exit;
        }

        r = sqrt(-log(r) );
        if ( r <= split2 ) {
            r = r - const2;
            value = dpolyValue(8, c, r) / dpolyValue(8, d, r);
        } else {
            r = r - split2;
            value = dpolyValue(8, e, r) / dpolyValue(8, f, r);
        }

        if ( q < 0.0 ) {
            value = -value;
        }
    }

    return value;
}


double Grid::dpolyValue(int n, double a[], double x)
{
    int i;
    double value;
    value = 0.0;
    for ( i = n - 1; 0 <= i; i-- ) {
        value = value * x + a [ i ];
    }

    return value;
}


#define SWAP(a, b) itemp = ( a ); ( a ) = ( b ); ( b ) = itemp;
#define M 7
#define NSTACK 50

void Grid::indexx(int n, double arr[], int indx[])
{
    int i, indxt, ir = n, itemp, j, k, l = 1;
    int jstack = 0, *istack;
    double a;

    istack = ivector(1, NSTACK);

    for ( j = 1; j <= n; j++ ) {
        indx [ j ] = j;
    }
    for ( ; ; ) {
        if ( ir - l < M ) {
            for ( j = l + 1; j <= ir; j++ ) {
                indxt = indx [ j ];
                a = arr [ indxt ];
                for ( i = j - 1; i >= l; i-- ) {
                    if ( arr [ indx [ i ] ] <= a ) {
                        break;
                    }
                    indx [ i + 1 ] = indx [ i ];
                }
                indx [ i + 1 ] = indxt;
            }
            if ( jstack == 0 ) {
                break;
            }
            ir = istack [ jstack-- ];
            l = istack [ jstack-- ];
        } else {
            k = ( l + ir ) >> 1;
            SWAP(indx [ k ], indx [ l + 1 ]);
            if ( arr [ indx [ l ] ] > arr [ indx [ ir ] ] ) {
                SWAP(indx [ l ], indx [ ir ])
            }
            if ( arr [ indx [ l + 1 ] ] > arr [ indx [ ir ] ] ) {
                SWAP(indx [ l + 1 ], indx [ ir ])
            }
            if ( arr [ indx [ l ] ] > arr [ indx [ l + 1 ] ] ) {
                SWAP(indx [ l ], indx [ l + 1 ])
            }
            i = l + 1;
            j = ir;
            indxt = indx [ l + 1 ];
            a = arr [ indxt ];
            for ( ; ; ) {
                do {
                    i++;
                } while ( arr [ indx [ i ] ] < a );
                do {
                    j--;
                } while ( arr [ indx [ j ] ] > a );
                if ( j < i ) {
                    break;
                }
                SWAP(indx [ i ], indx [ j ])
            }
            indx [ l + 1 ] = indx [ j ];
            indx [ j ] = indxt;
            jstack += 2;
            if ( jstack > NSTACK ) {
                printf("Error in indexx: NSTACK too small in indexx.");
                std::exit(1);
            }
            if ( ir - i + 1 >= j - l ) {
                istack [ jstack ] = ir;
                istack [ jstack - 1 ] = i;
                ir = j - 1;
            } else {
                istack [ jstack ] = j - 1;
                istack [ jstack - 1 ] = l;
                l = i;
            }
        }
    }
    free_ivector(istack, 1, NSTACK);
}



#define NR_END 1
#define FREE_ARG char *
int *Grid::ivector(int nl, int nh)
{
    int *v;

    v = ( int * ) malloc( ( size_t ) ( ( nh - nl + 1 + NR_END ) * sizeof( int ) ) );
    if ( !v ) {
        printf("allocation failure in ivector()");
        std::exit(1);
    }
    return v - nl + NR_END;
}

void Grid::free_ivector(int *v, int nl, int nh)
{
    free( ( FREE_ARG ) ( v + nl - NR_END ) );
}



void Grid::giveVtkOutput2(const std::string &fileName, int nb_of_mt)
{
    FILE *outputStream;


    if ( periodicityFlag.at(1) != 1 &&
         periodicityFlag.at(2) != 1 &&
         periodicityFlag.at(3) != 1 ) {
        const std::string fname = fileName + ".voronoicell.vtu";
        if ( FILE *f = converter::fopen_or_die(fname, "w") ) {
            giveVoronoiCellVTKOutput(f);
            std::fclose(f);
        }
    }


    for (int imat = 1; imat <= nb_of_mt; ++imat) {
        const std::string fname =
            fileName + ".delaunayelement.mat" + std::to_string(imat) + ".vtu";

        if ( FILE *f = converter::fopen_or_die(fname, "w") ) {
            giveDelaunayElementVTKOutput2(f, imat);
            std::fclose(f);
        }
    }

    // Fibre beam elements
    const std::string fname2 = fileName + ".fibre.beamElement.vtu";
    if ( FILE *f = converter::fopen_or_die(fname2, "w") ) {
        giveBeamElementVTKOutput(f);
        std::fclose(f);
    }

    // Fibre link elements
    const std::string fname3 = fileName + ".fibre.linkElement.vtu";
    if ( FILE *f = converter::fopen_or_die(fname3, "w") ) {
        giveLinkElementVTKOutput(f);
        std::fclose(f);
    }

    // Voronoi elements
    const std::string fname4 = fileName + ".voronoielement.vtu";
    if ( FILE *f = converter::fopen_or_die(fname4, "w") ) {
        giveVoronoiElementVTKOutput(f);
        std::fclose(f);
    }


    return;
}



void Grid::giveVtkOutputTetra(const std::string &fileName, int nb_of_mt)
{
    //    FILE *outputStream;

    const std::string fname = fileName + ".tetraElement.vtu";
    if ( FILE *f = converter::fopen_or_die(fname, "w") ) {
        giveTetraElementVTKOutput(f);
        std::fclose(f);
    }

    return;
}



void
Grid::giveDelaunayElementVTKOutput2(FILE *outputStream, int nb_mtx)
{
    //Here should go the code for the Voronoi nodes

    //This should be extended to write out also material properties.
    //In this way, it would be possible to show the mesh with projected properties without having to run any analyses.

    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    oofem::IntArray nodes;
    oofem::IntArray crossSectionNodes;
    Vertex *delaunayVertex;
    Line *delaunayLine;

    FILE *lengthFile;


    fprintf(outputStream, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(outputStream, "<UnstructuredGrid>\n");

    numberOfNodes = 0;

    //we eventually keep all nodes so as to display lines which cross boundaries


    for ( int inode = 0; inode < this->giveNumberOfDelaunayVertices(); inode++ ) {
        //if ( this->giveDelaunayVertex(inode + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(inode + 1)->giveOutsideFlag() == 2 ) {
        numberOfNodes++;
        // }
    }

    numberOfLines = 0;
    for ( int iline = 0; iline < this->giveNumberOfDelaunayLines(); iline++ ) {
        if ( this->giveDelaunayLine(iline + 1)->giveOutsideFlag() != 1 ) {
            if ( this->giveDelaunayLine(iline + 1)->giveMaterial() == nb_mtx ) {
                numberOfLines++;
            }
        }
    }


    fprintf(outputStream, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", numberOfNodes, numberOfLines);
    // export nodes in region as vtk vertices

    fprintf(outputStream, "<Points>\n <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"> \n");

    // new numeratation of nodes
    oofem::IntArray nodeConverter(this->giveNumberOfDelaunayVertices() );
    int newNodeNumber = 0;
    for ( int inode = 0; inode < this->giveNumberOfDelaunayVertices(); inode++ ) {
        //if ( this->giveDelaunayVertex(inode + 1)->giveOutsideFlag() == 0 || this->giveDelaunayVertex(inode + 1)->giveOutsideFlag() == 2 ) {
        newNodeNumber++;
        nodeConverter.at(inode + 1) = newNodeNumber;
        delaunayVertex = this->giveDelaunayVertex(inode + 1);
        delaunayVertex->giveCoordinates(coords);
        for ( int i = 1; i <= 3; i++ ) {
            fprintf(outputStream, "%e ", coords.at(i) );
        }
        fprintf(outputStream, "\n");
        //}
    }

    fprintf(outputStream, "</DataArray>\n</Points>\n");

    fprintf(outputStream, "<Cells>\n");
    // output the connectivity data
    fprintf(outputStream, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    //

    for ( int i = 0; i <  this->giveNumberOfDelaunayLines(); i++ ) {
        if ( this->giveDelaunayLine(i + 1)->giveOutsideFlag() != 1 ) {
            if ( ( this->giveDelaunayLine(i + 1)->giveMaterial() ) == nb_mtx ) {
                delaunayLine = this->giveDelaunayLine(i + 1);
                delaunayLine->giveLocalVertices(nodes);
                delaunayLine->giveCrossSectionVertices(crossSectionNodes);
                fprintf(outputStream, "%d %d ", nodeConverter.at(nodes.at(1) ) - 1, nodeConverter.at(nodes.at(2) ) - 1);
            }
        }
    }
    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");

    for ( int i = 0; i < numberOfLines; i++ ) {
        fprintf(outputStream, "%d ", 2 * ( i + 1 ) );
    }

    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"> ");
    for ( int cell = 0; cell < numberOfLines; cell++ ) {
        fprintf(outputStream, "%d ", 3);
    }
    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "</Cells>\n");
    fprintf(outputStream, "</Piece>\n");
    fprintf(outputStream, "</UnstructuredGrid>\n</VTKFile>");

    return;
}

// visualise lattice beams and lattice links

void
Grid::giveBeamElementVTKOutput(FILE *outputStream)
{ // give VTK output for fibre (part I : only beam elements)
    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    oofem::IntArray nodes;
    oofem::IntArray crossSectionNodes;
    Vertex *delaunayVertex;
    Line *delaunayLine;


    fprintf(outputStream, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(outputStream, "<UnstructuredGrid>\n");

    numberOfNodes = this->giveNumberOfReinforcementNode();

    //This needs to be changed for periodic cells. If elements cross the boundary then the nodes should be included as well.
    //Thus, it needs to be written differently to allow for the case of crossing nodes
    //Also, the structure should be changed so that info on material or cross-sections could be written which originate from underlying meso-structure.

    numberOfLines = 0;
    for ( int iline = 0; iline < this->giveNumberOfLatticeBeams(); iline++ ) {
        if ( this->giveLatticeBeam(iline + 1)->giveOutsideFlag() != 1 ) {
            numberOfLines++;
        }
    }


    fprintf(outputStream, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", numberOfNodes, numberOfLines);
    // export nodes in region as vtk vertices

    fprintf(outputStream, "<Points>\n <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"> \n");

    // writing of nodes

    for ( int inode = 0; inode < this->giveNumberOfReinforcementNode(); inode++ ) {
        delaunayVertex = this->giveReinforcementNode(inode + 1);
        delaunayVertex->giveCoordinates(coords);
        for ( int i = 1; i <= 3; i++ ) {
            fprintf(outputStream, "%e ", coords.at(i) );
        }
        fprintf(outputStream, "\n");
    }


    fprintf(outputStream, "</DataArray>\n</Points>\n");

    fprintf(outputStream, "<Cells>\n");
    // output the connectivity data
    fprintf(outputStream, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    //

    printf("\n nb of lattice beam= %d \n", this->giveNumberOfLatticeBeams() );////

    for ( int i = 0; i <  this->giveNumberOfLatticeBeams(); i++ ) {
        if ( this->giveLatticeBeam(i + 1)->giveOutsideFlag() != 1 ) {
            delaunayLine = this->giveLatticeBeam(i + 1);
            delaunayLine->giveLocalVertices(nodes);
            fprintf(outputStream, "%d %d ",
                    +nodes.at(1) - 1,
                    +nodes.at(2)  - 1);
        }
    }
    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");

    for ( int i = 0; i < numberOfLines; i++ ) {
        fprintf(outputStream, "%d ", 2 * ( i + 1 ) );
    }

    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"> ");
    for ( int cell = 0; cell < numberOfLines; cell++ ) {
        fprintf(outputStream, "%d ", 3);
    }
    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "</Cells>\n");
    fprintf(outputStream, "</Piece>\n");
    fprintf(outputStream, "</UnstructuredGrid>\n</VTKFile>");

    return;
}

void
Grid::giveLinkElementVTKOutput(FILE *outputStream)
{ // give VTK output for fibre (part I : only beam elements)
    int numberOfNodes, numberOfLines;
    oofem::FloatArray coords;
    oofem::IntArray nodes;
    oofem::IntArray crossSectionNodes;
    Vertex *delaunayVertex;
    Line *delaunayLine;


    fprintf(outputStream, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(outputStream, "<UnstructuredGrid>\n");

    numberOfNodes = this->giveNumberOfDelaunayVertices() + this->giveNumberOfReinforcementNode();

    //This needs to be changed for periodic cells. If elements cross the boundary then the nodes should be included as well.
    //Thus, it needs to be written differently to allow for the case of crossing nodes
    //Also, the structure should be changed so that info on material or cross-sections could be written which originate from underlying meso-structure.

    numberOfLines = 0;
    for ( int iline = 0; iline < this->giveNumberOfLatticeLinks(); iline++ ) {
        if ( this->giveLatticeLink(iline + 1)->giveOutsideFlag() != 1 ) {
            numberOfLines++;
        }
    }


    fprintf(outputStream, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", numberOfNodes, numberOfLines);
    // export nodes in region as vtk vertices

    fprintf(outputStream, "<Points>\n <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"> \n");

    // writing of nodes


    for ( int inode = 0; inode < this->giveNumberOfDelaunayVertices(); inode++ ) {
        delaunayVertex = this->giveDelaunayVertex(inode + 1);
        delaunayVertex->giveCoordinates(coords);
        for ( int i = 1; i <= 3; i++ ) {
            fprintf(outputStream, "%e ", coords.at(i) );
        }
        fprintf(outputStream, "\n");
    }


    for ( int inode = 0; inode < this->giveNumberOfReinforcementNode(); inode++ ) {
        delaunayVertex = this->giveReinforcementNode(inode + 1);
        delaunayVertex->giveCoordinates(coords);
        for ( int i = 1; i <= 3; i++ ) {
            fprintf(outputStream, "%e ", coords.at(i) );
        }
        fprintf(outputStream, "\n");
    }


    fprintf(outputStream, "</DataArray>\n</Points>\n");

    fprintf(outputStream, "<Cells>\n");
    // output the connectivity data
    fprintf(outputStream, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    //

    printf("\n nb of lattice links = %d \n", this->giveNumberOfLatticeLinks() );////

    for ( int i = 0; i <  this->giveNumberOfLatticeLinks(); i++ ) {
        if ( this->giveLatticeLink(i + 1)->giveOutsideFlag() != 1 ) {
            delaunayLine = this->giveLatticeLink(i + 1);
            delaunayLine->giveLocalVertices(nodes);
            fprintf(outputStream, "%d %d ",
                    this->giveNumberOfDelaunayVertices() +  nodes.at(1) - 1,
                    nodes.at(2)  - 1);
        }
    }
    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");

    for ( int i = 0; i < numberOfLines; i++ ) {
        fprintf(outputStream, "%d ", 2 * ( i + 1 ) );
    }

    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"> ");
    for ( int cell = 0; cell < numberOfLines; cell++ ) {
        fprintf(outputStream, "%d ", 3);
    }
    fprintf(outputStream, "</DataArray>\n");
    fprintf(outputStream, "</Cells>\n");
    fprintf(outputStream, "</Piece>\n");
    fprintf(outputStream, "</UnstructuredGrid>\n</VTKFile>");

    return;
}

oofem::IntArray Grid::findDelaunayNodesWithinBox(oofem::FloatArray coord, double TOL)
{ // function created to allow other objects to use the localizer
    nodeContainerType nodeSet;
    delaunayLocalizer->giveAllNodesWithinBox(nodeSet, coord, TOL, 0);

    oofem::IntArray nodeSetNumbers(nodeSet.size() );
    for (int i = 1; i < 1 + nodeSetNumbers.giveSize(); i++) {
        nodeSetNumbers.at(i) = * nodeSet.begin();
        nodeSet.pop_front();
    }

    return nodeSetNumbers;
}


void Grid::giveTetrahedronBarycentres(oofem::FloatArray &centres, oofem::FloatArray &tetraCoords, oofem::FloatArray &pointCoords)
{
    double *tetra = new double [ 12 ];

    for (int i = 0; i < 12; i++) {
        tetra [ i ] = tetraCoords.at(i + 1);
    }

    double *p = new double [ 3 ];
    for (int i = 0; i < 3; i++) {
        p [ i ] = pointCoords.at(i + 1);
    }

    double *c = new double [ 3 ];
    c = tetrahedron_barycentric(tetra, p);

    for (int i = 0; i < 3; i++) {
        centres.at(i + 1) = c [ i ];
    }

    delete[] tetra;
    delete[] p;
    delete[] c;

    return;
}

//****************************************************************************80

double * Grid::tetrahedron_barycentric(double tetra[ 3 * 4 ], double p[ 3 ])

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_BARYCENTRIC returns the barycentric coordinates of a point.
//
//  Discussion:
//
//    The barycentric coordinates of a point P with respect to
//    a tetrahedron are a set of four values C(1:4), each associated
//    with a vertex of the tetrahedron.  The values must sum to 1.
//    If all the values are between 0 and 1, the point is contained
//    within the tetrahedron.
//
//    The barycentric coordinate of point X related to vertex A can be
//    interpreted as the ratio of the volume of the tetrahedron with
//    vertex A replaced by vertex X to the volume of the original
//    tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Input, double P[3], the point to be checked.
//
//    Output, double C[4], the barycentric coordinates of the point with
//    respect to the tetrahedron.
//
{
# define N 3
# define RHS_NUM 1

    double a[ N * ( N + RHS_NUM ) ];
    double *c;
    int info;
    //
    //  Set up the linear system
    //
    //    ( X2-X1  X3-X1  X4-X1 ) C1    X - X1
    //    ( Y2-Y1  Y3-Y1  Y4-Y1 ) C2  = Y - Y1
    //    ( Z2-Z1  Z3-Z1  Z4-Z1 ) C3    Z - Z1
    //
    //  which is satisfied by the barycentric coordinates.
    //

    a [ 0 + 0 * N ] = tetra [ 0 + 1 * 3 ] - tetra [ 0 + 0 * 3 ];
    a [ 1 + 0 * N ] = tetra [ 1 + 1 * 3 ] - tetra [ 1 + 0 * 3 ];
    a [ 2 + 0 * N ] = tetra [ 2 + 1 * 3 ] - tetra [ 2 + 0 * 3 ];

    a [ 0 + 1 * N ] = tetra [ 0 + 2 * 3 ] - tetra [ 0 + 0 * 3 ];
    a [ 1 + 1 * N ] = tetra [ 1 + 2 * 3 ] - tetra [ 1 + 0 * 3 ];
    a [ 2 + 1 * N ] = tetra [ 2 + 2 * 3 ] - tetra [ 2 + 0 * 3 ];

    a [ 0 + 2 * N ] = tetra [ 0 + 3 * 3 ] - tetra [ 0 + 0 * 3 ];
    a [ 1 + 2 * N ] = tetra [ 1 + 3 * 3 ] - tetra [ 1 + 0 * 3 ];
    a [ 2 + 2 * N ] = tetra [ 2 + 3 * 3 ] - tetra [ 2 + 0 * 3 ];

    a [ 0 + 3 * N ] = p [ 0 ]         - tetra [ 0 + 0 * 3 ];
    a [ 1 + 3 * N ] = p [ 1 ]         - tetra [ 1 + 0 * 3 ];
    a [ 2 + 3 * N ] = p [ 2 ]         - tetra [ 2 + 0 * 3 ];
    //
    //  Solve the linear system.
    //
    info = r8mat_solve(N, RHS_NUM, a);

    if ( info != 0 ) {
        std::cout << "\n";
        std::cout << "TETRAHEDRON_BARYCENTRIC - Fatal error!\n";
        std::cout << "  The linear system is singular.\n";
        std::cout << "  The input data does not form a proper tetrahedron.\n";
        exit(1);
    }

    c = new double [ 4 ];

    c [ 1 ] = a [ 0 + 3 * N ];
    c [ 2 ] = a [ 1 + 3 * N ];
    c [ 3 ] = a [ 2 + 3 * N ];

    c [ 0 ] = 1.0 - c [ 1 ] - c [ 2 ] - c [ 3 ];

    return c;

# undef N
# undef RHS_NUM
}


int Grid::r8mat_solve(int n, int rhs_num, double a[])

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
//
//  Discussion:
//
//    A R8MAT is a doubly dimensioned array of double precision values, which
//    may be stored as a vector in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*N]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int RHS_NUM, the number of right hand sides.  RHS_NUM
//    must be at least 0.
//
//    Input/output, double A[N*(N+RHS_NUM)], contains in rows and columns 1
//    to N the coefficient matrix, and in columns N+1 through
//    N+RHS_NUM, the right hand sides.  On output, the coefficient matrix
//    area has been destroyed, while the right hand sides have
//    been overwritten with the corresponding solutions.
//
//    Output, int R8MAT_SOLVE, singularity flag.
//    0, the matrix was not singular, the solutions were computed;
//    J, factorization failed on step J, and the solutions could not
//    be computed.
//
{
    double apivot;
    double factor;
    int i;
    int ipivot;
    int j;
    int k;
    double temp;

    for ( j = 0; j < n; j++ ) {
        //
        //  Choose a pivot row.
        //
        ipivot = j;
        apivot = a [ j + j * n ];

        for ( i = j; i < n; i++ ) {
            if ( fabs(apivot) < fabs(a [ i + j * n ]) ) {
                apivot = a [ i + j * n ];
                ipivot = i;
            }
        }

        if ( apivot == 0.0 ) {
            return j;
        }
        //
        //  Interchange.
        //
        for ( i = 0; i < n + rhs_num; i++ ) {
            temp          = a [ ipivot + i * n ];
            a [ ipivot + i * n ] = a [ j + i * n ];
            a [ j + i * n ]      = temp;
        }
        //
        //  A(J,J) becomes 1.
        //
        a [ j + j * n ] = 1.0;
        for ( k = j; k < n + rhs_num; k++ ) {
            a [ j + k * n ] = a [ j + k * n ] / apivot;
        }
        //
        //  A(I,J) becomes 0.
        //
        for ( i = 0; i < n; i++ ) {
            if ( i != j ) {
                factor = a [ i + j * n ];
                a [ i + j * n ] = 0.0;
                for ( k = j; k < n + rhs_num; k++ ) {
                    a [ i + k * n ] = a [ i + k * n ] - factor * a [ j + k * n ];
                }
            }
        }
    }

    return 0;
}

//#endif
