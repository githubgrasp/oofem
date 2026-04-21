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
#include "octreegridlocalizer.h"
#include "error.h"
#include "floatarray.h"
#include "prism.h"
#include "convertererror.h"
#include <sstream>

#include <map>
#include <unordered_map>
#include <utility>


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


struct FaceKey {
    int a, b, c;

    FaceKey() : a(0), b(0), c(0) {}

    FaceKey(int i, int j, int k)
    {
        if ( i > j ) {
            std::swap(i, j);
        }
        if ( j > k ) {
            std::swap(j, k);
        }
        if ( i > j ) {
            std::swap(i, j);
        }
        a = i;
        b = j;
        c = k;
    }

    bool operator < ( const FaceKey &other ) const
    {
        if ( a != other.a ) {
            return a < other.a;
        }
        if ( b != other.b ) {
            return b < other.b;
        }
        return c < other.c;
    }
};


Grid::Grid(int i)
{
    delaunayLocalizer      = NULL;
    voronoiLocalizer       = NULL;
    reinforcementLocalizer = NULL;

    liveDir.resize(3);
    liveDir.at(1) = 0.0;
    liveDir.at(2) = 0.0;
    liveDir.at(3) = -1.0;
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
    } else {
        converter::errorf("Unknown grid type %s\n", name.c_str() );
    }
    return;
}



bool Grid::readT3d(const std::string &fn,
                   std::vector < Node > & nodes,
                   std::vector < Tri > & tris,
                   std::vector < Tet > & tets)
{
    std::ifstream in(fn);
    if ( !in ) {
        return false;
    }

    int meshType, deg, renum, outType;
    in >> meshType >> deg >> renum >> outType;

    this->t3dOutType = outType;

    int nNodes = 0, nEdges = 0, nTris = 0, nQuads = 0, nTets = 0, nPyr = 0, nWed = 0, nHex = 0;

    if ( meshType == 7 ) {
        in >> nNodes >> nEdges >> nTris >> nQuads >> nTets >> nPyr >> nWed >> nHex;
    } else if ( meshType == 3 ) {
        int nTetras;
        in >> nNodes >> nEdges >> nTris >> nTetras;
        nTets = nTetras;
    } else if ( meshType == 4 ) {
        int nHexas;
        in >> nNodes >> nEdges >> nQuads >> nHexas;
        nHex = nHexas;
    } else {
        return false;
    }

    entityNodes.clear();
    entityTris.clear();

    nodes.clear();
    nodes.reserve(nNodes);

    for (int i = 0; i < nNodes; i++) {
        Node n;
        int entProp;
        in >> n.id >> n.x >> n.y >> n.z >> n.entType >> n.entID >> entProp;
        nodes.push_back(n);
        entityNodes [ n.entType ] [ n.entID ].push_back(n.id);
    }

    // ---------- triangles ----------
    tris.clear();
    tris.reserve(nTris);

    for (int i = 0; i < nTris; i++) {
        Tri t {};
        in >> t.id >> t.n1 >> t.n2 >> t.n3 >> t.entType >> t.entID >> t.entProp;

        // optional triangle fields (order matters)

        // iso type (bit 64)
        if ( t3dOutType & 64 ) {
            int iso;
            in >> iso;
        }

        // neighbour element IDs (bit 32)
        if ( t3dOutType & 32 ) {
            int ng1, ng2, ng3;
            in >> ng1 >> ng2 >> ng3;
        }

        // boundary curve IDs + props (bit 8)
        if ( t3dOutType & 8 ) {
            in >> t.bndCurveId [ 0 ] >> t.bndCurveId [ 1 ] >> t.bndCurveId [ 2 ]
            >> t.bndCurveProp [ 0 ] >> t.bndCurveProp [ 1 ] >> t.bndCurveProp [ 2 ];
        }

        tris.push_back(t);
    }

    entityTris.clear();
    for (size_t i = 0; i < tris.size(); ++i) {
        entityTris [ tris [ i ].entType ] [ tris [ i ].entID ].push_back( ( int ) i);
    }

    int cnt = 0;
    for (const auto &t : tris) {
        if ( t.bndCurveId [ 0 ] || t.bndCurveId [ 1 ] || t.bndCurveId [ 2 ] ) {
            cnt++;
        }
    }
    printf("Triangles touching curves: %d / %d\n", cnt, nTris);


    // ---------- tetrahedra ----------
    tets.clear();
    tets.reserve(nTets);

    // consume remainder of current line before using getline
    std::string line;
    std::getline(in, line);

    for (int i = 0; i < nTets; ) {
        if ( !std::getline(in, line) ) {
            converter::error("Unexpected end of file while reading tetrahedra");
        }

        // skip empty lines
        if ( line.find_first_not_of(" \t\r\n") == std::string::npos ) {
            continue;
        }

        std::istringstream iss(line);

        Tet t {};

        if ( !( iss >> t.id >> t.n1 >> t.n2 >> t.n3 >> t.n4
                >> t.entType >> t.entID >> t.entProp ) ) {
            converter::error("Failed to parse tetrahedron line");
        }

        // For outType 8, tetra lines appear to contain:
        // 4 face entity IDs, 4 face entity types, 4 face props
        for (int k = 0; k < 4; ++k) {
            if ( !( iss >> t.faceEntID [ k ] ) ) {
                t.faceEntID [ k ] = 0;
            }
        }

        for (int k = 0; k < 4; ++k) {
            if ( !( iss >> t.faceEntType [ k ] ) ) {
                t.faceEntType [ k ] = 0;
            }
        }

        for (int k = 0; k < 4; ++k) {
            if ( !( iss >> t.faceEntProp [ k ] ) ) {
                t.faceEntProp [ k ] = 0;
            }
        }

        tets.push_back(t);
        ++i;
    }

    return true;
}


int Grid::entityTypeFromString(const std::string &s) const {
    if ( s == "vertex" ) {
        return 1;
    }
    if ( s == "curve" ) {
        return 2;
    }
    if ( s == "surface" ) {
        return 3;
    }
    if ( s == "patch" ) {
        return 5;
    }
    if ( s == "shell" ) {
        return 6;
    }
    return -1;
}


/*This function reconstructs the boundary triangle mesh from the tetrahedra and transfers the T3D face classification, so that the original T3D patch/surface numbering can be used for loads, BCs, sets, and other input control.*/
void Grid::buildBoundaryTrisFromTets()
{
    struct FaceInfo {
        int tetIndex;
        int localFace;
    };

    std::map < FaceKey, std::vector < FaceInfo >> faceMap;

    for (size_t ti = 0; ti < tets.size(); ++ti) {
        const Tet &t = tets [ ti ];

        int fn[ 4 ] [ 3 ] = {
            {
                t.n1, t.n2, t.n3
            },                // T3D face 1
            {
                t.n1, t.n2, t.n4
            },                // T3D face 2
            {
                t.n2, t.n3, t.n4
            },                // T3D face 3
            {
                t.n1, t.n3, t.n4
            }                 // T3D face 4
        };

        for (int k = 0; k < 4; ++k) {
            FaceKey key(fn [ k ] [ 0 ], fn [ k ] [ 1 ], fn [ k ] [ 2 ]);
            faceMap [ key ].push_back({ ( int ) ti, k });
        }
    }

    tris.clear();
    int triID = 1;

    for (const auto &kv : faceMap) {
        if ( kv.second.size() == 1 ) {
            const FaceInfo &fi = kv.second [ 0 ];
            const Tet &t = tets [ fi.tetIndex ];
            int lf = fi.localFace;

            Tri tr {};
            tr.id = triID++;
            tr.n1 = kv.first.a;
            tr.n2 = kv.first.b;
            tr.n3 = kv.first.c;

            tr.entType = t.faceEntType [ lf ];
            tr.entID   = t.faceEntID [ lf ];
            tr.entProp = t.faceEntProp [ lf ];

            tr.bndCurveId [ 0 ] = tr.bndCurveId [ 1 ] = tr.bndCurveId [ 2 ] = 0;
            tr.bndCurveProp [ 0 ] = tr.bndCurveProp [ 1 ] = tr.bndCurveProp [ 2 ] = 0;

            tris.push_back(tr);
        }
    }
}

void Grid::buildCurveSegsFromTris()
{
    curveSegs.clear();
    curveToSegIdx.clear();

    std::unordered_set < long long > seen;
    seen.reserve(tris.size() * 2);

    auto makeKey = [] ( int cid, int a, int b )->long long {
        if ( a > b ) {
            std::swap(a, b);
        }
        return ( ( long long ) cid << 42 ) ^ ( ( long long ) a << 21 ) ^ ( long long ) b;
    };

    for (const Tri &t : tris) {
        const int n[ 3 ] = {
            t.n1, t.n2, t.n3
        };
        const int ea[ 3 ] [ 2 ] = { {
                                        0, 1
                                    }, {
                                        1, 2
                                    }, {
                                        2, 0
                                    } };              // (n1-n2),(n2-n3),(n3-n1)

        for (int e = 0; e < 3; ++e) {
            const int cid = t.bndCurveId [ e ];
            if ( cid <= 0 ) {
                continue;
            }

            int a = n [ ea [ e ] [ 0 ] ];
            int b = n [ ea [ e ] [ 1 ] ];
            if ( a > b ) {
                std::swap(a, b);
            }

            long long key = makeKey(cid, a, b);
            if ( seen.find(key) != seen.end() ) {
                continue;
            }
            seen.insert(key);

            int idx = ( int ) curveSegs.size();
            curveSegs.push_back({ a, b, cid });
            curveToSegIdx [ cid ].push_back(idx);
        }
    }
}


void Grid::buildEdgeAdjacency3D()
{
    edgeToTets.assign(edges.size(), {});
    edgeToBoundaryTris.assign(edges.size(), {});

    std::map < std::pair < int, int >, int > edgeMap;
    for (size_t ei = 0; ei < edges.size(); ++ei) {
        int a = edges [ ei ].n1;
        int b = edges [ ei ].n2;
        if ( a > b ) {
            std::swap(a, b);
        }
        edgeMap [ { a, b } ] = ( int ) ei;
    }

    // tetra adjacency
    for (size_t ti = 0; ti < tets.size(); ++ti) {
        const Tet &t = tets [ ti ];

        int en[ 6 ] [ 2 ] = {
            {
                t.n1, t.n2
            }, {
                t.n1, t.n3
            }, {
                t.n1, t.n4
            },
            {
                t.n2, t.n3
            }, {
                t.n2, t.n4
            }, {
                t.n3, t.n4
            }
        };

        for (int k = 0; k < 6; ++k) {
            int a = en [ k ] [ 0 ];
            int b = en [ k ] [ 1 ];
            if ( a > b ) {
                std::swap(a, b);
            }

            auto it = edgeMap.find({ a, b });
            if ( it != edgeMap.end() ) {
                edgeToTets [ it->second ].push_back( ( int ) ti);
            }
        }
    }

    // boundary triangle adjacency
    for (size_t trii = 0; trii < tris.size(); ++trii) {
        const Tri &tr = tris [ trii ];

        int en[ 3 ] [ 2 ] = {
            {
                tr.n1, tr.n2
            },
            {
                tr.n2, tr.n3
            },
            {
                tr.n3, tr.n1
            }
        };

        for (int k = 0; k < 3; ++k) {
            int a = en [ k ] [ 0 ];
            int b = en [ k ] [ 1 ];
            if ( a > b ) {
                std::swap(a, b);
            }

            auto it = edgeMap.find({ a, b });
            if ( it != edgeMap.end() ) {
                edgeToBoundaryTris [ it->second ].push_back( ( int ) trii);
            }
        }
    }
}


oofem::FloatArray Grid::tetBarycentre(int tetIndex) const
{
    const Tet &t = tets [ tetIndex ];

    oofem::FloatArray x1 = getX(t.n1);
    oofem::FloatArray x2 = getX(t.n2);
    oofem::FloatArray x3 = getX(t.n3);
    oofem::FloatArray x4 = getX(t.n4);

    oofem::FloatArray c(3);
    c.zero();
    c.add(x1);
    c.add(x2);
    c.add(x3);
    c.add(x4);
    c.times(0.25);

    return c;
}


oofem::FloatArray Grid::faceBarycentre(int triIndex) const
{
    const Tri &t = tris [ triIndex ];

    oofem::FloatArray x1 = getX(t.n1);
    oofem::FloatArray x2 = getX(t.n2);
    oofem::FloatArray x3 = getX(t.n3);

    oofem::FloatArray c(3);
    c.zero();
    c.add(x1);
    c.add(x2);
    c.add(x3);
    c.times(1.0 / 3.0);

    return c;
}

void Grid::buildEdgePolygon3D(int edgeIndex, oofem::FloatArray &polycoords) const
{
    const Edge &e = edges [ edgeIndex ];

    oofem::FloatArray xi = getX(e.n1);
    oofem::FloatArray xj = getX(e.n2);

    // midpoint
    oofem::FloatArray xm(3);
    xm = xi;
    xm.add(xj);
    xm.times(0.5);

    // edge direction
    oofem::FloatArray u(3);
    u.beDifferenceOf(xj, xi);
    double L = u.computeNorm();
    if ( L <= 0.0 ) {
        converter::error("Zero length edge in buildEdgePolygon3D");
    }
    u.times(1.0 / L);

    // transverse basis r,s
    oofem::FloatArray ref(3), r(3), s(3);
    ref.zero();
    ref.at(3) = 1.0;

    if ( fabs(u.dotProduct(ref) ) > 0.9 ) {
        ref.zero();
        ref.at(2) = 1.0;
    }

    r.beVectorProductOf(ref, u);
    r.normalize();

    s.beVectorProductOf(u, r);
    s.normalize();

    struct PolarPoint {
        double ang;
        oofem::FloatArray p;
    };

    std::vector < PolarPoint > pts;

    auto addProjectedPoint = [ & ](const oofem::FloatArray & c)
    {
        oofem::FloatArray d(3), tmp(3), proj(3);
        d.beDifferenceOf(c, xm);

        tmp = u;
        tmp.times(d.dotProduct(u) );

        proj = d;
        proj.subtract(tmp);

        double pr = proj.dotProduct(r);
        double ps = proj.dotProduct(s);
        double ang = atan2(ps, pr);

        oofem::FloatArray xp(3);
        xp = xm;
        xp.add(proj);

        pts.push_back({ ang, xp });
    };

    const size_t ntets  = edgeToTets [ edgeIndex ].size();

    // --- special case: one tetra only ---
    if ( ntets == 1 ) {
        addProjectedPoint(tetBarycentre(edgeToTets [ edgeIndex ] [ 0 ]) );

        for (int triIdx : edgeToBoundaryTris [ edgeIndex ]) {
            addProjectedPoint(faceBarycentre(triIdx) );
        }

        // add midpoint itself as 4th stabilising point
        addProjectedPoint(xm);
    } else {
        // standard case
        for (int tetIdx : edgeToTets [ edgeIndex ]) {
            addProjectedPoint(tetBarycentre(tetIdx) );
        }

        for (int triIdx : edgeToBoundaryTris [ edgeIndex ]) {
            addProjectedPoint(faceBarycentre(triIdx) );
        }
    }

    if ( pts.size() < 3 ) {
        converter::error("Not enough points to build 3D edge polygon");
    }

    std::sort(pts.begin(), pts.end(),
              [] ( const PolarPoint &a, const PolarPoint &b ) {
        return a.ang < b.ang;
    });

    // remove duplicates
    const double tol = 1e-10;
    std::vector < oofem::FloatArray > uniquePts;

    for (const auto &pp : pts) {
        bool isNew = true;
        for (const auto &up : uniquePts) {
            oofem::FloatArray d(3);
            d.beDifferenceOf(pp.p, up);
            if ( d.computeNorm() < tol ) {
                isNew = false;
                break;
            }
        }
        if ( isNew ) {
            uniquePts.push_back(pp.p);
        }
    }

    if ( uniquePts.size() < 3 ) {
        converter::error("Polygon degenerates after duplicate removal");
    }

    polycoords.resize(3 * ( int ) uniquePts.size() );
    for (size_t k = 0; k < uniquePts.size(); ++k) {
        polycoords.at(3 * ( int ) k + 1) = uniquePts [ k ].at(1);
        polycoords.at(3 * ( int ) k + 2) = uniquePts [ k ].at(2);
        polycoords.at(3 * ( int ) k + 3) = uniquePts [ k ].at(3);
    }
}

void Grid::write3DEdgeSection(std::ostream &out, int &eid, const Edge &e, int edgeIndex)
{
    const EdgeSpec spec = resolveEdgeSpec(e, EdgeSpec{ "lattice3D", 1, 1 });

    oofem::FloatArray polycoords;
    buildEdgePolygon3D(edgeIndex, polycoords);

    out << spec.elementName << " " << eid++
        << " nodes 2 " << e.n1 << " " << e.n2
        << " crossSect " << spec.crossSect
        << " mat " << spec.material
        << " polycoords " << polycoords.giveSize() << " "
        << std::scientific;

    for (int k = 1; k <= polycoords.giveSize(); ++k) {
        out << polycoords.at(k);
        if ( k < polycoords.giveSize() ) {
            out << " ";
        }
    }
    out << "\n";
}


void Grid::readControlRecords()
{
    std::ifstream in(controlFileName);
    if ( !in ) {
        converter::errorf("Cannot open control file '%s'", controlFileName.c_str());
    }

    std::string line;

    while ( std::getline(in, line) ) {
        std::istringstream iss(line);

        std::string tag;
        if ( !( iss >> tag ) ) {
            continue;
        }

        // -----------------------
        // #@BC  (T3D pipeline)
        // -----------------------
        if ( tag == "#@BC" ) {
            BCRequest bc;

            std::string typeStr;
            iss >> typeStr >> bc.entID;

            bc.entType = entityTypeFromString(typeStr);
            if ( bc.entType < 0 ) {
                converter::error("Unknown BC entity type");
            }

            std::string token;
            while ( iss >> token ) {
                if ( token == "vertices" ) {
                    // read integers until next token is not an int
                    while ( true ) {
                        int v;
                        std::streampos p = iss.tellg();
                        if ( iss >> v ) {
                            bc.extraVertices.push_back(v);
                        } else {
                            iss.clear();
                            iss.seekg(p);
                            break;
                        }
                    }
                } else if ( token == "dofs" ) {
                    int n;
                    iss >> n;
                    bc.dofs.resize(n);
                    for (int i = 0; i < n; ++i) {
                        iss >> bc.dofs [ i ];
                    }
                } else if ( token == "values" ) {
                    int n;
                    iss >> n;
                    bc.values.resize(n);
                    for (int i = 0; i < n; ++i) {
                        iss >> bc.values [ i ];
                    }
                }
            }

            if ( bc.dofs.empty() || bc.values.empty() || bc.dofs.size() != bc.values.size() ) {
                converter::error("Invalid #@BC: dofs/values missing or size mismatch");
            }

            bcRequests.push_back(std::move(bc) );
        }
        // -----------------------
        // #@LOAD  (T3D pipeline)
        // -----------------------
        else if ( tag == "#@LOAD" ) {
            LoadRequest lr;

            std::string typeStr;
            iss >> typeStr >> lr.entID;            // e.g. "patch 1"
            lr.entType = entityTypeFromString(typeStr);
            if ( lr.entType < 0 ) {
                converter::error("Unknown LOAD entity type");
            }

            std::string qLabel;
            iss >> qLabel >> lr.q;                 // expects "q 3000"
            if ( qLabel != "q" ) {
                converter::error("Invalid #@LOAD: expected 'q <value>'");
            }

            loadRequests.push_back(std::move(lr) );
        }
        // ---- load direction ----  (T3D pipeline)
        else if ( tag == "#@DIR" ) {
            liveDir.resize(3);
            iss >> liveDir.at(1) >> liveDir.at(2) >> liveDir.at(3);
            liveDir.normalize();
        }
        // -----------------------
        // #@THICKNESS  (T3D pipeline)
        // -----------------------
        else if ( tag == "#@THICKNESS" ) {
            std::string next;
            if ( !( iss >> next ) ) {
                converter::error("Invalid #@THICKNESS line");
            }

            if ( std::isdigit(next [ 0 ]) || next [ 0 ] == '.' || next [ 0 ] == '-' ) {
                // Case 1: global thickness
                defaultThickness = std::stod(next);
            } else {
                // Case 2: entity-specific
                int entType = entityTypeFromString(next);
                if ( entType < 0 ) {
                    converter::error("Unknown entity type in #@THICKNESS");
                }
                int entID;
                double t;
                iss >> entID >> t;
                entityThickness [ entType ] [ entID ] = t;
            }
        } else if ( tag == "#@element" ) {
            // (T3D) #@element <entityKind> <entityID> <elementName> <crossSect> <mat>
            std::string kind;
            iss >> kind;
            int entType = entityTypeFromString(kind);
            if ( entType < 0 ) {
                converter::errorf("Unknown entity kind '%s' in #@element directive", kind.c_str());
            }
            int entID = 0;
            std::string elementName;
            int crossSect = 1, material = 1;
            if ( !( iss >> entID >> elementName >> crossSect >> material ) ) {
                converter::error("Malformed #@element directive — expected: <kind> <id> <name> <cs> <mat>");
            }
            elementSpecsByEntity[ { entType, entID } ] = EdgeSpec{ elementName, crossSect, material };
        } else if ( tag == "#@3DSECTION" ) {
            // (T3D)
            std::string mode;
            iss >> mode;
            use3DFrameSection = ( mode == "FRAME" );
        } else if ( tag == "#@SHELLWIDTHSCALE" ) {
            // (T3D)
            iss >> shellWidthScale;
        }
        // -----------------------
        // Qhull pipeline directives follow.
        // -----------------------
        else if ( tag == "#@grid" ) {
            std::string typeName;
            iss >> typeName;
            for (char &c : typeName) {
                c = std::tolower(static_cast< unsigned char >(c));
            }
            resolveGridType(typeName);
        } else if ( tag == "#@diam" ) {
            iss >> diameter;
            TOL = 1.e-6 * diameter;
        } else if ( tag == "#@perflag" ) {
            int n;
            iss >> n;
            periodicityFlag.resize(n);
            for (int i = 1; i <= n; ++i) {
                iss >> periodicityFlag.at(i);
            }
        } else if ( tag == "#@ranint" ) {
            iss >> randomInteger;
            if ( randomInteger >= 0 ) {
                randomInteger = -time(NULL);
            }
        } else if ( tag == "#@pov" ) {
            // Opt in to writing the POV-Ray rendering files alongside oofem.in.
            emitPovOutput = true;
        } else if ( tag == "#@vtk" ) {
            // Opt in to writing the ParaView .vtu files alongside oofem.in.
            emitVtkOutput = true;
        } else if ( tag == "#@prism" ) {
            int num;
            iss >> num;

            std::string boxKw;
            int boxSize;
            iss >> boxKw >> boxSize;

            oofem::FloatArray box(boxSize);
            for (int i = 1; i <= boxSize; ++i) {
                iss >> box.at(i);
            }

            auto *p = new Prism(num, this);
            p->setBox(box);
            regionList.resize(std::max(( int ) regionList.size(), num), nullptr);
            setRegion(num, p);
        } else if ( tag == "#@cylinder" ) {
            int num;
            iss >> num;

            std::string lineKw;
            int lineSize;
            iss >> lineKw >> lineSize;

            oofem::FloatArray lin(lineSize);
            for (int i = 1; i <= lineSize; ++i) {
                iss >> lin.at(i);
            }

            double rad = 0.0;
            std::string radKw;
            iss >> radKw >> rad;

            auto *c = new Cylinder(num, this);
            c->setLine(lin);
            c->setRadius(rad);
            regionList.resize(std::max(( int ) regionList.size(), num), nullptr);
            setRegion(num, c);
        } else if ( tag == "#@interfacecylinder" ) {
            int num;
            iss >> num;

            std::string lineKw;
            int lineSize;
            iss >> lineKw >> lineSize;

            oofem::FloatArray lin(lineSize);
            for (int i = 1; i <= lineSize; ++i) {
                iss >> lin.at(i);
            }

            double rad = 0.0;
            std::string radKw;
            iss >> radKw >> rad;

            double itz = diameter;
            std::string token;
            while ( iss >> token ) {
                if ( token == "itz" ) {
                    iss >> itz;
                }
            }

            auto *ic = new InterfaceCylinder(num, this);
            ic->setLine(lin);
            ic->setRadius(rad);
            ic->setITZThickness(itz);
            inclusionList.resize(std::max(( int ) inclusionList.size(), num), nullptr);
            setInclusion(num, ic);
        } else if ( tag == "#@fibre" ) {
            int num;
            iss >> num;

            std::string kw;
            int sz = 0;
            iss >> kw >> sz;
            oofem::FloatArray endpoints(sz);
            for ( int i = 1; i <= sz; ++i ) {
                iss >> endpoints.at(i);
            }

            iss >> kw;            // "diameter"
            double diam = 0.0;
            iss >> diam;

            auto *f = new Fibre(num, this);
            f->initializeFromCoords(endpoints, diam);
            fibreList.resize(std::max(( int ) fibreList.size(), num), nullptr);
            setFibre(num, f);
        } else if ( tag == "#@controlvertex" ) {
            // #@controlvertex <id> coords 3 x y z
            int id;
            iss >> id;
            std::string kw;
            int sz = 0;
            iss >> kw >> sz;
            if ( kw != "coords" || sz != 3 ) {
                converter::error("Malformed #@controlvertex — expected 'coords 3 x y z'");
            }
            oofem::FloatArray c(3);
            iss >> c.at(1) >> c.at(2) >> c.at(3);
            controlVertexDefinitions.emplace_back(id, c);
        } else if ( tag == "#@notch" ) {
            // #@notch <id> box 6 xmin ymin zmin xmax ymax zmax material <m>
            int num;
            iss >> num;
            std::string kw;
            int sz = 0;
            iss >> kw >> sz;              // "box" 6
            if ( kw != "box" || sz != 6 ) {
                converter::error("Malformed #@notch — expected 'box 6 xmin ymin zmin xmax ymax zmax material <m>'");
            }
            NotchSpec n;
            iss >> n.xmin >> n.ymin >> n.zmin >> n.xmax >> n.ymax >> n.zmax;
            iss >> kw;                    // "material"
            iss >> n.material;
            notchSpecs.push_back(n);
        } else if ( tag == "#@sphereinclusion" ) {
            // #@sphereinclusion <id> centre 3 x y z radius r itz t
            //   inside <mi> interface <mif>
            int num;
            iss >> num;
            std::string kw;
            int sz = 0;
            iss >> kw >> sz;              // "centre" 3
            if ( kw != "centre" || sz != 3 ) {
                converter::error("Malformed #@sphereinclusion — expected 'centre 3 x y z'");
            }
            SphereInclusionSpec s;
            iss >> s.cx >> s.cy >> s.cz;
            iss >> kw >> s.radius;        // "radius" r
            if ( kw != "radius" ) {
                converter::error("Malformed #@sphereinclusion — expected 'radius <r>'");
            }
            iss >> kw >> s.itz;           // "itz" t
            if ( kw != "itz" ) {
                converter::error("Malformed #@sphereinclusion — expected 'itz <t>'");
            }
            iss >> kw >> s.inside;        // "inside" mi
            if ( kw != "inside" ) {
                converter::error("Malformed #@sphereinclusion — expected 'inside <m>'");
            }
            iss >> kw >> s.interface_;    // "interface" mif
            if ( kw != "interface" ) {
                converter::error("Malformed #@sphereinclusion — expected 'interface <m>'");
            }
            sphereInclusionSpecs.push_back(s);
        } else if ( tag == "#@cylinderinclusion" ) {
            // #@cylinderinclusion <id> line 6 x1 y1 z1 x2 y2 z2
            //   radius r itz t inside <mi> interface <mif>
            int num;
            iss >> num;
            std::string kw;
            int sz = 0;
            iss >> kw >> sz;              // "line" 6
            if ( kw != "line" || sz != 6 ) {
                converter::error("Malformed #@cylinderinclusion — expected 'line 6 x1 y1 z1 x2 y2 z2'");
            }
            CylinderInclusionSpec c;
            iss >> c.x1 >> c.y1 >> c.z1 >> c.x2 >> c.y2 >> c.z2;
            iss >> kw >> c.radius;        // "radius" r
            if ( kw != "radius" ) {
                converter::error("Malformed #@cylinderinclusion — expected 'radius <r>'");
            }
            iss >> kw >> c.itz;           // "itz" t
            if ( kw != "itz" ) {
                converter::error("Malformed #@cylinderinclusion — expected 'itz <t>'");
            }
            iss >> kw >> c.inside;
            if ( kw != "inside" ) {
                converter::error("Malformed #@cylinderinclusion — expected 'inside <m>'");
            }
            iss >> kw >> c.interface_;
            if ( kw != "interface" ) {
                converter::error("Malformed #@cylinderinclusion — expected 'interface <m>'");
            }
            cylinderInclusionSpecs.push_back(c);
        } else if ( tag == "#@bodyload" ) {
            // #@bodyload <mat> <bc_id> — any element whose crossSect/mat
            // resolves to <mat> gets "bodyloads 1 <bc_id>" appended. Decouples
            // bodyload placement from inclusion directives: Wong uses
            // `#@bodyload 1 2` (matrix only, eigenstrain), the corrosion
            // cylinder test uses `#@bodyload 3 3` (interface only,
            // eigendisplacement).
            int mat = 0, bc = 0;
            if ( !( iss >> mat >> bc ) ) {
                converter::error("Malformed #@bodyload — expected '<mat> <bc_id>'");
            }
            bodyloadByMaterial[mat] = bc;
        } else if ( tag == "#@couplingflag" ) {
            // #@couplingflag — toggle. When present, emitters append
            // "couplingflag 1 couplingnumber N <ids>" to each element record.
            // Used by staggered SMTM analyses (e.g. Wong percolation) where
            // the SM and TM subproblems exchange data element-by-element.
            emitCouplingFlag = true;
        } else if ( tag == "#@rigidarm" ) {
            // #@rigidarm <master_ctl_id> face <axis> <side>
            //   mastermask 6 m1..m6 doftype 6 d1..d6
            // `<axis>` is 1|2|3 (x|y|z); `<side>` is min|max.
            RigidArmSpec ra;
            std::string kw, sideWord;
            int axis = 0;
            if ( !( iss >> ra.masterCtlId >> kw >> axis >> sideWord ) || kw != "face" ) {
                converter::error("Malformed #@rigidarm — expected '<master_ctl_id> face <axis> <side> mastermask 6 … doftype 6 …'");
            }
            if ( axis < 1 || axis > 3 ) {
                converter::error("#@rigidarm — axis must be 1, 2 or 3");
            }
            if ( sideWord != "min" && sideWord != "max" ) {
                converter::error("#@rigidarm — side must be 'min' or 'max'");
            }
            ra.axis = axis;
            ra.sideMax = ( sideWord == "max" );
            iss >> kw;                // "mastermask"
            if ( kw != "mastermask" ) {
                converter::error("#@rigidarm — expected 'mastermask 6 …'");
            }
            int sz = 0;
            iss >> sz;
            if ( sz != 6 ) {
                converter::error("#@rigidarm — mastermask must have size 6");
            }
            ra.mastermask.resize(6);
            for ( int i = 1; i <= 6; ++i ) iss >> ra.mastermask.at(i);
            iss >> kw;                // "doftype"
            if ( kw != "doftype" ) {
                converter::error("#@rigidarm — expected 'doftype 6 …' after mastermask");
            }
            iss >> sz;
            if ( sz != 6 ) {
                converter::error("#@rigidarm — doftype must have size 6");
            }
            ra.doftype.resize(6);
            for ( int i = 1; i <= 6; ++i ) iss >> ra.doftype.at(i);
            rigidArmSpecs.push_back(ra);
        } else if ( tag == "#@material_around" ) {
            // #@material_around <ctl_id> material <m>
            MaterialAroundSpec m;
            std::string kw;
            if ( !( iss >> m.ctlId >> kw >> m.material ) || kw != "material" ) {
                converter::error("Malformed #@material_around — expected '<ctl_id> material <m>'");
            }
            materialAroundSpecs.push_back(m);
        } else if ( tag == "#@inclusionfile" ) {
            // #@inclusionfile <path> itz <t> inside <m> interface <m>
            std::string path;
            iss >> path;
            double itz = 0.0;
            int insideMat = 2;
            int interfaceMat = 3;
            std::string sub;
            while ( iss >> sub ) {
                if ( sub == "itz" ) {
                    iss >> itz;
                } else if ( sub == "inside" ) {
                    iss >> insideMat;
                } else if ( sub == "interface" ) {
                    iss >> interfaceMat;
                } else {
                    converter::errorf("readControlRecords: '#@inclusionfile' unknown sub-keyword '%s'",
                                      sub.c_str());
                }
            }
            this->readInclusionFile(path, itz, insideMat, interfaceMat);
        }
    }
}


void Grid::readInclusionFile(const std::string &path,
                             double itz, int insideMaterial, int interfaceMaterial)
{
    std::ifstream in(path);
    if ( !in ) {
        converter::errorf("Grid::readInclusionFile: cannot open '%s'", path.c_str());
    }

    int sphereCount = 0;
    int ellipsoidCount = 0;
    int fibreCount = 0;
    std::string line;
    while ( std::getline(in, line) ) {
        if ( line.empty() || line[0] == '#' ) {
            continue;
        }
        std::istringstream iss(line);
        std::string keyword;
        if ( !( iss >> keyword ) ) {
            continue;
        }

        if ( keyword == "sphere" ) {
            int packingId;
            iss >> packingId;
            std::string kw;
            int sz = 0;
            iss >> kw >> sz;            // "centre" 3
            if ( kw != "centre" || sz != 3 ) {
                converter::errorf("Grid::readInclusionFile: sphere line malformed in '%s'", path.c_str());
            }
            SphereInclusionSpec s;
            iss >> s.cx >> s.cy >> s.cz;
            iss >> kw >> s.radius;       // "radius" r
            if ( kw != "radius" ) {
                converter::errorf("Grid::readInclusionFile: sphere line malformed (missing radius) in '%s'", path.c_str());
            }
            s.itz = itz;
            s.inside = insideMaterial;
            s.interface_ = interfaceMaterial;
            sphereInclusionSpecs.push_back(s);
            ++sphereCount;
        } else if ( keyword == "ellipsoid" ) {
            ++ellipsoidCount;
        } else if ( keyword == "fibre" ) {
            int packingId;
            iss >> packingId;
            std::string kw;
            int sz = 0;
            iss >> kw >> sz;            // "endpointcoords" 6
            oofem::FloatArray endpoints(sz);
            for ( int i = 1; i <= sz; ++i ) {
                iss >> endpoints.at(i);
            }
            iss >> kw;                  // "diameter"
            double diam = 0.0;
            iss >> diam;
            // Renumber to avoid collisions with `#@fibre`-declared fibres.
            const int newId = static_cast<int>( fibreList.size() ) + 1;
            auto *f = new Fibre(newId, this);
            f->initializeFromCoords(endpoints, diam);
            fibreList.resize(newId, nullptr);
            setFibre(newId, f);
            ++fibreCount;
        }
        // Unknown keywords are ignored to stay forward-compatible.
    }

    if ( ellipsoidCount > 0 ) {
        std::fprintf(stderr,
                     "Warning: Grid::readInclusionFile: %d ellipsoid line(s) in '%s' ignored "
                     "(converter only handles sphere inclusions for material assignment)\n",
                     ellipsoidCount, path.c_str());
    }
    std::printf("readInclusionFile('%s'): %d sphere(s), %d fibre(s) loaded\n",
                path.c_str(), sphereCount, fibreCount);
}


void Grid::rebuildEntityTris()
{
    entityTris.clear();
    for (size_t i = 0; i < tris.size(); ++i) {
        entityTris [ tris [ i ].entType ] [ tris [ i ].entID ].push_back( ( int ) i);
    }
}


double Grid::triArea(int triIndex) const
{
    const Tri &t = tris [ triIndex ];

    oofem::FloatArray x1 = getX(t.n1);
    oofem::FloatArray x2 = getX(t.n2);
    oofem::FloatArray x3 = getX(t.n3);

    oofem::FloatArray a(3), b(3), c(3);

    a.beDifferenceOf(x2, x1);
    b.beDifferenceOf(x3, x1);

    c.beVectorProductOf(a, b);

    return 0.5 * c.computeNorm();
}


double Grid::computePolygonAreaProjected(const oofem::FloatArray &polycoords,
                                         const oofem::FloatArray &xm,
                                         const oofem::FloatArray &r,
                                         const oofem::FloatArray &s) const
{
    int n = polycoords.giveSize() / 3;

    if ( n < 3 ) {
        return 0.0;
    }

    std::vector < double > x(n), y(n);

    for (int i = 0; i < n; ++i) {
        oofem::FloatArray p(3), d(3);
        p.at(1) = polycoords.at(3 * i + 1);
        p.at(2) = polycoords.at(3 * i + 2);
        p.at(3) = polycoords.at(3 * i + 3);

        d.beDifferenceOf(p, xm);

        x [ i ] = d.dotProduct(r);
        y [ i ] = d.dotProduct(s);
    }

    double A = 0.0;
    for (int i = 0; i < n; ++i) {
        int j = ( i + 1 ) % n;
        A += x [ i ] * y [ j ] - x [ j ] * y [ i ];
    }

    return 0.5 * fabs(A);
}


double Grid::tetVolume(int tetIndex) const
{
    const Tet &t = tets [ tetIndex ];

    oofem::FloatArray x1 = getX(t.n1);
    oofem::FloatArray x2 = getX(t.n2);
    oofem::FloatArray x3 = getX(t.n3);
    oofem::FloatArray x4 = getX(t.n4);

    oofem::FloatArray a(3), b(3), c(3), cross(3);

    a.beDifferenceOf(x2, x1);
    b.beDifferenceOf(x3, x1);
    c.beDifferenceOf(x4, x1);

    cross.beVectorProductOf(b, c);

    double vol = fabs(a.dotProduct(cross) ) / 6.0;

    return vol;
}


double Grid::computeTargetArea(int edgeIndex, double le) const
{
    double a = 0.0;

    for (int tetIdx : edgeToTets [ edgeIndex ]) {
        double v = tetVolume(tetIdx);
        a += v / ( 2.0 * le );
    }

    return a;
}


double Grid::giveThicknessForEntity(int entType, int entID) const
{
    auto itT = entityThickness.find(entType);
    if ( itT != entityThickness.end() ) {
        auto itID = itT->second.find(entID);
        if ( itID != itT->second.end() ) {
            return itID->second;
        }
    }
    return defaultThickness;
}


double Grid::edgeLength(const Edge &e) const
{
    oofem::FloatArray xi = getX(e.n1);
    oofem::FloatArray xj = getX(e.n2);

    oofem::FloatArray d(3);
    d.beDifferenceOf(xj, xi);

    return d.computeNorm();
}


void Grid::computeEdgeWidths(double thickness)
{
    edgeWidth.assign(edges.size(), 0.0);

    // PASS 1 — accumulate triangle area to edges
    for (size_t ti = 0; ti < tris.size(); ++ti) {
        double At = triArea(ti);

        const Tri &t = tris [ ti ];

        // edge lengths
        double L1 = ( getX(t.n1) - getX(t.n2) ).computeNorm();
        double L2 = ( getX(t.n2) - getX(t.n3) ).computeNorm();
        double L3 = ( getX(t.n3) - getX(t.n1) ).computeNorm();

        double Lsum = L1 + L2 + L3;

        // contributions
        double A1 = At * L1 / Lsum;
        double A2 = At * L2 / Lsum;
        double A3 = At * L3 / Lsum;

        // find edges and add
        auto addToEdge = [ & ](int nA, int nB, double Ae)
        {
            if ( nA > nB ) {
                std::swap(nA, nB);
            }

            for (size_t ei = 0; ei < edges.size(); ++ei) {
                if ( edges [ ei ].n1 == nA && edges [ ei ].n2 == nB ) {
                    edgeWidth [ ei ] += Ae;
                    break;
                }
            }
        };

        addToEdge(t.n1, t.n2, A1);
        addToEdge(t.n2, t.n3, A2);
        addToEdge(t.n3, t.n1, A3);
    }

    // PASS 2 — convert area to width
    for (size_t ei = 0; ei < edges.size(); ++ei) {
        double Le = edgeLength(edges [ ei ]);
        edgeWidth [ ei ] /= ( Le );
        edgeWidth [ ei ] *= 2.;
    }
}


void Grid::buildEdges(const std::vector < Tri > & tris,
                      const std::vector < Tet > & tets,
                      std::vector < Edge > & edges)
{
    std::map < std::pair < int, int >, int > edgeMap;
    edges.clear();

    auto addEdge = [ & ](int a, int b, int triIdx, int tetIdx)
    {
        if ( a > b ) {
            std::swap(a, b);
        }
        auto key = std::make_pair(a, b);

        auto it = edgeMap.find(key);
        if ( it == edgeMap.end() ) {
            Edge e;
            e.n1 = a;
            e.n2 = b;
            if ( triIdx >= 0 ) {
                e.tri1 = triIdx;
            }
            if ( tetIdx >= 0 ) {
                e.tet1 = tetIdx;
            }

            edges.push_back(e);
            edgeMap [ key ] = ( int ) edges.size() - 1;
        } else {
            Edge &e = edges [ it->second ];
            if ( triIdx >= 0 ) {
                if ( e.tri1 < 0 ) {
                    e.tri1 = triIdx;
                } else {
                    e.tri2 = triIdx;
                }
            }
            if ( tetIdx >= 0 ) {
                if ( e.tet1 < 0 ) {
                    e.tet1 = tetIdx;
                } else {
                    e.tet2 = tetIdx;
                }
            }
        }
    };

    // --- triangles ---
    for (int i = 0; i < ( int ) tris.size(); ++i) {
        const auto &t = tris [ i ];
        addEdge(t.n1, t.n2, i, -1);
        addEdge(t.n2, t.n3, i, -1);
        addEdge(t.n3, t.n1, i, -1);
    }

    // --- tetrahedra ---
    for (int i = 0; i < ( int ) tets.size(); ++i) {
        const auto &t = tets [ i ];
        addEdge(t.n1, t.n2, -1, i);
        addEdge(t.n1, t.n3, -1, i);
        addEdge(t.n1, t.n4, -1, i);
        addEdge(t.n2, t.n3, -1, i);
        addEdge(t.n2, t.n4, -1, i);
        addEdge(t.n3, t.n4, -1, i);
    }
}

void Grid::writeLiveLoads(std::ostream &out, int &bcID)
{
    const size_t n = nodes.size();

    for (size_t i = 0; i < n; ++i) {
        if ( loadNodeSetID [ i ] < 0 ) {
            continue;
        }

        out << "NodalLoad " << bcID++
            << " loadTimeFunction 2"
            << " dofs 3 1 2 3"
            << " components 3 "
            << std::scientific
            << loadFx [ i ] << " "
            << loadFy [ i ] << " "
            << loadFz [ i ]
            << " set " << loadNodeSetID [ i ]
            << "\n";
    }
}

void Grid::computeNodalAreasOnTriEntity(int entType, int entID, std::vector < double > & A) const
{
    A.assign(nodes.size(), 0.0);

    auto itT = entityTris.find(entType);
    if ( itT == entityTris.end() ) {
        return;
    }

    auto itID = itT->second.find(entID);
    if ( itID == itT->second.end() ) {
        return;
    }

    for (int ti : itID->second) {
        const Tri &tr = tris [ ti ];
        double At = triArea(ti);

        A [ nodeIndex.at(tr.n1) ] += At / 3.0;
        A [ nodeIndex.at(tr.n2) ] += At / 3.0;
        A [ nodeIndex.at(tr.n3) ] += At / 3.0;
    }
}


int Grid::instanciateYourselfFromT3d(const std::string &t3d, const std::string &control)
{
    this->controlFileName = control;

    // 1) read mesh
    if ( !readT3d(t3d, nodes, tris, tets) ) {
        converter::error("Failed to read T3d file");
    }

    bool is3D = !tets.empty();

    // build node index (only once needed)
    nodeIndex.clear();
    nodeIndex.reserve(nodes.size() );
    for (size_t i = 0; i < nodes.size(); ++i) {
        nodeIndex [ nodes [ i ].id ] = ( int ) i;
    }

    //Preprocess based on mesh type
    if ( is3D ) {
        // construct boundary triangles if missing
        if ( tris.empty() ) {
            buildBoundaryTrisFromTets();
            rebuildEntityTris();
        }
    } else {
        // existing surface logic
        buildCurveSegsFromTris();
    }

    // Use edges of triangles/tretras for lattice elemenst
    buildEdges(tris, tets, edges);

    //Read control commands to generate oofem input
    readControlRecords();


    // Carry out modifications to the geometry
    if ( is3D ) {
        buildEdgeAdjacency3D();
    }

    // -------------------------------------------------
    // 5) OUTPUT INFO
    // -------------------------------------------------

    int nNodes = ( int ) nodes.size();
    int nElems = ( int ) edges.size();

    int bnd = 0, interior = 0;
    for ( const auto &e : edges ) {
        ( e.tri2 < 0 ? bnd : interior )++;
    }

    printf("T3D mesh: %d nodes, %zu triangles, %zu tetrahedra, %d edges (%d boundary, %d interior)\n",
           nNodes, tris.size(), tets.size(), nElems, bnd, interior);
    return 1;
}


void Grid::computeNodalLengthsOnCurve(int curveID, std::vector < double > & L) const
{
    L.assign(nodes.size(), 0.0);

    auto it = curveToSegIdx.find(curveID);
    if ( it == curveToSegIdx.end() ) {
        return;
    }

    for (int si : it->second) {
        const CurveSeg &s = curveSegs [ si ];
        const double len = segLength(s.n1, s.n2);

        auto it1 = nodeIndex.find(s.n1);
        auto it2 = nodeIndex.find(s.n2);
        if ( it1 == nodeIndex.end() || it2 == nodeIndex.end() ) {
            converter::errorf("Missing nodeIndex for segment (%d, %d) on curve %d",
                              s.n1, s.n2, curveID);
        }

        L [ it1->second ] += 0.5 * len;
        L [ it2->second ] += 0.5 * len;
    }
}


std::set < int > Grid::collectBCNodes(const BCRequest &bc) const
{
    std::set < int > result;

    // interior nodes
    auto tIt = entityNodes.find(bc.entType);
    if ( tIt != entityNodes.end() ) {
        auto idIt = tIt->second.find(bc.entID);
        if ( idIt != tIt->second.end() ) {
            result.insert(idIt->second.begin(), idIt->second.end() );
        }
    }

    // add endpoint vertices if curve
    for (int vid : bc.extraVertices) {
        auto vIt = entityNodes.find(1);
        if ( vIt != entityNodes.end() ) {
            auto idIt = vIt->second.find(vid);
            if ( idIt != vIt->second.end() ) {
                result.insert(idIt->second.begin(), idIt->second.end() );
            }
        }
    }

    return result;
}

double Grid::segLength(int n1, int n2) const
{
    const Node &a = nodes.at(nodeIndex.at(n1) );
    const Node &b = nodes.at(nodeIndex.at(n2) );
    const double dx = b.x - a.x, dy = b.y - a.y, dz = b.z - a.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}


oofem::FloatArray Grid::getX(int nodeID) const
{
    auto it = nodeIndex.find(nodeID);
    if ( it == nodeIndex.end() ) {
        converter::errorf("getX: nodeID %d not found in nodeIndex", nodeID);
    }

    const Node &n = nodes [ it->second ];
    oofem::FloatArray x(3);
    x.at(1) = n.x;
    x.at(2) = n.y;
    x.at(3) = n.z;
    return x;
}

oofem::FloatArray Grid::triNormal(int triIndex) const
{
    const Tri &t = tris [ triIndex ];

    oofem::FloatArray x1 = getX(t.n1);
    oofem::FloatArray x2 = getX(t.n2);
    oofem::FloatArray x3 = getX(t.n3);

    oofem::FloatArray a(3), b(3), n(3);

    a.beDifferenceOf(x2, x1);
    b.beDifferenceOf(x3, x1);
    n.beVectorProductOf(a, b);
    n.normalize();

    return n;
}


int Grid::instanciateYourselfFromQhull(const std::string &controlFile,
                                       const char *nodeFileName,
                                       const char *voronoiFileName)
{
    controlFileName = controlFile;
    randomInteger   = -time(NULL);
    periodicityFlag.resize(3);
    periodicityFlag.zero();
    gridType = _3dSM; // default; overridden by #@grid if present

    readControlRecords();

    // #@perflag may have resized periodicityFlag to a size ≠ 3; the qhull
    // writers index it as at(1..3), so pad back out to 3 zeros if needed.
    if ( periodicityFlag.giveSize() != 3 ) {
        periodicityFlag.resize(3);
        periodicityFlag.zero();
    }

    if ( delaunayLocalizer == nullptr ) {
        delaunayLocalizer = new OctreeGridLocalizer(1, this, 0);
    }
    if ( voronoiLocalizer == nullptr ) {
        voronoiLocalizer = new OctreeGridLocalizer(1, this, 1);
    }
    if ( reinforcementLocalizer == nullptr ) {
        reinforcementLocalizer = new OctreeGridLocalizer(1, this, 2);
    }

    // read Delaunay vertices
    std::ifstream vertexField(nodeFileName);
    if ( !vertexField.is_open() ) {
        converter::errorf("instanciateYourselfFromQhull: Unable to open node file %s", nodeFileName);
    }
    vertexField.precision(16);

    int junk, nDelaunayVertices;
    oofem::FloatArray coords(3);
    vertexField >> junk >> nDelaunayVertices;

    delaunayVertexList.resize(nDelaunayVertices, nullptr);
    for (int i = 0; i < nDelaunayVertices; ++i) {
        if ( !( vertexField >> coords.at(1) >> coords.at(2) >> coords.at(3) ) ) {
            converter::errorf("instanciateYourselfFromQhull: failed to read coordinates for vertex %d", i + 1);
        }
        auto *v = new Vertex(i + 1, this);
        v->setCoordinates(coords);
        setDelaunayVertex(i + 1, v);
    }
    delaunayLocalizer->init(true);
    resolveControlVertices();
    printf("Finished Delaunay vertices (%d)\n", nDelaunayVertices);

    // read Voronoi vertices and Delaunay lines
    std::ifstream voronoiField(voronoiFileName);
    if ( !voronoiField.is_open() ) {
        converter::errorf("instanciateYourselfFromQhull: Unable to open voronoi file %s", voronoiFileName);
    }
    voronoiField.precision(16);

    int nVoronoiVertices = 0;
    voronoiField >> junk >> nVoronoiVertices;

    voronoiVertexList.resize(nVoronoiVertices, nullptr);
    for (int i = 0; i < nVoronoiVertices; ++i) {
        if ( !( voronoiField >> coords.at(1) >> coords.at(2) >> coords.at(3) ) ) {
            converter::errorf("Voronoi file: unexpected EOF reading vertex %d/%d", i + 1, nVoronoiVertices);
        }
        auto *v = new Vertex(i + 1, this);
        v->setCoordinates(coords);
        setVoronoiVertex(i + 1, v);
    }
    if ( voronoiLocalizer ) {
        // Two-arg init with nodeType=1 so the localiser actually indexes Voronoi
        // vertices; the single-arg overload defaults to nodeType=0 (Delaunay),
        // which silently miscategorises the localiser and causes downstream
        // periodic-partner lookups to return 0 matches.
        static_cast< OctreeGridLocalizer * >(voronoiLocalizer)->init(true, /*nodeType=*/ 1);
    }
    printf("Finished Voronoi vertices (%d)\n", nVoronoiVertices);

    int nDelaunayLines;
    voronoiField >> nDelaunayLines;

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

        auto *delaunayLine = new Line(i + 1, this);
        delaunayLine->setVertices(delaunayNodes);
        this->giveDelaunayVertex(delaunayNodes.at(1))->setLocalLine(i + 1);
        this->giveDelaunayVertex(delaunayNodes.at(2))->setLocalLine(i + 1);

        const int nVorNodes = size - 2;
        oofem::IntArray voronoiNodes(nVorNodes);
        for (int k = 0; k < nVorNodes; ++k) {
            voronoiField >> voronoiNodes.at(k + 1);
        }
        delaunayLine->updateCrossSectionVertices(voronoiNodes);

        oofem::IntArray nodesA(2);
        for (int m = 0; m < nVorNodes; ++m) {
            nodesA.at(1) = voronoiNodes.at(m + 1);
            nodesA.at(2) = ( m < nVorNodes - 1 ) ? voronoiNodes.at(m + 2) : voronoiNodes.at(1);

            oofem::IntArray localVoronoiLines;
            if ( nodesA.at(1) != 0 ) {
                this->giveVoronoiVertex(nodesA.at(1))->giveLocalLines(localVoronoiLines);
            } else if ( nodesA.at(2) != 0 ) {
                this->giveVoronoiVertex(nodesA.at(2))->giveLocalLines(localVoronoiLines);
            } else {
                std::fprintf(stderr, "error: cannot have two zero Voronoi nodes\n");
                std::exit(1);
            }

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
                    this->giveVoronoiLine(lid)->updateCrossSectionVertices(delaunayNodes);
                    this->giveVoronoiLine(lid)->updateCrossSectionElement(i + 1);
                    delaunayLine->updateCrossSectionElement(lid);
                    break;
                }
            }

            if ( !exists ) {
                const int newId = ++voronoiLineCounter;
                auto *vorLine = new Line(newId, this);
                vorLine->setVertices(nodesA);
                vorLine->updateCrossSectionVertices(delaunayNodes);
                vorLine->updateCrossSectionElement(i + 1);
                delaunayLine->updateCrossSectionElement(newId);
                converter::put1_replace(voronoiLineList, newId, vorLine);
                if ( nodesA.at(1) != 0 ) {
                    this->giveVoronoiVertex(nodesA.at(1))->setLocalLine(newId);
                }
                if ( nodesA.at(2) != 0 ) {
                    this->giveVoronoiVertex(nodesA.at(2))->setLocalLine(newId);
                }
            }
        }

        converter::put1_replace(delaunayLineList, i + 1, delaunayLine);
    }

    printf("Finished Delaunay and Voronoi lines\n");

    // build cell vertex lists (needed for VTK output)
    oofem::IntArray localLines, crossSectionNodes;
    for (int i = 0; i < nDelaunayVertices; i++) {
        this->giveDelaunayVertex(i + 1)->giveLocalLines(localLines);
        for (int m = 0; m < localLines.giveSize(); m++) {
            this->giveDelaunayLine(localLines.at(m + 1))->giveCrossSectionVertices(crossSectionNodes);
            this->giveDelaunayVertex(i + 1)->updateCellVertices(crossSectionNodes);
        }
    }

    // Discretise any fibres declared via #@fibre directives.
    this->discretizeFibres();

    return 1;
}




void
Grid::discretizeFibres()
{
    const int nfibre = giveNumberOfFibres();
    if ( nfibre == 0 ) return;

    Line *beamLine, *linkLine;
    oofem::IntArray beamNodes(2), linkNodes(2);
    oofem::FloatArray coordP1, coordP2;
    int beamElementCounter = converter::size1(latticeBeamList);
    int linkElementCounter = converter::size1(latticeLinkList);

    printf("Generating beam and link elements for fibres\n");

    for ( int i = 1; i <= nfibre; i++ ) {
        const double fibreDiameter = giveFibre(i)->giveDiameter();
        oofem::FloatArray fibreDirection = giveFibre(i)->giveDirectionVector();

        // Reinforcement-node placement: nodes sit at the centre of the cross-section
        // where the fibre intersects each Voronoi cell along its length.
        giveFibre(i)->discretize();

        const int numberOfReinforcementNodes = giveFibre(i)->giveNumberOfReinforcementNodes();
        const int numberOfBeams = numberOfReinforcementNodes - 1;
        const int numberOfLinks = numberOfReinforcementNodes;

        const int indexOfLinkElements = converter::size1(latticeLinkList);
        const int indexOfBeamElements = converter::size1(latticeBeamList);
        latticeLinkList.resize(indexOfLinkElements + numberOfLinks, nullptr);
        latticeBeamList.resize(indexOfBeamElements + numberOfBeams, nullptr);

        for ( int j = 1; j <= numberOfBeams; j++ ) {
            beamElementCounter++;
            beamLine = ( Line * ) ( Line(beamElementCounter + 1, this).ofType() );
            setLatticeBeam(beamElementCounter, beamLine);
            beamNodes.at(1) = giveFibre(i)->giveReinforcementNodeNumber(j);
            beamNodes.at(2) = giveFibre(i)->giveReinforcementNodeNumber(j + 1);
            beamLine->setVertices(beamNodes);
            giveReinforcementNode(beamNodes.at(1))->setLocalLine(beamElementCounter);
            giveReinforcementNode(beamNodes.at(2))->setLocalLine(beamElementCounter);
            beamLine->setDiameter(fibreDiameter);
            beamLine->setDirectionVector(fibreDirection);
        }

        for ( int j = 1; j <= numberOfLinks; j++ ) {
            linkElementCounter++;
            linkLine = ( Line * ) ( Line(linkElementCounter + 1, this).ofType() );
            setLatticeLink(linkElementCounter, linkLine);
            linkNodes.at(1) = giveFibre(i)->giveReinforcementNodeNumber(j);
            linkNodes.at(2) = giveFibre(i)->giveDelaunayNodeNumber(j);
            linkLine->setVertices(linkNodes);

            giveInterNode(giveFibre(i)->giveIntersectionPointNumber(j))->giveCoordinates(coordP1);
            giveInterNode(giveFibre(i)->giveIntersectionPointNumber(j + 1))->giveCoordinates(coordP2);
            linkLine->setAssociatedLength(Fibre::computeDistance(coordP1, coordP2));

            giveReinforcementNode(linkNodes.at(1))->setLocalLink(linkElementCounter);
            giveDelaunayVertex(linkNodes.at(2))->setLocalLink(linkElementCounter);

            linkLine->setDiameter(fibreDiameter);
            linkLine->setDirectionVector(fibreDirection);
            linkLine->setEndLength(giveFibre(i)->giveEndLength(j));
        }
    }
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
        std::fprintf(stderr, "warning: degenerate normal — orientation unchanged\n");
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
        } else if ( ( z < 0 && count > 0 ) || ( z > 0 && count < 0 ) ) {         //detected problem
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
    printf("Writing outputs\n");
    giveOofemOutput(fileName);

    if ( emitVtkOutput ) {
        giveVtkOutput2(fileName, 3);
    }
    if ( emitPovOutput ) {
        givePOVOutput(fileName);
    }
}

void Grid::giveOofemOutput(const std::string &fileName)
{
    //Start with oofem output
    printf("Writing OOFEM input file\n");
    if ( gridType == _3dSM ) { //Base implementation
        give3DSMOutput(fileName);
    } else if ( gridType == _3dTM ) { //Base implementation
        give3DTMOutput(fileName);
    } else {
        converter::error("Unknown grid type\n");
    }
    return;
};

void Grid::giveOutputT3d(const std::string &fileName)
{
    int bcID = 1;

    // Prepare BC and load sets first
    prepareBCSets();
    prepareLiveLoadSets();

    int nLoads = 0;
    for (int sid : loadNodeSetID) {
        if ( sid >= 0 ) {
            ++nLoads;
        }
    }

    int nBCs = 0;
    for (const auto &bc : bcRequests) {
        if ( bc.setID >= 0 ) {
            ++nBCs;
        }
    }

    int totalNBC = nLoads + nBCs;
    int nSets = ( int ) generatedNodeSets.size();

    printf("Writing OOFEM input file\n");

    std::ifstream ctrl(controlFileName);
    std::ofstream out(fileName);

    if ( !ctrl ) {
        converter::error("Cannot open control file");
    }
    if ( !out ) {
        converter::error("Cannot open output file");
    }

    std::string line;
    bool injected = false;

    const int nNodes = ( int ) nodes.size();
    const int nElems = ( int ) edges.size();

    while ( std::getline(ctrl, line) ) {
        std::string t = line;
        size_t pos = t.find_first_not_of(" \t");
        if ( pos != std::string::npos ) {
            t.erase(0, pos);
        } else {
            t.clear();
        }

        if ( !injected && t.rfind("ncrosssect", 0) == 0 ) {
            std::istringstream iss(t);
            std::string token;

            out << "ndofman " << nNodes
                << " nelem " << nElems << " ";

            while ( iss >> token ) {
                if ( token == "nbc" ) {
                    out << "nbc " << totalNBC << " ";
                    iss >> token; // skip old value
                } else if ( token == "nset" ) {
                    out << "nset " << nSets << " ";
                    iss >> token; // skip old value
                } else {
                    out << token << " ";
                }
            }

            out << "\n";

            writeT3dNodesOofem(out);
            writeT3dElemsOofem(out);

            injected = true;
            continue;
        }


        if ( t.rfind("#@INSERT_CROSSSECTION", 0) == 0 ) {
            if ( use3DFrameSection ) {
                out << "latticecs 1 material 1 shape 1 radius "
                    << std::scientific << frameRadius
                    << "\n";
            } else {
                converter::error("Encountered #@INSERT_CROSSSECTION but use3DFrameSection is false");
            }
            continue;
        }



        if ( t.rfind("#@INSERT_LIVELOADS", 0) == 0 ) {
            writeLiveLoads(out, bcID);
            writeBCRecords(out, bcID);
            continue;
        }

        if ( t.rfind("#@INSERT_SETS", 0) == 0 ) {
            writeGeneratedSets(out);
            continue;
        }

        if ( !isConverterDirective(t) ) {
            out << line << "\n";
        }
    }
}

void Grid::prepareBCSets()
{
    generatedNodeSets.clear();

    int setID = 1;

    for (size_t i = 0; i < bcRequests.size(); ++i) {
        std::set < int > ns = collectBCNodes(bcRequests [ i ]);
        if ( ns.empty() ) {
            continue;
        }

        SetDef sd;
        sd.setID = setID++;
        sd.nodeIDs.assign(ns.begin(), ns.end() );

        generatedNodeSets.push_back(sd);

        bcRequests [ i ].setID = sd.setID;
    }
}

void Grid::prepareLiveLoadSets()
{
    const size_t n = nodes.size();

    loadFx.assign(n, 0.0);
    loadFy.assign(n, 0.0);
    loadFz.assign(n, 0.0);
    loadNodeSetID.assign(n, -1);

    // --- accumulate forces ---
    for (const auto &req : loadRequests) {
        std::vector < double > w(n, 0.0);

        if ( req.entType == 1 ) {
            for (size_t i = 0; i < n; ++i) {
                if ( nodes [ i ].entType == 1 && nodes [ i ].entID == req.entID ) {
                    w [ i ] = 1.0;
                }
            }
        } else if ( req.entType == 2 ) {
            computeNodalLengthsOnCurve(req.entID, w);
        } else {
            computeNodalAreasOnTriEntity(req.entType, req.entID, w);
        }

        double sumW = 0.0;
        int nLoaded = 0;
        for (size_t i = 0; i < n; ++i) {
            sumW += w [ i ];
            if ( w [ i ] > 0.0 ) {
                ++nLoaded;
            }
        }

        for (size_t i = 0; i < n; ++i) {
            if ( w [ i ] <= 0.0 ) {
                continue;
            }

            double Fi = req.q * w [ i ];

            loadFx [ i ] += Fi * liveDir.at(1);
            loadFy [ i ] += Fi * liveDir.at(2);
            loadFz [ i ] += Fi * liveDir.at(3);
        }
    }

    // --- create sets ---
    int nextSetID = 1;
    for (const auto &sd : generatedNodeSets) {
        if ( sd.setID >= nextSetID ) {
            nextSetID = sd.setID + 1;
        }
    }

    for (size_t i = 0; i < n; ++i) {
        const double tol = 1e-15;
        if ( std::abs(loadFx [ i ]) < tol &&
             std::abs(loadFy [ i ]) < tol &&
             std::abs(loadFz [ i ]) < tol ) {
            continue;
        }

        SetDef sd;
        sd.setID = nextSetID++;
        sd.nodeIDs.push_back(nodes [ i ].id);

        generatedNodeSets.push_back(sd);

        loadNodeSetID [ i ] = sd.setID;
    }
}


void Grid::writeT3dNodesOofem(std::ostream &out)
{
    for (const auto &n : nodes) {
        out << "node " << n.id
            << " coords 3 "
            << std::scientific
            << n.x << " " << n.y << " " << n.z
            << "\n";
    }
}

void Grid::writeT3dElemsOofem(std::ostream &out)
{
    int eid = 1;
    const bool is3D = !tets.empty();
    const bool useFrame = use3DFrameSection;

    for (size_t i = 0; i < edges.size(); ++i) {
        const Edge &e = edges [ i ];

        // =================================================
        // 3D tetrahedral case: barycentric polygon section
        // =================================================

        if ( is3D ) {
            if ( useFrame ) {
                const EdgeSpec spec = resolveEdgeSpec(e, EdgeSpec{ "lattice3d", 1, 1 });
                out << spec.elementName << " " << eid++
                    << " nodes 2 " << e.n1 << " " << e.n2
                    << " crossSect " << spec.crossSect
                    << " mat " << spec.material
                    << "\n";
                continue;
            }

            oofem::FloatArray polycoords;
            buildEdgePolygon3D( ( int ) i, polycoords);

            // edge endpoints
            oofem::FloatArray xi = getX(e.n1);
            oofem::FloatArray xj = getX(e.n2);

            // midpoint
            oofem::FloatArray xm(3);
            xm = xi;
            xm.add(xj);
            xm.times(0.5);

            // edge direction
            oofem::FloatArray u(3);
            u.beDifferenceOf(xj, xi);
            double L = u.computeNorm();
            if ( L <= 0.0 ) {
                converter::error("Zero edge length in writeT3dElemsOofem");
            }
            u.times(1.0 / L);

            // local transverse basis r,s
            oofem::FloatArray ref(3), r(3), s(3);
            ref.zero();
            ref.at(3) = 1.0;

            if ( fabs(u.dotProduct(ref) ) > 0.9 ) {
                ref.zero();
                ref.at(2) = 1.0;
            }

            r.beVectorProductOf(ref, u);
            r.normalize();

            s.beVectorProductOf(u, r);
            s.normalize();

            // area of current polygon
            double A_geom = computePolygonAreaProjected(polycoords, xm, r, s);
            if ( A_geom <= 0.0 ) {
                converter::errorf("Zero or negative geometric area for edge %zu (%d, %d): A_geom = %g",
                                  i, e.n1, e.n2, A_geom);
            }

            // target area based on V/(3*l)
            double A_target = computeTargetArea( ( int ) i, L);
            if ( A_target <= 0.0 ) {
                converter::errorf("Zero or negative target area for edge %zu (%d, %d): A_target = %g",
                                  i, e.n1, e.n2, A_target);
            }

            double scale = sqrt(A_target / A_geom);

            // scale polygon uniformly about midpoint
            int npts = polycoords.giveSize() / 3;
            for (int k = 0; k < npts; ++k) {
                oofem::FloatArray p(3), d(3);

                p.at(1) = polycoords.at(3 * k + 1);
                p.at(2) = polycoords.at(3 * k + 2);
                p.at(3) = polycoords.at(3 * k + 3);

                d.beDifferenceOf(p, xm);
                d.times(scale);

                p = xm;
                p.add(d);

                polycoords.at(3 * k + 1) = p.at(1);
                polycoords.at(3 * k + 2) = p.at(2);
                polycoords.at(3 * k + 3) = p.at(3);
            }

            // write 3D element
            const EdgeSpec spec = resolveEdgeSpec(e, EdgeSpec{ "lattice3D", 1, 1 });
            out << spec.elementName << " " << eid++
                << " nodes 2 " << e.n1 << " " << e.n2
                << " crossSect " << spec.crossSect
                << " mat " << spec.material
                << " polycoords " << polycoords.giveSize() << " "
                << std::scientific;

            for (int k = 1; k <= polycoords.giveSize(); ++k) {
                out << polycoords.at(k);
                if ( k < polycoords.giveSize() ) {
                    out << " ";
                }
            }
            out << "\n";

            continue;
        }

        // ===========
        // Shell case:
        // ===========

        if ( useFrame ) {
            const EdgeSpec spec = resolveEdgeSpec(e, EdgeSpec{ "lattice3d", 1, 1 });
            out << spec.elementName << " " << eid++
                << " nodes 2 " << e.n1 << " " << e.n2
                << " crossSect " << spec.crossSect
                << " mat " << spec.material
                << "\n";
            continue;
        }

        // thickness from adjacent triangle(s)
        double t = 0.0;

        if ( e.tri1 >= 0 && e.tri2 >= 0 ) {
            const Tri &tr1 = tris [ e.tri1 ];
            const Tri &tr2 = tris [ e.tri2 ];

            double t1 = giveThicknessForEntity(tr1.entType, tr1.entID);
            double t2 = giveThicknessForEntity(tr2.entType, tr2.entID);

            t = 0.5 * ( t1 + t2 );
        } else if ( e.tri1 >= 0 ) {
            const Tri &tr1 = tris [ e.tri1 ];
            t = giveThicknessForEntity(tr1.entType, tr1.entID);
        } else if ( e.tri2 >= 0 ) {
            const Tri &tr2 = tris [ e.tri2 ];
            t = giveThicknessForEntity(tr2.entType, tr2.entID);
        } else {
            converter::error("Edge has no adjacent triangle; cannot determine shell thickness");
        }

        const double b = 0.5 * t;

        // endpoints
        oofem::FloatArray xi = getX(e.n1);
        oofem::FloatArray xj = getX(e.n2);

        // axis direction
        oofem::FloatArray u(3);
        u.beDifferenceOf(xj, xi);
        u.normalize();

        // midpoint
        oofem::FloatArray xm = xi;
        xm.add(xj);
        xm.times(0.5);

        // averaged normal
        oofem::FloatArray n = triNormal(e.tri1);

        if ( e.tri2 >= 0 ) {
            oofem::FloatArray n2 = triNormal(e.tri2);
            n.add(n2);
            n.normalize();
        }

        // width direction
        oofem::FloatArray r(3);
        r.beVectorProductOf(n, u);
        r.normalize();

        // --- projected side distances from barycentres ---
        double aPlus  = 0.0;
        double aMinus = 0.0;

        if ( e.tri1 >= 0 ) {
            oofem::FloatArray bc1 = triBarycentre(e.tri1);
            oofem::FloatArray d1(3);
            d1.beDifferenceOf(bc1, xm);

            double s1 = d1.dotProduct(r);

            if ( s1 >= 0.0 ) {
                aPlus = std::max(aPlus, s1);
            } else {
                aMinus = std::max(aMinus, -s1);
            }
        }

        if ( e.tri2 >= 0 ) {
            oofem::FloatArray bc2 = triBarycentre(e.tri2);
            oofem::FloatArray d2(3);
            d2.beDifferenceOf(bc2, xm);

            double s2 = d2.dotProduct(r);

            if ( s2 >= 0.0 ) {
                aPlus = std::max(aPlus, s2);
            } else {
                aMinus = std::max(aMinus, -s2);
            }
        }

        // fallback safety
        if ( aPlus <= 0.0 && aMinus <= 0.0 ) {
            converter::error("Both projected side widths are zero in writeT3dElemsOofem");
        }

        double wGeom = aPlus + aMinus;
        if ( wGeom <= 0.0 ) {
            converter::error("Geometric width is zero in writeT3dElemsOofem");
        }


        double s = 1.0;

        if ( shellWidthScale > 0.0 ) {
            s = shellWidthScale;
        } else {
            double Le = edgeLength(e);
            double wTarget = 0.0;

            if ( e.tri1 >= 0 ) {
                wTarget += 2.0 * triArea(e.tri1) / ( 3.0 * Le );
            }
            if ( e.tri2 >= 0 ) {
                wTarget += 2.0 * triArea(e.tri2) / ( 3.0 * Le );
            }

            s = wTarget / wGeom;
        }

        aPlus  *= s;
        aMinus *= s;

        // compute 4 corners
        oofem::FloatArray x1(3), x2(3), x3(3), x4(3), tmp(3);
        x1 = xm;
        x2 = xm;
        x3 = xm;
        x4 = xm;

        // x1 = xm + aPlus*r + b*n
        tmp = r;
        tmp.times(+aPlus);
        x1.add(tmp);
        tmp = n;
        tmp.times(+b);
        x1.add(tmp);

        // x2 = xm - aMinus*r + b*n
        tmp = r;
        tmp.times(-aMinus);
        x2.add(tmp);
        tmp = n;
        tmp.times(+b);
        x2.add(tmp);

        // x3 = xm - aMinus*r - b*n
        tmp = r;
        tmp.times(-aMinus);
        x3.add(tmp);
        tmp = n;
        tmp.times(-b);
        x3.add(tmp);

        // x4 = xm + aPlus*r - b*n
        tmp = r;
        tmp.times(+aPlus);
        x4.add(tmp);
        tmp = n;
        tmp.times(-b);
        x4.add(tmp);

        // write shell element
        const EdgeSpec spec = resolveEdgeSpec(e, EdgeSpec{ "lattice3Dnl", 1, 1 });
        out << spec.elementName << " " << eid++
            << " nodes 2 " << e.n1 << " " << e.n2
            << " crossSect " << spec.crossSect
            << " mat " << spec.material
            << " thickness " << t
            << " polycoords 12 "
            << std::scientific
            << x1.at(1) << " " << x1.at(2) << " " << x1.at(3) << " "
            << x2.at(1) << " " << x2.at(2) << " " << x2.at(3) << " "
            << x3.at(1) << " " << x3.at(2) << " " << x3.at(3) << " "
            << x4.at(1) << " " << x4.at(2) << " " << x4.at(3)
            << "\n";
    }
}

std::pair< int, int >
Grid::entityForEdge(const Edge &e) const
{
    auto idxA = nodeIndex.find(e.n1);
    auto idxB = nodeIndex.find(e.n2);
    if ( idxA != nodeIndex.end() && idxB != nodeIndex.end() ) {
        const Node &na = nodes [ idxA->second ];
        const Node &nb = nodes [ idxB->second ];
        // Both endpoints classified to the same low-dim entity (curve/vertex).
        // Lowest entType wins by virtue of T3D's own classification.
        if ( na.entType > 0 && na.entType == nb.entType && na.entID == nb.entID ) {
            return { na.entType, na.entID };
        }
    }

    // Adjacent triangles (shell case) — surface entity.
    if ( e.tri1 >= 0 ) {
        const Tri &t1 = tris [ e.tri1 ];
        if ( e.tri2 >= 0 ) {
            const Tri &t2 = tris [ e.tri2 ];
            if ( t1.entType == t2.entType && t1.entID == t2.entID ) {
                return { t1.entType, t1.entID };
            }
        }
        return { t1.entType, t1.entID };          // boundary edge of this surface
    }

    // Adjacent tetrahedra (3D case) — region entity.
    if ( e.tet1 >= 0 && ( int ) tets.size() > e.tet1 ) {
        const Tet &t1 = tets [ e.tet1 ];
        if ( e.tet2 >= 0 ) {
            const Tet &t2 = tets [ e.tet2 ];
            if ( t1.entType == t2.entType && t1.entID == t2.entID ) {
                return { t1.entType, t1.entID };
            }
        }
        return { t1.entType, t1.entID };
    }

    return { 0, 0 };
}


Grid::EdgeSpec
Grid::resolveEdgeSpec(const Edge &e, const EdgeSpec &defaultSpec) const
{
    auto key = entityForEdge(e);
    auto it  = elementSpecsByEntity.find(key);
    if ( it != elementSpecsByEntity.end() ) {
        return it->second;
    }
    return defaultSpec;
}


void
Grid::resolveControlVertices()
{
    if ( controlVertexDefinitions.empty() ) return;
    const double tol = 0.1 * diameter;          // loose — user-declared point must land on the nearest mesh node
    oofem::FloatArray c;
    for ( const auto &def : controlVertexDefinitions ) {
        const int id = def.first;
        const oofem::FloatArray &target = def.second;
        int bestIdx = -1;
        double bestD2 = tol * tol;
        for ( int i = 1; i <= giveNumberOfDelaunayVertices(); i++ ) {
            giveDelaunayVertex(i)->giveCoordinates(c);
            double d2 = 0.;
            for ( int k = 1; k <= 3; k++ ) {
                double dd = c.at(k) - target.at(k);
                d2 += dd * dd;
            }
            if ( d2 < bestD2 ) {
                bestD2  = d2;
                bestIdx = i;
            }
        }
        if ( bestIdx <= 0 ) {
            converter::errorf("No Delaunay vertex within tol of #@controlvertex %d at (%g, %g, %g)",
                              id, target.at(1), target.at(2), target.at(3));
        }
        controlNodeIds[ id ] = bestIdx;
    }
}


int
Grid::resolveNotchMaterial(const oofem::FloatArray &A, const oofem::FloatArray &B,
                           int defaultMat) const
{
    if ( notchSpecs.empty() ) return defaultMat;
    const double mx = 0.5 * ( A.at(1) + B.at(1) );
    const double my = 0.5 * ( A.at(2) + B.at(2) );
    const double mz = 0.5 * ( A.at(3) + B.at(3) );
    for ( const auto &n : notchSpecs ) {
        if ( mx >= n.xmin && mx <= n.xmax &&
             my >= n.ymin && my <= n.ymax &&
             mz >= n.zmin && mz <= n.zmax ) {
            return n.material;
        }
    }
    return defaultMat;
}


int
Grid::resolveInclusionMaterial(const oofem::FloatArray &A, const oofem::FloatArray &B,
                               int defaultMat) const
{
    for ( const auto &s : sphereInclusionSpecs ) {
        const double effR = s.radius + 0.5 * s.itz;
        const double d1 = std::sqrt(( A.at(1) - s.cx ) * ( A.at(1) - s.cx ) +
                                    ( A.at(2) - s.cy ) * ( A.at(2) - s.cy ) +
                                    ( A.at(3) - s.cz ) * ( A.at(3) - s.cz ) );
        const double d2 = std::sqrt(( B.at(1) - s.cx ) * ( B.at(1) - s.cx ) +
                                    ( B.at(2) - s.cy ) * ( B.at(2) - s.cy ) +
                                    ( B.at(3) - s.cz ) * ( B.at(3) - s.cz ) );
        const bool in1 = d1 < effR;
        const bool in2 = d2 < effR;
        if ( in1 && in2 ) {
            return s.inside;
        } else if ( in1 != in2 ) {
            return s.interface_;
        }
    }
    for ( const auto &c : cylinderInclusionSpecs ) {
        const double effR = c.radius + 0.5 * c.itz;
        const double ax = c.x2 - c.x1, ay = c.y2 - c.y1, az = c.z2 - c.z1;
        const double aLen2 = ax * ax + ay * ay + az * az;
        if ( aLen2 <= 0. ) continue;
        auto perp = [ & ](const oofem::FloatArray &P) {
            const double dx = P.at(1) - c.x1, dy = P.at(2) - c.y1, dz = P.at(3) - c.z1;
            const double cx = dy * az - dz * ay;
            const double cy = dz * ax - dx * az;
            const double cz = dx * ay - dy * ax;
            return std::sqrt(( cx * cx + cy * cy + cz * cz ) / aLen2);
        };
        const bool in1 = perp(A) < effR;
        const bool in2 = perp(B) < effR;
        if ( in1 && in2 ) {
            return c.inside;
        } else if ( in1 != in2 ) {
            return c.interface_;
        }
    }
    return defaultMat;
}


void Grid::writeBCRecords(std::ostream &out, int &bcID) const
{
    for (const auto &bc : bcRequests) {
        if ( bc.setID < 0 ) {
            continue;
        }

        out << "BoundaryCondition " << bcID++
            << " loadTimeFunction 1"
            << " dofs " << bc.dofs.size();
        for (int d : bc.dofs) {
            out << " " << d;
        }

        out << " values " << bc.values.size();
        for (double v : bc.values) {
            out << " " << std::scientific << v;
        }

        out << " set " << bc.setID << "\n";
    }
}


void Grid::writeGeneratedSets(std::ostream &out) const
{
    for (const auto &sd : generatedNodeSets) {
        out << "set " << sd.setID
            << " nodes " << sd.nodeIDs.size();

        for (int nid : sd.nodeIDs) {
            out << " " << nid;
        }

        out << "\n";
    }
}

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
    // Unified SM writer: handles both non-periodic (plain 3DSM) and periodic
    // (formerly 3DFPZ) cases. Periodicity is driven by `perflag`. With at
    // least one periodic axis the writer emits the periodic control node,
    // pins the first inside Delaunay vertex, and uses lattice3Dboundary +
    // periodic-image node references for boundary-crossing lines.
    const bool periodic = ( periodicityFlag.at(1) == 1 ||
                            periodicityFlag.at(2) == 1 ||
                            periodicityFlag.at(3) == 1 );

    int numberOfNodes = 0;
    int numberOfLines = 0;
    oofem::FloatArray coords(3);
    oofem::IntArray lineNodes;
    oofem::IntArray crossSectionNodes;
    oofem::IntArray location(2);

    for ( int i = 0; i < this->giveNumberOfDelaunayVertices(); i++ ) {
        int flag = this->giveDelaunayVertex(i + 1)->giveOutsideFlag();
        if ( flag == 0 || flag == 2 ) {
            numberOfNodes++;
        }
    }

    // Reinforcement (fibre) nodes are inside-only.
    for ( int i = 0; i < this->giveNumberOfReinforcementNode(); i++ ) {
        if ( this->giveReinforcementNode(i + 1)->giveOutsideFlag() == 0 ) {
            numberOfNodes++;
        }
    }

    for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
        int flag = this->giveDelaunayLine(i + 1)->giveOutsideFlag();
        if ( periodic ) {
            if ( flag == 0 || flag == 2 || flag == 3 ) numberOfLines++;
        } else {
            if ( ( flag == 0 || flag == 3 ) && this->giveDelaunayLine(i + 1)->delaunayAreaCheck() == 1 ) {
                numberOfLines++;
            }
        }
    }

    // Fibre beam segments (lattice3D / lattice3Dboundary, circular cross-section).
    for ( int i = 0; i < this->giveNumberOfLatticeBeams(); i++ ) {
        int flag = this->giveLatticeBeam(i + 1)->giveOutsideFlag();
        if ( flag == 0 || ( periodic && flag == 2 ) ) numberOfLines++;
    }

    // Fibre/matrix coupling links (latticelink3D / latticelink3Dboundary).
    for ( int i = 0; i < this->giveNumberOfLatticeLinks(); i++ ) {
        int flag = this->giveLatticeLink(i + 1)->giveOutsideFlag();
        if ( flag == 0 || ( periodic && flag == 2 ) ) numberOfLines++;
    }

    // Compact node map for non-periodic mode (renumbers inside vertices to 1..N).
    // Periodic mode uses raw vertex indices because periodic-partner references
    // returned by Vertex::givePeriodicNode() are raw indices.
    const int nDelV = this->giveNumberOfDelaunayVertices();
    std::vector< int >nodeMap(nDelV + 1, 0);
    if ( !periodic ) {
        int nodeCounter = 0;
        for ( int i = 0; i < nDelV; i++ ) {
            int flag = this->giveDelaunayVertex(i + 1)->giveOutsideFlag();
            if ( flag == 0 || flag == 2 ) {
                nodeMap [ i + 1 ] = ++nodeCounter;
            }
        }
    }
    auto mapId = [ & ](int rawId) {
        return periodic ? rawId : nodeMap [ rawId ];
    };

    // Periodic control-node id (only meaningful when periodic). Includes
    // reinforcement-node count so the fibre case still reaches a unique id.
    const int nReinf = this->giveNumberOfReinforcementNode();
    const int ctlNode = nDelV + nReinf + 1;
    const std::string ctlPlaceholder = "#@CTLNODE";

    // #@CTL<id> → resolved Delaunay-vertex id from #@controlvertex definitions.
    // Replace longer ids first so e.g. `#@CTL12` isn't consumed by the `#@CTL1` rule.
    std::vector< std::pair< std::string, int > >sortedCtlTokens;
    sortedCtlTokens.reserve(controlNodeIds.size());
    for ( const auto &kv : controlNodeIds ) {
        // Translate raw Delaunay index through the compact nodeMap in non-periodic mode;
        // periodic mode uses raw ids so mapId is identity there.
        sortedCtlTokens.emplace_back("#@CTL" + std::to_string(kv.first), mapId(kv.second));
    }
    std::sort(sortedCtlTokens.begin(), sortedCtlTokens.end(),
              [](const auto &a, const auto &b) { return a.first.size() > b.first.size(); });

    auto replaceAll = [](std::string &s, const std::string &needle, const std::string &replacement) {
        if ( needle.empty() ) return;
        size_t pos = 0;
        while ( ( pos = s.find(needle, pos) ) != std::string::npos ) {
            s.replace(pos, needle.size(), replacement);
            pos += replacement.size();
        }
    };

    auto substitute = [ & ](std::string &s) {
        // #@CTLNODE (periodic control node) — periodic mode only.
        if ( periodic ) {
            replaceAll(s, ctlPlaceholder, std::to_string(ctlNode));
        }
        // #@CTL<id> (mesh control vertices from #@controlvertex).
        for ( const auto &tok : sortedCtlTokens ) {
            replaceAll(s, tok.first, std::to_string(tok.second));
        }
    };

    // Specimen geometry for the control-node coordinates and (non-periodic) midpoint closure.
    oofem::FloatArray bounds;
    this->giveRegion(1)->defineBoundaries(bounds);
    oofem::FloatArray specimenDim(3);
    specimenDim.at(1) = bounds.at(2) - bounds.at(1);
    specimenDim.at(2) = bounds.at(4) - bounds.at(3);
    specimenDim.at(3) = bounds.at(6) - bounds.at(5);
    const double tol = this->giveTol();

    auto faceMask = [ & ](const oofem::FloatArray &c) {
        int m = 0;
        if ( std::abs(c.at(1) - bounds.at(1)) < tol ) m |= 1;
        if ( std::abs(c.at(1) - bounds.at(2)) < tol ) m |= 2;
        if ( std::abs(c.at(2) - bounds.at(3)) < tol ) m |= 4;
        if ( std::abs(c.at(2) - bounds.at(4)) < tol ) m |= 8;
        if ( std::abs(c.at(3) - bounds.at(5)) < tol ) m |= 16;
        if ( std::abs(c.at(3) - bounds.at(6)) < tol ) m |= 32;
        return m;
    };

    std::ifstream ctrl(controlFileName);
    std::ofstream out(fileName);
    if ( !ctrl ) {
        converter::error("give3DSMOutput: Cannot open control file");
    }
    if ( !out ) {
        converter::error("give3DSMOutput: Cannot open output file");
    }

    std::string line;
    bool injected = false;

    while ( std::getline(ctrl, line) ) {
        std::string t = line;
        size_t pos = t.find_first_not_of(" \t");
        if ( pos != std::string::npos ) {
            t.erase(0, pos);
        } else {
            t.clear();
        }

        if ( !injected && t.rfind("ncrosssect", 0) == 0 ) {
            std::istringstream iss(t);
            std::string token;
            out << "ndofman " << ( numberOfNodes + ( periodic ? 1 : 0 ) )
                << " nelem " << numberOfLines << " ";
            while ( iss >> token ) {
                out << token << " ";
            }
            out << "\n";

            // Delaunay nodes (inside or on boundary). In periodic mode, the
            // first inside vertex is pinned and raw vertex IDs are used.
            // Vertices matching a #@rigidarm face are emitted as
            // rigidarmnode slaved to the master controlvertex (except the
            // master itself, which is kept as a regular node).
            int compactCounter = 0;
            int firstFlag = 0;
            for ( int i = 0; i < nDelV; i++ ) {
                int flag = this->giveDelaunayVertex(i + 1)->giveOutsideFlag();
                if ( flag != 0 && flag != 2 ) continue;

                this->giveDelaunayVertex(i + 1)->giveCoordinates(coords);
                int nid = periodic ? ( i + 1 ) : ( ++compactCounter );

                // Check whether this vertex lands on any #@rigidarm face.
                const RigidArmSpec *raMatch = nullptr;
                for ( const auto &ra : rigidArmSpecs ) {
                    double faceCoord = ra.sideMax ? bounds.at(2 * ra.axis) : bounds.at(2 * ra.axis - 1);
                    if ( std::abs(coords.at(ra.axis) - faceCoord) >= tol ) continue;
                    auto it = controlNodeIds.find(ra.masterCtlId);
                    if ( it == controlNodeIds.end() ) continue;
                    int masterRaw = it->second;
                    if ( masterRaw == i + 1 ) continue;  // the master itself stays a regular node
                    raMatch = &ra;
                    break;
                }

                if ( raMatch ) {
                    auto it = controlNodeIds.find(raMatch->masterCtlId);
                    int masterNid = mapId(it->second);
                    out << "rigidarmnode " << nid << " coords 3 " << std::scientific
                        << coords.at(1) << " " << coords.at(2) << " " << coords.at(3)
                        << " master " << masterNid
                        << " mastermask 6";
                    for ( int k = 1; k <= 6; k++ ) out << " " << raMatch->mastermask.at(k);
                    out << " doftype 6";
                    for ( int k = 1; k <= 6; k++ ) out << " " << raMatch->doftype.at(k);
                    out << "\n";
                    continue;
                }

                out << "node " << nid << " coords 3 " << std::scientific
                    << coords.at(1) << " " << coords.at(2) << " " << coords.at(3);
                if ( periodic && firstFlag == 0 ) {
                    firstFlag = 1;
                    out << " bc 6 1 1 1 1 1 1";
                }
                out << "\n";
            }

            // Reinforcement (fibre) nodes — id = nDelV + reinforcement-index.
            for ( int i = 0; i < nReinf; i++ ) {
                if ( this->giveReinforcementNode(i + 1)->giveOutsideFlag() != 0 ) continue;
                this->giveReinforcementNode(i + 1)->giveCoordinates(coords);
                out << "node " << ( nDelV + i + 1 ) << " coords 3 " << std::scientific
                    << coords.at(1) << " " << coords.at(2) << " " << coords.at(3) << "\n";
            }

            // Periodic control node.
            if ( periodic ) {
                out << "node " << ctlNode << " coords 3 " << std::scientific
                    << specimenDim.at(1) << " " << specimenDim.at(2) << " " << specimenDim.at(3)
                    << " load 1 2\n";
            }

            int elemCounter = 0;

            // Helper: starting from the default material 1, apply the
            // #@material_around overrides (element touches a named
            // controlvertex), then #@notch (midpoint inside a box), then
            // #@sphere/cylinderinclusion. Later stages can override earlier
            // ones.
            auto resolveLineMaterial = [ & ](int endpoint1, int endpoint2,
                                             const oofem::FloatArray &Acoord,
                                             const oofem::FloatArray &Bcoord) -> int {
                int mat = 1;
                for ( const auto &ma : materialAroundSpecs ) {
                    auto it = controlNodeIds.find(ma.ctlId);
                    if ( it == controlNodeIds.end() ) continue;
                    if ( endpoint1 == it->second || endpoint2 == it->second ) {
                        mat = ma.material;
                        break;
                    }
                }
                mat = resolveNotchMaterial(Acoord, Bcoord, mat);
                mat = resolveInclusionMaterial(Acoord, Bcoord, mat);
                return mat;
            };

            // Delaunay lines → lattice3D (inside) or lattice3Dboundary (periodic crossing).
            oofem::FloatArray A(3), B(3), M(3), cCurr(3), cNext(3);
            for ( int i = 0; i < this->giveNumberOfDelaunayLines(); i++ ) {
                int flag = this->giveDelaunayLine(i + 1)->giveOutsideFlag();

                const bool emitInside = ( flag == 0 || flag == 3 ) &&
                    ( periodic || this->giveDelaunayLine(i + 1)->delaunayAreaCheck() == 1 );
                const bool emitBoundary = periodic && flag == 2;
                if ( !emitInside && !emitBoundary ) continue;

                this->giveDelaunayLine(i + 1)->giveLocalVertices(lineNodes);
                this->giveDelaunayLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

                // Midpoint-based notch check uses the raw (pre-substitution) endpoints.
                oofem::FloatArray boundaryA(3), boundaryB(3);
                this->giveDelaunayVertex(lineNodes.at(1))->giveCoordinates(boundaryA);
                this->giveDelaunayVertex(lineNodes.at(2))->giveCoordinates(boundaryB);
                int matBoundary = resolveLineMaterial(lineNodes.at(1), lineNodes.at(2),
                                                      boundaryA, boundaryB);
                auto boundaryBodyloadIt = bodyloadByMaterial.find(matBoundary);

                if ( emitBoundary ) {
                    location.zero();
                    for ( int m = 0; m < 2; m++ ) {
                        if ( this->giveDelaunayVertex(lineNodes.at(m + 1))->giveOutsideFlag() == 1 ) {
                            location.at(m + 1) = this->giveDelaunayVertex(lineNodes.at(m + 1))->giveLocation();
                            lineNodes.at(m + 1) = this->giveDelaunayVertex(lineNodes.at(m + 1))->givePeriodicNode();
                        }
                    }

                    out << "lattice3Dboundary " << ++elemCounter
                        << " nodes 3 " << lineNodes.at(1) << " " << lineNodes.at(2) << " " << ctlNode
                        << " crossSect " << matBoundary << " mat " << matBoundary
                        << " polycoords " << ( 3 * crossSectionNodes.giveSize() );
                    for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                        this->giveVoronoiVertex(crossSectionNodes.at(m + 1))->giveCoordinates(coords);
                        out << " " << coords.at(1) << " " << coords.at(2) << " " << coords.at(3);
                    }
                    if ( emitCouplingFlag ) {
                        oofem::IntArray crossSectionElements;
                        this->giveDelaunayLine(i + 1)->giveCrossSectionElements(crossSectionElements);
                        out << " couplingflag 1 couplingnumber " << crossSectionElements.giveSize();
                        for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                            int vid = crossSectionElements.at(m + 1);
                            int mapped = ( this->giveVoronoiLine(vid)->giveOutsideFlag() == 1 )
                                ? this->giveVoronoiLine(vid)->givePeriodicElement()
                                : vid;
                            out << " " << mapped;
                        }
                    }
                    if ( boundaryBodyloadIt != bodyloadByMaterial.end() ) {
                        out << " bodyloads 1 " << boundaryBodyloadIt->second;
                    }
                    out << " location 2 " << location.at(1) << " " << location.at(2) << "\n";
                    continue;
                }

                // Inside element. In non-periodic mode, apply midpoint closure when
                // the cross-section polygon wraps around a specimen edge.
                this->giveDelaunayVertex(lineNodes.at(1))->giveCoordinates(A);
                this->giveDelaunayVertex(lineNodes.at(2))->giveCoordinates(B);
                M.at(1) = 0.5 * ( A.at(1) + B.at(1) );
                M.at(2) = 0.5 * ( A.at(2) + B.at(2) );
                M.at(3) = 0.5 * ( A.at(3) + B.at(3) );

                const int nPoly = crossSectionNodes.giveSize();
                std::vector< oofem::FloatArray >polyOut;
                polyOut.reserve(nPoly + 3);
                for ( int m = 0; m < nPoly; m++ ) {
                    this->giveVoronoiVertex(crossSectionNodes.at(m + 1))->giveCoordinates(cCurr);
                    polyOut.push_back(cCurr);

                    if ( !periodic ) {
                        int next = ( m + 1 ) % nPoly;
                        this->giveVoronoiVertex(crossSectionNodes.at(next + 1))->giveCoordinates(cNext);
                        int mCurr = faceMask(cCurr);
                        int mNext = faceMask(cNext);
                        if ( mCurr != 0 && mNext != 0 && ( mCurr & mNext ) == 0 ) {
                            polyOut.push_back(M);
                        }
                    }
                }

                int matInside = resolveLineMaterial(lineNodes.at(1), lineNodes.at(2), A, B);
                auto insideBodyloadIt = bodyloadByMaterial.find(matInside);
                out << "lattice3D " << ++elemCounter
                    << " nodes 2 " << mapId(lineNodes.at(1)) << " " << mapId(lineNodes.at(2))
                    << " crossSect " << matInside << " mat " << matInside
                    << " polycoords " << 3 * ( int ) polyOut.size();
                for ( const auto &c : polyOut ) {
                    out << " " << c.at(1) << " " << c.at(2) << " " << c.at(3);
                }
                if ( emitCouplingFlag ) {
                    oofem::IntArray crossSectionElements;
                    this->giveDelaunayLine(i + 1)->giveCrossSectionElements(crossSectionElements);
                    out << " couplingflag 1 couplingnumber " << crossSectionElements.giveSize();
                    for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                        int vid = crossSectionElements.at(m + 1);
                        int mapped = ( this->giveVoronoiLine(vid)->giveOutsideFlag() == 1 )
                            ? this->giveVoronoiLine(vid)->givePeriodicElement()
                            : vid;
                        out << " " << mapped;
                    }
                }
                if ( insideBodyloadIt != bodyloadByMaterial.end() ) {
                    out << " bodyloads 1 " << insideBodyloadIt->second;
                }
                out << "\n";
            }

            // Fibre beam segments: lattice3D / lattice3Dboundary with circular
            // cross-section (defined by latticecs 2 in the control file).
            for ( int i = 0; i < this->giveNumberOfLatticeBeams(); i++ ) {
                int flag = this->giveLatticeBeam(i + 1)->giveOutsideFlag();
                const bool emitInside = ( flag == 0 );
                const bool emitBoundary = periodic && flag == 2;
                if ( !emitInside && !emitBoundary ) continue;

                this->giveLatticeBeam(i + 1)->giveLocalVertices(lineNodes);

                if ( emitBoundary ) {
                    location.zero();
                    for ( int m = 0; m < 2; m++ ) {
                        if ( this->giveReinforcementNode(lineNodes.at(m + 1))->giveOutsideFlag() == 1 ) {
                            location.at(m + 1) = this->giveReinforcementNode(lineNodes.at(m + 1))->giveLocation();
                            lineNodes.at(m + 1) = this->giveReinforcementNode(lineNodes.at(m + 1))->givePeriodicNode();
                        }
                    }
                    out << "lattice3Dboundary " << ++elemCounter
                        << " nodes 3 " << ( nDelV + lineNodes.at(1) ) << " " << ( nDelV + lineNodes.at(2) )
                        << " " << ctlNode << " crossSect 2 mat 2"
                        << " location 2 " << location.at(1) << " " << location.at(2) << "\n";
                } else {
                    out << "lattice3D " << ++elemCounter
                        << " nodes 2 " << ( nDelV + lineNodes.at(1) ) << " " << ( nDelV + lineNodes.at(2) )
                        << " crossSect 2 mat 2\n";
                }
            }

            // Fibre/matrix coupling links: latticelink3D / latticelink3Dboundary.
            // Endpoint 1 is a reinforcement node (offset by nDelV);
            // endpoint 2 is a matrix Delaunay vertex (raw id).
            for ( int i = 0; i < this->giveNumberOfLatticeLinks(); i++ ) {
                int flag = this->giveLatticeLink(i + 1)->giveOutsideFlag();
                const bool emitInside = ( flag == 0 );
                const bool emitBoundary = periodic && flag == 2;
                if ( !emitInside && !emitBoundary ) continue;

                this->giveLatticeLink(i + 1)->giveLocalVertices(lineNodes);
                const double linkLen = this->giveLatticeLink(i + 1)->giveAssociatedLength();
                const double linkDiam = this->giveLatticeLink(i + 1)->giveDiameter();
                const double linkLend = this->giveLatticeLink(i + 1)->giveEndLength();
                oofem::FloatArray linkDir = this->giveLatticeLink(i + 1)->giveDirectionVector();

                if ( emitBoundary ) {
                    location.zero();
                    if ( this->giveReinforcementNode(lineNodes.at(1))->giveOutsideFlag() == 1 ) {
                        location.at(1) = this->giveReinforcementNode(lineNodes.at(1))->giveLocation();
                        lineNodes.at(1) = this->giveReinforcementNode(lineNodes.at(1))->givePeriodicNode();
                    }
                    if ( this->giveDelaunayVertex(lineNodes.at(2))->giveOutsideFlag() == 1 ) {
                        location.at(2) = this->giveDelaunayVertex(lineNodes.at(2))->giveLocation();
                        lineNodes.at(2) = this->giveDelaunayVertex(lineNodes.at(2))->givePeriodicNode();
                    }
                    out << "latticelink3Dboundary " << ++elemCounter
                        << " nodes 3 " << ( nDelV + lineNodes.at(1) ) << " " << lineNodes.at(2)
                        << " " << ctlNode << " crossSect 3 mat 3"
                        << " length " << linkLen << " diameter " << linkDiam
                        << " dirvector 3 " << linkDir.at(1) << " " << linkDir.at(2) << " " << linkDir.at(3)
                        << " L_end " << linkLend
                        << " location 2 " << location.at(1) << " " << location.at(2) << "\n";
                } else {
                    out << "latticelink3D " << ++elemCounter
                        << " nodes 2 " << ( nDelV + lineNodes.at(1) ) << " " << lineNodes.at(2)
                        << " crossSect 3 mat 3"
                        << " length " << linkLen << " diameter " << linkDiam
                        << " dirvector 3 " << linkDir.at(1) << " " << linkDir.at(2) << " " << linkDir.at(3)
                        << " L_end " << linkLend << "\n";
                }
            }

            injected = true;
            continue;
        }

        if ( !isConverterDirective(t) ) {
            std::string written = line;
            substitute(written);
            out << written << "\n";
        }
    }
}


void
Grid::give3DTMOutput(const std::string &fileName)
{
    // Unified TM writer: handles both non-periodic (the original 3DTMQhull
    // case) and periodic cases. Mirrors give3DSMOutput's periodic structure.
    // In non-periodic mode compact Voronoi-vertex ids (1..N) are used; in
    // periodic mode raw ids are used because periodic-partner references
    // returned by Vertex::givePeriodicNode() are raw indices.
    const bool periodic = ( periodicityFlag.at(1) == 1 ||
                            periodicityFlag.at(2) == 1 ||
                            periodicityFlag.at(3) == 1 );

    oofem::FloatArray coords(3);
    oofem::IntArray nodes;
    oofem::IntArray crossSectionNodes;
    oofem::IntArray location(2);

    // Compact Voronoi-vertex id map for non-periodic mode (1..N over
    // flag==0 | flag==2 vertices). Periodic mode uses raw ids.
    const int nVorV = this->giveNumberOfVoronoiVertices();
    std::unordered_map< int, int >nodeIdMap;
    int numberOfNodes = 0;
    for ( int i = 0; i < nVorV; i++ ) {
        int flag = this->giveVoronoiVertex(i + 1)->giveOutsideFlag();
        if ( flag == 0 || flag == 2 ) {
            if ( !periodic ) {
                nodeIdMap [ i + 1 ] = ++numberOfNodes;
            } else {
                ++numberOfNodes;
            }
        }
    }
    auto mapId = [ & ](int rawId) {
        return periodic ? rawId : nodeIdMap.at(rawId);
    };

    int numberOfLines = 0;
    for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
        int flag = this->giveVoronoiLine(i + 1)->giveOutsideFlag();
        if ( flag == 0 ) {
            numberOfLines++;
        } else if ( periodic && flag == 2 ) {
            numberOfLines++;
        }
    }

    // Periodic control node id: one past the last raw Voronoi vertex.
    const int ctlNode = nVorV + 1;
    const std::string ctlPlaceholder = "#@CTLNODE";

    auto replaceAll = [](std::string &s, const std::string &needle, const std::string &replacement) {
        if ( needle.empty() ) return;
        size_t pos = 0;
        while ( ( pos = s.find(needle, pos) ) != std::string::npos ) {
            s.replace(pos, needle.size(), replacement);
            pos += replacement.size();
        }
    };

    auto substitute = [ & ](std::string &s) {
        if ( periodic ) {
            replaceAll(s, ctlPlaceholder, std::to_string(ctlNode));
        }
    };

    std::ifstream ctrl(controlFileName);
    std::ofstream out(fileName);
    if ( !ctrl ) {
        converter::error("give3DTMOutput: Cannot open control file");
    }
    if ( !out ) {
        converter::error("give3DTMOutput: Cannot open output file");
    }

    std::string line_s;
    bool injected = false;

    while ( std::getline(ctrl, line_s) ) {
        std::string t = line_s;
        size_t pos = t.find_first_not_of(" \t");
        if ( pos != std::string::npos ) {
            t.erase(0, pos);
        } else {
            t.clear();
        }

        if ( !injected && t.rfind("ncrosssect", 0) == 0 ) {
            std::istringstream iss(t);
            std::string token;
            out << "ndofman " << ( numberOfNodes + ( periodic ? 1 : 0 ) )
                << " nelem " << numberOfLines << " ";
            while ( iss >> token ) {
                out << token << " ";
            }
            out << "\n";

            // write Voronoi nodes. In periodic mode the first emitted vertex
            // is pinned (bc 1 1) to remove the transport rigid-body mode.
            int firstFlag = 0;
            for ( int i = 0; i < nVorV; i++ ) {
                int flag = this->giveVoronoiVertex(i + 1)->giveOutsideFlag();
                if ( flag != 0 && flag != 2 ) continue;
                this->giveVoronoiVertex(i + 1)->giveCoordinates(coords);
                int nid = periodic ? ( i + 1 ) : nodeIdMap.at(i + 1);
                out << "node " << nid << " coords 3 " << std::scientific
                    << coords.at(1) << " " << coords.at(2) << " " << coords.at(3);
                if ( periodic && firstFlag == 0 ) {
                    firstFlag = 1;
                    out << " bc 1 1";
                }
                out << "\n";
            }

            // Periodic control node — 3 DOFs (IDs 1,2,3) with the Wong-style
            // DOF-BC mapping "bc 3 1 1 2" (first two DOFs constrained by BC 1,
            // the last by BC 2 — tuned to the percolation loading setup).
            oofem::FloatArray bounds;
            this->giveRegion(1)->defineBoundaries(bounds);
            oofem::FloatArray specimenDim(3);
            specimenDim.at(1) = bounds.at(2) - bounds.at(1);
            specimenDim.at(2) = bounds.at(4) - bounds.at(3);
            specimenDim.at(3) = bounds.at(6) - bounds.at(5);
            if ( periodic ) {
                out << "node " << ctlNode << " coords 3 " << std::scientific
                    << specimenDim.at(1) << " " << specimenDim.at(2) << " " << specimenDim.at(3)
                    << " ndofs 3 dofIDmask 3 1 2 3 bc 3 1 1 2\n";
            }

            // Voronoi elements: latticemt3D (inside) or latticemt3Dboundary
            // (periodic crossing). Cross-section vertices are Delaunay
            // vertices (dual mesh). Material is resolved via notch +
            // inclusion directives using the endpoint coords.
            int elemCounter = 0;
            oofem::FloatArray A(3), B(3);
            for ( int i = 0; i < this->giveNumberOfVoronoiLines(); i++ ) {
                int flag = this->giveVoronoiLine(i + 1)->giveOutsideFlag();
                const bool emitInside = ( flag == 0 );
                const bool emitBoundary = periodic && flag == 2;
                if ( !emitInside && !emitBoundary ) continue;

                this->giveVoronoiLine(i + 1)->giveLocalVertices(nodes);
                this->giveVoronoiLine(i + 1)->giveCrossSectionVertices(crossSectionNodes);

                // Pre-substitution endpoint coords for material resolution.
                oofem::FloatArray cA(3), cB(3);
                this->giveVoronoiVertex(nodes.at(1))->giveCoordinates(cA);
                this->giveVoronoiVertex(nodes.at(2))->giveCoordinates(cB);
                int mat = resolveNotchMaterial(cA, cB, 1);
                mat = resolveInclusionMaterial(cA, cB, mat);
                auto bodyloadIt = bodyloadByMaterial.find(mat);

                if ( emitBoundary ) {
                    location.zero();
                    for ( int m = 0; m < 2; m++ ) {
                        if ( this->giveVoronoiVertex(nodes.at(m + 1))->giveOutsideFlag() == 1 ) {
                            location.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1))->giveLocation();
                            nodes.at(m + 1) = this->giveVoronoiVertex(nodes.at(m + 1))->givePeriodicNode();
                        }
                    }

                    out << "latticemt3Dboundary " << ++elemCounter
                        << " nodes 3 " << nodes.at(1) << " " << nodes.at(2) << " " << ctlNode
                        << " crossSect " << mat << " mat " << mat
                        << " polycoords " << 3 * crossSectionNodes.giveSize();
                    for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                        this->giveDelaunayVertex(crossSectionNodes.at(m + 1))->giveCoordinates(coords);
                        out << " " << coords.at(1) << " " << coords.at(2) << " " << coords.at(3);
                    }
                    if ( emitCouplingFlag ) {
                        oofem::IntArray crossSectionElements;
                        this->giveVoronoiLine(i + 1)->giveCrossSectionElements(crossSectionElements);
                        out << " couplingflag 1 couplingnumber " << crossSectionElements.giveSize();
                        for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                            int did = crossSectionElements.at(m + 1);
                            int mapped = ( this->giveDelaunayLine(did)->giveOutsideFlag() == 1 )
                                ? this->giveDelaunayLine(did)->givePeriodicElement()
                                : did;
                            out << " " << mapped;
                        }
                    }
                    if ( bodyloadIt != bodyloadByMaterial.end() ) {
                        out << " bodyloads 1 " << bodyloadIt->second;
                    }
                    out << " location 2 " << location.at(1) << " " << location.at(2) << "\n";
                    continue;
                }

                out << "latticemt3D " << ++elemCounter
                    << " nodes 2 " << mapId(nodes.at(1)) << " " << mapId(nodes.at(2))
                    << " crossSect " << mat << " mat " << mat
                    << " polycoords " << 3 * crossSectionNodes.giveSize();
                for ( int m = 0; m < crossSectionNodes.giveSize(); m++ ) {
                    this->giveDelaunayVertex(crossSectionNodes.at(m + 1))->giveCoordinates(coords);
                    out << " " << coords.at(1) << " " << coords.at(2) << " " << coords.at(3);
                }
                if ( emitCouplingFlag ) {
                    oofem::IntArray crossSectionElements;
                    this->giveVoronoiLine(i + 1)->giveCrossSectionElements(crossSectionElements);
                    out << " couplingflag 1 couplingnumber " << crossSectionElements.giveSize();
                    for ( int m = 0; m < crossSectionElements.giveSize(); m++ ) {
                        int did = crossSectionElements.at(m + 1);
                        int mapped = ( this->giveDelaunayLine(did)->giveOutsideFlag() == 1 )
                            ? this->giveDelaunayLine(did)->givePeriodicElement()
                            : did;
                        out << " " << mapped;
                    }
                }
                if ( bodyloadIt != bodyloadByMaterial.end() ) {
                    out << " bodyloads 1 " << bodyloadIt->second;
                }
                out << "\n";
            }

            injected = true;
            continue;
        }

        if ( !isConverterDirective(t) ) {
            std::string written = line_s;
            substitute(written);
            out << written << "\n";
        }
    }
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

    printf("Lattice beams: %d\n", this->giveNumberOfLatticeBeams());

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

oofem::FloatArray Grid::triBarycentre(int triIndex) const
{
    const Tri &t = tris [ triIndex ];

    oofem::FloatArray x1 = getX(t.n1);
    oofem::FloatArray x2 = getX(t.n2);
    oofem::FloatArray x3 = getX(t.n3);

    oofem::FloatArray b(3);
    b.at(1) = ( x1.at(1) + x2.at(1) + x3.at(1) ) / 3.0;
    b.at(2) = ( x1.at(2) + x2.at(2) + x3.at(2) ) / 3.0;
    b.at(3) = ( x1.at(3) + x2.at(3) + x3.at(3) ) / 3.0;

    return b;
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

    printf("Lattice links: %d\n", this->giveNumberOfLatticeLinks());

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

//****************************************************************************80



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
