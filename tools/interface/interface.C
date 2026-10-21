// interface.C
//
// Post-processor that inserts node-to-node contact (latticecontact3d) along the
// interfaces between bodies of an existing OOFEM continuum mesh.
//
// Bodies are distinguished by their cross-section id (region number in OOFEM:
// Element::giveRegionNumber() == giveCrossSection()->giveNumber()). An interface
// is the set of tetrahedron faces shared by two tets of different bodies. Along
// the interfaces the tool:
//   1. duplicates the shared interface nodes (as real nodes, not hanging nodes),
//      making one copy per incident body so that a node where k bodies meet is
//      split k ways (2 for an ordinary bar/concrete face, 3 along the triple line
//      where two tied bars and the concrete meet, ...),
//   2. reassigns each body's tets to that body's copy (the "unzip"); the body
//      listed first keeps the original node id,
//   3. inserts one latticecontact3d per adjacent body-pair per node, with the
//      tributary area and area-weighted normal taken from that pair's interface
//      faces. Every internal face separates exactly one pair, so the tributary
//      area at a triple line is partitioned among the pairs with no double count.
//
// Preconditions:
//   * The mesh is conforming across the interfaces (adjacent bodies share the
//     interface nodes) so node-to-node pairing is well defined.
//   * The contact cross-section (LatticeCS) and material (e.g. latticebondplast)
//     already exist in the input mesh; pass their ids. The tool does not add
//     cross-sections or materials, only nodes, contact elements and one set.
//
// Usage (internal, k bodies):
//   interface in.in out.in --bodies c1,c2,...  --contact-cs C --contact-mat M
//             [--pair P:Q[:cs:mat] ...] [--no-set]
//     Without --pair, every inter-body pair gets a contact with the default
//     --contact-cs/--contact-mat. With one or more --pair, only the listed pairs
//     get contacts (nodes are still split by full incidence); P:Q:cs:mat overrides
//     the cross-section/material for that pair.
//
// Usage (internal, legacy 2-body form, equivalent to --bodies A,B):
//   interface in.in out.in --sideA-cs A --sideB-cs B --contact-cs C --contact-mat M [--no-set]

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <array>
#include <algorithm>
#include <cmath>
#include <cctype>

using namespace std;

struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3() {}
    Vec3(double a, double b, double c) : x(a), y(b), z(c) {}
    Vec3 operator-(const Vec3 &o) const { return { x - o.x, y - o.y, z - o.z }; }
    Vec3 operator+(const Vec3 &o) const { return { x + o.x, y + o.y, z + o.z }; }
    Vec3 operator*(double s) const { return { x * s, y * s, z * s }; }
    Vec3 operator/(double s) const { return { x / s, y / s, z / s }; }
    double dot(const Vec3 &o) const { return x * o.x + y * o.y + z * o.z; }
    Vec3 cross(const Vec3 &o) const { return { y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x }; }
    double norm() const { return sqrt(x * x + y * y + z * z); }
    Vec3 normalized() const { double n = norm(); return n == 0 ? Vec3() : (*this) / n; }
};

static string lower(string s) { for ( char &c : s ) c = (char) tolower((unsigned char) c); return s; }

static vector<string> tokenize(const string &line)
{
    vector<string> t;
    istringstream iss(line);
    string w;
    while ( iss >> w ) t.push_back(w);
    return t;
}

static string join(const vector<string> &t)
{
    string s;
    for ( size_t i = 0; i < t.size(); ++i ) { if ( i ) s += " "; s += t[i]; }
    return s;
}

// Value following a keyword in a token list (case-insensitive key match).
static int tokenValueAfter(const vector<string> &t, const string &key, int fallback = -1)
{
    for ( size_t i = 0; i + 1 < t.size(); ++i )
        if ( lower(t[i]) == lower(key) ) return stoi(t[i + 1]);
    return fallback;
}

// Contiguous element-id range formatting: {(1 12) (20 24)}
static string formatRanges(vector<int> ids)
{
    if ( ids.empty() ) return "";
    sort(ids.begin(), ids.end());
    ostringstream ss;
    ss << "elementranges {";
    int start = ids[0], prev = ids[0];
    for ( size_t i = 1; i < ids.size(); ++i ) {
        if ( ids[i] == prev + 1 ) prev = ids[i];
        else { ss << "(" << start << " " << prev << ") "; start = prev = ids[i]; }
    }
    ss << "(" << start << " " << prev << ")}";
    return ss.str();
}

static double coordAxis(const Vec3 &p, int axis) { return axis == 1 ? p.x : ( axis == 2 ? p.y : p.z ); }

// Remove a "bc <ndof> <v1..vndof>" clause from a node token list, if present.
static vector<string> stripBC(const vector<string> &t)
{
    vector<string> r;
    for ( size_t i = 0; i < t.size(); ++i ) {
        if ( lower(t[i]) == "bc" && i + 1 < t.size() ) {
            int nd = stoi(t[i + 1]);
            i += 1 + nd;               // skip "bc", the count, and nd flags
            continue;
        }
        r.push_back(t[i]);
    }
    return r;
}

// Node ids from a Set line: supports "nodes <n> <ids>" and "noderanges { (a b) c ... }".
static vector<int> parseSetNodes(const string &line)
{
    vector<int> ids;
    string low = lower(line);
    size_t p = low.find(" nodes ");
    if ( p != string::npos ) {
        istringstream iss(line.substr(p + 7));
        int cnt;
        if ( iss >> cnt ) { int v; for ( int j = 0; j < cnt && ( iss >> v ); ++j ) ids.push_back(v); }
        return ids;
    }
    p = low.find("noderanges");
    if ( p != string::npos ) {
        size_t b1 = line.find('{', p), b2 = line.find('}', p);
        if ( b1 != string::npos && b2 != string::npos && b2 > b1 ) {
            string blob = line.substr(b1 + 1, b2 - b1 - 1);
            for ( size_t i = 0; i < blob.size(); ) {
                if ( blob[i] == '(' ) {
                    size_t close = blob.find(')', i);
                    if ( close == string::npos ) break;
                    istringstream iss(blob.substr(i + 1, close - i - 1));
                    int a, b; if ( iss >> a >> b ) for ( int v = a; v <= b; ++v ) ids.push_back(v);
                    i = close + 1;
                } else if ( isdigit((unsigned char) blob[i]) ) {
                    size_t j = i; while ( j < blob.size() && isdigit((unsigned char) blob[j]) ) ++j;
                    ids.push_back(stoi(blob.substr(i, j - i)));
                    i = j;
                } else ++i;
            }
        }
    }
    return ids;
}

struct TetRec {
    int id;
    vector<string> tokens;   // full element line, tokens; node-id tokens rewritten on unzip
    int nodeTokIdx[4];       // token indices of the four node ids
    int n[4];                // current node ids
    int cs;                  // cross-section id
};

int main(int argc, char **argv)
{
    string inFile, outFile;
    int csA = -1, csB = -1, contactCS = -1, contactMat = -1;
    bool addSet = true;
    // free-surface mode
    bool freeMode = false, transferBC = true;
    bool bcToDup = false;   // internal: at a split node keep the prescribed BC on the
                            // duplicate(s) (the unzipped, non-rank-0 bodies) and strip it
                            // from the kept original -- e.g. keep a loaded/supported bar's
                            // BC on the bar, not on the concrete it was unzipped from
    int  planeAxis = 0, setId = -1;   // planeAxis 1/2/3 = x/y/z; 0 = unset
    double planeVal = 0., planeTol = 1e-6;
    // internal k-way mode
    vector<int> bodies;                                  // participating cross-sections, first = keeps original
    struct PairSpec { int p, q, cs, mat; };              // p,q as given (order-normalized later)
    vector<PairSpec> pairSpecs;                          // explicit --pair list (empty = all pairs)
    double splitAngleDeg = -1.0;                         // <0 = one contact per node (blended normal);
                                                         // >=0 = split a node's faces into clusters whose
                                                         // normals differ by more than this angle (e.g. a
                                                         // mantle/end-cap rim -> a radial + an axial link)
    // geometry-based contact-material zones: a contact whose interface-node coord
    // on `axis` (1/2/3 = x/y/z) lies in [lo,hi] takes cross-section `cs`, material
    // `mat` instead of the pair default -- e.g. a frictionless debonded sleeve over
    // a z-band. First matching zone wins; mesh-size independent (keyed on coords).
    struct ZoneMat { int axis; double lo, hi; int cs, mat; };
    vector<ZoneMat> zoneMats;
    vector<string> pos;

    // "3,7,9" or "3 7 9" -> {3,7,9}
    auto parseIntList = [](string s) {
        for ( char &c : s ) if ( c == ',' ) c = ' ';
        vector<int> v; istringstream iss(s); int x; while ( iss >> x ) v.push_back(x); return v;
    };

    for ( int i = 1; i < argc; ++i ) {
        string a = argv[i];
        auto need = [&](const char *nm) -> string {
            if ( i + 1 >= argc ) { cerr << "Missing value for " << nm << "\n"; exit(1); }
            return argv[++i];
        };
        if      ( a == "--sideA-cs" )        csA        = stoi(need("--sideA-cs"));
        else if ( a == "--sideB-cs" )        csB        = stoi(need("--sideB-cs"));
        else if ( a == "--bodies" )          bodies     = parseIntList(need("--bodies"));
        else if ( a == "--contact-cs" )      contactCS  = stoi(need("--contact-cs"));
        else if ( a == "--contact-mat" )     contactMat = stoi(need("--contact-mat"));
        else if ( a == "--pair" ) {
            // P:Q  or  P:Q:cs:mat
            string s = need("--pair");
            for ( char &c : s ) if ( c == ':' ) c = ' ';
            istringstream iss(s);
            PairSpec ps { -1, -1, -1, -1 };
            iss >> ps.p >> ps.q;
            iss >> ps.cs >> ps.mat;          // optional; leave -1 to fall back to defaults
            if ( ps.p < 0 || ps.q < 0 ) { cerr << "--pair needs P:Q[:cs:mat]\n"; return 1; }
            pairSpecs.push_back(ps);
        }
        else if ( a == "--split-angle" )     splitAngleDeg = stod(need("--split-angle"));
        else if ( a == "--zone-mat" ) {
            string ax = lower(need("--zone-mat"));
            int axis = ( ax == "x" || ax == "1" ) ? 1 : ( ax == "y" || ax == "2" ) ? 2 : ( ax == "z" || ax == "3" ) ? 3 : 0;
            if ( !axis ) { cerr << "--zone-mat axis must be x/y/z\n"; return 1; }
            double lo = stod(need("--zone-mat")), hi = stod(need("--zone-mat"));
            int zcs = stoi(need("--zone-mat")), zmat = stoi(need("--zone-mat"));
            zoneMats.push_back({ axis, lo, hi, zcs, zmat });
        }
        else if ( a == "--bc-to-dup" )       bcToDup    = true;
        else if ( a == "--no-set" )          addSet     = false;
        else if ( a == "--free" )            freeMode   = true;
        else if ( a == "--no-transfer-bc" )  transferBC = false;
        else if ( a == "--node-set" )        setId      = stoi(need("--node-set"));
        else if ( a == "--plane-tol" )       planeTol   = stod(need("--plane-tol"));
        else if ( a == "--plane" ) {
            string ax = lower(need("--plane"));
            planeAxis = ( ax == "x" || ax == "1" ) ? 1 : ( ax == "y" || ax == "2" ) ? 2 : ( ax == "z" || ax == "3" ) ? 3 : 0;
            if ( !planeAxis ) { cerr << "--plane axis must be x/y/z\n"; return 1; }
            planeVal = stod(need("--plane"));
        }
        else                                 pos.push_back(a);
    }
    // legacy 2-body form folds into the bodies list (A first so A keeps the original id)
    if ( !freeMode && bodies.empty() && csA >= 0 && csB >= 0 ) bodies = { csA, csB };

    if ( pos.size() < 2 || contactCS < 0 || contactMat < 0 ) {
        cerr << "Usage:\n"
             << "  internal: " << argv[0] << " in.in out.in --bodies c1,c2,... --contact-cs C --contact-mat M [--pair P:Q[:cs:mat] ...] [--no-set]\n"
             << "  legacy:   " << argv[0] << " in.in out.in --sideA-cs A --sideB-cs B --contact-cs C --contact-mat M [--no-set]\n"
             << "  free:     " << argv[0] << " in.in out.in --free (--plane <x|y|z> <coord> [--plane-tol t] | --node-set <id>)"
             << " --contact-cs C --contact-mat M [--no-transfer-bc] [--no-set]\n";
        return 1;
    }
    if ( freeMode ) {
        if ( planeAxis == 0 && setId < 0 ) { cerr << "--free needs --plane <axis> <coord> or --node-set <id>\n"; return 1; }
    } else {
        if ( bodies.size() < 2 ) { cerr << "internal mode needs --bodies with >=2 cross-sections (or --sideA-cs/--sideB-cs)\n"; return 1; }
    }
    inFile = pos[0];
    outFile = pos[1];

    ifstream in(inFile);
    if ( !in ) { cerr << "Cannot open " << inFile << "\n"; return 1; }

    // ---- Read preamble up to and including the counts line ---------------------
    vector<string> preamble;
    string countsLine;
    string line;
    while ( getline(in, line) ) {
        preamble.push_back(line);
        if ( line.find("ndofman") != string::npos && line.find("nelem") != string::npos ) {
            countsLine = line;
            break;
        }
    }
    if ( countsLine.empty() ) { cerr << "No counts (ndofman/nelem) line found\n"; return 1; }

    vector<string> ct = tokenize(countsLine);
    auto countOf = [&](const string &k) { return tokenValueAfter(ct, k, 0); };
    const int nNode = countOf("ndofman");
    const int nElem = countOf("nelem");
    const int nCS   = countOf("ncrosssect");
    const int nMat  = countOf("nmat");
    const int nBC   = countOf("nbc");
    const int nIC   = countOf("nic");
    const int nLTF  = countOf("nltf");
    const int nSet  = countOf("nset");

    // Read the next non-blank record line.
    auto readRec = [&](string &out) -> bool {
        while ( getline(in, out) ) {
            size_t p = out.find_first_not_of(" \t\r\n");
            if ( p != string::npos ) return true;
        }
        return false;
    };

    // ---- Nodes -----------------------------------------------------------------
    // Keep raw lines (emitted unchanged); parse coords for geometry. Duplicates are
    // made by cloning the raw line and swapping the id, so BCs on the original
    // interface node carry over to both split nodes.
    vector<string>          nodeRawById;          // index by insertion, but we need id->raw
    unordered_map<int, string> nodeRaw;           // id -> raw line
    unordered_map<int, Vec3>   nodeCoord;         // id -> coords
    vector<int>             nodeOrder;            // preserve original order
    int maxNodeId = 0;
    for ( int k = 0; k < nNode; ++k ) {
        if ( !readRec(line) ) { cerr << "Unexpected EOF reading nodes\n"; return 1; }
        vector<string> t = tokenize(line);
        int id = stoi(t[1]);
        // tokens: node id coords ncoord x y z ...
        int ci = -1;
        for ( size_t j = 0; j < t.size(); ++j ) if ( lower(t[j]) == "coords" ) { ci = (int) j; break; }
        if ( ci < 0 || ci + 4 >= (int) t.size() ) { cerr << "Bad node line: " << line << "\n"; return 1; }
        double x = stod(t[ci + 2]), y = stod(t[ci + 3]), z = stod(t[ci + 4]);
        nodeRaw[id]   = line;
        nodeCoord[id] = { x, y, z };
        nodeOrder.push_back(id);
        maxNodeId = max(maxNodeId, id);
    }

    // ---- Elements --------------------------------------------------------------
    // LTRSpace tets are parsed structurally (for unzip); everything else passes
    // through unchanged.
    vector<TetRec>  tets;
    vector<string>  otherElems;
    int maxElemId = 0;
    for ( int k = 0; k < nElem; ++k ) {
        if ( !readRec(line) ) { cerr << "Unexpected EOF reading elements\n"; return 1; }
        vector<string> t = tokenize(line);
        int id = stoi(t[1]);
        maxElemId = max(maxElemId, id);
        if ( lower(t[0]) == "ltrspace" ) {
            // tokens: LTRSpace id nodes 4 n1 n2 n3 n4 crossSect cs mat m ...
            int ni = -1;
            for ( size_t j = 0; j < t.size(); ++j ) if ( lower(t[j]) == "nodes" ) { ni = (int) j; break; }
            if ( ni < 0 ) { cerr << "LTRSpace without 'nodes': " << line << "\n"; return 1; }
            int nn = stoi(t[ni + 1]);
            if ( nn != 4 ) { cerr << "Only 4-node LTRSpace supported: " << line << "\n"; return 1; }
            TetRec tr;
            tr.id = id;
            tr.tokens = t;
            tr.cs = tokenValueAfter(t, "crossSect");
            for ( int a = 0; a < 4; ++a ) {
                tr.nodeTokIdx[a] = ni + 2 + a;
                tr.n[a] = stoi(t[tr.nodeTokIdx[a]]);
            }
            tets.push_back(tr);
        } else {
            otherElems.push_back(line);
        }
    }

    // ---- Remaining sections: keep raw, preserve order --------------------------
    auto slurp = [&](int n, vector<string> &dst) {
        for ( int k = 0; k < n; ++k ) { if ( readRec(line) ) dst.push_back(line); }
    };
    vector<string> csLines, matLines, bcLines, icLines, ltfLines, setLines;
    slurp(nCS,  csLines);
    slurp(nMat, matLines);
    slurp(nBC,  bcLines);
    slurp(nIC,  icLines);
    slurp(nLTF, ltfLines);
    slurp(nSet, setLines);
    vector<string> trailing;                 // extractor block etc.
    while ( getline(in, line) ) trailing.push_back(line);

    int maxSetId = 0;
    for ( const string &s : setLines ) {
        vector<string> t = tokenize(s);
        if ( t.size() >= 2 ) maxSetId = max(maxSetId, stoi(t[1]));
    }

    // ---- Build tet-face adjacency ---------------------------------------------
    // key = sorted node triple; value = list of (tetIndex, oppositeNode).
    struct FaceUse { int tet; int opp; };
    map<array<int, 3>, vector<FaceUse>> faceMap;
    for ( size_t ti = 0; ti < tets.size(); ++ti ) {
        const TetRec &tr = tets[ti];
        for ( int excl = 0; excl < 4; ++excl ) {
            array<int, 3> f;
            int c = 0;
            for ( int a = 0; a < 4; ++a ) if ( a != excl ) f[c++] = tr.n[a];
            sort(f.begin(), f.end());
            faceMap[f].push_back({ (int) ti, tr.n[excl] });
        }
    }

    // Shared emission outputs, filled by whichever mode runs below.
    vector<int>    dupNodeIds;                 // new duplicate node ids, in emission order
    vector<string> contactLines;
    vector<int>    contactIds;
    int    nFaces = 0;
    double totalArea = 0.;
    int    eid = maxElemId;

    if ( freeMode ) {
        // ==== FREE-SURFACE MODE (unchanged: one support duplicate per node) ======
        map<int, double> nodeArea;             // node -> tributary area
        map<int, Vec3>   nodeNrm;              // node -> area-weighted outward normal

        unordered_set<int> supportNodes;       // node-set variant: gather support ids
        if ( setId >= 0 ) {
            for ( const string &s : setLines ) {
                vector<string> t = tokenize(s);
                if ( t.size() >= 2 && stoi(t[1]) == setId )
                    for ( int id : parseSetNodes(s) ) supportNodes.insert(id);
            }
            if ( supportNodes.empty() ) { cerr << "Node set " << setId << " not found or has no nodes\n"; return 1; }
        }

        for ( const auto &kv : faceMap ) {
            const auto &uses = kv.second;
            const array<int, 3> &f = kv.first;
            if ( uses.size() != 1 ) continue;                        // boundary face only
            Vec3 x1 = nodeCoord[f[0]], x2 = nodeCoord[f[1]], x3 = nodeCoord[f[2]];
            Vec3 fn = ( x2 - x1 ).cross(x3 - x1);
            double area = 0.5 * fn.norm();
            Vec3 n = fn.normalized();
            Vec3 centroid = ( x1 + x2 + x3 ) / 3.0;
            bool onSupport;
            if ( setId >= 0 )
                onSupport = supportNodes.count(f[0]) && supportNodes.count(f[1]) && supportNodes.count(f[2]);
            else
                onSupport = fabs(coordAxis(x1, planeAxis) - planeVal) < planeTol
                         && fabs(coordAxis(x2, planeAxis) - planeVal) < planeTol
                         && fabs(coordAxis(x3, planeAxis) - planeVal) < planeTol;
            if ( !onSupport ) continue;
            int apex = uses[0].opp;                                  // orient outward (away from interior)
            if ( n.dot(nodeCoord[apex] - centroid) > 0 ) n = n * ( -1.0 );
            for ( int a = 0; a < 3; ++a ) {
                nodeArea[f[a]] += area / 3.0;
                nodeNrm[f[a]]   = nodeNrm[f[a]] + n * area;
            }
            ++nFaces;
            totalArea += area;
        }
        if ( nodeArea.empty() ) { cerr << "No support faces found on the requested surface. Nothing to do.\n"; return 1; }

        for ( const auto &kv : nodeArea ) {                          // duplicate + free specimen
            int oldId = kv.first, newId = ++maxNodeId;
            vector<string> t = tokenize(nodeRaw[oldId]);
            t[1] = to_string(newId);
            nodeRaw[newId] = join(t);                                // duplicate keeps the BC
            if ( transferBC ) nodeRaw[oldId] = join(stripBC(tokenize(nodeRaw[oldId])));
            dupNodeIds.push_back(newId);
            Vec3 n = nodeNrm[oldId].normalized();
            ostringstream ss;
            ss << "latticecontact3d " << ++eid
               << " nodes 2 " << oldId << " " << newId
               << " crossSect " << contactCS << " mat " << contactMat
               << " area " << scientific << kv.second
               << " normal 3 " << n.x << " " << n.y << " " << n.z;
            contactLines.push_back(ss.str());
            contactIds.push_back(eid);
        }
    } else {
        // ==== INTERNAL k-WAY MODE ================================================
        // Bodies participate by cross-section id; rank = position in --bodies (the
        // rank-0 body keeps the original node id at every split).
        unordered_map<int, int> rank;                                // cs -> rank
        for ( size_t i = 0; i < bodies.size(); ++i ) rank[bodies[i]] = (int) i;
        auto isBody = [&](int cs) { return rank.count(cs) > 0; };

        // Which pairs get a contact, and with what cross-section/material.
        // Pair key is order-normalized by rank (lo rank first). Empty pairSpecs =>
        // every inter-body pair, using the default contact cs/mat.
        auto keyOf = [&](int a, int b) {                             // -> {loRankCs, hiRankCs}
            return ( rank[a] < rank[b] ) ? make_pair(a, b) : make_pair(b, a);
        };
        map<pair<int,int>, pair<int,int>> pairMat;                   // key -> {cs, mat}, present == emit
        bool allPairs = pairSpecs.empty();
        for ( const PairSpec &ps : pairSpecs ) {
            if ( !isBody(ps.p) || !isBody(ps.q) ) { cerr << "--pair " << ps.p << ":" << ps.q << " names a cross-section not in --bodies\n"; return 1; }
            pairMat[keyOf(ps.p, ps.q)] = { ps.cs < 0 ? contactCS : ps.cs, ps.mat < 0 ? contactMat : ps.mat };
        }
        auto pairEmits = [&](const pair<int,int> &k) { return allPairs || pairMat.count(k) > 0; };
        auto pairCS    = [&](const pair<int,int> &k) { return allPairs ? contactCS  : pairMat[k].first;  };
        auto pairMatId = [&](const pair<int,int> &k) { return allPairs ? contactMat : pairMat[k].second; };

        // Node incidence: which bodies touch each node (from all body tets).
        unordered_map<int, set<int>> nodeBodies;
        for ( const TetRec &tr : tets )
            if ( isBody(tr.cs) )
                for ( int a = 0; a < 4; ++a ) nodeBodies[tr.n[a]].insert(tr.cs);

        // Per-pair, per-node list of face contributions (unit normal lo->hi + tributary
        // area). Kept per-face so a rim node's faces can be split by normal at emit time.
        struct Contrib { Vec3 n; double a; };
        map<pair<int,int>, map<int, vector<Contrib>>> pairContribs;
        for ( const auto &kv : faceMap ) {
            const auto &uses = kv.second;
            if ( uses.size() != 2 ) continue;                        // internal face only
            int cs0 = tets[uses[0].tet].cs, cs1 = tets[uses[1].tet].cs;
            if ( !isBody(cs0) || !isBody(cs1) || cs0 == cs1 ) continue;
            pair<int,int> key = keyOf(cs0, cs1);
            const array<int, 3> &f = kv.first;
            Vec3 x1 = nodeCoord[f[0]], x2 = nodeCoord[f[1]], x3 = nodeCoord[f[2]];
            Vec3 fn = ( x2 - x1 ).cross(x3 - x1);
            double area = 0.5 * fn.norm();
            Vec3 n = fn.normalized();
            Vec3 centroid = ( x1 + x2 + x3 ) / 3.0;
            // orient lo-rank -> hi-rank (toward the hi-rank body's apex)
            int hiUse   = ( rank[cs0] > rank[cs1] ) ? 0 : 1;
            int hiApex  = uses[hiUse].opp;
            if ( n.dot(nodeCoord[hiApex] - centroid) < 0 ) n = n * ( -1.0 );
            for ( int a = 0; a < 3; ++a )
                pairContribs[key][f[a]].push_back({ n, area / 3.0 });
            ++nFaces;
            totalArea += area;
        }
        if ( pairContribs.empty() ) { cerr << "No interface faces found between the requested bodies. Nothing to do.\n"; return 1; }

        // Split every node incident to >=2 bodies: rank-0-most body keeps the
        // original id, each other incident body gets a fresh duplicate.
        unordered_map<int, map<int, int>> nodeCopy;                  // node -> (cs -> node id used by that body)
        for ( const auto &kv : nodeBodies ) {
            int node = kv.first;
            const set<int> &bs = kv.second;
            if ( bs.size() < 2 ) continue;                           // interior to one body
            int keep = *min_element(bs.begin(), bs.end(),
                                    [&](int a, int b) { return rank[a] < rank[b]; });
            for ( int cs : bs ) {
                if ( cs == keep ) { nodeCopy[node][cs] = node; continue; }
                int newId = ++maxNodeId;
                vector<string> t = tokenize(nodeRaw[node]);
                t[1] = to_string(newId);
                nodeRaw[newId] = join(t);                            // clone (BC carries to the duplicate)
                nodeCopy[node][cs] = newId;
                dupNodeIds.push_back(newId);
            }
            // move the prescribed BC onto the duplicate(s): strip it from the kept
            // original so a loaded/supported bar drives the bar, not the concrete.
            if ( bcToDup )
                nodeRaw[node] = join(stripBC(tokenize(nodeRaw[node])));
        }

        // Unzip: every body tet references its own copy of each split node.
        for ( TetRec &tr : tets ) {
            if ( !isBody(tr.cs) ) continue;
            auto nc = nodeCopy.end();
            for ( int a = 0; a < 4; ++a ) {
                nc = nodeCopy.find(tr.n[a]);
                if ( nc == nodeCopy.end() ) continue;
                auto jt = nc->second.find(tr.cs);
                if ( jt == nc->second.end() ) continue;
                tr.n[a] = jt->second;
                tr.tokens[tr.nodeTokIdx[a]] = to_string(jt->second);
            }
        }

        // Emit contacts per requested pair per node (node1 = lo-rank copy). With
        // --split-angle, a node's faces are greedily clustered by normal so a rim
        // (mantle + end cap) yields separate radial and axial links; otherwise all
        // faces at the node merge into one area-weighted normal.
        const double cosThr = ( splitAngleDeg >= 0 ) ? cos(splitAngleDeg * M_PI / 180.0) : -2.0;
        for ( const auto &pk : pairContribs ) {
            const pair<int,int> &key = pk.first;
            if ( !pairEmits(key) ) continue;
            int ccs = pairCS(key), cmat = pairMatId(key);
            for ( const auto &na : pk.second ) {
                int node = na.first;
                int n1 = nodeCopy[node][key.first];                  // lo-rank body copy
                int n2 = nodeCopy[node][key.second];                 // hi-rank body copy
                // geometry-based material override (first matching zone wins)
                int eccs = ccs, ecmat = cmat;
                for ( const ZoneMat &z : zoneMats ) {
                    double c = coordAxis(nodeCoord[node], z.axis);
                    if ( c >= z.lo && c <= z.hi ) { eccs = z.cs; ecmat = z.mat; break; }
                }
                // greedy normal clustering: sum(weighted normal), sum(area) per cluster
                vector<Vec3>   cnrm;
                vector<double> carea;
                for ( const Contrib &c : na.second ) {
                    int best = -1;
                    if ( splitAngleDeg >= 0 ) {
                        double bestDot = cosThr;
                        for ( size_t j = 0; j < cnrm.size(); ++j ) {
                            double d = cnrm[j].normalized().dot(c.n);
                            if ( d > bestDot ) { bestDot = d; best = (int) j; }
                        }
                    } else if ( !cnrm.empty() ) best = 0;            // disabled: single cluster
                    if ( best < 0 ) { cnrm.push_back(c.n * c.a); carea.push_back(c.a); }
                    else            { cnrm[best] = cnrm[best] + c.n * c.a; carea[best] += c.a; }
                }
                for ( size_t j = 0; j < cnrm.size(); ++j ) {
                    Vec3 n = cnrm[j].normalized();
                    ostringstream ss;
                    ss << "latticecontact3d " << ++eid
                       << " nodes 2 " << n1 << " " << n2
                       << " crossSect " << eccs << " mat " << ecmat
                       << " area " << scientific << carea[j]
                       << " normal 3 " << n.x << " " << n.y << " " << n.z;
                    contactLines.push_back(ss.str());
                    contactIds.push_back(eid);
                }
            }
        }
    }

    // ---- Write -----------------------------------------------------------------
    ofstream out(outFile);
    if ( !out ) { cerr << "Cannot open " << outFile << " for writing\n"; return 1; }

    const int newNNode  = (int) nodeRaw.size();
    const int newNElem  = (int) tets.size() + (int) otherElems.size() + (int) contactLines.size();
    const int addedSets = addSet ? ( freeMode ? 2 : 1 ) : 0;   // free adds contact set + support-node set
    const int newNSet   = nSet + addedSets;

    // preamble, rewriting only the counts line
    for ( const string &pl : preamble ) {
        if ( &pl == &preamble.back() ) {
            // rewrite counts: ndofman / nelem / nset
            vector<string> t = tokenize(countsLine);
            ostringstream oss;
            bool sawSet = false;
            for ( size_t j = 0; j < t.size(); ++j ) {
                string k = lower(t[j]);
                if ( k == "ndofman" )      { oss << t[j] << " " << newNNode << " "; ++j; }
                else if ( k == "nelem" )   { oss << t[j] << " " << newNElem << " "; ++j; }
                else if ( k == "nset" )    { oss << t[j] << " " << newNSet  << " "; ++j; sawSet = true; }
                else                        { oss << t[j] << " "; }
            }
            if ( !sawSet && addSet ) oss << "nset " << newNSet << " ";
            out << oss.str() << "\n";
        } else {
            out << pl << "\n";
        }
    }

    // nodes: originals (in order) then duplicates (in order)
    for ( int id : nodeOrder )   out << nodeRaw[id] << "\n";
    for ( int id : dupNodeIds )  out << nodeRaw[id] << "\n";

    // elements: tets (possibly unzipped) + other elements + contacts
    for ( const TetRec &tr : tets ) out << join(tr.tokens) << "\n";
    for ( const string &e : otherElems ) out << e << "\n";
    for ( const string &c : contactLines ) out << c << "\n";

    // remaining sections
    for ( const string &s : csLines )  out << s << "\n";
    for ( const string &s : matLines ) out << s << "\n";
    for ( const string &s : bcLines )  out << s << "\n";
    for ( const string &s : icLines )  out << s << "\n";
    for ( const string &s : ltfLines ) out << s << "\n";
    for ( const string &s : setLines ) out << s << "\n";
    if ( addSet ) {
        out << "Set " << ( maxSetId + 1 ) << " " << formatRanges(contactIds) << "\n";
        if ( freeMode ) {
            // support-node set (the duplicates) so the support BC can be assigned/adjusted
            ostringstream sn;
            sn << "Set " << ( maxSetId + 2 ) << " nodes " << dupNodeIds.size();
            for ( int id : dupNodeIds ) sn << " " << id;
            out << sn.str() << "\n";
        }
    }

    for ( const string &s : trailing ) out << s << "\n";

    cerr << "interface: mode=" << ( freeMode ? "free" : "internal" )
         << ", " << nFaces << " faces, " << dupNodeIds.size() << " duplicate nodes, "
         << contactIds.size() << " contacts, total area " << totalArea << "\n";
    return 0;
}
