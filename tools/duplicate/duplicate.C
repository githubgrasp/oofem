#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <regex>
#include <algorithm>
#include <map>  

using namespace std;


#include <cmath>
#include <map>


// Helper for 3D vector
struct Vec3 {
    double x, y, z;

    Vec3(double x=0, double y=0, double z=0) : x(x), y(y), z(z) {}

    Vec3 operator-(const Vec3& other) const {
        return Vec3(x - other.x, y - other.y, z - other.z);
    }

    Vec3 operator+(const Vec3& other) const {
        return Vec3(x + other.x, y + other.y, z + other.z);
    }

    Vec3 operator/(double val) const {
        return Vec3(x / val, y / val, z / val);
    }

    double norm() const {
        return std::sqrt(x*x + y*y + z*z);
    }

    Vec3 normalized() const {
        double n = norm();
        return (n == 0) ? Vec3(0, 0, 0) : (*this) / n;
    }
};




struct LinkInfo {
    int id;
    int originalNodeId;
    int newNodeId;
    double diameter;      // user-defined
    double length;        // computed from neighboring elements
    Vec3 dirVector;       // averaged normal vector
    int crossSect = 4;    // default bond section; overridable per zone (--zone-mat)
    int mat = 4;          // default bond material; overridable per zone (--zone-mat)
};



struct NodeEntry {
    int id;
    string line;
};

struct ElementEntry {
    int id;
    string type;
    string line;
};

string formatRanges(const vector<int> &ids) {
    if (ids.empty()) return "";

    vector<int> sorted = ids;
    sort(sorted.begin(), sorted.end());

    stringstream ss;
    ss << "elementranges {";
    int start = sorted[0], prev = sorted[0];

    for (size_t i = 1; i < sorted.size(); ++i) {
        if (sorted[i] == prev + 1) {
            prev = sorted[i];
        } else {
            ss << "(" << start << " " << prev << ") ";
            start = prev = sorted[i];
        }
    }
    ss << "(" << start << " " << prev << ")}";
    return ss.str();
}


// A geometric band along one axis that overrides the bond link's cross-section /
// material. Mirrors the interface tool's --zone-mat: used to model a debonded
// (zero-bond) sleeve over a portion of the bar while keeping the link element so
// the bar stays laterally constrained. First matching zone wins.
struct ZoneMat {
    char axis;       // 'x','y','z'
    double lo, hi;   // inclusive coordinate band
    int crossSect;
    int mat;
};

// ---------------------------------------------------------------------------
// 2D mode
//
// Input: an OOFEM record file in which the bars are truss2d elements on
// ordinary nodes, e.g. straight from t3d2oofem. The file is expected in the
// usual order: header record with ndofman/nelem, then ndofman node records,
// then nelem element records, then everything else (cross sections,
// materials, BCs, time functions, extractor block). That order is kept.
//
// For every truss2d node the tool
//   - turns the bar node into a hanging node with u free and v slaved to the
//     host continuum element (doftype 2 0 2), so the bar can slip along its
//     axis but follows the concrete transversally,
//   - adds a duplicated, fully hanging node at the same position (doftype
//     2 2 2), which carries the concrete displacement, and
//   - joins the pair with a bondlink1d element.
//
// bondlink1d couples the u DOF only, so all bars must be parallel to x.
// ---------------------------------------------------------------------------

static string lowerCase(string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}

static bool isRecord(const string &line) {
    size_t p = line.find_first_not_of(" \t\r");
    return p != string::npos && line[p] != '#';
}

// 2 for a plane stress/strain domain record, 3 otherwise.
static int detectDimension(const vector<string> &lines) {
    for (const string &l : lines) {
        istringstream iss(l);
        string kw, type;
        iss >> kw >> type;
        if (lowerCase(kw) == "domain") {
            type = lowerCase(type);
            return (type.rfind("2d", 0) == 0 || type == "planestrain") ? 2 : 3;
        }
    }
    return 3;
}

static int run2D(const vector<string> &lines, ofstream &outfile, const vector<ZoneMat> &zones,
                 double diameter, int bondCS, int bondMat) {
    // header: locate ndofman / nelem
    int header = -1, nDofMan = 0, nElem = 0;
    for (size_t i = 0; i < lines.size(); ++i) {
        if (lines[i].find("ndofman") != string::npos && lines[i].find("nelem") != string::npos) {
            istringstream iss(lines[i]);
            string tok;
            while (iss >> tok) {
                if (tok == "ndofman") iss >> nDofMan;
                else if (tok == "nelem") iss >> nElem;
            }
            header = (int)i;
            break;
        }
    }
    if (header < 0) {
        cerr << "Error: no record with ndofman and nelem found." << endl;
        return 1;
    }

    // classify the records that follow the header
    struct NodeRec { int line; int id; string kw; vector<string> coords; };
    vector<NodeRec> nodeRecs;
    vector<int> elemLines;
    unordered_map<int, size_t> nodeIndex;   // node id -> index in nodeRecs
    int maxNodeId = 0, maxElemId = 0, lastNodeLine = -1, lastElemLine = -1;
    int count = 0;
    for (size_t i = header + 1; i < lines.size() && count < nDofMan + nElem; ++i) {
        if (!isRecord(lines[i])) continue;
        istringstream iss(lines[i]);
        string kw;
        int id;
        iss >> kw >> id;
        if (count < nDofMan) {
            NodeRec r{(int)i, id, lowerCase(kw), {}};
            string tok;
            while (iss >> tok) {
                if (lowerCase(tok) == "coords") {
                    int n;
                    iss >> n;
                    for (int k = 0; k < n; ++k) { iss >> tok; r.coords.push_back(tok); }
                    break;
                }
            }
            if (r.coords.size() < 2) {
                cerr << "Error: node " << id << " has no coords." << endl;
                return 1;
            }
            nodeIndex[id] = nodeRecs.size();
            nodeRecs.push_back(r);
            maxNodeId = max(maxNodeId, id);
            lastNodeLine = (int)i;
        } else {
            elemLines.push_back((int)i);
            maxElemId = max(maxElemId, id);
            lastElemLine = (int)i;
        }
        ++count;
    }
    if (count < nDofMan + nElem) {
        cerr << "Error: expected " << nDofMan << " nodes and " << nElem << " elements, found "
             << count << " records." << endl;
        return 1;
    }

    auto coordOf = [&](int id) {
        const auto &c = nodeRecs[nodeIndex.at(id)].coords;
        return Vec3(std::stod(c[0]), std::stod(c[1]), c.size() > 2 ? std::stod(c[2]) : 0.0);
    };

    // bar elements: truss2d. Collect node order, directions and lengths.
    vector<int> barNodes;   // in order of first appearance
    unordered_map<int, vector<Vec3>> nodeDirs;
    unordered_map<int, vector<double>> nodeLengths;
    for (int li : elemLines) {
        istringstream iss(lines[li]);
        string kw, tok;
        int id;
        iss >> kw >> id;
        if (lowerCase(kw) != "truss2d") continue;
        int n1 = 0, n2 = 0;
        while (iss >> tok) {
            if (lowerCase(tok) == "nodes") { int n; iss >> n >> n1 >> n2; break; }
        }
        if (!nodeIndex.count(n1) || !nodeIndex.count(n2)) {
            cerr << "Error: truss2d " << id << " refers to an unknown node." << endl;
            return 1;
        }
        Vec3 d = coordOf(n2) - coordOf(n1);
        for (int n : {n1, n2}) {
            if (!nodeDirs.count(n)) barNodes.push_back(n);
            nodeDirs[n].push_back(d.normalized());
            nodeLengths[n].push_back(d.norm());
        }
    }
    if (barNodes.empty()) {
        cerr << "Error: no truss2d elements found." << endl;
        return 1;
    }

    // one duplicated node and one link per bar node
    vector<LinkInfo> links;
    unordered_set<int> isBarNode(barNodes.begin(), barNodes.end());
    int nextNode = maxNodeId, nextElem = maxElemId;
    for (int n : barNodes) {
        LinkInfo link;
        link.id = ++nextElem;
        link.originalNodeId = n;
        link.newNodeId = ++nextNode;
        link.diameter = diameter;
        link.crossSect = bondCS;
        link.mat = bondMat;

        Vec3 dir;
        double len = 0.;
        for (const Vec3 &d : nodeDirs[n]) dir = dir + d;
        for (double l : nodeLengths[n]) len += l;
        size_t k = nodeDirs[n].size();
        link.dirVector = (dir / (double)k).normalized();
        link.length = len / (double)k;
        if (k == 1) link.length *= 0.5;   // bar end: half the adjacent element

        if (std::fabs(std::fabs(link.dirVector.x) - 1.0) > 1.e-6) {
            cerr << "Error: bar node " << n << " is not on a bar parallel to x. "
                 << "bondlink1d couples u only." << endl;
            return 1;
        }

        Vec3 c = coordOf(n);
        for (const auto &z : zones) {
            double v = (z.axis == 'x') ? c.x : (z.axis == 'y') ? c.y : c.z;
            if (v >= z.lo && v <= z.hi) {
                link.crossSect = z.crossSect;
                link.mat = z.mat;
                break;
            }
        }
        links.push_back(link);
    }

    // write
    for (int i = 0; i < (int)lines.size(); ++i) {
        const string &line = lines[i];
        if (i == header) {
            istringstream iss(line);
            ostringstream oss;
            string tok;
            while (iss >> tok) {
                if (tok == "ndofman") { int v; iss >> v; oss << tok << " " << nDofMan + links.size() << " "; }
                else if (tok == "nelem") { int v; iss >> v; oss << tok << " " << nElem + links.size() << " "; }
                else oss << tok << " ";
            }
            outfile << oss.str() << "\n";
            continue;
        }

        bool written = false;
        if (i <= lastNodeLine && isRecord(line)) {
            istringstream iss(line);
            string kw;
            int id;
            iss >> kw >> id;
            if (isBarNode.count(id) && lowerCase(kw) == "node") {
                // bar node: u free (slip), v from the host element
                string rest = line.substr(line.find(kw) + kw.size());
                outfile << "hangingnode" << rest;
                if (lowerCase(line).find("doftype") == string::npos) outfile << " doftype 2 0 2";
                outfile << "\n";
                written = true;
            }
        }
        if (!written) outfile << line << "\n";

        if (i == lastNodeLine) {
            for (const auto &link : links) {
                const auto &c = nodeRecs[nodeIndex.at(link.originalNodeId)].coords;
                outfile << "hangingnode " << link.newNodeId << " coords " << c.size();
                for (const string &v : c) outfile << " " << v;
                outfile << " doftype 2 2 2\n";
            }
        }
        if (i == lastElemLine) {
            for (const auto &link : links) {
                outfile << "bondlink1d " << link.id
                        << " nodes 2 " << link.originalNodeId << " " << link.newNodeId
                        << " crossSect " << link.crossSect << " mat " << link.mat
                        << " length " << link.length
                        << " diameter " << link.diameter
                        << " dirvector 3 " << link.dirVector.x << " " << link.dirVector.y << " "
                        << link.dirVector.z << "\n";
            }
        }
    }
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0]
             << " inputfile outputfile [--2d | --3d] [--diameter <d>] [--bond <cs> <mat>]"
             << " [--zone-mat <axis> <lo> <hi> <cs> <mat>]..." << endl;
        return 1;
    }

    ifstream infile(argv[1]);
    ofstream outfile(argv[2]);

    if (!infile.is_open() || !outfile.is_open()) {
        cerr << "Failed to open input or output file." << endl;
        return 1;
    }

    // --- options ----------------------------------------------------------------
    std::vector<ZoneMat> zones;
    double diameter = 0.016;   // bond diameter written to every link
    int bondCS = 4, bondMat = 4;  // default bond section / material
    int forceDim = 0;          // 0 = detect from the domain record, 2 or 3 = forced
    for (int a = 3; a < argc; ++a) {
        std::string opt = argv[a];
        if (opt == "--diameter") {
            if (a + 1 >= argc) {
                cerr << "Error: --diameter needs <d>" << endl;
                return 1;
            }
            diameter = std::stod(argv[++a]);
        } else if (opt == "--bond") {
            if (a + 2 >= argc) {
                cerr << "Error: --bond needs <cs> <mat>" << endl;
                return 1;
            }
            bondCS = std::stoi(argv[++a]);
            bondMat = std::stoi(argv[++a]);
        } else if (opt == "--2d") {
            forceDim = 2;
        } else if (opt == "--3d") {
            forceDim = 3;
        } else if (opt == "--zone-mat") {
            if (a + 5 >= argc) {
                cerr << "Error: --zone-mat needs <axis> <lo> <hi> <cs> <mat>" << endl;
                return 1;
            }
            ZoneMat z;
            std::string axis = argv[++a];
            z.axis = axis.empty() ? 'z' : axis[0];
            z.lo  = std::stod(argv[++a]);
            z.hi  = std::stod(argv[++a]);
            z.crossSect = std::stoi(argv[++a]);
            z.mat = std::stoi(argv[++a]);
            zones.push_back(z);
        } else {
            cerr << "Unknown option: " << opt << endl;
            return 1;
        }
    }

    // --- 2D mode (plane stress / plane strain): see run2D() -------------------
    {
        vector<string> all;
        string l;
        while (getline(infile, l)) all.push_back(l);
        int dim = forceDim;
        if (dim == 0) dim = detectDimension(all);
        if (dim == 2) return run2D(all, outfile, zones, diameter, bondCS, bondMat);
        infile.clear();
        infile.seekg(0);
    }

    string line;
 unordered_map<int, NodeEntry> nodes;
 unordered_map<int, NodeEntry> hangingNodes;
 unordered_map<int, ElementEntry> elements;

// This stores info for each link between original and duplicated nodes
std::vector<LinkInfo> linkInfos;

    vector<string> preambleLines;

    vector<string> crossSectionLines;
    vector<string> materialLines;
    vector<string> bcLoadLines;
    vector<string> timeFunctionLines;
    vector<string> setLines;
    vector<string> otherLines;

    vector<string> linkLines;

    int nextNodeId = 0;
    int nextElementId = 0;
    bool headerParsed = false;

    while (getline(infile, line)) {
        if (!headerParsed) {
            preambleLines.push_back(line);
            if (line.find("ndofman") != string::npos && line.find("nelem") != string::npos) {
                headerParsed = true;
                continue;
            }
            continue;
        }

        istringstream iss(line);
        string type;
        iss >> type;

        if (type == "node" || type == "hangingnode") {
            int id;
            iss >> id;
            nodes[id] = {id, line};
            nextNodeId = max(nextNodeId, id);
        } else if (type == "ltrspace" || type == "truss3d" || type == "latticeframe3d") {
            int id;
            iss >> id;
            elements[id] = {id, type, line};
            nextElementId = max(nextElementId, id);
        } else if (type == "SimpleCS" || type == "LatticeCS" || type == "latticecs") {
            crossSectionLines.push_back(line);
        } else if (type == "IsotropicMaterial" || type == "con2dpm" || type == "LatticeFrameSteelPlastic" || type == "latticeslip") {
            materialLines.push_back(line);
        } else if (type == "BoundaryCondition" || type == "NodalLoad") {
            bcLoadLines.push_back(line);
        } else if (type == "ConstantFunction" || type == "PiecewiseLinFunction") {
            timeFunctionLines.push_back(line);
        } else if (type == "Set") {
            setLines.push_back(line);
        } else {
            otherLines.push_back(line);
        }
    }

    std::unordered_map<int, Vec3> nodeCoords;
std::unordered_map<int, std::vector<std::pair<int, Vec3>>> nodeToElementDirs;

// Store node coordinates
for (const auto& [id, node] : nodes) {
    std::istringstream iss(node.line);
    std::string dummy;
    int nodeId;
    double x, y, z;
    iss >> dummy >> nodeId >> dummy >> dummy >> x >> y >> z;
    nodeCoords[nodeId] = Vec3(x, y, z);
}

// Build map: node -> list of (element ID, dir vector)
for (const auto& [id, elem] : elements) {
    if (elem.type != "truss3d" && elem.type != "latticeframe3d") continue;

    std::istringstream iss(elem.line);
    std::string dummy;
    int eid, n1, n2;
    iss >> dummy >> eid >> dummy >> dummy >> n1 >> n2;

    Vec3 dir = (nodeCoords[n2] - nodeCoords[n1]).normalized();
    nodeToElementDirs[n1].emplace_back(eid, dir);
    nodeToElementDirs[n2].emplace_back(eid, dir);
}


    unordered_map<int, int> duplicatedNodeMap;

// Loop over all truss3d and latticeframe3d elements
    for (const auto& [eid, elem] : elements) {
        if (elem.type != "truss3d" && elem.type != "latticeframe3d") continue;

        std::istringstream iss(elem.line);
        std::string dummy;
        int id, mat, cs, n1, n2;
        iss >> dummy >> id >> dummy >> dummy >> n1 >> n2;

        for (int origId : {n1, n2}) {
            if (duplicatedNodeMap.count(origId)) continue; // already duplicated

            int newId = ++nextNodeId;
            duplicatedNodeMap[origId] = newId;

            // Create hanging node based on original node coordinates
            auto coordIt = nodeCoords.find(origId);
            if (coordIt != nodeCoords.end()) {
                Vec3 coord = coordIt->second;

                std::ostringstream newLine;
                newLine << "hangingnode " << newId
                        << " coords 3 " << coord.x << " " << coord.y << " " << coord.z
                        << " dofType 6 2 2 2 2 2 2"
                        << " dofidmask 6 1 2 3 4 5 6";

                nodes[newId] = {newId, newLine.str()};
            }

            // Add new link element
            LinkInfo link;
            link.id = ++nextElementId;
            link.originalNodeId = origId;
            link.newNodeId = newId;
            link.diameter = diameter; // --diameter, default 0.016
            link.crossSect = bondCS;  // --bond, default 4
            link.mat = bondMat;       // --bond, default 4
            linkInfos.push_back(link);
        }
    }


    for (auto& link : linkInfos) {
        int n1 = link.originalNodeId;

        Vec3 avgDir(0, 0, 0);
        double avgLength = 0.0;
        int count = 0;

        auto it = nodeToElementDirs.find(n1);
        if (it != nodeToElementDirs.end()) {
            for (const auto& [eid, dir] : it->second) {
                const auto& element = elements.at(eid);

                std::istringstream iss(element.line);
                std::string dummy;
                int id, mat, cs, nodeA, nodeB;
                iss >> dummy >> id >> dummy >> dummy >> nodeA >> nodeB;

                if (nodeCoords.count(nodeA) && nodeCoords.count(nodeB)) {
                    Vec3 d = nodeCoords[nodeB] - nodeCoords[nodeA];
                    double len = d.norm();
                    avgDir = avgDir + d.normalized();
                    avgLength += len;
                    count++;
                }
            }
        }

        if (count > 0) {
            avgDir = avgDir / count;
            avgDir = avgDir.normalized();
            avgLength = avgLength / count;

	    // Correction for boundary node (only one adjacent element)
	    if (count == 1) {
	      avgLength *= 0.5;
	    }
	}

        link.dirVector = avgDir;
        link.length = avgLength;
    }



for (auto& link : linkInfos) {
    // Resolve the link's section/material from the original bar-node coordinate.
    // First matching --zone-mat band wins; outside all bands keep the --bond default.
    auto cIt = nodeCoords.find(link.originalNodeId);
    if (cIt != nodeCoords.end()) {
        const Vec3& c = cIt->second;
        for (const auto& z : zones) {
            double v = (z.axis == 'x') ? c.x : (z.axis == 'y') ? c.y : c.z;
            if (v >= z.lo && v <= z.hi) {
                link.crossSect = z.crossSect;
                link.mat = z.mat;
                break;
            }
        }
    }

    stringstream ss;
    ss << "latticelink3d " << link.id
       << " nodes 2 " << link.originalNodeId << " " << link.newNodeId
       << " crossSect " << link.crossSect << " mat " << link.mat
       << " length " << link.length
       << " diameter " << link.diameter
       << " dirvector 3 " << link.dirVector.x << " " << link.dirVector.y << " " << link.dirVector.z;
    linkLines.push_back(ss.str());
}


    for (string &line : preambleLines) {
        if (line.find("ndofman") != string::npos && line.find("nelem") != string::npos) {
            istringstream iss(line);
            ostringstream oss;
            string token;
            while (iss >> token) {
                if (token == "ndofman") {
                    oss << token << " " << nodes.size() << " ";
                    iss >> token; // skip old number
                } else if (token == "nelem") {
                    oss << token << " " << elements.size() + linkLines.size() << " ";
                    iss >> token; // skip old number
                } else {
                    oss << token << " ";
                }
            }
            line = oss.str();  // update the line in preambleLines
            break;
        }
    }

    for (const string &line : preambleLines) {
        outfile << line << "\n";
    }

    // Write node and hanging node entries in order
    map<int, NodeEntry> orderedNodes;
    for (auto &[id, entry] : nodes) orderedNodes[id] = entry;
    for (auto &[id, entry] : hangingNodes) orderedNodes[id] = entry;
    for (auto &[id, entry] : orderedNodes) outfile << entry.line << "\n";

    // Write original elements
    map<int, ElementEntry> orderedElements;
    for (auto &[id, entry] : elements) orderedElements[id] = entry;
    for (auto &[id, entry] : orderedElements) outfile << entry.line << "\n";

    // Write latticelink3d elements
    for (const string &line : linkLines) outfile << line << "\n";

    for (const string &line : crossSectionLines) outfile << line << "\n";
    for (const string &line : materialLines) outfile << line << "\n";
    for (const string &line : bcLoadLines) outfile << line << "\n";
    for (const string &line : timeFunctionLines) outfile << line << "\n";

    // --- Generate element sets ---
    vector<int> trussOrFrameElems, linkElems, otherElems;

    for (const auto& [id, elem] : elements) {
        if (elem.type == "truss3d" || elem.type == "latticeframe3d") {
            trussOrFrameElems.push_back(id);
        } else {
            otherElems.push_back(id);
        }
    }
    for (const auto& link : linkInfos) {
        linkElems.push_back(link.id);
    }

/* // --- Format set definitions --- */
/*     auto generateSetLine = [](int setId, const string& name, const vector<int>& elems) { */
/*         stringstream ss; */
/*         ss << "Set " << setId << " elementranges " << elems.size(); */
/*         for (int id : elems) ss << " " << id; */
/*         ss << " name " << name; */
/*         return ss.str(); */
/*     }; */

/*     setLines.push_back(generateSetLine(1, "truss_or_frame", trussOrFrameElems)); */
/*     setLines.push_back(generateSetLine(2, "latticelink", linkElems)); */
/*     setLines.push_back(generateSetLine(3, "others", otherElems)); */


    if (!trussOrFrameElems.empty())
        outfile << "Set 1 " << formatRanges(trussOrFrameElems) << "\n";
    if (!linkElems.empty())
        outfile << "Set 2 " << formatRanges(linkElems) << "\n";
    if (!otherElems.empty())
        outfile << "Set 3 " << formatRanges(otherElems) << "\n";

    outfile << "Set 4 nodes 1 9 \n";
    outfile << "Set 5 nodes 1 12 \n";

    for (const string &line : otherLines) outfile << line << "\n";

    return 0;
}
