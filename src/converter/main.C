#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <math.h>

#include "converterdatareader.h"
#include "convertertxtdatareader.h"
#include "convertertxtinputrecord.h"
#include "converterinputrecord.h"
#include "convertererror.h"

#include "vertex.h"
#include "grid.h"


const char *giveClassName() { return "Main"; }

namespace {
void printUsage(const char *prog)
{
    std::fprintf(stderr,
                 "Usage:\n"
                 "  %s [--mesher t3d|qhull] <control.in> <mesh files…>\n"
                 "\n"
                 "Mesher is required and may be supplied either via the --mesher\n"
                 "(or -m) flag or via a `#@mesher t3d` / `#@mesher qhull` directive\n"
                 "in the control file. CLI flag wins if both are present.\n"
                 "\n"
                 "Mesh-file expectations per mesher:\n"
                 "  t3d   : <mesh.t3d>\n"
                 "  qhull : <mesh.nodes> <mesh.voronoi>\n"
                 "        : <mesh.nodes> <mesh.delaunay> <mesh.voronoi>\n",
                 prog);
}

/// Scan the control file for a `#@mesher` directive.
/// Returns "t3d", "qhull", or "" if none found.
std::string peekMesherFromControl(const std::string &controlFile)
{
    std::ifstream in(controlFile);
    if ( !in.is_open() ) return "";
    std::string line;
    while ( std::getline(in, line) ) {
        std::istringstream iss(line);
        std::string tag;
        if ( !( iss >> tag ) ) continue;
        if ( tag != "#@mesher" ) continue;
        std::string value;
        iss >> value;
        for ( char &c : value ) c = std::tolower(static_cast< unsigned char >(c));
        return value;
    }
    return "";
}
} // namespace

int main(int argc, char *argv[])
{
    // Pull --mesher / -m off argv first, then process remaining positional args.
    std::string mesher;
    std::vector< const char * >positional;
    positional.reserve(argc);
    for ( int i = 0; i < argc; ++i ) {
        std::string a = argv[i];
        if ( ( a == "--mesher" || a == "-m" ) && i + 1 < argc ) {
            mesher = argv[++i];
            for ( char &c : mesher ) c = std::tolower(static_cast< unsigned char >(c));
        } else {
            positional.push_back(argv[i]);
        }
    }

    const int npos = ( int ) positional.size();
    if ( npos < 3 || npos > 5 ) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    const char *controlFile = positional[1];

    // CLI takes priority; fall back to the directive in control.in.
    if ( mesher.empty() ) {
        mesher = peekMesherFromControl(controlFile);
    }
    if ( mesher != "t3d" && mesher != "qhull" ) {
        std::fprintf(stderr,
                     "error: mesher not specified or unknown ('%s').\n"
                     "       Pass --mesher t3d|qhull or add `#@mesher t3d|qhull` to %s.\n",
                     mesher.c_str(), controlFile);
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    auto *grid = new Grid(1);

    if ( mesher == "t3d" ) {
        if ( npos != 3 ) {
            std::fprintf(stderr, "error: t3d mesher expects exactly one mesh file (mesh.t3d), got %d.\n", npos - 2);
            printUsage(argv[0]);
            return EXIT_FAILURE;
        }
        const char *t3dFile = positional[2];
        grid->instanciateYourselfFromT3d(t3dFile, controlFile);
        std::string outName = "oofem.in";
        grid->giveOutputT3d(outName);
    } else { // qhull
        if ( npos == 4 ) {
            const char *nodeFile = positional[2];
            const char *vorFile  = positional[3];
            grid->instanciateYourselfFromQhull(controlFile, nodeFile, vorFile);
            grid->generateOutput();
            grid->giveOutput("oofem.in");
        } else if ( npos == 5 ) {
            // Legacy 3-file qhull form for tetrahedral-element analyses
            // (_3dTetraSM, _3dRCSM, etc.). The mesh.delaunay file comes from a
            // separate `qdelaunay i < mesh.nodes` call. Still goes through the
            // old ConverterTXTDataReader path; should eventually be folded into
            // instanciateYourselfFromQhull (auto-detect Delaunay from positional
            // arg or a `#@delaunay` directive) and the tet output writers
            // migrated to the `#@`-template form like the lattice ones.
            const char *nodeFile = positional[2];
            const char *delFile  = positional[3];
            const char *vorFile  = positional[4];
            ConverterTXTDataReader dr(controlFile);
            const std::string outName = dr.giveOutputFileName();
            if ( outName.empty() ) {
                converter::error("Output filename missing in first line of input file");
            }
            grid->instanciateYourself(&dr, nodeFile, delFile, vorFile);
            grid->generateOutput();
            grid->giveOutput(outName);
        } else {
            std::fprintf(stderr,
                         "error: qhull mesher expects 2 mesh files (nodes + voronoi)\n"
                         "       or 3 (nodes + delaunay + voronoi), got %d.\n",
                         npos - 2);
            printUsage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    std::printf("Conversion complete\n");
    delete grid;
    return 0;
}
