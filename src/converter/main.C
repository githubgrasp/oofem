#include <stdio.h>
#include <fstream>
#include <iostream>
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

int main(int argc, char *argv[]) {
    if ( !( argc == 2 || argc == 4 || argc == 5 ) ) {
        std::fprintf(stderr,
                     "Usage:\n"
                     "  %s <control>\n"
                     "  %s <control> <nodes> <voronoi>\n"
                     "  %s <control> <nodes> <delaunay> <voronoi>\n",
                     argv [ 0 ], argv [ 0 ], argv [ 0 ]);
        return EXIT_FAILURE;
    }

    const char *controlFile = argv [ 1 ];
    const char *nodeFile    = ( argc >= 4 ) ? argv [ 2 ] : nullptr;
    const char *delFile     = ( argc == 5 ) ? argv [ 3 ] : nullptr;
    const char *vorFile     = ( argc >= 4 ) ? argv [ argc - 1 ] : nullptr;

    ConverterTXTDataReader dr(controlFile);

    const std::string outName = dr.giveOutputFileName();
    if ( outName.empty() ) {
        converter::error("Output filename missing in first line of input file");
    }

    auto * grid = new Grid(1);

    grid->instanciateYourself(& dr, nodeFile, delFile, vorFile);

    grid->generateOutput();

    grid->giveOutput(outName);

    std::printf("Conversion complete\n");

    delete grid;
    return 0;
}
