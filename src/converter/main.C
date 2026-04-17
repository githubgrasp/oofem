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

  if ( !( argc == 2 || argc == 3 || argc == 4 || argc == 5 ) ) {
std::fprintf(stderr,
    "Usage:\n"
    "  %s <control>\n"
    "  %s <control> <t3dmesh>\n"
    "  %s <control> <nodes> <voronoi>\n"
    "  %s <control> <nodes> <delaunay> <voronoi>\n",
    argv[0], argv[0], argv[0], argv[0]);
        return EXIT_FAILURE;
    }

    const char *controlFile = argv[1];
    const char *nodeFile    = nullptr;
    const char *delFile     = nullptr;
    const char *vorFile     = nullptr;
    const char *t3dFile     = nullptr;

    if (argc == 3) {
      t3dFile = argv[2];
    }
    else if (argc == 4) {
      nodeFile = argv[2];
      vorFile  = argv[3];
    }
    else if (argc == 5) {
      nodeFile = argv[2];
      delFile  = argv[3];
      vorFile  = argv[4];
    }

    auto * grid = new Grid(1);

    if (t3dFile) {
      grid->instanciateYourselfFromT3d(t3dFile,controlFile);
      std::string outName = "oofem.in";
      grid->giveOutputT3d(outName);
    } else if (nodeFile && !delFile) {
      // qhull path (nodes + voronoi only): use template-style control file
      grid->instanciateYourselfFromQhull(controlFile, nodeFile, vorFile);
      grid->generateOutput();
      grid->giveOutput("oofem.in");
    } else {
      ConverterTXTDataReader dr(controlFile);
      const std::string outName = dr.giveOutputFileName();
      if ( outName.empty() ) {
        converter::error("Output filename missing in first line of input file");
      }
      grid->instanciateYourself(&dr, nodeFile, delFile, vorFile);
      grid->generateOutput();
      grid->giveOutput(outName);
    }
    
    std::printf("Conversion complete\n");

    delete grid;
    return 0;
}
