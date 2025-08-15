#include <stdio.h>

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include "converterdatareader.h"
#include "convertertxtdatareader.h"
#include "convertertxtinputrecord.h"
#include "converterinputrecord.h"
#include "generatorinputrecord.h"
#include "generatorerror.h"

#include "vertex.h"
#include "grid.h"

const char *giveClassName() { return "Main"; }

int main(int argc, char *argv[]) {

  if (!(argc == 2 || argc == 4 || argc == 5)) {
        std::fprintf(stderr,
            "Usage:\n"
            "  %s <control>\n"
            "  %s <control> <nodes> <voronoi>\n"
            "  %s <control> <nodes> <delaunay> <voronoi>\n",
            argv[0], argv[0], argv[0]);
        return EXIT_FAILURE;
    }

    const char* controlFile = argv[1];
    const char* nodeFile    = (argc >= 4) ? argv[2] : nullptr;
    const char* delFile     = (argc == 5) ? argv[3] : nullptr;
    const char* vorFile     = (argc >= 4) ? argv[argc - 1] : nullptr;

    ConverterTXTDataReader dr(controlFile);

    // Set output file name from first line of control file
    const std::string outName = dr.giveOutputFileName();
    if (outName.empty()) {
        generator::error("Output filename missing in first line of input file");
    }

    clock_t start = std::clock();

    auto* grid = new Grid(1);

    // Safe: pass nullptrs when the files aren’t provided
    grid->instanciateYourself(&dr, nodeFile, delFile, vorFile);

    grid->generateOutput();

    // If your Grid::giveOutput expects a C string filename:
    grid->giveOutput(const_cast<char*>(outName.c_str()));

    std::printf("Conversion complete\n");

    delete grid;
    return 0;
}


/*   if ( ! ( argc == 2 ||argc == 4 || argc == 5 ) ) { */
/*     printf("wrong number of arguments. Either\n 1) name of control and t3d file \n or \n 2) name of control file, the names of nodes, delaunay and voronoi file, if continuum \n or \n 3) name of control file, the names of nodes, and voronoi file\n"); */
/*         std :: exit(1); */
/*     } */

/*   const char *inputFileName = argv [ 1 ]; */
  
/*   Grid grid(1); */
      
/*     //These two things are required by the macro */
/*     const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro */

/*     IRResultType result;              // Required by IR_GIVE_FIELD macro */
/*     char inputFileName [ MAX_FILENAME_LENGTH + 10 ], buff [ MAX_FILENAME_LENGTH ]; */
/*     char outFileName [ MAX_FILENAME_LENGTH + 10 ]; */

/*     strcpy(inputFileName, argv [ 1 ]); */

/*     char nodeFileName [ MAX_FILENAME_LENGTH + 10 ]; */
/*     strcpy(nodeFileName, argv [ 2 ]); */

/*     char delaunayFileName [ MAX_FILENAME_LENGTH + 10 ]; */
/*     char voronoiFileName [ MAX_FILENAME_LENGTH + 10 ];     */
/*     if(argc == 5){//continuum mesh */
/*       strcpy(delaunayFileName, argv [ 3 ]); */
/*       strcpy(voronoiFileName, argv [ 4 ]); */
/*     } */
/*     else if(argc == 4){//lattice mesh */
/*       strcpy(voronoiFileName, argv [ 3 ]); */
/*     } */
    
    
/*     ConverterTXTDataReader dr(inputFileName); */

/*     dr->giveOutputFileName(outputFileName); */

/*     // measure time consumed by converter */
/*     clock_t start; */
/*     clock_t sc; */
/*     clock_t ec; */
/*     long nsec; */
    
/*     start = :: clock();     */
/*     sc = :: clock(); */

/*     grid->instanciateYourself(& dr, nodeFileName, delaunayFileName, voronoiFileName); */

/*     ec = :: clock(); */
/*     nsec = ( ec - sc ) / CLOCKS_PER_SEC; */
/*     OOFEM_LOG_INFO("instanciateYourself completed in %lds \n", nsec); */

/*     sc = :: clock(); */
    
/*     grid->generateOutput(); */

/*     ec = :: clock(); */
/*     nsec = ( ec - sc ) / CLOCKS_PER_SEC; */

/*     OOFEM_LOG_INFO("generateOutput completed in %lds \n", nsec); */

/*     sc = :: clock(); */
        
/*     grid->giveOutput(outputFileName); */

/*     ec = :: clock(); */
/*     nsec = ( ec - sc ) / CLOCKS_PER_SEC; */
/*     OOFEM_LOG_INFO("giveOutput completed in %lds \n", nsec); */

/*     nsec = ( ec - start ) / CLOCKS_PER_SEC; */
/*     OOFEM_LOG_INFO("conversion completed in %lds \n", nsec); */
    
/*     printf("Conversion complete \n"); */

/*     delete grid; */
    
    

/*     return 1; */
/* } */
