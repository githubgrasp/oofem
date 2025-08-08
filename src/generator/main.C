/* #include "logger.h" */


/* #include <stdio.h> */
/* #include <fstream> */
/* #include <iostream> */
/* #include <cstdlib> */
/* #include <math.h> */
/* #include "floatarray.h" */
/* #include "intarray.h" */
/* #include "datareader.h" */
/* #include "oofemtxtdatareader.h" */
/* #include "oofemtxtinputrecord.h" */
/* #include "vertex.h" */
/* #include "grid.h" */


/* /// Returns class name of the receiver. */
/*  char *giveClassName() { return "Main"; } */

/* int main(int argc, char* argv[]){ */
/*   if(argc != 2){ */
/*     printf("wrong number of arguments: name of input file is needed\n"); */
/*     std::exit(1); */
/*   } */
  
/*   Grid *grid = new Grid(1); */

  
/*   char inputFileName [ MAX_FILENAME_LENGTH + 10 ], buff [ MAX_FILENAME_LENGTH ]; */
/*   char dataOutputFileName [ MAX_FILENAME_LENGTH + 10 ]; */
/*   strcpy(inputFileName, argv[1]); */
/*   OOFEMTXTDataReader dr(inputFileName); */

/*   //Read the output file name */
/*   InputRecord *ir = dr.giveInputRecord(DataReader :: IR_outFileRec, 1); */
/*   //  __keyword = NULL; */
/*   result = ir->giveField(dataOutputFileName, MAX_FILENAME_LENGTH, IFT_outfile, NULL); */
/*   if ( result != IRRT_OK ) { */
/*     IR_IOERR("", __proc, IFT_outfile, "Output file record", ir, result); */
/*   } */
 
/*   FILE *outputStream;   */
/*   if ( ( outputStream = fopen(dataOutputFileName, "w") ) == NULL ) { */
/*     printf("Can't open output file %s", dataOutputFileName); */
/*     exit(1); */
/*   } */
  
/*   grid->instanciateYourself(&dr); */
  
/*   grid->generatePoints(); */
  
/*   grid->giveOutput(outputStream); */

/*   printf("numberOfVertices = %d\n", grid->giveNumberOfVertices()); */
/*   printf("Point generation completed\n"); */
/*   return 1; */
/* } */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include "floatarray.h"
#include "intarray.h"

// your generator I/O
#include "generatordatareader.h"
#include "generatortxtdatareader.h"
#include "generatorinputrecord.h"
#include "generatorerror.h"
// your domain
#include "grid.h"

// If your IR_* macros ever need this:
const char* giveClassName() { return "Main"; }

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: generator <input-file>\n";
        return EXIT_FAILURE;
    }

    const char* inputFileName = argv[1];

    // Concrete reader (your copied/renamed version of OOFEMTXTDataReader)
    GeneratorTXTDataReader dr(inputFileName);


    // --- read output file record
    auto &irOut = dr.giveInputRecord(GeneratorDataReader::GIR_outManRec, 1);
    std::string outFileName;
    irOut.giveRawLine(outFileName);   // no keyword, just the line

    if (outFileName.empty()) {
      generator::error("Output filename missing in first line of input file");
    }

    irOut.finish();
    
    /* auto &irOut = dr.giveInputRecord(GeneratorDataReader::GIR_outManRec, 1); */

    /* std::string outFileName; */
    /* IR_GIVE_FIELD(irOut, outFileName, _IFT_outfile); */
    /* irOut.finish(); */

    FILE* outputStream = std::fopen(outFileName.c_str(), "w");
    if (!outputStream) {
      std::cerr << "Can't open output file: " << outFileName << "\n";
      return EXIT_FAILURE;
    }


    /* if (!outputStream) { */
    /*     std::perror(("Can't open output file: " + outFileName).c_str()); */
    /*     return EXIT_FAILURE; */
    /* } */

    // --- build grid
    Grid grid(1);

    // NOTE: instanciateYourself signature should take GeneratorDataReader*
    // If yours still takes DataReader*, this still compiles because of inheritance.
    grid.instanciateYourself(&dr);

    grid.generatePoints();
    grid.giveOutput(outputStream);

    std::fclose(outputStream);

    std::cout << "numberOfVertices = " << grid.giveNumberOfVertices() << "\n";
    std::cout << "Point generation completed\n";
    return EXIT_SUCCESS;
}
