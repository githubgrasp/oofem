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
#include <fstream>
#include <string>
#include <cctype>


#include "floatarray.h"
#include "intarray.h"

// your generator I/O
#include "generatordatareader.h"
#include "generatortxtdatareader.h"
#include "generatorinputrecord.h"
#include "generatorerror.h"
// your domain
#include "grid.h"

const char * giveClassName() { return "Main"; }

static inline void rtrim(std::string &s) {
    while ( !s.empty() && ( s.back() == '\r' || s.back() == '\n' || s.back() == ' ' || s.back() == '\t' ) ) {
        s.pop_back();
    }
}
static inline bool is_comment_or_empty(const std::string &s) {
    for (unsigned char ch : s) {
        if ( !std::isspace(ch) ) {
            return ch == '#';
        }
    }
    return true;
}
static std::string read_first_line_filename(const char *path) {
    std::ifstream in(path);
    if ( !in ) {
        generator::errorf("Can't open input file (%s)", path);
    }

    std::string line;
    while ( std::getline(in, line) ) {
        rtrim(line);
        if ( !is_comment_or_empty(line) ) {
            return line;
        }
    }
    return {};
}



int main(int argc, char *argv[])
{
    if ( argc != 2 ) {
        std::cerr << "Usage: generator <input-file>\n";
        return EXIT_FAILURE;
    }

    const char *inputFileName = argv [ 1 ];


    // read the output file name directly from the first non-empty, non-comment /* line */
    /* std::string outFileName = read_first_line_filename(inputFileName); */
    /* if (outFileName.empty()) { */
    /*   generator::error("Output filename missing in first line of input file"); */
    /* } */

    /* FILE* outputStream = std::fopen(outFileName.c_str(), "w"); */
    /* if (!outputStream) { */
    /*   std::perror(("Can't open output file: " + outFileName).c_str()); */
    /*   return EXIT_FAILURE; */
    /* } */



    // Concrete reader (your copied/renamed version of OOFEMTXTDataReader)
    GeneratorTXTDataReader dr(inputFileName);

    /*     auto &irOut = dr.giveInputRecord(GeneratorDataReader::GIR_outFileRec, 1); */

    /*     std::string outFileName; */
    /*     irOut.giveRawLine(outFileName);   // fill it first */

    /*     // debug */
    /*     std::cout << "first line = [" << outFileName << "] len=" << outFileName.size() << "\n"; */
    /*     for (unsigned char c : outFileName) std::cout << std::hex << (int)c << ' '; */
    /*     std::cout << std::dec << "\n"; */

    /* // trim trailing CR/LF/spaces */
    /* while (!outFileName.empty() && (outFileName.back()=='\r' || outFileName.back()=='\n' || outFileName.back()==' ')) */
    /*     outFileName.pop_back(); */

    /* if (outFileName.empty()) { */
    /*     generator::error("Output filename missing in first line of input file"); */
    /* } */

    /* FILE* outputStream = std::fopen(outFileName.c_str(), "w"); */
    /* if (!outputStream) { */
    /*   std::perror(("Can't open output file: " + outFileName).c_str()); */
    /*   return EXIT_FAILURE; */
    /* } */

    // irOut.finish();


    /* // --- read output file record */
    /* auto &irOut = dr.giveInputRecord(GeneratorDataReader::GIR_outManRec, 1); */
    /* std::string outFileName; */
    /* printf("line = %s\n", outFileName.c_str()); */
    /* irOut.giveRawLine(outFileName);   // no keyword, just the line */

    /* if (outFileName.empty()) { */
    /*   generator::error("Output filename missing in first line of input file"); */
    /* } */



    /* auto &irOut = dr.giveInputRecord(GeneratorDataReader::GIR_outManRec, 1); */

    /* std::string outFileName; */
    /* IR_GIVE_FIELD(irOut, outFileName, _IFT_outfile); */
    /* irOut.finish(); */

    /* FILE* outputStream = std::fopen(outFileName.c_str(), "w"); */
    /* if (!outputStream) { */
    /*   std::cerr << "Can't open output file: " << outFileName << "\n"; */
    /*   return EXIT_FAILURE; */
    /* } */


    /* if (!outputStream) { */
    /*     std::perror(("Can't open output file: " + outFileName).c_str()); */
    /*     return EXIT_FAILURE; */
    /* } */

    // --- build grid
    Grid grid(1);

    // NOTE: instanciateYourself signature should take GeneratorDataReader*
    // If yours still takes DataReader*, this still compiles because of inheritance.
    grid.instanciateYourself(& dr);



    grid.generatePoints();
    grid.exportVTK("points.vtk");

    FILE *outputStream = std::fopen(dr.giveOutputFileName().c_str(), "w");
    if ( !outputStream ) {
        std::perror( ( "Can't open output file: " + dr.giveOutputFileName() ).c_str() );
        return EXIT_FAILURE;
    }

    grid.giveOutput(outputStream);

    std::fclose(outputStream);

    std::cout << "numberOfVertices = " << grid.giveNumberOfVertices() << "\n";
    std::cout << "Point generation completed\n";
    return EXIT_SUCCESS;
}
