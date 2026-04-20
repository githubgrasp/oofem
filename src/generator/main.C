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

#include "generatordatareader.h"
#include "generatortxtdatareader.h"
#include "generatorinputrecord.h"
#include "generatorerror.h"

#include "grid.h"

const char * giveClassName() { return "Main"; }

static inline void rtrim(std::string &s) {
    while ( !s.empty() && ( s.back() == '\r' || s.back() == '\n' || s.back() == ' ' || s.back() == '\t' ) ) {
        s.pop_back();
    }
}
/// Return true if the input file contains any `#@` directive line.
/// Old-format files never contain `#@`; `#@`-format files start with them.
static bool uses_directive_format(const char *path) {
    std::ifstream in(path);
    if ( !in ) {
        generator::errorf("Can't open input file (%s)", path);
    }
    std::string line;
    while ( std::getline(in, line) ) {
        rtrim(line);
        if ( line.rfind("#@", 0) == 0 ) {
            return true;
        }
    }
    return false;
}



int main(int argc, char *argv[])
{
    if ( argc != 2 ) {
        std::cerr << "Usage: generator <input-file>\n";
        return EXIT_FAILURE;
    }

    const char *inputFileName = argv [ 1 ];

    Grid grid(1);
    std::string outputPath;

    if ( uses_directive_format(inputFileName) ) {
        grid.readControlRecords(inputFileName);
        outputPath = grid.giveOutputFileName();
    } else {
        GeneratorTXTDataReader dr(inputFileName);
        grid.instanciateYourself(& dr);
        outputPath = dr.giveOutputFileName();
    }

    grid.generatePoints();
    if ( grid.giveVtkFlag() ) {
        grid.exportVTK("points.vtk");
    }

    FILE *outputStream = std::fopen(outputPath.c_str(), "w");
    if ( !outputStream ) {
        std::perror( ( "Can't open output file: " + outputPath ).c_str() );
        return EXIT_FAILURE;
    }

    grid.giveOutput(outputStream);

    std::fclose(outputStream);

    std::cout << "numberOfVertices = " << grid.giveNumberOfVertices() << "\n";
    std::cout << "Point generation completed\n";
    return EXIT_SUCCESS;
}
