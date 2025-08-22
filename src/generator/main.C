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

    GeneratorTXTDataReader dr(inputFileName);
    Grid grid(1);

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
