#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>


#include "floatarray.h"
#include "intarray.h"

#include "generatorerror.h"

#include "grid.h"

const char * giveClassName() { return "Main"; }


int main(int argc, char *argv[])
{
    if ( argc != 2 ) {
        std::cerr << "Usage: generator <input-file>\n";
        return EXIT_FAILURE;
    }

    const char *inputFileName = argv [ 1 ];

    Grid grid(1);

    grid.readControlRecords(inputFileName);
    const std::string &outputPath = grid.giveOutputFileName();

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
