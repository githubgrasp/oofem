#include <cstdlib>
#include <iostream>

#include "box.h"

int main(int argc, char *argv[])
{
    if ( argc != 2 ) {
        std::cerr << "Usage: aggregate <control-file>\n";
        return EXIT_FAILURE;
    }

    aggregate::Box box;
    box.readControlRecords(argv[1]);
    box.writePackingFile();

    std::cout << "Packing written to " << box.giveOutputFileName() << "\n";
    return EXIT_SUCCESS;
}
