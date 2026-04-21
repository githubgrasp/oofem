#include "box.h"

#include <fstream>
#include <sstream>

#include "aggregateerror.h"

namespace aggregate {

Box::Box()
{
    dimensions.resize(3);
    dimensions.zero();
    periodicityFlag.resize(3);
    periodicityFlag.zero();
}

void Box::readControlRecords(const std::string &fileName)
{
    std::ifstream stream(fileName);
    if ( !stream ) {
        errorf("Box::readControlRecords: cannot open '%s'", fileName.c_str());
    }

    std::string line;
    while ( std::getline(stream, line) ) {
        if ( line.size() < 2 || line[0] != '#' || line[1] != '@' ) {
            continue;
        }
        applyDirective(line.substr(2));
    }

    if ( outputFileName.empty() ) {
        error("Box::readControlRecords: missing '#@output' directive");
    }
}

void Box::applyDirective(const std::string &line)
{
    std::istringstream iss(line);
    std::string keyword;
    if ( !( iss >> keyword ) ) {
        return;
    }

    if ( keyword == "output" ) {
        iss >> outputFileName;
    } else if ( keyword == "box" ) {
        int n;
        iss >> n;
        if ( n != 3 ) {
            errorf("Box::applyDirective: '#@box' expects 3 components, got %d", n);
        }
        dimensions.resize(3);
        for ( int i = 1; i <= 3; ++i ) {
            iss >> dimensions.at(i);
        }
    } else if ( keyword == "periodicity" ) {
        int n;
        iss >> n;
        if ( n != 3 ) {
            errorf("Box::applyDirective: '#@periodicity' expects 3 flags, got %d", n);
        }
        periodicityFlag.resize(3);
        for ( int i = 1; i <= 3; ++i ) {
            iss >> periodicityFlag.at(i);
        }
    } else if ( keyword == "seed" ) {
        iss >> randomSeed;
    } else if ( keyword == "maxiter" ) {
        iss >> maximumIterations;
    } else {
        errorf("Box::applyDirective: unknown keyword '%s'", keyword.c_str());
    }
}

void Box::writePackingFile() const
{
    std::ofstream out(outputFileName);
    if ( !out ) {
        errorf("Box::writePackingFile: cannot open '%s' for writing", outputFileName.c_str());
    }
    out << "# packing v1\n";
    out << "# box " << dimensions.at(1) << ' '
        << dimensions.at(2) << ' '
        << dimensions.at(3) << '\n';
    out << "# periodicity " << periodicityFlag.at(1) << ' '
        << periodicityFlag.at(2) << ' '
        << periodicityFlag.at(3) << '\n';
}

} // namespace aggregate
