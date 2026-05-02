#include "box.h"

#include <fstream>
#include <ios>
#include <iomanip>
#include <limits>
#include <sstream>

#include "aggregateerror.h"

namespace aggregate {

Box::Box()
    : dimensions(Eigen::Vector3d::Zero()),
      periodicityFlag(Eigen::Vector3i::Zero())
{}

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
        if ( n != 2 && n != 3 ) {
            errorf("Box::applyDirective: '#@box' expects 2 (2D) or 3 (3D) components, got %d", n);
        }
        dim = n;
        if ( n == 2 ) {
            iss >> dimensions(0) >> dimensions(1);
            dimensions(2) = 0.;
        } else {
            iss >> dimensions(0) >> dimensions(1) >> dimensions(2);
        }
    } else if ( keyword == "periodicity" ) {
        int n;
        iss >> n;
        if ( n != 2 && n != 3 ) {
            errorf("Box::applyDirective: '#@periodicity' expects 2 (2D) or 3 (3D) flags, got %d", n);
        }
        if ( n == 2 ) {
            iss >> periodicityFlag(0) >> periodicityFlag(1);
            periodicityFlag(2) = 0;
        } else {
            iss >> periodicityFlag(0) >> periodicityFlag(1) >> periodicityFlag(2);
        }
    } else if ( keyword == "seed" ) {
        iss >> randomSeed;
    } else if ( keyword == "maxiter" ) {
        iss >> maximumIterations;
    } else if ( keyword == "vtu" ) {
        iss >> vtuFileName;
    } else if ( keyword == "grading" ) {
        std::string sub;
        while ( iss >> sub ) {
            if ( sub == "dmin" ) {
                iss >> grading.dmin;
            } else if ( sub == "dmax" ) {
                iss >> grading.dmax;
            } else if ( sub == "fraction" ) {
                iss >> grading.aggregateFraction;
            } else if ( sub == "shape" ) {
                std::string shapeName;
                iss >> shapeName;
                if ( shapeName == "sphere" ) {
                    grading.shape = GradingShape::Sphere;
                } else if ( shapeName == "ellipsoid" ) {
                    grading.shape = GradingShape::Ellipsoid;
                } else {
                    errorf("Box::applyDirective: '#@grading shape' must be 'sphere' or 'ellipsoid', got '%s'",
                           shapeName.c_str());
                }
            } else {
                errorf("Box::applyDirective: '#@grading' unknown sub-keyword '%s'", sub.c_str());
            }
        }
        grading.present = true;
    } else if ( keyword == "fibres" ) {
        std::string sub;
        while ( iss >> sub ) {
            if ( sub == "fraction" ) {
                iss >> fibres.fibreFraction;
            } else if ( sub == "length" ) {
                iss >> fibres.fibreLength;
            } else if ( sub == "diameter" ) {
                iss >> fibres.fibreDiameter;
            } else {
                errorf("Box::applyDirective: '#@fibres' unknown sub-keyword '%s'", sub.c_str());
            }
        }
        fibres.present = true;
    } else {
        errorf("Box::applyDirective: unknown keyword '%s'", keyword.c_str());
    }
}

void Box::addReal(std::unique_ptr<Inclusion> inclusion)
{
    realInclusions.push_back(std::move(inclusion));
}

void Box::addGhost(std::unique_ptr<Inclusion> inclusion)
{
    ghostInclusions.push_back(std::move(inclusion));
}

void Box::writePackingFile() const
{
    std::ofstream out(outputFileName);
    if ( !out ) {
        errorf("Box::writePackingFile: cannot open '%s' for writing", outputFileName.c_str());
    }
    // Lossless decimal formatting — guarantees exact-bit reproducibility
    // of doubles for diff-based regression tests, regardless of locale.
    out << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << "# packing v1\n";
    if ( dim == 2 ) {
        out << "# box " << dimensions(0) << ' ' << dimensions(1) << '\n';
        out << "# periodicity " << periodicityFlag(0) << ' '
            << periodicityFlag(1) << '\n';
    } else {
        out << "# box " << dimensions(0) << ' '
            << dimensions(1) << ' '
            << dimensions(2) << '\n';
        out << "# periodicity " << periodicityFlag(0) << ' '
            << periodicityFlag(1) << ' '
            << periodicityFlag(2) << '\n';
    }
    for ( const auto &inc : realInclusions ) {
        inc->writeTo(out);
    }
}

} // namespace aggregate
