#include "vtu.h"

#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <ios>
#include <limits>
#include <vector>

#include "aggregateerror.h"
#include "box.h"
#include "ellipsoid.h"
#include "fibre.h"

namespace aggregate {

namespace {

constexpr int kKindEllipsoid = 0;
constexpr int kKindFibre = 1;

constexpr unsigned char kVtkQuad = 9;
constexpr unsigned char kVtkLine = 3;

/// Append the (meshSize+1)² parametric points of one ellipsoid to `points`.
void appendEllipsoidPoints(const Ellipsoid &e, int meshSize,
                           std::vector<Eigen::Vector3d> &points)
{
    const Eigen::Matrix3d R = e.giveRotation();
    const Eigen::Vector3d &c = e.giveCentre();
    const Eigen::Vector3d &s = e.giveSemiAxes();
    for ( int j = 0; j <= meshSize; ++j ) {
        const double t2 = static_cast<double>( j ) * M_PI / meshSize;
        const double sinT2 = std::sin(t2);
        const double cosT2 = std::cos(t2);
        for ( int i = 0; i <= meshSize; ++i ) {
            const double t1 = static_cast<double>( i ) * 2.0 * M_PI / meshSize;
            const Eigen::Vector3d local(s(0) * std::cos(t1) * sinT2,
                                        s(1) * std::sin(t1) * sinT2,
                                        s(2) * cosT2);
            points.push_back(c + R * local);
        }
    }
}

} // namespace

void writeVtu(const Box &box, const std::string &fileName, int meshSize)
{
    if ( meshSize < 1 ) {
        errorf("writeVtu: meshSize must be ≥ 1, got %d", meshSize);
    }

    std::vector<Eigen::Vector3d> points;
    std::vector<int> connectivity;
    std::vector<int> offsets;
    std::vector<unsigned char> cellTypes;
    std::vector<int> cellKind;
    std::vector<int> cellId;

    const int patchPoints = ( meshSize + 1 ) * ( meshSize + 1 );

    for ( const auto &incPtr : box.giveRealInclusions() ) {
        if ( const auto *e = dynamic_cast<const Ellipsoid *>( incPtr.get() ) ) {
            const int basePoint = static_cast<int>( points.size() );
            appendEllipsoidPoints(*e, meshSize, points);
            for ( int j = 0; j < meshSize; ++j ) {
                for ( int i = 0; i < meshSize; ++i ) {
                    const int p00 = basePoint + j * ( meshSize + 1 ) + i;
                    const int p10 = basePoint + j * ( meshSize + 1 ) + ( i + 1 );
                    const int p11 = basePoint + ( j + 1 ) * ( meshSize + 1 ) + ( i + 1 );
                    const int p01 = basePoint + ( j + 1 ) * ( meshSize + 1 ) + i;
                    connectivity.push_back(p00);
                    connectivity.push_back(p10);
                    connectivity.push_back(p11);
                    connectivity.push_back(p01);
                    offsets.push_back(static_cast<int>( connectivity.size() ));
                    cellTypes.push_back(kVtkQuad);
                    cellKind.push_back(kKindEllipsoid);
                    cellId.push_back(e->giveNumber());
                }
            }
            (void)patchPoints;
        } else if ( const auto *f = dynamic_cast<const Fibre *>( incPtr.get() ) ) {
            const int basePoint = static_cast<int>( points.size() );
            points.push_back(f->giveEndpointA());
            points.push_back(f->giveEndpointB());
            connectivity.push_back(basePoint);
            connectivity.push_back(basePoint + 1);
            offsets.push_back(static_cast<int>( connectivity.size() ));
            cellTypes.push_back(kVtkLine);
            cellKind.push_back(kKindFibre);
            cellId.push_back(f->giveNumber());
        }
    }

    std::ofstream out(fileName);
    if ( !out ) {
        errorf("writeVtu: cannot open '%s' for writing", fileName.c_str());
    }
    out << std::scientific << std::setprecision(std::numeric_limits<float>::max_digits10);

    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "  <UnstructuredGrid>\n"
        << "    <Piece NumberOfPoints=\"" << points.size()
        << "\" NumberOfCells=\"" << cellTypes.size() << "\">\n";

    out << "      <Points>\n"
        << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for ( const auto &p : points ) {
        out << "          " << static_cast<float>( p(0) ) << ' '
            << static_cast<float>( p(1) ) << ' '
            << static_cast<float>( p(2) ) << '\n';
    }
    out << "        </DataArray>\n"
        << "      </Points>\n";

    out << "      <Cells>\n"
        << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for ( int v : connectivity ) {
        out << ' ' << v;
    }
    out << "\n        </DataArray>\n"
        << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for ( int v : offsets ) {
        out << ' ' << v;
    }
    out << "\n        </DataArray>\n"
        << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for ( unsigned char v : cellTypes ) {
        out << ' ' << static_cast<int>( v );
    }
    out << "\n        </DataArray>\n"
        << "      </Cells>\n";

    out << "      <CellData>\n"
        << "        <DataArray type=\"Int32\" Name=\"kind\" format=\"ascii\">\n";
    for ( int v : cellKind ) {
        out << ' ' << v;
    }
    out << "\n        </DataArray>\n"
        << "        <DataArray type=\"Int32\" Name=\"id\" format=\"ascii\">\n";
    for ( int v : cellId ) {
        out << ' ' << v;
    }
    out << "\n        </DataArray>\n"
        << "      </CellData>\n";

    out << "    </Piece>\n"
        << "  </UnstructuredGrid>\n"
        << "</VTKFile>\n";
}

} // namespace aggregate
