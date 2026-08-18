#include "disk.h"

#include <ostream>

namespace aggregate {

Disk::Disk(int n, const Eigen::Vector2d &centre_, double radius_)
    : Inclusion(n), centre(centre_), radius(radius_)
{}

void Disk::writeTo(std::ostream &os) const
{
    os << "disk " << number
       << " centre 2 " << centre(0) << ' ' << centre(1)
       << " radius " << radius
       << '\n';
}

} // namespace aggregate
