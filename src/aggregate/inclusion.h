#pragma once

#include <iosfwd>
#include <string>

namespace aggregate {

/**
 * Abstract base class for any geometric inclusion that can be placed inside
 * the periodic box (ellipsoid, fibre, sphere, ...).
 *
 * Concrete subclasses own their own geometric state and know how to write
 * themselves to the shared packing-file format. Pairwise intersection tests
 * are implemented as free functions in `intersection.h` rather than as
 * virtuals to avoid double dispatch.
 */
class Inclusion
{
public:
    explicit Inclusion(int n) : number(n) {}
    virtual ~Inclusion() = default;

    /// Sequential id, written to the packing file.
    int giveNumber() const { return number; }

    /// Token written as the first field of a packing-file line ("ellipsoid", "fibre", ...).
    virtual std::string typeName() const = 0;

    /// Append one line to the stream in the shared packing-file format.
    virtual void writeTo(std::ostream &os) const = 0;

protected:
    int number;
};

} // namespace aggregate
