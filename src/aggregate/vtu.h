#pragma once

#include <string>

namespace aggregate {

class Box;

/**
 * Write the box's real inclusions to a ParaView-compatible VTU file.
 *
 * Each ellipsoid is tessellated into a (`meshSize`+1)×(`meshSize`+1) grid
 * of points and `meshSize`² quad cells (matching the Matlab visualisation
 * routine). Each fibre becomes 2 points and 1 line cell. Cell-data arrays
 * `kind` (0 for ellipsoid surface, 1 for fibre) and `id` (the inclusion's
 * sequence number) let ParaView filter and colour by inclusion type.
 */
void writeVtu(const Box &box, const std::string &fileName, int meshSize = 10);

} // namespace aggregate
