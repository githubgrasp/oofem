#pragma once

#include <string>

namespace aggregate {

class Box;

/**
 * Write the box's inclusions to a ParaView-compatible VTU file.
 *
 * Each ellipsoid is tessellated into a (`meshSize`+1)×(`meshSize`+1) grid
 * of points and `meshSize`² quad cells (matching the Matlab visualisation
 * routine). Each fibre becomes 2 points and 1 line cell. Cell-data arrays
 * `kind` (0 for ellipsoid surface, 1 for fibre), `id` (the inclusion's
 * sequence number), and `ghost` (0 for a real inclusion, 1 for a periodic
 * image) let ParaView filter and colour by type, identity, and real/ghost.
 *
 * When `includeGhosts` is true (the default) the periodic ghost images are
 * drawn as well, so the VTU shows the packing tiling exactly as `packing.dat`
 * describes it. Passing false draws only the real inclusions. This is purely
 * a visualisation choice: the ghosts always exist on the box (for overlap
 * tests) and are always written to `packing.dat` regardless.
 */
void writeVtu(const Box &box, const std::string &fileName, int meshSize = 10,
              bool includeGhosts = true);

} // namespace aggregate
