#ifndef NBT_IO_HXX
#define NBT_IO_HXX

#include "LoopBuilder.hxx" // for Loop, OrientedCurve
#include "RegionBuilder.hxx" // for Region
#include "Nurbs.hxx"
#include <string>
#include <vector>

namespace NbtIO {
// Write topology (.nbt) referencing a .nbs file.
// Format (v1):
// * XF-stream NBT topology
// nbs_file <relative_or_name>
// loops <L>               # number of loops
// loop <i> curves <c0> <c1> ... <ck>   # use negative id for reversed orientation
// regions <R>
// region <i> outer <loop_index> inners <l0> <l1> ...
// end
// curves: geometry curves whose ids are referenced by loops.
// Each Nurbs must have id() >= 0.
bool writeFile(const std::string& path,
               const std::string& nbsFileName,
               const std::vector<Nurbs>& curves,
               const std::vector<Loop>& loops,
               const std::vector<Region>& regions,
               std::string* errorMessage = nullptr);

// Read topology file; does not load geometry. Returns loops & regions.
// Reconstruct loops as sequences of OrientedCurve with curve indices looked up
// by matching ids from provided curves vector.
bool readFile(const std::string& path,
              const std::vector<Nurbs>& curves,
              std::string& outNbsFileName,
              std::vector<Loop>& outLoops,
              std::vector<Region>& outRegions,
              std::string* errorMessage = nullptr);
}

#endif // NBT_IO_HXX
