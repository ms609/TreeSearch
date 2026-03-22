#ifndef TS_COLLAPSED_H
#define TS_COLLAPSED_H

#include "ts_tree.h"
#include "ts_data.h"
#include <vector>
#include <cstdint>

namespace ts {

/// Compute collapsed flags for all nodes.
///
/// collapsed[c] == 1 means the edge from c to parent[c] is zero-length
/// and clipping c cannot improve the score (safe to skip as TBR/SPR clip).
///
/// Requires valid state arrays — call after full_rescore / score_tree.
void compute_collapsed_flags(
    const TreeState& tree,
    const DataSet& ds,
    std::vector<uint8_t>& collapsed);

} // namespace ts

#endif // TS_COLLAPSED_H
