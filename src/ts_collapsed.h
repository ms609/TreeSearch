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

/// Collapsed-region information for regraft merging.
///
/// After compute_collapsed_regions(), every node that lies on a collapsed
/// edge has region_id[node] >= 0 identifying its connected component.
/// Nodes on non-collapsed edges have region_id[node] == -1.
///
/// A "collapsed region" is a maximal connected set of nodes linked by
/// collapsed edges. Within a region, all regraft positions produce the
/// same indirect-evaluation score, so only one representative per region
/// need be evaluated.
struct CollapsedRegions {
  std::vector<uint8_t> collapsed;  ///< per-node collapsed flag
  std::vector<int> region_id;     ///< per-node region assignment (-1 = not collapsed)
  int n_regions = 0;              ///< total number of collapsed regions
  int n_collapsed = 0;            ///< total collapsed edges
};

/// Compute collapsed flags and group connected collapsed edges into regions.
///
/// Two nodes share a region if they are connected by a parent-child edge
/// where the child's collapsed flag is set (the edge from child to parent
/// is zero-length). Both parent and child receive the same region_id.
///
/// Requires valid state arrays — call after full_rescore / score_tree.
void compute_collapsed_regions(
    const TreeState& tree,
    const DataSet& ds,
    CollapsedRegions& info);

} // namespace ts

#endif // TS_COLLAPSED_H
