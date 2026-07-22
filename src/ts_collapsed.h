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

/// Aggressive (TNT `collapse 3`) flags: collapsed[c] == 1 means the edge from c
/// to parent[c] has MINIMUM POSSIBLE length 0 — there exists a most-parsimonious
/// reconstruction with no state change along it (a weaker, more permissive
/// criterion than compute_collapsed_flags, which requires provable score-identity
/// for exact regraft merging).  Using these flags to skip clips/regrafts is a
/// HEURISTIC neighbourhood reduction (Goloboff's asymmetric reachability): scoring
/// stays exact, but some improving moves may be skipped.  For the standard
/// (no-inapplicable) Fitch path the criterion is final_[p] & final_[c] != 0 for
/// every character (validated bit-for-bit against a brute-force MPR oracle,
/// dev/benchmarks/b2_minlength_oracle.R).  For datasets with inapplicable
/// characters it falls back to the conservative compute_collapsed_flags (NA
/// soft-collapse is not yet derived), so NA datasets are unaffected.
///
/// Requires valid state arrays — call after full_rescore / score_tree.
void compute_collapsed_flags_aggressive(
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
