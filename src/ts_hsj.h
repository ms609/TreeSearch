#ifndef TS_HSJ_H
#define TS_HSJ_H

// Hopkins & St. John (2021) scoring for hierarchical characters.
//
// Each hierarchy block consists of a controlling primary character
// (absent/present) and m secondary characters that are applicable only
// when the primary is present. The contribution of each block to the
// tree score is computed via a modified Fitch traversal tracking:
//   a(n) = min score at node n if controlling primary is ABSENT
//   p(n) = min score at node n if controlling primary is PRESENT
//
// Secondary character labels at internal nodes are obtained via a
// standard Fitch first-pass (with inapplicable treated as a separate
// state).
//
// Non-hierarchy characters are scored via standard Fitch.

#include "ts_data.h"
#include "ts_tree.h"
#include <vector>

namespace ts {

// Compute adjusted pattern weights that exclude hierarchy characters.
//
// Given a hierarchy specification and the phyDat index (mapping original
// characters to pattern indices), returns a copy of weight_r with hierarchy
// characters' contributions subtracted.  Patterns that only appear in
// hierarchy characters will have weight 0 and be dropped by build_dataset().
//
// index_r:    n_orig_chars vector; index_r[c] = pattern index (0-based) for
//             original character c.
// weight_r:   n_patterns vector; pattern frequencies (original weights).
// hierarchy_chars: vector of original character indices (0-based) belonging
//                  to any hierarchy block.
// n_patterns: number of unique patterns.
//
// Returns: adjusted weight vector of length n_patterns.
std::vector<int> partition_weights(
    const int* index_r, int n_orig_chars,
    const int* weight_r, int n_patterns,
    const std::vector<int>& hierarchy_chars);

// Score a tree under the HSJ dissimilarity-metric criterion.
//
// Characters referenced by hierarchy_blocks are scored via the HSJ
// algorithm; all other characters are scored via standard Fitch.
//
// tree: must be in valid postorder; states will be modified.
// ds: dataset (used for non-hierarchy characters AND hierarchy data)
// hierarchy_blocks: hierarchy specification
// alpha: HSJ scaling parameter in [0, 1]
// tip_labels: per-tip, per-original-char state labels (0-based token index).
//   Layout: tip_labels[tip * n_orig_chars + char]. This is the full
//   (uncompressed) original matrix needed for secondary character matching.
// n_orig_chars: number of original characters (before compression)
//
// Returns total tree score (Fitch for non-hierarchy + HSJ for hierarchy).
double hsj_score(
    TreeState& tree,
    const DataSet& ds,
    const std::vector<HierarchyBlock>& hierarchy_blocks,
    double alpha,
    const std::vector<int>& tip_labels,
    int n_orig_chars);

// Convenience overload using hierarchy data stored in DataSet.
// Requires ds.scoring_mode == HSJ with hierarchy fields populated.
double hsj_score(TreeState& tree, const DataSet& ds);

} // namespace ts

#endif // TS_HSJ_H
