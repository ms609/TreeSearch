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
#include <algorithm>
#include <vector>

namespace ts {

// A hierarchy block describes one controlling primary + its secondaries.
struct HierarchyBlock {
  int primary_char;               // original character index (0-based)
  std::vector<int> secondary_chars; // original character indices (0-based)
  int n_secondaries;              // = secondary_chars.size()
};

// Score a tree under the HSJ dissimilarity-metric criterion.
//
// Characters referenced by hierarchy_blocks are scored via the HSJ
// algorithm; all other characters are scored via standard Fitch.
//
// tree: must be in valid postorder; states will be modified.
// ds: dataset (used for non-hierarchy characters)
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

} // namespace ts

#endif // TS_HSJ_H
