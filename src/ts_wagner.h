#ifndef TS_WAGNER_H
#define TS_WAGNER_H

// Wagner tree construction: greedy taxon addition for parsimony.
//
// Builds a tree by adding taxa one at a time at the position that
// minimizes the parsimony score. Used to produce starting trees for
// heuristic search (TBR, ratchet, etc.).

#include "ts_data.h"
#include "ts_tree.h"
#include <vector>

namespace ts {

struct WagnerResult {
  double score;
};

// Build a Wagner tree by greedy addition.
// `tree` is populated in place (overwritten).
// `addition_order`: tip indices in insertion order (size n_tips).
//   If empty, uses sequential order 0..n_tips-1.
WagnerResult wagner_tree(TreeState& tree, const DataSet& ds,
                         const std::vector<int>& addition_order);

// Build a random-addition-sequence Wagner tree.
// Uses R's RNG (respects set.seed).
WagnerResult random_wagner_tree(TreeState& tree, const DataSet& ds);

} // namespace ts

#endif // TS_WAGNER_H
