#ifndef TS_SPLITS_H
#define TS_SPLITS_H

// Split (bipartition) computation and hashing for phylogenetic trees.
//
// A split is the bipartition of tips induced by removing an internal edge.
// Represented as bitsets: for T tips, each split uses ceil(T/64) uint64_t words.
// Canonical form: tip 0 always in the "0" partition (bit 0 unset).
//
// An unrooted tree with T tips has T-3 non-trivial splits.

#include "ts_tree.h"
#include <vector>
#include <cstdint>

namespace ts {

struct SplitSet {
  int n_tips;
  int words_per_split;             // ceil(n_tips / 64)
  std::vector<uint64_t> splits;    // n_splits * words_per_split, contiguous
  int n_splits;                    // typically n_tips - 3

  const uint64_t* split(int i) const {
    return &splits[static_cast<size_t>(i) * words_per_split];
  }

  uint64_t* split(int i) {
    return &splits[static_cast<size_t>(i) * words_per_split];
  }
};

// Compute the set of non-trivial splits from a rooted TreeState.
// The tree is treated as unrooted (root edge splits are excluded).
SplitSet compute_splits(const TreeState& tree);

// Order-independent hash of a SplitSet.
// Two trees with identical split sets will produce the same hash regardless
// of split ordering.
uint64_t hash_splits(const SplitSet& ss);

// Exact equality check: two SplitSets represent the same unrooted topology
// iff they contain the same set of splits.
bool splits_equal(const SplitSet& a, const SplitSet& b);

} // namespace ts

#endif // TS_SPLITS_H
