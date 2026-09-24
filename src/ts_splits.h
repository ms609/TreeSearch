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

// Compute the set of non-trivial splits, skipping collapsed edges.
// Edges where collapsed[node] == 1 are excluded.  Two binary trees that
// differ only in zero-length resolutions produce the same collapsed split
// set.  If collapsed is empty, falls back to compute_splits().
SplitSet compute_collapsed_splits(const TreeState& tree,
                                  const std::vector<uint8_t>& collapsed);

// Order-independent hash of a SplitSet.
// Two trees with identical split sets will produce the same hash regardless
// of split ordering.
uint64_t hash_splits(const SplitSet& ss);

// Exact equality check: two SplitSets represent the same unrooted topology
// iff they contain the same set of splits.
bool splits_equal(const SplitSet& a, const SplitSet& b);

// Lightweight topology hash computed in a single postorder pass.
// Equivalent to hash_splits(compute_splits(tree)) but avoids allocating
// a SplitSet. Requires tree.postorder to be valid.
uint64_t hash_tree(const TreeState& tree);

// FNV-1a hash of a single canonicalized split bitset.
// Used by the pool for per-split frequency counting and consensus hashing.
inline uint64_t hash_single_split(const uint64_t* s, int wps) {
  uint64_t h = 0xcbf29ce484222325ULL; // FNV offset basis
  for (int w = 0; w < wps; ++w) {
    h ^= s[w];
    h *= 0x100000001b3ULL; // FNV prime
  }
  return h;
}

} // namespace ts

#endif // TS_SPLITS_H
