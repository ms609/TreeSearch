#ifndef TS_POOL_H
#define TS_POOL_H

// TreePool: stores unique trees found during search.
//
// Deduplication via split hashing with full equality fallback.
// Eviction policy: discard trees worse than best + suboptimal threshold.
// Tracks hits_to_best for convergence detection.

#include "ts_tree.h"
#include "ts_data.h"
#include "ts_splits.h"
#include "ts_collapsed.h"
#include <vector>
#include <cstdint>
#include <unordered_set>
#include <unordered_map>

namespace ts {

// Per-split frequency table for conflict-guided sector selection.
// Maps per-split hash → count across best-score pool trees.
struct SplitFrequencyTable {
  std::unordered_map<uint64_t, int> freq;  // split hash → occurrence count
  int n_trees = 0;  // number of best-score trees used to build the table
};

struct PoolEntry {
  TreeState tree;
  double score;
  uint64_t split_hash;
  SplitSet splits;
};

class TreePool {
public:
  int max_size;        // maximum number of entries (e.g. 100)
  double suboptimal;   // keep trees within this many steps of best

  TreePool(int max_sz = 100, double subopt = 0.0)
    : max_size(max_sz), suboptimal(subopt),
      best_score_(1e18), hits_to_best_(0),
      consensus_hash_(0), consensus_unchanged_(0) {}

  // Add a tree if it's not a duplicate and meets score threshold.
  // Returns true if the tree was actually added.
  bool add(const TreeState& tree, double score);

  // Add with collapsed-topology dedup: two trees that collapse to the
  // same polytomy are treated as duplicates. Uses collapsed flags to
  // filter out zero-length-edge splits before hashing.
  bool add_collapsed(const TreeState& tree, double score,
                     const std::vector<uint8_t>& collapsed);

  // Get the best (lowest-scoring) entry.
  const PoolEntry& best() const;

  // Get all entries.
  const std::vector<PoolEntry>& all() const { return entries_; }

  // Evict entries worse than best_score + suboptimal.
  void evict();

  // Number of times the current best score has been independently found.
  int hits_to_best() const { return hits_to_best_; }

  // Override hits_to_best (used to undo inflation from non-independent hits).
  void set_hits_to_best(int n) { hits_to_best_ = n; }

  // Current best score.
  double best_score() const { return best_score_; }

  // Number of entries in the pool.
  int size() const { return static_cast<int>(entries_.size()); }

  // Count entries at exactly the best score.
  int count_at_best() const {
    int count = 0;
    for (const auto& e : entries_) {
      if (e.score == best_score_) ++count;
    }
    return count;
  }

  // Clear the pool.
  void clear();

  // Compute strict-consensus hash of all best-score trees in the pool.
  // A split is "in the consensus" if it appears in ALL best-score trees.
  // Returns an order-independent hash of the consensus split set.
  uint64_t compute_consensus_hash() const;

  // Update consensus stability tracker. Call after each replicate.
  // Returns the number of consecutive replicates where the consensus
  // hash has been unchanged.
  int update_consensus_stability();

  // Current consecutive-unchanged count.
  int consensus_unchanged() const { return consensus_unchanged_; }

  // Compute per-split frequency table across best-score pool trees.
  // For each non-trivial split in any best-score tree, records how many
  // best-score trees contain it. Used by conflict-guided sector selection.
  SplitFrequencyTable compute_split_frequencies() const;

  // Extract splits that appear in ALL best-score pool trees (strict consensus).
  // Returns contiguous bitset data: n_splits * words_per_split uint64_t values.
  // Sets n_unanimous and words_per_split. Returns empty if pool has <2 trees.
  std::vector<uint64_t> extract_consensus_splits(
      int& n_unanimous, int& words_per_split) const;

private:
  std::vector<PoolEntry> entries_;
  double best_score_;
  int hits_to_best_;

  // Consensus stability tracking
  uint64_t consensus_hash_;
  int consensus_unchanged_;

  bool is_duplicate(uint64_t hash, const SplitSet& ss) const;
};

} // namespace ts

#endif // TS_POOL_H
