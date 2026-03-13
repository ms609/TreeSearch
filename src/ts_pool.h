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
#include <vector>
#include <cstdint>

namespace ts {

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
      best_score_(1e18), hits_to_best_(0) {}

  // Add a tree if it's not a duplicate and meets score threshold.
  // Returns true if the tree was actually added.
  bool add(const TreeState& tree, double score);

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

  // Clear the pool.
  void clear();

private:
  std::vector<PoolEntry> entries_;
  double best_score_;
  int hits_to_best_;

  bool is_duplicate(uint64_t hash, const SplitSet& ss) const;
};

} // namespace ts

#endif // TS_POOL_H
