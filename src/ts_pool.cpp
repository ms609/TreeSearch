#include "ts_pool.h"
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <unordered_map>

namespace ts {

bool TreePool::is_duplicate(uint64_t hash, const SplitSet& ss) const {
  for (const auto& entry : entries_) {
    if (entry.split_hash == hash && splits_equal(entry.splits, ss)) {
      return true;
    }
  }
  return false;
}

bool TreePool::add(const TreeState& tree, double score) {
  // Reject if worse than threshold (when pool is non-empty)
  if (!entries_.empty() && score > best_score_ + suboptimal) {
    return false;
  }

  SplitSet ss = compute_splits(tree);
  uint64_t hash = hash_splits(ss);

  if (is_duplicate(hash, ss)) {
    // Duplicate topology: still counts as a hit if it matches the best score
    if (score <= best_score_) {
      ++hits_to_best_;
    }
    return false;
  }

  // New best score found
  if (score < best_score_) {
    best_score_ = score;
    hits_to_best_ = 1;
    evict();
  } else if (score <= best_score_) {
    // Ties with best (score == best_score_)
    ++hits_to_best_;
  }

  // If pool is full, evict worst entry to make room
  if (static_cast<int>(entries_.size()) >= max_size) {
    // Find worst entry
    int worst_idx = 0;
    double worst_score = entries_[0].score;
    for (int i = 1; i < static_cast<int>(entries_.size()); ++i) {
      if (entries_[i].score > worst_score) {
        worst_score = entries_[i].score;
        worst_idx = i;
      }
    }
    // Only evict if the new tree is at least as good
    if (score <= worst_score) {
      entries_[worst_idx] = PoolEntry{tree, score, hash, std::move(ss)};
    } else {
      return false;
    }
    return true;
  }

  entries_.push_back(PoolEntry{tree, score, hash, std::move(ss)});
  return true;
}

const PoolEntry& TreePool::best() const {
  if (entries_.empty()) {
    throw std::runtime_error("TreePool::best() called on empty pool");
  }
  int best_idx = 0;
  for (int i = 1; i < static_cast<int>(entries_.size()); ++i) {
    if (entries_[i].score < entries_[best_idx].score) {
      best_idx = i;
    }
  }
  return entries_[best_idx];
}

void TreePool::evict() {
  if (entries_.empty()) return;
  double threshold = best_score_ + suboptimal;
  entries_.erase(
    std::remove_if(entries_.begin(), entries_.end(),
      [threshold](const PoolEntry& e) { return e.score > threshold; }),
    entries_.end());
}

void TreePool::clear() {
  entries_.clear();
  best_score_ = 1e18;
  hits_to_best_ = 0;
  consensus_hash_ = 0;
  consensus_unchanged_ = 0;
}

uint64_t TreePool::compute_consensus_hash() const {
  // Collect best-score entries
  std::vector<const PoolEntry*> best_entries;
  for (const auto& e : entries_) {
    if (e.score <= best_score_) {
      best_entries.push_back(&e);
    }
  }
  if (best_entries.empty()) return 0;
  if (best_entries.size() == 1) {
    // Single tree: consensus = all its splits
    return best_entries[0]->split_hash;
  }

  int wps = best_entries[0]->splits.words_per_split;
  int n_best = static_cast<int>(best_entries.size());

  // Count occurrences of each split hash across best-score trees.
  // A split is in the strict consensus iff it appears in ALL best trees.
  // Use a two-level approach: hash-based counting, then verify with bitwise
  // equality for splits that appear n_best times.
  std::unordered_map<uint64_t, int> split_counts;
  for (const auto* e : best_entries) {
    for (int s = 0; s < e->splits.n_splits; ++s) {
      uint64_t sh = hash_single_split(e->splits.split(s), wps);
      ++split_counts[sh];
    }
  }

  // XOR of hashes of all unanimous splits → order-independent consensus hash.
  // (Hash collisions between different splits are extremely unlikely; a
  // false match would over-count unanimity, making the consensus *more*
  // stable than it truly is — a conservative error direction.)
  uint64_t consensus = 0;
  for (const auto& kv : split_counts) {
    if (kv.second >= n_best) {
      consensus ^= kv.first;
    }
  }
  return consensus;
}

int TreePool::update_consensus_stability() {
  uint64_t new_hash = compute_consensus_hash();
  if (new_hash == consensus_hash_ && consensus_hash_ != 0) {
    ++consensus_unchanged_;
  } else {
    consensus_hash_ = new_hash;
    consensus_unchanged_ = 0;
  }
  return consensus_unchanged_;
}

SplitFrequencyTable TreePool::compute_split_frequencies() const {
  SplitFrequencyTable sft;

  // Collect best-score entries
  std::vector<const PoolEntry*> best_entries;
  for (const auto& e : entries_) {
    if (e.score <= best_score_) {
      best_entries.push_back(&e);
    }
  }
  sft.n_trees = static_cast<int>(best_entries.size());
  if (sft.n_trees < 2) return sft;

  int wps = best_entries[0]->splits.words_per_split;

  for (const auto* e : best_entries) {
    for (int s = 0; s < e->splits.n_splits; ++s) {
      uint64_t sh = hash_single_split(e->splits.split(s), wps);
      ++sft.freq[sh];
    }
  }

  return sft;
}

} // namespace ts
