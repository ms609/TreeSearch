#include "ts_pool.h"
#include <algorithm>
#include <stdexcept>

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
}

} // namespace ts
