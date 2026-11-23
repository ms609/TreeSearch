#include "ts_pool.h"
#include "ts_constraint.h"  // displays_forbidden_clade (negative-constraint filter)
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <unordered_map>
#include <unordered_set>

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

  // Reject any tree displaying a forbidden clade (converse-constraint search).
  if (forbidden_ && displays_forbidden_clade(tree, *forbidden_)) {
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

  // If pool is full, evict an entry to make room
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
    if (score > worst_score) return false;  // new tree worse than all

    if (score < worst_score) {
      // Strictly better than worst — just evict worst
      entries_[worst_idx] = PoolEntry{tree, score, hash, std::move(ss)};
    } else {
      // Tied with worst (common when pool has converged to best score).
      // Diversity-aware eviction: among entries at worst_score, evict the
      // one most similar to the new tree (most shared splits).
      int wps = ss.words_per_split;
      int evict_idx = worst_idx;
      if (wps > 0 && ss.n_splits > 0) {
        std::unordered_set<uint64_t> new_hashes;
        new_hashes.reserve(ss.n_splits);
        for (int s = 0; s < ss.n_splits; ++s) {
          new_hashes.insert(hash_single_split(ss.split(s), wps));
        }

        int max_shared = -1;
        for (int i = 0; i < static_cast<int>(entries_.size()); ++i) {
          if (entries_[i].score != worst_score) continue;
          int shared = 0;
          for (int s = 0; s < entries_[i].splits.n_splits; ++s) {
            uint64_t sh = hash_single_split(
                entries_[i].splits.split(s), wps);
            if (new_hashes.count(sh)) ++shared;
          }
          if (shared > max_shared) {
            max_shared = shared;
            evict_idx = i;
          }
        }
      }
      entries_[evict_idx] = PoolEntry{tree, score, hash, std::move(ss)};
    }
    return true;
  }

  entries_.push_back(PoolEntry{tree, score, hash, std::move(ss)});
  return true;
}

bool TreePool::add_collapsed(const TreeState& tree, double score,
                             const std::vector<uint8_t>& collapsed) {
  if (collapsed.empty()) return add(tree, score);

  // Reject if worse than threshold
  if (!entries_.empty() && score > best_score_ + suboptimal) {
    return false;
  }

  // Reject any tree displaying a forbidden clade (converse-constraint search).
  if (forbidden_ && displays_forbidden_clade(tree, *forbidden_)) {
    return false;
  }

  // Compute collapsed splits (skip zero-length edges) for dedup
  SplitSet css = compute_collapsed_splits(tree, collapsed);
  uint64_t chash = hash_splits(css);

  if (is_duplicate(chash, css)) {
    if (score <= best_score_) ++hits_to_best_;
    return false;
  }

  // New best score found
  if (score < best_score_) {
    best_score_ = score;
    hits_to_best_ = 1;
    evict();
  } else if (score <= best_score_) {
    ++hits_to_best_;
  }

  // If pool is full, evict an entry to make room
  if (static_cast<int>(entries_.size()) >= max_size) {
    int worst_idx = 0;
    double worst_score = entries_[0].score;
    for (int i = 1; i < static_cast<int>(entries_.size()); ++i) {
      if (entries_[i].score > worst_score) {
        worst_score = entries_[i].score;
        worst_idx = i;
      }
    }
    if (score > worst_score) return false;

    if (score < worst_score) {
      entries_[worst_idx] =
          PoolEntry{tree, score, chash, std::move(css)};
    } else {
      // Diversity-aware eviction using collapsed splits
      int wps = css.words_per_split;
      int evict_idx = worst_idx;
      if (wps > 0 && css.n_splits > 0) {
        std::unordered_set<uint64_t> new_hashes;
        new_hashes.reserve(css.n_splits);
        for (int s = 0; s < css.n_splits; ++s) {
          new_hashes.insert(hash_single_split(css.split(s), wps));
        }
        int max_shared = -1;
        for (int i = 0; i < static_cast<int>(entries_.size()); ++i) {
          if (entries_[i].score != worst_score) continue;
          int shared = 0;
          for (int s = 0; s < entries_[i].splits.n_splits; ++s) {
            uint64_t sh = hash_single_split(
                entries_[i].splits.split(s), wps);
            if (new_hashes.count(sh)) ++shared;
          }
          if (shared > max_shared) {
            max_shared = shared;
            evict_idx = i;
          }
        }
      }
      entries_[evict_idx] =
          PoolEntry{tree, score, chash, std::move(css)};
    }
    return true;
  }

  entries_.push_back(PoolEntry{tree, score, chash, std::move(css)});
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

std::vector<uint64_t> TreePool::extract_consensus_splits(
    int& n_unanimous, int& words_per_split) const {
  n_unanimous = 0;
  words_per_split = 0;

  // Collect best-score entries
  std::vector<const PoolEntry*> best_entries;
  for (const auto& e : entries_) {
    if (e.score <= best_score_) {
      best_entries.push_back(&e);
    }
  }
  int n_best = static_cast<int>(best_entries.size());
  if (n_best < 2) return {};

  int wps = best_entries[0]->splits.words_per_split;
  words_per_split = wps;

  // Count per-split-hash occurrences, keeping one representative bitset
  struct SplitRecord {
    const uint64_t* bits;
    int count;
  };
  std::unordered_map<uint64_t, SplitRecord> seen;

  for (const auto* e : best_entries) {
    for (int s = 0; s < e->splits.n_splits; ++s) {
      const uint64_t* bits = e->splits.split(s);
      uint64_t sh = hash_single_split(bits, wps);
      auto it = seen.find(sh);
      if (it == seen.end()) {
        seen[sh] = SplitRecord{bits, 1};
      } else {
        ++(it->second.count);
      }
    }
  }

  // Collect splits present in all best-score trees
  std::vector<uint64_t> result;
  for (const auto& kv : seen) {
    if (kv.second.count >= n_best) {
      result.insert(result.end(), kv.second.bits,
                    kv.second.bits + wps);
      ++n_unanimous;
    }
  }

  return result;
}

} // namespace ts
