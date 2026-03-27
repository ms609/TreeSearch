#ifndef TS_SPIC_H
#define TS_SPIC_H

// SPIC (Splitwise Phylogenetic Information Content) scoring for consensus search.
//
// Scores a candidate tree by summing the phylogenetic information content of
// each split, weighted by its frequency in the input tree set (interpreted as
// the probability that the split is correct).
//
// Reference: Smith (2022) Systematic Biology syab099.
//
// Key advantage over MCI (Mutual Clustering Information):
//   O(splits * n_bins) per evaluation vs O(splits^2) for MCI (no LAP solver).
//   Split-additive: each split contributes independently to the total score.

#include "ts_data.h"
#include "ts_tree.h"
#include "ts_cid.h"    // reuse CidSplitSet, compute_splits_cid
#include <vector>
#include <cstdint>
#include <cmath>
#include <unordered_map>

namespace ts {

// --------------------------------------------------------------------------
// SpicData: precomputed data for SPIC scoring
// --------------------------------------------------------------------------
struct SpicData {
  int n_tips;
  int n_bins;      // ceil(n_tips / 64)
  int n_trees;     // number of input trees

  // Split frequency table: maps FNV-1a hash of split data to
  // (count, in_split) for O(1) lookup per candidate split.
  // Collisions handled by storing a list of entries per hash bucket.
  struct SplitEntry {
    int count;        // number of input trees containing this split
    int in_split;     // popcount (partition size a)
    // For collision resolution, store the split data
    std::vector<uint64_t> data;  // n_bins words
  };
  std::unordered_multimap<uint64_t, SplitEntry> freq_table;

  // Precomputed log2 tables for the IC formula.
  // log2_rooted[k] = Log2Rooted(k) = log2((2k-3)!!) for k >= 2; 0 for k <= 1.
  // Indexed from 0..n_tips.
  std::vector<double> log2_rooted;

  // log2_unrooted_n = Log2Unrooted(n_tips)
  double log2_unrooted_n;

  // Precomputed log2(P_consistent(a)) for each possible partition size a.
  // = Log2Rooted(a) + Log2Rooted(n-a) - Log2Unrooted(n)
  // Indexed from 0..n_tips.
  std::vector<double> log2_p_consistent;

  // Precomputed log2(1 - P_consistent(a)) for each possible partition size a.
  // Indexed from 0..n_tips.  May be -Inf for trivial splits where P = 1.
  std::vector<double> log2_p_inconsistent;

  // Persistent candidate buffers (reused across spic_score calls)
  mutable std::vector<uint64_t> cand_tip_bits;
  mutable CidSplitSet cand_buf;
};

// --------------------------------------------------------------------------
// Public API
// --------------------------------------------------------------------------

// Build SpicData from input tree split sets (already stored in CidData).
// Builds the frequency table and precomputes log2 tables.
SpicData build_spic_data(const CidData& cd);

// Build SpicData directly from a list of CidSplitSets.
// When n_bins_override > 0, use that for the bit-width of split data
// instead of computing from n_tips.  This is needed when the split data
// was produced by bit-masking (e.g. rogue prescreen) and retains the
// original wider representation.
SpicData build_spic_data_from_splits(const std::vector<CidSplitSet>& tree_splits,
                                      int n_tips, int n_trees,
                                      int n_bins_override = 0);

// Score a candidate tree using SPIC.
// Returns negated SPIC sum (lower = better consensus), consistent with
// the cid_score() sign convention.
double spic_score(TreeState& tree, const SpicData& sd);

// Per-split IC contribution (for debugging / diagnostics).
// Returns the IC of a split with partition size `a` (b = n - a),
// support frequency `freq` (in [0, 1]), using precomputed tables.
double split_ic(int a, int n, double freq, const SpicData& sd);

} // namespace ts

#endif // TS_SPIC_H
