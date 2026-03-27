// SPIC (Splitwise Phylogenetic Information Content) scoring.
//
// See ts_spic.h for overview and API.

#include "ts_spic.h"
#include <cstring>      // memcmp

namespace ts {

// ==========================================================================
// Precomputed log2 tables
// ==========================================================================

// Log2Rooted(k) = log2((2k-3)!!) = sum(log2(2i-1), i=1..k-1) for k >= 2.
// Log2Rooted(0) = Log2Rooted(1) = 0.
static std::vector<double> build_log2_rooted(int n) {
  std::vector<double> tbl(n + 1, 0.0);
  for (int k = 2; k <= n; ++k) {
    tbl[k] = tbl[k - 1] + std::log2(static_cast<double>(2 * k - 3));
  }
  return tbl;
}

// Log2Unrooted(n) = log2((2n-5)!!) for n >= 3.
// = Log2Rooted(n) - log2(2n - 3).
static double compute_log2_unrooted(int n, const std::vector<double>& log2_rooted) {
  if (n < 3) return 0.0;
  return log2_rooted[n] - std::log2(static_cast<double>(2 * n - 3));
}


// ==========================================================================
// build_spic_data
// ==========================================================================

SpicData build_spic_data_from_splits(const std::vector<CidSplitSet>& tree_splits,
                                      int n_tips, int n_trees,
                                      int n_bins_override) {
  SpicData sd;
  sd.n_tips = n_tips;
  sd.n_bins = (n_bins_override > 0) ? n_bins_override : (n_tips + 63) / 64;
  sd.n_trees = n_trees;

  // Build log2 tables
  sd.log2_rooted = build_log2_rooted(n_tips);
  sd.log2_unrooted_n = compute_log2_unrooted(n_tips, sd.log2_rooted);

  // Precompute P_consistent and P_inconsistent for each partition size a
  sd.log2_p_consistent.resize(n_tips + 1, 0.0);
  sd.log2_p_inconsistent.resize(n_tips + 1, 0.0);
  for (int a = 1; a <= n_tips; ++a) {
    int b = n_tips - a;
    if (a < 2 || b < 2) {
      // Trivial split (1|n-1 or 0|n): present in all trees
      sd.log2_p_consistent[a] = 0.0;  // log2(1) = 0
      sd.log2_p_inconsistent[a] = -std::numeric_limits<double>::infinity();
    } else {
      double l2_n_consistent = sd.log2_rooted[a] + sd.log2_rooted[b];
      sd.log2_p_consistent[a] = l2_n_consistent - sd.log2_unrooted_n;
      double p_consistent = std::exp2(sd.log2_p_consistent[a]);
      if (p_consistent >= 1.0) {
        sd.log2_p_inconsistent[a] = -std::numeric_limits<double>::infinity();
      } else {
        sd.log2_p_inconsistent[a] = std::log2(1.0 - p_consistent);
      }
    }
  }

  // Build split frequency table from all input tree splits
  int n_bins = sd.n_bins;
  for (int t = 0; t < n_trees; ++t) {
    const CidSplitSet& ss = tree_splits[t];
    for (int s = 0; s < ss.n_splits; ++s) {
      const uint64_t* sp = ss.split(s);
      uint64_t h = hash_split_key(sp, n_bins);

      // Check for existing entry with same hash and same data
      bool found = false;
      auto range = sd.freq_table.equal_range(h);
      for (auto it = range.first; it != range.second; ++it) {
        if (it->second.in_split == ss.in_split[s] &&
            std::memcmp(sp, it->second.data.data(),
                        sizeof(uint64_t) * n_bins) == 0) {
          it->second.count++;
          found = true;
          break;
        }
      }
      if (!found) {
        SpicData::SplitEntry entry;
        entry.count = 1;
        entry.in_split = ss.in_split[s];
        entry.data.assign(sp, sp + n_bins);
        sd.freq_table.emplace(h, std::move(entry));
      }
    }
  }

  // Allocate candidate buffers
  sd.cand_tip_bits.resize(static_cast<size_t>(n_tips) * n_bins, 0);
  sd.cand_buf.n_bins = n_bins;

  return sd;
}


SpicData build_spic_data(const CidData& cd) {
  return build_spic_data_from_splits(cd.tree_splits, cd.n_tips, cd.n_trees);
}


// ==========================================================================
// split_ic: information content of one split
// ==========================================================================
//
// Per-split IC formula (from SplitwiseInfo with probability weighting):
//
//   IC = log2(N_unrooted) + p*(log2(p) - log2(N_consistent))
//        + (1-p)*(log2(1-p) - log2(N_inconsistent))
//
// where:
//   p = split frequency in input trees
//   N_consistent = TreesMatchingSplit(a, b) = Rooted(a) * Rooted(b)
//   N_inconsistent = N_unrooted - N_consistent
//   log2(N_consistent) = Log2Rooted(a) + Log2Rooted(b)
//   log2(N_inconsistent) = log2(P_inconsistent) + log2(N_unrooted)
//
// Rewritten using precomputed tables:
//   IC = l2n + p*(log2(p) - l2n - l2pC) + q*(log2(q) - l2n - l2pI)
//      = l2n * (1 - p - q) + p*(log2(p) - l2pC) + q*(log2(q) - l2pI)
//   But 1-p-q = 0, so:
//   IC = p*(log2(p) - log2(N_consistent)) + q*(log2(q) - log2(N_inconsistent))
//   Wait, the formula uses l2n + p*(log2(p) - l2nConsistent):
//   l2nConsistent = log2(N_consistent) = Log2Rooted(a) + Log2Rooted(b)
//   l2nInconsistent = log2(N_inconsistent) = log2(P_inconsistent) + l2n
//   IC = l2n + p*(log2(p) - l2nConsistent) + q*(log2(q) - l2nInconsistent)

double split_ic(int a, int n, double freq, const SpicData& sd) {
  int b = n - a;
  // Canonical: ensure a <= b for partition size
  if (a > b) { int tmp = a; a = b; b = tmp; }

  // Trivial splits (1|n-1 or smaller): always present, IC = 0
  if (a < 2) return 0.0;

  double p = freq;
  double q = 1.0 - p;

  double l2n = sd.log2_unrooted_n;
  double l2n_consistent = sd.log2_rooted[a] + sd.log2_rooted[b];
  double l2p_inconsistent = sd.log2_p_inconsistent[a];
  double l2n_inconsistent = l2p_inconsistent + l2n;

  double ic = l2n;

  // p * (log2(p) - log2(N_consistent))
  if (p > 0.0) {
    ic += p * (std::log2(p) - l2n_consistent);
  }
  // q * (log2(q) - log2(N_inconsistent))
  if (q > 0.0 && std::isfinite(l2n_inconsistent)) {
    ic += q * (std::log2(q) - l2n_inconsistent);
  }

  return ic;
}


// ==========================================================================
// spic_score: score candidate tree using SPIC
// ==========================================================================

double spic_score(TreeState& tree, const SpicData& sd) {
  // Extract candidate tree splits
  compute_splits_cid(tree, sd.cand_tip_bits, sd.cand_buf);
  const CidSplitSet& cand = sd.cand_buf;

  double total_ic = 0.0;
  int n_bins = sd.n_bins;

  for (int s = 0; s < cand.n_splits; ++s) {
    const uint64_t* sp = cand.split(s);
    int a = cand.in_split[s];
    int b = sd.n_tips - a;

    // Skip trivial splits
    if (a < 2 || b < 2) continue;

    // Look up frequency in input tree set
    uint64_t h = hash_split_key(sp, n_bins);
    int count = 0;

    auto range = sd.freq_table.equal_range(h);
    for (auto it = range.first; it != range.second; ++it) {
      if (it->second.in_split == a &&
          std::memcmp(sp, it->second.data.data(),
                      sizeof(uint64_t) * n_bins) == 0) {
        count = it->second.count;
        break;
      }
    }

    double freq = static_cast<double>(count) / sd.n_trees;
    total_ic += split_ic(a, sd.n_tips, freq, sd);
  }

  // Return negated (lower = better consensus), consistent with cid_score()
  return -total_ic;
}


} // namespace ts
