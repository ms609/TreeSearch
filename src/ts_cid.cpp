// CID (Clustering Information Distance) scoring for the driven search.
//
// Contains:
//   - LAP solver (Jonker-Volgenant, adapted from TreeDist)
//   - Mutual clustering information computation
//   - Clustering entropy computation
//   - cid_score(): scores a candidate tree against input trees
//   - build_mrp_dataset(): constructs MRP characters for Fitch screening
//   - prepare_cid_data(): one-time initialization of hash indices and logs
//
// When TreeDist exposes LinkingTo headers, the LAP solver and MCI
// computation should be replaced with the TreeDist implementations.

#include "ts_cid.h"
#include "ts_data.h"
#include <cassert>
#include <cstring>
#include <limits>
#include <numeric>

namespace ts {

// ==========================================================================
// LAP solver (Jonker-Volgenant shortest augmenting path algorithm)
// Adapted from TreeDist's lap.cpp (Smith 2020, after Volgenant 1996)
// ==========================================================================

lap_cost lap_solve(int dim, LapMatrix& cost,
                   int* rowsol, int* colsol, LapScratch& scratch) {
  if (dim <= 0) return 0;
  if (dim == 1) {
    rowsol[0] = 0;
    colsol[0] = 0;
    return cost(0, 0);
  }

  scratch.ensure(dim);
  lap_cost* v = scratch.v.data();
  int* matches = scratch.matches.data();
  int* free_rows = scratch.free_rows.data();
  int* col_list = scratch.col_list.data();
  lap_cost* d = scratch.d.data();
  int* pred = scratch.pred.data();

  std::fill(matches, matches + dim, 0);
  std::fill(rowsol, rowsol + dim, -1);
  std::fill(colsol, colsol + dim, -1);

  for (int j = dim; j--; ) {
    lap_cost min_val = cost(0, j);
    int min_row = 0;
    for (int i = 1; i < dim; ++i) {
      if (cost(i, j) < min_val) { min_val = cost(i, j); min_row = i; }
    }
    v[j] = min_val;
    matches[min_row]++;
    if (matches[min_row] == 1) { rowsol[min_row] = j; colsol[j] = min_row; }
    else { colsol[j] = -1; }
  }

  int n_free = 0;
  for (int i = 0; i < dim; ++i) {
    if (matches[i] == 0) {
      free_rows[n_free++] = i;
    } else if (matches[i] == 1) {
      int j1 = rowsol[i];
      lap_cost min2 = LAP_BIG;
      for (int j = 0; j < dim; ++j) {
        if (j != j1) {
          lap_cost reduced = cost(i, j) - v[j];
          if (reduced < min2) min2 = reduced;
        }
      }
      v[j1] -= min2;
    }
  }

  for (int pass = 0; pass < 2; ++pass) {
    int k = 0; int prev_free = n_free; n_free = 0;
    while (k < prev_free) {
      int i = free_rows[k++];
      lap_cost u1 = cost(i, 0) - v[0]; int j1 = 0; lap_cost u2 = LAP_BIG;
      for (int j = 1; j < dim; ++j) {
        lap_cost reduced = cost(i, j) - v[j];
        if (reduced < u2) {
          if (reduced >= u1) { u2 = reduced; }
          else { u2 = u1; u1 = reduced; j1 = j; }
        }
      }
      int i0 = colsol[j1];
      if (u1 < u2) { v[j1] += u1 - u2; }
      else if (i0 >= 0) {
        int j2 = -1; lap_cost u2_find = LAP_BIG;
        for (int j = 0; j < dim; ++j) {
          if (j != j1) { lap_cost r2 = cost(i, j) - v[j]; if (r2 < u2_find) { u2_find = r2; j2 = j; } }
        }
        if (j2 >= 0) { j1 = j2; i0 = colsol[j1]; }
      }
      rowsol[i] = j1; colsol[j1] = i;
      if (i0 >= 0) { if (u1 < u2) free_rows[--k] = i0; else free_rows[n_free++] = i0; }
    }
  }

  for (int f = 0; f < n_free; ++f) {
    int i1 = free_rows[f]; int j1 = -1;
    for (int j = 0; j < dim; ++j) { d[j] = cost(i1, j) - v[j]; pred[j] = i1; col_list[j] = j; }
    int low = 0, up = 0;
    while (true) {
      if (up == low) {
        lap_cost last_min = LAP_BIG;
        for (int k = up; k < dim; ++k) {
          int j = col_list[k];
          if (d[j] <= last_min) { if (d[j] < last_min) { last_min = d[j]; up = low; } col_list[k] = col_list[up]; col_list[up] = j; ++up; }
        }
        for (int k = low; k < up; ++k) { if (colsol[col_list[k]] < 0) { j1 = col_list[k]; goto augment; } }
      }
      int j_scan = col_list[low++]; int i_scan = colsol[j_scan];
      lap_cost scan_cost = cost(i_scan, j_scan) - v[j_scan] - d[j_scan];
      for (int k = up; k < dim; ++k) {
        int j = col_list[k]; lap_cost new_d = cost(i_scan, j) - v[j] - scan_cost;
        if (new_d < d[j]) { d[j] = new_d; pred[j] = i_scan;
          if (new_d == d[col_list[low]]) { if (colsol[j] < 0) { j1 = j; goto augment; } col_list[k] = col_list[up]; col_list[up++] = j; }
        }
      }
    }
    augment:
    for (int k = 0; k < low; ++k) { int j = col_list[k]; v[j] += d[j] - d[j1]; }
    while (true) { int i_aug = pred[j1]; colsol[j1] = i_aug; int j_prev = rowsol[i_aug]; rowsol[i_aug] = j1; j1 = j_prev; if (i_aug == i1) break; }
  }

  lap_cost total = 0;
  for (int i = 0; i < dim; ++i) total += cost(i, rowsol[i]);
  return total;
}


// ==========================================================================
// Clustering entropy
// ==========================================================================

double clustering_entropy(const CidSplitSet& ss, int n_tips) {
  if (n_tips <= 1 || ss.n_splits == 0) return 0.0;
  double lg2n = std::log2(static_cast<double>(n_tips));
  double ce = 0.0;
  for (int i = 0; i < ss.n_splits; ++i) {
    int a = ss.in_split[i];
    int b = n_tips - a;
    if (a > 1 && b > 1) {
      ce += (a * std::log2(static_cast<double>(a)) +
             b * std::log2(static_cast<double>(b)))
            / n_tips - lg2n;
    }
  }
  return -ce;
}

double clustering_entropy_fast(const CidSplitSet& ss, int n_tips,
                               double lg2_n) {
  if (n_tips <= 1 || ss.n_splits == 0) return 0.0;
  double ce = 0.0;
  for (int i = 0; i < ss.n_splits; ++i) {
    int a = ss.in_split[i];
    int b = n_tips - a;
    if (a > 1 && b > 1) {
      ce += (a * ss.lg2_in[i] + b * ss.lg2_out[i]) / n_tips - lg2_n;
    }
  }
  return -ce;
}


// ==========================================================================
// Split hashing for O(1) exact-match lookup
// ==========================================================================

static inline uint64_t hash_split_key(const uint64_t* data, int n_bins) {
  if (n_bins == 1) return data[0];
  uint64_t h = 14695981039346656037ULL;
  for (int i = 0; i < n_bins; ++i) {
    h ^= data[i];
    h *= 1099511628211ULL;
  }
  return h;
}


// ==========================================================================
// Mutual clustering information (adapted from TreeDist)
// ==========================================================================

double mutual_clustering_info(const CidSplitSet& a, const CidSplitSet& b,
                              int n_tips, LapScratch& scratch) {
  int a_n = a.n_splits;
  int b_n = b.n_splits;
  int n_bins = a.n_bins;

  if (a_n == 0 || b_n == 0 || n_tips == 0) return 0.0;

  double n_recip = 1.0 / n_tips;
  double lg2_n = std::log2(static_cast<double>(n_tips));
  bool a_has_logs = (static_cast<int>(a.lg2_in.size()) == a_n);
  bool b_has_logs = (static_cast<int>(b.lg2_in.size()) == b_n);

  // Phase 1: find exact matches
  double exact_score = 0.0;
  int n_exact = 0;
  std::vector<bool> a_matched(a_n, false);
  std::vector<bool> b_matched(b_n, false);

  if (!b.hash_index.empty()) {
    // Hash-based O(n) exact match
    for (int ai = 0; ai < a_n; ++ai) {
      uint64_t h = hash_split_key(a.split(ai), n_bins);
      auto it = b.hash_index.find(h);
      if (it != b.hash_index.end()) {
        int bi = it->second;
        if (!b_matched[bi] &&
            std::memcmp(a.split(ai), b.split(bi),
                        sizeof(uint64_t) * n_bins) == 0) {
          int na = a.in_split[ai];
          int nb_val = n_tips - na;
          if (a_has_logs) {
            exact_score += n_tips * lg2_n - na * a.lg2_in[ai]
                           - nb_val * a.lg2_out[ai];
          } else {
            exact_score += n_tips * lg2_n
                - na * std::log2(static_cast<double>(na))
                - nb_val * std::log2(static_cast<double>(nb_val));
          }
          a_matched[ai] = true;
          b_matched[bi] = true;
          ++n_exact;
        }
      }
    }
  } else {
    // Fallback: O(n^2) scan with complement check
    int unset = (n_tips % 64) ? 64 - (n_tips % 64) : 0;
    uint64_t last_mask = unset ? (~uint64_t(0)) >> unset : ~uint64_t(0);
    for (int ai = 0; ai < a_n; ++ai) {
      const uint64_t* a_row = a.split(ai);
      for (int bi = 0; bi < b_n; ++bi) {
        if (b_matched[bi]) continue;
        const uint64_t* b_row = b.split(bi);
        bool match = true;
        for (int bin = 0; bin < n_bins; ++bin) {
          if (a_row[bin] != b_row[bin]) { match = false; break; }
        }
        if (!match) {
          match = true;
          for (int bin = 0; bin < n_bins - 1; ++bin) {
            if (a_row[bin] != ~b_row[bin]) { match = false; break; }
          }
          if (match && n_bins > 0) {
            match = (a_row[n_bins - 1] == (b_row[n_bins - 1] ^ last_mask));
          }
        }
        if (match) {
          int na = a.in_split[ai];
          int nb_val = n_tips - na;
          if (a_has_logs) {
            exact_score += n_tips * lg2_n - na * a.lg2_in[ai]
                           - nb_val * a.lg2_out[ai];
          } else {
            exact_score += n_tips * lg2_n
                - na * std::log2(static_cast<double>(na))
                - nb_val * std::log2(static_cast<double>(nb_val));
          }
          a_matched[ai] = true;
          b_matched[bi] = true;
          ++n_exact;
          break;
        }
      }
    }
  }

  if (n_exact == a_n || n_exact == b_n) {
    return exact_score * n_recip;
  }

  // Phase 2: LAP for remaining unmatched splits
  int a_unmatched = a_n - n_exact;
  int b_unmatched = b_n - n_exact;
  int lap_dim = std::max(a_unmatched, b_unmatched);

  scratch.ensure(lap_dim);
  scratch.score_pool.resize(lap_dim);
  LapMatrix& score = scratch.score_pool;

  constexpr lap_cost max_score = LAP_BIG / 4096;
  double max_over_tips = static_cast<double>(max_score) * n_recip;

  int a_pos = 0;
  for (int ai = 0; ai < a_n; ++ai) {
    if (a_matched[ai]) continue;
    const uint64_t* a_row = a.split(ai);
    int na = a.in_split[ai];
    int nA = n_tips - na;
    double offset_a = a_has_logs
        ? (lg2_n - a.lg2_in[ai])
        : (lg2_n - std::log2(static_cast<double>(na)));
    double offset_A = a_has_logs
        ? (lg2_n - a.lg2_out[ai])
        : (lg2_n - std::log2(static_cast<double>(nA)));

    int b_pos = 0;
    for (int bi = 0; bi < b_n; ++bi) {
      if (b_matched[bi]) continue;
      const uint64_t* b_row = b.split(bi);
      int nb = b.in_split[bi];
      int nB = n_tips - nb;

      int a_and_b = 0;
      for (int bin = 0; bin < n_bins; ++bin) {
        a_and_b += popcount64(a_row[bin] & b_row[bin]);
      }
      int a_and_B = na - a_and_b;
      int A_and_b = nb - a_and_b;
      int A_and_B = nA - A_and_b;

      if (a_and_b == A_and_b && a_and_b == a_and_B && a_and_b == A_and_B) {
        score(a_pos, b_pos) = max_score;
      } else {
        double lg2_nb = b_has_logs ? b.lg2_in[bi]
            : ((nb > 0) ? std::log2(static_cast<double>(nb)) : 0);
        double lg2_nB = b_has_logs ? b.lg2_out[bi]
            : ((nB > 0) ? std::log2(static_cast<double>(nB)) : 0);
        double ic_sum = 0.0;
        if (a_and_b > 0)
          ic_sum += a_and_b * (std::log2(static_cast<double>(a_and_b))
                    + offset_a - lg2_nb);
        if (a_and_B > 0)
          ic_sum += a_and_B * (std::log2(static_cast<double>(a_and_B))
                    + offset_a - lg2_nB);
        if (A_and_b > 0)
          ic_sum += A_and_b * (std::log2(static_cast<double>(A_and_b))
                    + offset_A - lg2_nb);
        if (A_and_B > 0)
          ic_sum += A_and_B * (std::log2(static_cast<double>(A_and_B))
                    + offset_A - lg2_nB);
        score(a_pos, b_pos) = max_score -
            static_cast<lap_cost>(ic_sum * max_over_tips);
      }
      ++b_pos;
    }
    for (int j = b_unmatched; j < lap_dim; ++j) score(a_pos, j) = max_score;
    ++a_pos;
  }
  for (int i = a_unmatched; i < lap_dim; ++i) {
    for (int j = 0; j < lap_dim; ++j) score(i, j) = max_score;
  }

  scratch.ensure(lap_dim);
  lap_cost lap_total = lap_solve(lap_dim, score,
                                  scratch.rowsol.data(),
                                  scratch.colsol.data(), scratch);

  double lap_score = static_cast<double>(
      static_cast<lap_cost>(max_score) * lap_dim - lap_total)
      / static_cast<double>(max_score);

  return exact_score * n_recip + lap_score;
}


// ==========================================================================
// Convert ts::SplitSet to CidSplitSet
// ==========================================================================

CidSplitSet splitset_to_cid(const SplitSet& ss, int n_tips) {
  CidSplitSet cs;
  cs.n_splits = ss.n_splits;
  cs.n_bins = ss.words_per_split;
  cs.data.assign(ss.splits.begin(), ss.splits.end());
  cs.in_split.resize(ss.n_splits);
  for (int i = 0; i < ss.n_splits; ++i) {
    int cnt = 0;
    const uint64_t* sp = ss.split(i);
    for (int w = 0; w < ss.words_per_split; ++w) {
      cnt += popcount64(sp[w]);
    }
    cs.in_split[i] = cnt;
  }
  return cs;
}


// ==========================================================================
// Compute splits directly into a preallocated CidSplitSet
// ==========================================================================

void compute_splits_cid(const TreeState& tree,
                        std::vector<uint64_t>& tip_bits_work,
                        CidSplitSet& out) {
  int n_tip = tree.n_tip;
  int wps = (n_tip + 63) / 64;

  size_t total = static_cast<size_t>(tree.n_node) * wps;
  if (tip_bits_work.size() < total) tip_bits_work.resize(total);
  std::fill(tip_bits_work.begin(), tip_bits_work.begin() + total, 0);

  for (int t = 0; t < n_tip; ++t) {
    tip_bits_work[static_cast<size_t>(t) * wps + t / 64] = 1ULL << (t % 64);
  }

  for (int pi = 0; pi < static_cast<int>(tree.postorder.size()); ++pi) {
    int node = tree.postorder[pi];
    int ni = node - n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    uint64_t* dst = &tip_bits_work[static_cast<size_t>(node) * wps];
    const uint64_t* lbits = &tip_bits_work[static_cast<size_t>(lc) * wps];
    const uint64_t* rbits = &tip_bits_work[static_cast<size_t>(rc) * wps];
    for (int w = 0; w < wps; ++w) dst[w] = lbits[w] | rbits[w];
  }

  int root = n_tip;
  int root_right = tree.right[0];
  int trailing = n_tip % 64;
  uint64_t trail_mask = (trailing != 0) ? ((1ULL << trailing) - 1) : ~0ULL;

  int n_splits = 0;
  for (int pi = 0; pi < static_cast<int>(tree.postorder.size()); ++pi) {
    int node = tree.postorder[pi];
    if (node == root || node == root_right) continue;
    const uint64_t* bits = &tip_bits_work[static_cast<size_t>(node) * wps];
    int count = 0;
    for (int w = 0; w < wps; ++w) count += popcount64(bits[w]);
    if (count <= 1 || count >= n_tip - 1) continue;
    ++n_splits;
  }

  out.n_splits = n_splits;
  out.n_bins = wps;
  size_t data_size = static_cast<size_t>(n_splits) * wps;
  if (out.data.size() < data_size) out.data.resize(data_size);
  if (static_cast<int>(out.in_split.size()) < n_splits)
    out.in_split.resize(n_splits);
  if (static_cast<int>(out.lg2_in.size()) < n_splits)
    out.lg2_in.resize(n_splits);
  if (static_cast<int>(out.lg2_out.size()) < n_splits)
    out.lg2_out.resize(n_splits);

  int idx = 0;
  for (int pi = 0; pi < static_cast<int>(tree.postorder.size()); ++pi) {
    int node = tree.postorder[pi];
    if (node == root || node == root_right) continue;
    const uint64_t* bits = &tip_bits_work[static_cast<size_t>(node) * wps];
    int count = 0;
    for (int w = 0; w < wps; ++w) count += popcount64(bits[w]);
    if (count <= 1 || count >= n_tip - 1) continue;

    uint64_t* dst = &out.data[static_cast<size_t>(idx) * wps];
    std::memcpy(dst, bits, sizeof(uint64_t) * wps);
    if (dst[0] & 1ULL) {
      for (int w = 0; w < wps; ++w) dst[w] = ~dst[w];
      if (trailing != 0) dst[wps - 1] &= trail_mask;
      count = n_tip - count;
    }
    out.in_split[idx] = count;
    out.lg2_in[idx] = (count > 0)
        ? std::log2(static_cast<double>(count)) : 0.0;
    out.lg2_out[idx] = (n_tip - count > 0)
        ? std::log2(static_cast<double>(n_tip - count)) : 0.0;
    ++idx;
  }
}


// ==========================================================================
// prepare_cid_data: one-time initialization
// ==========================================================================

void prepare_cid_data(CidData& cd) {
  cd.lg2_n = std::log2(static_cast<double>(cd.n_tips));
  cd.max_splits = 0;

  for (int t = 0; t < cd.n_trees; ++t) {
    CidSplitSet& ss = cd.tree_splits[t];
    ss.lg2_in.resize(ss.n_splits);
    ss.lg2_out.resize(ss.n_splits);
    for (int i = 0; i < ss.n_splits; ++i) {
      int a = ss.in_split[i];
      int b = cd.n_tips - a;
      ss.lg2_in[i] = (a > 0) ? std::log2(static_cast<double>(a)) : 0.0;
      ss.lg2_out[i] = (b > 0) ? std::log2(static_cast<double>(b)) : 0.0;
    }
    ss.hash_index.clear();
    ss.hash_index.reserve(ss.n_splits);
    for (int i = 0; i < ss.n_splits; ++i) {
      uint64_t h = hash_split_key(ss.split(i), ss.n_bins);
      ss.hash_index.emplace(h, i);
    }
    if (ss.n_splits > cd.max_splits) cd.max_splits = ss.n_splits;
  }

  if (cd.max_splits > 0) {
    cd.lap_scratch.ensure(cd.max_splits);
    cd.lap_scratch.score_pool.resize(cd.max_splits);
  }

  int n_node = 2 * cd.n_tips - 1;
  cd.cand_tip_bits.resize(static_cast<size_t>(n_node) * cd.n_bins, 0);
}


// ==========================================================================
// cid_score: score candidate tree against all input trees.
// Returns negated weighted mean MCI (lower = better consensus).
// ==========================================================================

double cid_score(TreeState& tree, const CidData& cd) {
  compute_splits_cid(tree, cd.cand_tip_bits, cd.cand_buf);
  CidSplitSet& cand = cd.cand_buf;
  double budget = cd.score_budget;

  // Upper bound on per-tree MCI (used for early termination)
  double cand_ce = (budget < HUGE_VAL)
      ? clustering_entropy_fast(cand, cd.n_tips, cd.lg2_n)
      : 0.0;

  double mci_sum = 0.0;
  double weight_done = 0.0;
  for (int i = 0; i < cd.n_trees; ++i) {
    if (cd.tree_weights[i] <= 0.0) continue;
    double mci = mutual_clustering_info(cand, cd.tree_splits[i],
                                         cd.n_tips, cd.lap_scratch);
    mci_sum += cd.tree_weights[i] * mci;
    weight_done += cd.tree_weights[i];
    // Early termination: even with perfect MCI for remaining trees,
    // the negated mean can't beat budget
    if (budget < HUGE_VAL) {
      double remaining = cd.weight_sum - weight_done;
      double best_possible = -(mci_sum + remaining * cand_ce) / cd.weight_sum;
      if (best_possible > budget) return best_possible;
    }
  }
  return -mci_sum / cd.weight_sum;
}


// ==========================================================================
// build_mrp_dataset: construct MRP binary characters for Fitch screening
// ==========================================================================

DataSet build_mrp_dataset(CidData& cd) {
  int n_tips = cd.n_tips;
  int n_bins = cd.n_bins;

  int total_chars = 0;
  for (int t = 0; t < cd.n_trees; ++t) {
    total_chars += cd.tree_splits[t].n_splits;
  }

  int n_blocks = (total_chars + MAX_CHARS_PER_BLOCK - 1) / MAX_CHARS_PER_BLOCK;
  if (n_blocks == 0) n_blocks = 1;

  DataSet ds;
  ds.n_tips = n_tips;
  ds.scoring_mode = ScoringMode::CID;
  ds.cid_data = &cd;
  ds.concavity = cd.mrp_concavity;

  static constexpr int N_STATES = 2;

  ds.blocks.clear();
  cd.mrp_tree_block_start.clear();
  cd.mrp_tree_block_start.reserve(cd.n_trees + 1);

  int char_idx = 0;
  int block_idx = 0;

  for (int t = 0; t < cd.n_trees; ++t) {
    cd.mrp_tree_block_start.push_back(block_idx);
    int n_splits = cd.tree_splits[t].n_splits;
    for (int s = 0; s < n_splits; ++s) {
      if (char_idx % MAX_CHARS_PER_BLOCK == 0) {
        CharBlock blk;
        blk.n_chars = 0;
        blk.n_states = N_STATES;
        blk.weight = 1;
        blk.has_inapplicable = false;
        blk.active_mask = 0;
        blk.upweight_mask = 0;
        std::memset(blk.pattern_index, 0, sizeof(blk.pattern_index));
        ds.blocks.push_back(blk);
        ++block_idx;
      }
      int local_idx = char_idx % MAX_CHARS_PER_BLOCK;
      ds.blocks.back().n_chars = local_idx + 1;
      ds.blocks.back().active_mask |= (uint64_t(1) << local_idx);
      ds.blocks.back().pattern_index[local_idx] = char_idx;
      ++char_idx;
    }
  }
  cd.mrp_tree_block_start.push_back(block_idx);

  ds.n_blocks = static_cast<int>(ds.blocks.size());
  ds.total_words = ds.n_blocks * N_STATES;
  ds.block_word_offset.resize(ds.n_blocks);
  for (int b = 0; b < ds.n_blocks; ++b) {
    ds.block_word_offset[b] = b * N_STATES;
  }

  size_t tip_state_size = static_cast<size_t>(n_tips) * ds.total_words;
  ds.tip_states.assign(tip_state_size, 0);

  char_idx = 0;
  int blk_i = 0;
  for (int t = 0; t < cd.n_trees; ++t) {
    const CidSplitSet& ss = cd.tree_splits[t];
    for (int s = 0; s < ss.n_splits; ++s) {
      int local_idx = char_idx % MAX_CHARS_PER_BLOCK;
      uint64_t bit = uint64_t(1) << local_idx;
      if (local_idx == 0 && char_idx > 0) ++blk_i;
      int word_off = ds.block_word_offset[blk_i];
      for (int tip = 0; tip < n_tips; ++tip) {
        size_t base = static_cast<size_t>(tip) * ds.total_words + word_off;
        int word_in_split = tip / 64;
        int bit_in_word = tip % 64;
        bool in_split = (word_in_split < n_bins) &&
            ((ss.data[static_cast<size_t>(s) * n_bins + word_in_split]
              >> bit_in_word) & 1);
        if (in_split) {
          ds.tip_states[base + 1] |= bit;
        } else {
          ds.tip_states[base + 0] |= bit;
        }
      }
      ++char_idx;
    }
  }

  ds.n_patterns = total_chars;
  ds.ew_offset = 0;
  ds.precomputed_steps.assign(total_chars, 0);
  ds.min_steps.assign(total_chars, 1);
  ds.pattern_freq.assign(total_chars, 1);

  // Populate IW weight arrays required by compute_iw() / compute_weighted_score().
  // Binary MRP characters have min_steps = 1, so eff_k and phi use the
  // standard IW formula: eff_k = concavity, phi = 1.0 (no XPIWE correction).
  // If concavity is not finite (EW mode), fill with 0 / 1 as a safe default.
  double k = std::isfinite(ds.concavity) ? ds.concavity : 0.0;
  ds.eff_k.assign(total_chars, k);
  ds.phi.assign(total_chars, 1.0);

  return ds;
}


// ==========================================================================
// CID weight management for ratchet
// ==========================================================================

void save_cid_weights(CidData& cd) {
  cd.saved_tree_weights = cd.tree_weights;
  cd.saved_weight_sum = cd.weight_sum;
}

void restore_cid_weights(CidData& cd) {
  cd.tree_weights = cd.saved_tree_weights;
  cd.weight_sum = cd.saved_weight_sum;
}

void sync_cid_weights_from_mrp(CidData& cd, const DataSet& ds) {
  cd.weight_sum = 0.0;
  for (int t = 0; t < cd.n_trees; ++t) {
    int b_start = cd.mrp_tree_block_start[t];
    int b_end = cd.mrp_tree_block_start[t + 1];
    int total_chars = 0;
    int active_chars = 0;
    for (int b = b_start; b < b_end; ++b) {
      total_chars += ds.blocks[b].n_chars;
      active_chars += popcount64(ds.blocks[b].active_mask);
      active_chars += popcount64(ds.blocks[b].upweight_mask);
    }
    cd.tree_weights[t] = (total_chars > 0)
        ? static_cast<double>(active_chars) / total_chars
        : 0.0;
    cd.weight_sum += cd.tree_weights[t];
  }
}

} // namespace ts
