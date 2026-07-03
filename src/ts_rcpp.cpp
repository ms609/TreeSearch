#include <Rcpp.h>
#include <chrono>
#include <cstdio>
#include <random>
#include "ts_data.h"
#include "ts_tree.h"
#include "ts_fitch.h"
#include "ts_search.h"
#include "ts_tbr.h"
#include "ts_drift.h"
#include "ts_ratchet.h"
#include "ts_splits.h"
#include "ts_pool.h"
#include "ts_wagner.h"
#include "ts_sector.h"
#include "ts_fuse.h"
#include "ts_driven.h"
#include "ts_constraint.h"
#include "ts_resample.h"
#include "ts_rng.h"
#include "ts_parallel.h"
#include "ts_simplify.h"
#include "ts_hsj.h"
#include "ts_ls.h"
#include "ts_collapsed.h"
// ts_temper.h removed — parallel tempering lives on feature/parallel-temper
#include "ts_strategy.h"

using namespace Rcpp;

namespace {

// Sentinel: concavity = -1 means equal weights (Inf).
// Rcpp can't auto-generate R_PosInf as an R default, so we use -1
// and convert here at the single gateway into the C++ engine.
ts::DataSet make_dataset(
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0,
    Nullable<NumericMatrix> infoAmounts = R_NilValue,
    bool xpiwe = false,
    double xpiwe_r = 0.5,
    double xpiwe_max_f = 5.0,
    IntegerVector obs_count = IntegerVector())
{
  if (concavity < 0) concavity = HUGE_VAL;
  int n_tips = tip_data.nrow();
  int n_patterns = tip_data.ncol();
  int n_tokens = contrast.nrow();
  int n_states = contrast.ncol();

  // Validate parallel-vector lengths at the boundary: build_dataset() indexes
  // these by raw pointer over [0, n_patterns)/[0, n_states), so a short vector
  // is an out-of-bounds read. Public wrappers always size them correctly; this
  // turns a direct internal call with mismatched lengths into a clean error
  // rather than undefined behaviour.
  if (weight.size() != n_patterns) {
    Rcpp::stop("`weight` length (%d) must equal the number of characters (%d)",
               static_cast<int>(weight.size()), n_patterns);
  }
  if (levels.size() != n_states) {
    Rcpp::stop("`levels` length (%d) must equal ncol(contrast) (%d)",
               static_cast<int>(levels.size()), n_states);
  }
  if (min_steps.size() > 0 && min_steps.size() != n_patterns) {
    Rcpp::stop("`min_steps` length (%d) must equal the number of characters "
               "(%d)", static_cast<int>(min_steps.size()), n_patterns);
  }
  if (xpiwe && obs_count.size() > 0 && obs_count.size() != n_patterns) {
    Rcpp::stop("`obs_count` length (%d) must equal the number of characters "
               "(%d)", static_cast<int>(obs_count.size()), n_patterns);
  }

  std::vector<std::string> level_strs(n_states);
  std::vector<const char*> level_ptrs(n_states);
  for (int i = 0; i < n_states; ++i) {
    level_strs[i] = as<std::string>(levels[i]);
    level_ptrs[i] = level_strs[i].c_str();
  }

  const int* min_steps_ptr = (min_steps.size() > 0) ? INTEGER(min_steps)
                                                     : nullptr;

  // Profile parsimony: extract info_amounts table if provided
  const double* info_amounts_ptr = nullptr;
  int info_max_steps = 0;
  if (infoAmounts.isNotNull()) {
    NumericMatrix ia(infoAmounts.get());
    info_amounts_ptr = REAL(ia);
    info_max_steps = ia.nrow();
  }

  // XPIWE: per-pattern observed-taxa count
  const int* obs_count_ptr = (xpiwe && obs_count.size() > 0)
                                 ? INTEGER(obs_count)
                                 : nullptr;

  return ts::build_dataset(
      REAL(contrast), n_tokens, n_states,
      INTEGER(tip_data), n_tips, n_patterns,
      INTEGER(weight),
      level_ptrs.data(),
      min_steps_ptr,
      concavity,
      info_amounts_ptr,
      info_max_steps,
      xpiwe,
      xpiwe_r,
      xpiwe_max_f,
      obs_count_ptr);
}

// Convert TreeState topology back to R edge matrix (2-column, 1-based)
IntegerMatrix tree_to_edge(const ts::TreeState& tree) {
  int n_edge = 2 * (tree.n_tip - 1);
  IntegerMatrix edge(n_edge, 2);
  int row = 0;
  for (int node = tree.n_tip; node < tree.n_node; ++node) {
    int ni = node - tree.n_tip;
    // +1 to convert to R's 1-based indexing
    edge(row, 0) = node + 1;
    edge(row, 1) = tree.left[ni] + 1;
    ++row;
    edge(row, 0) = node + 1;
    edge(row, 1) = tree.right[ni] + 1;
    ++row;
  }
  return edge;
}

// Emit a (possibly multifurcating) R edge matrix (2-col, 1-based) from a binary
// TreeState, contracting every internal edge flagged in `collapsed`.
//
// collapsed[c] == 1 (c internal) means the edge c -> parent[c] is zero-length:
// c is removed and its children reattach to c's nearest retained ancestor,
// producing a polytomy.  Tips and the root pseudo-node are always retained;
// compute_collapsed_flags[_aggressive] never flags the root or its children, so
// when the tree is rooted on a tip (the precondition for a rooting-invariant
// collapse) no informative edge is root-adjacent and the basal split survives.
IntegerMatrix tree_to_collapsed_edge(const ts::TreeState& tree,
                                     const std::vector<uint8_t>& collapsed) {
  const int n_tip = tree.n_tip;
  const int n_node = tree.n_node;
  const int root = n_tip;

  // Which nodes survive contraction.
  std::vector<char> kept(n_node, 1);
  for (int v = n_tip + 1; v < n_node; ++v) {
    if (v < static_cast<int>(collapsed.size()) && collapsed[v]) kept[v] = 0;
  }
  // Tips (0..n_tip-1) and root are kept by construction.

  // New 1-based ids: tips keep v+1; retained internals numbered after the tips,
  // root first (ape convention: root id == n_tip + 1).
  std::vector<int> new_id(n_node, 0);
  for (int v = 0; v < n_tip; ++v) new_id[v] = v + 1;
  int next = n_tip + 1;
  new_id[root] = next++;
  for (int v = n_tip + 1; v < n_node; ++v) {
    if (kept[v]) new_id[v] = next++;
  }
  const int n_kept = next - 1;                 // total retained nodes (incl. tips)
  const int n_edge = n_kept - 1;               // rooted tree: edges = nodes - 1

  IntegerMatrix edge(n_edge, 2);
  int row = 0;
  for (int v = 0; v < n_node; ++v) {
    if (v == root || !kept[v]) continue;
    int p = tree.parent[v];                    // nearest retained ancestor
    while (!kept[p]) p = tree.parent[p];        // root is retained -> terminates
    edge(row, 0) = new_id[p];
    edge(row, 1) = new_id[v];
    ++row;
  }
  return edge;
}

// Build a topology-only TreeState (no character data) from a 1-based edge
// matrix.  Used by the least-squares path, which reads only the tree shape.
// Assumes the standard TreeSearch numbering: tips 1..n_tip, root = n_tip + 1,
// first-encountered child of each node goes left.
ts::TreeState build_topology_tree(const IntegerMatrix& edge) {
  int n_edge = edge.nrow();
  int n_tip = n_edge / 2 + 1;

  ts::TreeState tree;
  tree.n_tip = n_tip;
  tree.n_internal = n_tip - 1;
  tree.n_node = 2 * n_tip - 1;
  tree.total_words = 0;
  tree.n_blocks = 0;
  tree.parent.assign(tree.n_node, -1);
  tree.left.assign(tree.n_internal, -1);
  tree.right.assign(tree.n_internal, -1);

  std::vector<int> child_count(tree.n_internal, 0);
  for (int i = 0; i < n_edge; ++i) {
    int par = edge(i, 0) - 1;
    int child = edge(i, 1) - 1;
    tree.parent[child] = par;
    int ni = par - n_tip;
    if (child_count[ni] == 0) {
      tree.left[ni] = child;
    } else {
      tree.right[ni] = child;
    }
    ++child_count[ni];
  }
  tree.parent[n_tip] = n_tip;
  tree.build_postorder();
  return tree;
}

} // anonymous namespace

// Forward declarations for helpers defined later in this file
static ts::ConstraintData build_constraint_from_r(
    int n_tips,
    Nullable<IntegerMatrix> consSplitMatrix,
    Nullable<NumericMatrix> consContrast,
    Nullable<IntegerMatrix> consTipData,
    Nullable<IntegerVector> consWeight,
    Nullable<CharacterVector> consLevels,
    int consExpectedScore);

// [[Rcpp::export]]
double ts_fitch_score(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0,
    Nullable<NumericMatrix> infoAmounts = R_NilValue,
    bool xpiwe = false,
    double xpiwe_r = 0.5,
    double xpiwe_max_f = 5.0,
    IntegerVector obs_count = IntegerVector())
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity, infoAmounts,
                                xpiwe, xpiwe_r, xpiwe_max_f, obs_count);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  return ts::score_tree(tree, ds);
}

// ---------------------------------------------------------------------------
//  Least-squares distance fitting (Lapointe & Cucumel average consensus)
// ---------------------------------------------------------------------------

namespace {

// Map per-branch fitted lengths back onto edge-matrix rows.  Each unrooted
// branch is keyed by its child node; the root's "other" child edge (which has
// no design column) is left at zero so that unrooting sums the pair correctly.
NumericVector ls_edge_lengths(const ts::LSFit& fit, const IntegerMatrix& edge,
                              int n_node) {
  std::vector<double> node_len(n_node, 0.0);
  for (int b = 0; b < fit.n_branch; ++b) {
    node_len[fit.branch_node[b]] = fit.branch_length[b];
  }
  int n_edge = edge.nrow();
  NumericVector elen(n_edge);
  for (int i = 0; i < n_edge; ++i) {
    elen[i] = node_len[edge(i, 1) - 1];
  }
  return elen;
}

const double* ls_weight_ptr(Nullable<NumericMatrix> weight, NumericMatrix& hold) {
  if (weight.isNull()) return nullptr;
  hold = NumericMatrix(weight.get());
  return REAL(hold);
}

} // anonymous namespace

// Fit branch lengths on a fixed topology to a target distance matrix.
// `method`: 0 = OLS (closed form), 1 = NNLS (non-negative, matches phangorn).
// Returns the per-edge lengths (aligned to `edge`'s rows) and the RSS.
// [[Rcpp::export]]
List ts_ls_fit(
    IntegerMatrix edge,
    NumericMatrix dist,
    Nullable<NumericMatrix> weight = R_NilValue,
    int method = 1)
{
  ts::TreeState tree = build_topology_tree(edge);
  NumericMatrix whold;
  const double* wptr = ls_weight_ptr(weight, whold);
  ts::LSData ls = ts::build_ls_data(REAL(dist), tree.n_tip, wptr);
  ts::LSMethod m = (method == 0) ? ts::LSMethod::OLS : ts::LSMethod::NNLS;

  ts::LSFit fit = ts::ls_fit(tree, ls, m);

  return List::create(
    Named("edge_length") = ls_edge_lengths(fit, edge, tree.n_node),
    Named("rss") = fit.rss,
    Named("ok") = fit.ok
  );
}

// Heuristic topology search minimising the least-squares fit to `dist`.
// Alternates NNI and (optionally) SPR hill-climbing until neither improves,
// then fits final branch lengths on the best topology.
// [[Rcpp::export]]
List ts_ls_search(
    IntegerMatrix edge,
    NumericMatrix dist,
    Nullable<NumericMatrix> weight = R_NilValue,
    int method = 1,
    int maxHits = 1,
    bool doSpr = true)
{
  ts::TreeState tree = build_topology_tree(edge);
  NumericMatrix whold;
  const double* wptr = ls_weight_ptr(weight, whold);
  ts::LSData ls = ts::build_ls_data(REAL(dist), tree.n_tip, wptr);
  ts::LSMethod m = (method == 0) ? ts::LSMethod::OLS : ts::LSMethod::NNLS;

  int total_moves = 0, total_iters = 0;
  ts::LSSearchResult r = ts::ls_nni_search(tree, ls, m, maxHits);
  double rss = r.rss;
  total_moves += r.n_moves;
  total_iters += r.n_iterations;

  if (doSpr) {
    bool improving = true;
    while (improving) {
      double before = rss;
      ts::LSSearchResult rs = ts::ls_spr_search(tree, ls, m, maxHits);
      ts::LSSearchResult rn = ts::ls_nni_search(tree, ls, m, maxHits);
      rss = rn.rss;
      total_moves += rs.n_moves + rn.n_moves;
      total_iters += rs.n_iterations + rn.n_iterations;
      improving = (rss < before - 1e-9);
    }
  }

  ts::LSFit fit = ts::ls_fit(tree, ls, m);
  IntegerMatrix out_edge = tree_to_edge(tree);

  return List::create(
    Named("edge") = out_edge,
    Named("edge_length") = ls_edge_lengths(fit, out_edge, tree.n_node),
    Named("rss") = fit.rss,
    Named("n_moves") = total_moves,
    Named("n_iterations") = total_iters
  );
}

// Diagnostic: per-node collapsed flags for a scored tree, for validating the
// aggressive (min-length-0) collapse criterion against a brute-force MPR oracle
// (dev/benchmarks/b2_minlength_oracle.R + b2_collapsed_kernel_validate.R).
// Returns 0-based parent[] (root = n_tip, self-parent), n_tip, and the flag
// vector.  Explicit down+up pass guarantees prelim/final_ are fresh.
// [[Rcpp::export]]
List ts_collapsed_flags_debug(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    bool aggressive)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);
  ts::TreeState tree;
  tree.init_from_edge(&edge(0, 0), &edge(0, 1), edge.nrow(), ds);
  ts::fitch_downpass(tree, ds);
  ts::fitch_uppass(tree, ds);              // populate final_ (MPR sets)

  std::vector<uint8_t> collapsed;
  if (aggressive) ts::compute_collapsed_flags_aggressive(tree, ds, collapsed);
  else            ts::compute_collapsed_flags(tree, ds, collapsed);

  IntegerVector parent(tree.n_node), flag(tree.n_node);
  for (int i = 0; i < tree.n_node; ++i) {
    parent[i] = tree.parent[i];            // 0-based; root (n_tip) is self-parent
    flag[i] = collapsed[i];
  }
  return List::create(
    Named("parent") = parent,
    Named("collapsed") = flag,
    Named("n_tip") = tree.n_tip,
    Named("n_node") = tree.n_node
  );
}

// [[Rcpp::export]]
List ts_na_debug_char(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int target_pattern)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);
  ts::TreeState tree;
  tree.init_from_edge(&edge(0, 0), &edge(0, 1), edge.nrow(), ds);
  int total_score = ts::fitch_na_score(tree, ds);

  // Find which block/bit contains target_pattern (0-based)
  int tgt = target_pattern - 1;  // R 1-based to C 0-based
  int tgt_block = -1, tgt_bit = -1;
  for (int b = 0; b < ds.n_blocks; ++b) {
    for (int c = 0; c < ds.blocks[b].n_chars; ++c) {
      if (ds.blocks[b].pattern_index[c] == tgt) {
        tgt_block = b;
        tgt_bit = c;
      }
    }
  }
  if (tgt_block < 0) return List::create(Named("error") = "pattern not found");

  const ts::CharBlock& blk = ds.blocks[tgt_block];
  int off = ds.block_word_offset[tgt_block];
  int k = blk.n_states;
  uint64_t mask = 1ULL << tgt_bit;

  // Per-node info
  int n_node = tree.n_node;
  IntegerVector node_id(n_node);
  CharacterVector prelim_str(n_node), final_str(n_node), down2_str(n_node);
  IntegerVector is_step(n_node);

  auto state_str = [&](const uint64_t* base) -> std::string {
    std::string s;
    for (int st = 0; st < k; ++st) {
      if (base[st] & mask) {
        if (!s.empty()) s += "/";
        if (blk.has_inapplicable && st == 0) s += "-";
        else s += std::to_string(st - (blk.has_inapplicable ? 1 : 0));
      }
    }
    return s.empty() ? "." : s;
  };

  for (int nd = 0; nd < n_node; ++nd) {
    node_id[nd] = nd + 1;
    size_t base = static_cast<size_t>(nd) * tree.total_words + off;
    prelim_str[nd] = state_str(&tree.prelim[base]);
    final_str[nd] = state_str(&tree.final_[base]);
    down2_str[nd] = state_str(&tree.down2[base]);
    is_step[nd] = 0;
  }

  // Determine which nodes had steps in Pass 3
  // Must match fitch_na_score: l_act & r_act & ~(ss_app & any_d2_isect)
  for (int node : tree.postorder) {
    if (!blk.has_inapplicable) continue;
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    size_t nb = static_cast<size_t>(node) * tree.total_words + off;
    size_t lb = static_cast<size_t>(lc) * tree.total_words + off;
    size_t rb = static_cast<size_t>(rc) * tree.total_words + off;
    const uint64_t* F = &tree.final_[nb];
    const uint64_t* L2 = &tree.down2[lb];
    const uint64_t* R2 = &tree.down2[rb];

    uint64_t ss_app = 0;
    for (int s = 1; s < k; ++s) ss_app |= F[s];

    uint64_t any_d2_isect = 0;
    for (int s = 0; s < k; ++s) any_d2_isect |= (L2[s] & R2[s]);

    const uint64_t* la = &tree.subtree_actives[lb];
    const uint64_t* ra = &tree.subtree_actives[rb];
    uint64_t l_act = 0, r_act = 0;
    for (int s = 1; s < k; ++s) { l_act |= la[s]; r_act |= ra[s]; }

    uint64_t needs_step = l_act & r_act
                        & ~(ss_app & any_d2_isect) & blk.active_mask;
    if (needs_step & mask) is_step[node] = 1;
  }

  // Parent info
  IntegerVector parent_id(n_node);
  IntegerVector left_child(tree.n_internal), right_child(tree.n_internal);
  for (int nd = 0; nd < n_node; ++nd) parent_id[nd] = tree.parent[nd] + 1;
  for (int ni = 0; ni < tree.n_internal; ++ni) {
    left_child[ni] = tree.left[ni] + 1;
    right_child[ni] = tree.right[ni] + 1;
  }

  return List::create(
    Named("total_score") = total_score,
    Named("block") = tgt_block + 1,
    Named("bit") = tgt_bit,
    Named("has_inapp") = blk.has_inapplicable,
    Named("n_states") = k,
    Named("node") = node_id,
    Named("parent") = parent_id,
    Named("left") = left_child,
    Named("right") = right_child,
    Named("prelim") = prelim_str,
    Named("final_state") = final_str,
    Named("down2") = down2_str,
    Named("is_step") = is_step,
    Named("n_tip") = tree.n_tip
  );
}

// [[Rcpp::export]]
List ts_na_char_steps(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);

  ts::TreeState tree;
  tree.init_from_edge(&edge(0, 0), &edge(0, 1), edge.nrow(), ds);

  // Run the three-pass scoring
  int total_score = ts::fitch_na_score(tree, ds);

  int n_pat = tip_data.ncol();
  IntegerVector steps(n_pat, 0);
  int debug_std_total = 0, debug_na_total = 0;

  // Standard blocks: count from local_cost
  for (int node : tree.postorder) {
    for (int b = 0; b < ds.n_blocks; ++b) {
      const ts::CharBlock& blk = ds.blocks[b];
      if (blk.has_inapplicable) continue;
      uint64_t mask = tree.local_cost[static_cast<size_t>(node) * tree.n_blocks + b];
      debug_std_total += blk.weight * ts::popcount64(mask);
      for (int c = 0; c < blk.n_chars; ++c) {
        if (mask & (1ULL << c)) {
          steps[blk.pattern_index[c]] += 1;
        }
      }
    }
  }

  // NA blocks: use same formula as Pass 3 (needs subtree_actives)
  for (int node : tree.postorder) {
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    size_t nb = static_cast<size_t>(node) * tree.total_words;
    size_t lb = static_cast<size_t>(lc) * tree.total_words;
    size_t rb = static_cast<size_t>(rc) * tree.total_words;

    for (int b = 0; b < ds.n_blocks; ++b) {
      if (!ds.blocks[b].has_inapplicable) continue;
      const ts::CharBlock& blk = ds.blocks[b];
      int off = ds.block_word_offset[b];
      int k = blk.n_states;

      const uint64_t* L2 = &tree.down2[lb + off];
      const uint64_t* R2 = &tree.down2[rb + off];
      const uint64_t* D2 = &tree.down2[nb + off];

      uint64_t ss_app = 0;
      for (int s = 1; s < k; ++s) ss_app |= D2[s];

      uint64_t any_d2_isect = 0;
      for (int s = 0; s < k; ++s) any_d2_isect |= (L2[s] & R2[s]);

      const uint64_t* la = &tree.subtree_actives[lb + off];
      const uint64_t* ra = &tree.subtree_actives[rb + off];
      uint64_t l_act = 0, r_act = 0;
      for (int s = 1; s < k; ++s) { l_act |= la[s]; r_act |= ra[s]; }

      uint64_t needs_step = l_act & r_act
                          & ~(ss_app & any_d2_isect) & blk.active_mask;
      debug_na_total += blk.weight * ts::popcount64(needs_step);

      for (int c = 0; c < blk.n_chars; ++c) {
        if (needs_step & (1ULL << c)) {
          steps[blk.pattern_index[c]] += 1;
        }
      }
    }
  }

  // Block info
  int n_std = 0, n_na = 0;
  for (int b = 0; b < ds.n_blocks; ++b) {
    if (ds.blocks[b].has_inapplicable) ++n_na; else ++n_std;
  }

  return List::create(
    Named("steps") = steps,
    Named("total_score") = total_score,
    Named("std_total") = debug_std_total,
    Named("na_total") = debug_na_total,
    Named("debug_sum") = debug_std_total + debug_na_total,
    Named("n_std_blocks") = n_std,
    Named("n_na_blocks") = n_na
  );
}

// [[Rcpp::export]]
IntegerVector ts_char_steps(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);

  ts::TreeState tree;
  tree.init_from_edge(&edge(0, 0), &edge(0, 1), edge.nrow(), ds);

  // Score tree (handles both standard and NA-aware scoring)
  ts::score_tree(tree, ds);

  // Extract per-pattern step counts
  int n_pat = tip_data.ncol();
  std::vector<int> char_steps(n_pat, 0);
  ts::extract_char_steps(tree, ds, char_steps);

  // Add precomputed steps from simplification
  IntegerVector result(n_pat);
  for (int p = 0; p < n_pat; ++p) {
    result[p] = char_steps[p];
    if (!ds.precomputed_steps.empty()) {
      result[p] += ds.precomputed_steps[p];
    }
  }
  return result;
}

// [[Rcpp::export]]
List ts_debug_clip(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int clip_node_1based)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);
  ts::TreeState tree;
  tree.init_from_edge(&edge(0, 0), &edge(0, 1), edge.nrow(), ds);

  int whole_score = ts::fitch_na_score(tree, ds);
  int clip_node = clip_node_1based - 1;  // 0-based

  tree.spr_clip(clip_node);
  tree.build_postorder();

  int divided_score = ts::fitch_score(tree, ds);

  // Also check postorder size
  int postorder_size = static_cast<int>(tree.postorder.size());

  // Unclip
  tree.spr_unclip();
  tree.build_postorder();
  int restored_score = ts::fitch_score(tree, ds);

  return List::create(
    Named("whole_score") = whole_score,
    Named("divided_score") = divided_score,
    Named("restored_score") = restored_score,
    Named("postorder_size") = postorder_size,
    Named("n_internal") = tree.n_internal,
    Named("clip_parent_0based") = tree.clip_state.clip_parent,
    Named("clip_sibling_0based") = tree.clip_state.clip_sibling,
    Named("clip_grandpar_0based") = tree.clip_state.clip_grandpar
  );
}

// [[Rcpp::export]]
List ts_test_indirect(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int clip_node_1based,
    int above_1based,
    int below_1based)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);
  ts::TreeState tree;
  tree.init_from_edge(&edge(0, 0), &edge(0, 1), edge.nrow(), ds);

  int whole_score = ts::fitch_na_score(tree, ds);
  int clip_node = clip_node_1based - 1;
  int above = above_1based - 1;
  int below = below_1based - 1;

  // Clip
  tree.spr_clip(clip_node);
  tree.build_postorder();
  tree.reset_states(ds);
  int main_score = ts::fitch_score(tree, ds);

  // Score clipped subtree
  int clip_score = 0;
  {
    std::vector<int> clip_stack;
    clip_stack.push_back(clip_node);
    std::vector<int> clip_preorder;
    while (!clip_stack.empty()) {
      int nd = clip_stack.back();
      clip_stack.pop_back();
      if (nd < tree.n_tip) continue;
      clip_preorder.push_back(nd);
      int ni = nd - tree.n_tip;
      clip_stack.push_back(tree.right[ni]);
      clip_stack.push_back(tree.left[ni]);
    }
    for (int j = static_cast<int>(clip_preorder.size()) - 1; j >= 0; --j) {
      int nd = clip_preorder[j];
      int ni = nd - tree.n_tip;
      int lc = tree.left[ni];
      int rc = tree.right[ni];
      for (int b = 0; b < ds.n_blocks; ++b) {
        const ts::CharBlock& blk = ds.blocks[b];
        int offset = ds.block_word_offset[b];
        clip_score += blk.weight * ts::fitch_downpass_node(
            &tree.prelim[static_cast<size_t>(lc) * tree.total_words + offset],
            &tree.prelim[static_cast<size_t>(rc) * tree.total_words + offset],
            &tree.prelim[static_cast<size_t>(nd) * tree.total_words + offset],
            blk.n_states, blk.active_mask);
      }
    }
  }

  int divided_length = main_score + clip_score;
  const uint64_t* clip_prelim =
      &tree.prelim[static_cast<size_t>(clip_node) * tree.total_words];
  int extra = ts::fitch_indirect_length(clip_prelim, tree, ds, above, below);
  int candidate = divided_length + extra;

  // Now regraft and get actual score
  tree.spr_regraft(above, below);
  tree.build_postorder();
  tree.reset_states(ds);
  int actual = ts::fitch_na_score(tree, ds);

  return List::create(
    Named("whole_score") = whole_score,
    Named("main_score") = main_score,
    Named("clip_score") = clip_score,
    Named("divided_length") = divided_length,
    Named("extra") = extra,
    Named("candidate") = candidate,
    Named("actual") = actual,
    Named("match") = (candidate == actual)
  );
}

// [[Rcpp::export]]
List ts_spr_search(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int maxHits = 20,
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  ts::SearchResult result = ts::spr_search(tree, ds, maxHits);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.score,
    Named("n_moves") = result.n_moves,
    Named("n_iterations") = result.n_iterations
  );
}

// [[Rcpp::export]]
List ts_tbr_search(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int maxHits = 1,
    bool acceptEqual = false,
    int maxChanges = 0,
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  ts::TBRParams params;
  params.max_hits = maxHits;
  params.accept_equal = acceptEqual;
  params.max_accepted_changes = maxChanges;

  ts::TBRResult result = ts::tbr_search(tree, ds, params);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.best_score,
    Named("n_accepted") = result.n_accepted,
    Named("n_evaluated") = result.n_evaluated,
    Named("n_zero_skipped") = result.n_zero_skipped,
    Named("converged") = result.converged
  );
}

// [[Rcpp::export]]
List ts_ratchet_search(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int nCycles = 10,
    double perturbProb = 0.04,
    int maxHits = 1,
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0,
    int perturbMode = 0,
    int perturbMaxMoves = 0,
    bool adaptive = false,
    double targetEscapeRate = 0.3)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  ts::RatchetParams params;
  params.n_cycles = nCycles;
  params.perturb_prob = perturbProb;
  params.max_hits = maxHits;
  params.perturb_mode = static_cast<ts::PerturbMode>(perturbMode);
  params.perturb_max_moves = perturbMaxMoves;
  params.adaptive = adaptive;
  params.target_escape_rate = targetEscapeRate;

  ts::RatchetResult result = ts::ratchet_search(tree, ds, params);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.best_score,
    Named("n_cycles") = result.n_cycles_completed,
    Named("total_tbr_moves") = result.total_tbr_moves,
    Named("n_escapes") = result.n_escapes,
    Named("final_perturb_prob") = result.final_perturb_prob
  );
}

// [[Rcpp::export]]
List ts_drift_search(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int nCycles = 10,
    int afdLimit = 3,
    double rfdLimit = 0.1,
    int maxHits = 1,
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  ts::DriftParams params;
  params.n_cycles = nCycles;
  params.afd_limit = afdLimit;
  params.rfd_limit = rfdLimit;
  params.max_hits = maxHits;

  ts::DriftResult result = ts::drift_search(tree, ds, params);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.best_score,
    Named("n_cycles_completed") = result.n_cycles_completed,
    Named("total_drift_moves") = result.total_drift_moves,
    Named("total_tbr_moves") = result.total_tbr_moves
  );
}

// [[Rcpp::export]]
List ts_wagner_tree(
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    IntegerVector addition_order = IntegerVector(),
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0,
    Nullable<NumericMatrix> infoAmounts = R_NilValue,
    Nullable<IntegerMatrix> consSplitMatrix = R_NilValue,
    Nullable<NumericMatrix> consContrast = R_NilValue,
    Nullable<IntegerMatrix> consTipData = R_NilValue,
    Nullable<IntegerVector> consWeight = R_NilValue,
    Nullable<CharacterVector> consLevels = R_NilValue,
    int consExpectedScore = 0)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity, infoAmounts);

  int n_tips = tip_data.nrow();
  ts::ConstraintData cd = build_constraint_from_r(
      n_tips, consSplitMatrix, consContrast, consTipData,
      consWeight, consLevels, consExpectedScore);
  ts::ConstraintData* cd_ptr = cd.active ? &cd : nullptr;

  std::vector<int> order;
  if (addition_order.size() > 0) {
    // Convert from 1-based R indices to 0-based
    order.resize(addition_order.size());
    for (int i = 0; i < addition_order.size(); ++i) {
      order[i] = addition_order[i] - 1;
    }
  }

  ts::TreeState tree;
  ts::WagnerResult result = ts::wagner_tree(tree, ds, order, cd_ptr);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.score
  );
}

// [[Rcpp::export]]
List ts_random_wagner_tree(
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0,
    Nullable<NumericMatrix> infoAmounts = R_NilValue,
    Nullable<IntegerMatrix> consSplitMatrix = R_NilValue,
    Nullable<NumericMatrix> consContrast = R_NilValue,
    Nullable<IntegerMatrix> consTipData = R_NilValue,
    Nullable<IntegerVector> consWeight = R_NilValue,
    Nullable<CharacterVector> consLevels = R_NilValue,
    int consExpectedScore = 0)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity, infoAmounts);

  int n_tips = tip_data.nrow();
  ts::ConstraintData cd = build_constraint_from_r(
      n_tips, consSplitMatrix, consContrast, consTipData,
      consWeight, consLevels, consExpectedScore);
  ts::ConstraintData* cd_ptr = cd.active ? &cd : nullptr;

  ts::TreeState tree;
  ts::WagnerResult result = ts::random_wagner_tree(tree, ds, cd_ptr);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.score
  );
}

// [[Rcpp::export]]
List ts_compute_splits(
    IntegerMatrix edge,
    int n_tip)
{
  // Build a minimal TreeState from the edge matrix (no character data needed).
  ts::TreeState tree;
  // We need a dummy DataSet just for init_from_edge dimensions.
  // Instead, manually initialize the tree topology.
  tree.n_tip = n_tip;
  tree.n_internal = n_tip - 1;
  tree.n_node = 2 * n_tip - 1;
  tree.total_words = 0;
  tree.n_blocks = 0;

  tree.parent.assign(tree.n_node, 0);
  tree.left.assign(tree.n_internal, 0);
  tree.right.assign(tree.n_internal, 0);

  // Parse edge matrix (R is 1-based)
  int n_edge = edge.nrow();
  for (int i = 0; i < n_edge; ++i) {
    int par = edge(i, 0) - 1;
    int child = edge(i, 1) - 1;
    tree.parent[child] = par;
    int ni = par - n_tip;
    if (tree.left[ni] == 0 && par != child) {
      // Check if left is unset (0 could be tip 0, but root's self-reference
      // won't use this path). We need a cleaner approach.
      // Actually: first child encountered goes left, second goes right.
      // Use a flag: if left == 0 and we haven't set it yet...
      // Simpler: track how many children we've seen per node.
      tree.left[ni] = child;
    } else {
      tree.right[ni] = child;
    }
  }
  // Fix: root's parent is itself
  tree.parent[n_tip] = n_tip;

  // Fix left/right assignment: re-parse properly
  // Reset and use a counter
  std::vector<int> child_count(tree.n_internal, 0);
  tree.left.assign(tree.n_internal, -1);
  tree.right.assign(tree.n_internal, -1);
  for (int i = 0; i < n_edge; ++i) {
    int par = edge(i, 0) - 1;
    int child = edge(i, 1) - 1;
    tree.parent[child] = par;
    int ni = par - n_tip;
    if (child_count[ni] == 0) {
      tree.left[ni] = child;
    } else {
      tree.right[ni] = child;
    }
    ++child_count[ni];
  }
  tree.parent[n_tip] = n_tip;

  tree.build_postorder();

  ts::SplitSet ss = ts::compute_splits(tree);

  // Return as a list of raw vectors (one per split)
  List result(ss.n_splits);
  for (int i = 0; i < ss.n_splits; ++i) {
    const uint64_t* s = ss.split(i);
    // Convert to an integer vector of tip indices in the split
    IntegerVector tips;
    for (int t = 0; t < n_tip; ++t) {
      int w = t / 64;
      int b = t % 64;
      if (s[w] & (1ULL << b)) {
        tips.push_back(t + 1);  // 1-based for R
      }
    }
    result[i] = tips;
  }

  return result;
}

// [[Rcpp::export]]
bool ts_trees_equal(
    IntegerMatrix edge1,
    IntegerMatrix edge2,
    int n_tip)
{
  auto build_tree = [&](IntegerMatrix edge) -> ts::TreeState {
    ts::TreeState tree;
    tree.n_tip = n_tip;
    tree.n_internal = n_tip - 1;
    tree.n_node = 2 * n_tip - 1;
    tree.total_words = 0;
    tree.n_blocks = 0;
    tree.parent.assign(tree.n_node, 0);
    tree.left.assign(tree.n_internal, -1);
    tree.right.assign(tree.n_internal, -1);

    int n_edge = edge.nrow();
    std::vector<int> child_count(tree.n_internal, 0);
    for (int i = 0; i < n_edge; ++i) {
      int par = edge(i, 0) - 1;
      int child = edge(i, 1) - 1;
      tree.parent[child] = par;
      int ni = par - n_tip;
      if (child_count[ni] == 0) {
        tree.left[ni] = child;
      } else {
        tree.right[ni] = child;
      }
      ++child_count[ni];
    }
    tree.parent[n_tip] = n_tip;
    tree.build_postorder();
    return tree;
  };

  ts::TreeState t1 = build_tree(edge1);
  ts::TreeState t2 = build_tree(edge2);

  ts::SplitSet s1 = ts::compute_splits(t1);
  ts::SplitSet s2 = ts::compute_splits(t2);

  return ts::splits_equal(s1, s2);
}

// [[Rcpp::export]]
List ts_pool_test(
    List edges,        // list of edge matrices
    NumericVector scores,
    int n_tip,
    int max_size = 100,
    double suboptimal = 0.0)
{
  auto build_tree = [&](IntegerMatrix edge) -> ts::TreeState {
    ts::TreeState tree;
    tree.n_tip = n_tip;
    tree.n_internal = n_tip - 1;
    tree.n_node = 2 * n_tip - 1;
    tree.total_words = 0;
    tree.n_blocks = 0;
    tree.parent.assign(tree.n_node, 0);
    tree.left.assign(tree.n_internal, -1);
    tree.right.assign(tree.n_internal, -1);

    int n_edge = edge.nrow();
    std::vector<int> child_count(tree.n_internal, 0);
    for (int i = 0; i < n_edge; ++i) {
      int par = edge(i, 0) - 1;
      int child = edge(i, 1) - 1;
      tree.parent[child] = par;
      int ni = par - n_tip;
      if (child_count[ni] == 0) {
        tree.left[ni] = child;
      } else {
        tree.right[ni] = child;
      }
      ++child_count[ni];
    }
    tree.parent[n_tip] = n_tip;
    tree.build_postorder();
    return tree;
  };

  ts::TreePool pool(max_size, suboptimal);

  int n_trees = edges.size();
  LogicalVector added(n_trees);
  for (int i = 0; i < n_trees; ++i) {
    IntegerMatrix e = as<IntegerMatrix>(edges[i]);
    ts::TreeState t = build_tree(e);
    added[i] = pool.add(t, scores[i]);
  }

  return List::create(
    Named("added") = added,
    Named("pool_size") = pool.size(),
    Named("best_score") = pool.best_score(),
    Named("hits_to_best") = pool.hits_to_best()
  );
}

// [[Rcpp::export]]
List ts_nni_search(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int maxHits = 20,
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  ts::SearchResult result = ts::nni_search(tree, ds, maxHits);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.score,
    Named("n_moves") = result.n_moves,
    Named("n_iterations") = result.n_iterations
  );
}

// [[Rcpp::export]]
List ts_tree_fuse(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    List pool_edges,
    NumericVector pool_scores,
    bool accept_equal = false,
    int max_rounds = 10,
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  // Build recipient tree
  ts::TreeState recipient;
  recipient.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  // Build pool from R-side inputs
  ts::TreePool pool(static_cast<int>(pool_edges.size()), 1e18);
  for (int i = 0; i < pool_edges.size(); ++i) {
    IntegerMatrix pe = as<IntegerMatrix>(pool_edges[i]);
    ts::TreeState t;
    t.init_from_edge(&pe(0, 0), &pe(0, 1), pe.nrow(), ds);
    double sc = pool_scores[i];
    pool.add(t, sc);
  }

  ts::FuseParams params;
  params.accept_equal = accept_equal;
  params.max_rounds = max_rounds;

  ts::FuseResult result = ts::tree_fuse(recipient, ds, pool, params);

  return List::create(
    Named("edge") = tree_to_edge(recipient),
    Named("score") = result.best_score,
    Named("n_exchanges") = result.n_exchanges,
    Named("n_rounds") = result.n_rounds
  );
}

// [[Rcpp::export]]
List ts_sector_diag(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int sector_root_1based)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  int full_score = ts::fitch_na_score(tree, ds);

  int sector_root = sector_root_1based - 1;
  int clade_size = ts::count_clade_tips(tree, sector_root);

  ts::ReducedDataset rd = ts::build_reduced_dataset(tree, ds, sector_root);

  int sector_score = ts::fitch_score(rd.subtree, rd.data);

  return List::create(
    Named("full_score") = full_score,
    Named("sector_root") = sector_root,
    Named("clade_size") = clade_size,
    Named("n_sector_tips") = rd.subtree.n_tip,
    Named("n_sector_nodes") = rd.subtree.n_node,
    Named("sector_score") = sector_score,
    Named("sector_total_words") = rd.data.total_words,
    Named("sector_n_blocks") = rd.data.n_blocks
  );
}

// [[Rcpp::export]]
List ts_rss_search(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int minSectorSize = 6,
    int maxSectorSize = 50,
    bool acceptEqual = false,
    int rssPicks = 0,
    int ratchetCycles = 6,
    int maxHits = 1,
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  ts::SectorParams params;
  params.min_sector_size = minSectorSize;
  params.max_sector_size = maxSectorSize;
  params.accept_equal = acceptEqual;
  params.rss_picks_per_round = rssPicks;
  params.internal_ratchet_cycles = ratchetCycles;
  params.internal_max_hits = maxHits;

  ts::SectorResult result = ts::rss_search(tree, ds, params);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.best_score,
    Named("n_sectors_searched") = result.n_sectors_searched,
    Named("n_sectors_improved") = result.n_sectors_improved,
    Named("total_steps_saved") = result.total_steps_saved
  );
}

// [[Rcpp::export]]
List ts_xss_search(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int nPartitions = 4,
    int xssRounds = 3,
    bool acceptEqual = false,
    int ratchetCycles = 6,
    int maxHits = 1,
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  ts::SectorParams params;
  params.n_partitions = nPartitions;
  params.xss_rounds = xssRounds;
  params.accept_equal = acceptEqual;
  params.internal_ratchet_cycles = ratchetCycles;
  params.internal_max_hits = maxHits;

  ts::SectorResult result = ts::xss_search(tree, ds, params);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.best_score,
    Named("n_sectors_searched") = result.n_sectors_searched,
    Named("n_sectors_improved") = result.n_sectors_improved,
    Named("total_steps_saved") = result.total_steps_saved
  );
}

// Helper: build ConstraintData from Nullable R objects.
// Returns cd with cd.active == false if no constraint provided.
static ts::ConstraintData build_constraint_from_r(
    int n_tips,
    Nullable<IntegerMatrix> consSplitMatrix,
    Nullable<NumericMatrix> consContrast,
    Nullable<IntegerMatrix> consTipData,
    Nullable<IntegerVector> consWeight,
    Nullable<CharacterVector> consLevels,
    int consExpectedScore)
{
  ts::ConstraintData cd;
  if (consSplitMatrix.isNotNull()) {
    IntegerMatrix csm(consSplitMatrix.get());
    int n_cons_splits = csm.nrow();
    if (n_cons_splits > 0 && csm.ncol() == n_tips) {
      cd = ts::build_constraint(INTEGER(csm), n_cons_splits, n_tips);

      if (consContrast.isNotNull() && consTipData.isNotNull() &&
          consWeight.isNotNull() && consLevels.isNotNull()) {
        NumericMatrix cc(consContrast.get());
        IntegerMatrix ctd(consTipData.get());
        IntegerVector cw(consWeight.get());
        CharacterVector cl(consLevels.get());

        int n_cons_tokens = cc.nrow();
        int n_cons_states = cc.ncol();
        int n_cons_patterns = ctd.ncol();

        std::vector<std::string> cons_level_strs(n_cons_states);
        std::vector<const char*> cons_level_ptrs(n_cons_states);
        for (int i = 0; i < n_cons_states; ++i) {
          cons_level_strs[i] = as<std::string>(cl[i]);
          cons_level_ptrs[i] = cons_level_strs[i].c_str();
        }

        ts::build_constraint_posthoc(
            cd,
            REAL(cc), n_cons_tokens, n_cons_states,
            INTEGER(ctd), n_tips, n_cons_patterns,
            INTEGER(cw),
            cons_level_ptrs.data(),
            consExpectedScore);
      }
    }
  }
  return cd;
}

// Helper: extract profile info_amounts pointer and max_steps from Nullable.
static void extract_info_amounts(
    Nullable<NumericMatrix> infoAmounts,
    const double*& info_amounts_ptr,
    int& info_max_steps)
{
  info_amounts_ptr = nullptr;
  info_max_steps = 0;
  if (infoAmounts.isNotNull()) {
    NumericMatrix ia(infoAmounts.get());
    info_amounts_ptr = REAL(ia);
    info_max_steps = ia.nrow();
  }
}

// --- Helpers to unpack R lists into C++ structs ---

// Unpack a searchControl list (from R's SearchControl object) into DrivenParams.
// Anneal fields are now part of searchControl (no longer a separate list).
static void unpack_search_control(List ctrl, ts::DrivenParams& params) {
  // TBR / initial descent
  params.tbr_max_hits       = as<int>(ctrl["tbrMaxHits"]);
  params.clip_order         = as<int>(ctrl["clipOrder"]);
  params.nni_first          = as<bool>(ctrl["nniFirst"]);
  params.spr_first          = as<bool>(ctrl["sprFirst"]);
  params.tabu_size          = as<int>(ctrl["tabuSize"]);
  params.wagner_starts      = as<int>(ctrl["wagnerStarts"]);
  params.wagner_bias        = as<int>(ctrl["wagnerBias"]);
  params.wagner_bias_temp   = as<double>(ctrl["wagnerBiasTemp"]);
  params.outer_cycles       = as<int>(ctrl["outerCycles"]);
  params.max_outer_resets   = as<int>(ctrl["maxOuterResets"]);

  // Ratchet
  params.ratchet_cycles           = as<int>(ctrl["ratchetCycles"]);
  params.ratchet_perturb_prob     = as<double>(ctrl["ratchetPerturbProb"]);
  params.ratchet_perturb_mode     = as<int>(ctrl["ratchetPerturbMode"]);
  params.ratchet_perturb_max_moves = as<int>(ctrl["ratchetPerturbMaxMoves"]);
  params.ratchet_adaptive         = as<bool>(ctrl["ratchetAdaptive"]);
  if (ctrl.containsElementNamed("ratchetTaper"))
    params.ratchet_taper          = as<bool>(ctrl["ratchetTaper"]);
  if (ctrl.containsElementNamed("stallEscalateFactor"))
    params.stall_escalate_factor  = as<double>(ctrl["stallEscalateFactor"]);

  // NNI perturbation
  params.nni_perturb_cycles   = as<int>(ctrl["nniPerturbCycles"]);
  params.nni_perturb_fraction = as<double>(ctrl["nniPerturbFraction"]);

  // Drift
  params.drift_cycles    = as<int>(ctrl["driftCycles"]);
  params.drift_afd_limit = as<int>(ctrl["driftAfdLimit"]);
  params.drift_rfd_limit = as<double>(ctrl["driftRfdLimit"]);

  // Sectorial search
  params.xss_rounds      = as<int>(ctrl["xssRounds"]);
  params.xss_partitions  = as<int>(ctrl["xssPartitions"]);
  params.rss_rounds      = as<int>(ctrl["rssRounds"]);
  params.css_rounds      = as<int>(ctrl["cssRounds"]);
  params.css_partitions  = as<int>(ctrl["cssPartitions"]);
  params.sector_min_size = as<int>(ctrl["sectorMinSize"]);
  params.sector_max_size = as<int>(ctrl["sectorMaxSize"]);
  params.post_ratchet_sectorial = as<bool>(ctrl["postRatchetSectorial"]);
  if (ctrl.containsElementNamed("rasStarts"))
    params.ras_starts = as<int>(ctrl["rasStarts"]);
  if (ctrl.containsElementNamed("sectorAcceptEqual"))
    params.sector_accept_equal = as<bool>(ctrl["sectorAcceptEqual"]);
  if (ctrl.containsElementNamed("sectorMaxHits"))
    params.sector_max_hits = as<int>(ctrl["sectorMaxHits"]);
  if (ctrl.containsElementNamed("sectorCollapseTarget"))
    params.sector_collapse_target = as<int>(ctrl["sectorCollapseTarget"]);

  // Fuse / pool
  params.fuse_interval      = as<int>(ctrl["fuseInterval"]);
  params.fuse_accept_equal  = as<bool>(ctrl["fuseAcceptEqual"]);
  params.intra_fuse         = as<bool>(ctrl["intraFuse"]);
  params.pool_max_size      = as<int>(ctrl["poolMaxSize"]);
  params.pool_suboptimal    = as<double>(ctrl["poolSuboptimal"]);

  // Stopping / adaptive
  params.consensus_stable_reps = as<int>(ctrl["consensusStableReps"]);
  params.perturb_stop_factor   = as<int>(ctrl["perturbStopFactor"]);
  params.adaptive_level        = as<bool>(ctrl["adaptiveLevel"]);
  params.consensus_constrain   = as<bool>(ctrl["consensusConstrain"]);
  params.adaptive_start        = as<bool>(ctrl["adaptiveStart"]);
  params.enum_time_fraction    = as<double>(ctrl["enumTimeFraction"]);

  // Taxon pruning-reinsertion (T-266)
  params.prune_reinsert_cycles    = as<int>(ctrl["pruneReinsertCycles"]);
  params.prune_reinsert_drop      = as<double>(ctrl["pruneReinsertDrop"]);
  params.prune_reinsert_selection = as<int>(ctrl["pruneReinsertSelection"]);
  params.prune_reinsert_tbr_moves = as<int>(ctrl["pruneReinsertTbrMoves"]);
  params.prune_reinsert_full_moves = as<int>(ctrl["pruneReinsertFullMoves"]);
  params.prune_reinsert_nni       = as<int>(ctrl["pruneReinsertNni"]);

  // Simulated annealing perturbation (PCSA)
  params.anneal_cycles          = as<int>(ctrl["annealCycles"]);
  params.anneal_phases          = as<int>(ctrl["annealPhases"]);
  params.anneal_t_start         = as<double>(ctrl["annealTStart"]);
  params.anneal_t_end           = as<double>(ctrl["annealTEnd"]);
  params.anneal_moves_per_phase = as<int>(ctrl["annealMovesPerPhase"]);
}

// Unpack a runtimeConfig list into DrivenParams fields + return nThreads.
// startEdge and progressCallback are handled separately (need special types).
static int unpack_runtime(List rt, ts::DrivenParams& params) {
  params.max_replicates = as<int>(rt["maxReplicates"]);
  params.target_hits    = as<int>(rt["targetHits"]);
  params.max_seconds    = as<double>(rt["maxSeconds"]);
  params.verbosity      = as<int>(rt["verbosity"]);
  int nThreads          = as<int>(rt["nThreads"]);

  // Starting tree edge matrix (optional)
  if (rt.containsElementNamed("startEdge") &&
      !Rf_isNull(rt["startEdge"])) {
    IntegerMatrix se = as<IntegerMatrix>(rt["startEdge"]);
    int n_edge = se.nrow();
    params.start_n_edge = n_edge;
    params.start_edge.resize(2 * n_edge);
    for (int i = 0; i < n_edge; ++i) {
      params.start_edge[i] = se(i, 0);
      params.start_edge[n_edge + i] = se(i, 1);
    }
  }

  // Progress callback (optional)
  if (rt.containsElementNamed("progressCallback") &&
      !Rf_isNull(rt["progressCallback"])) {
    Rcpp::Function r_cb = as<Function>(rt["progressCallback"]);
    params.progress_callback = [r_cb](const ts::ProgressInfo& pi) {
      r_cb(Rcpp::List::create(
        Rcpp::Named("replicate") = pi.replicate,
        Rcpp::Named("max_replicates") = pi.max_replicates,
        Rcpp::Named("best_score") = pi.best_score,
        Rcpp::Named("hits_to_best") = pi.hits_to_best,
        Rcpp::Named("target_hits") = pi.target_hits,
        Rcpp::Named("pool_size") = pi.pool_size,
        Rcpp::Named("phase") = std::string(pi.phase),
        Rcpp::Named("elapsed") = pi.elapsed_seconds,
        Rcpp::Named("phase_score") = pi.phase_score
      ));
    };
  }

  return nThreads;
}

// Unpack scoringConfig list and build DataSet.
static ts::DataSet unpack_scoring(
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    List scoring)
{
  IntegerVector min_steps;
  if (scoring.containsElementNamed("min_steps") &&
      !Rf_isNull(scoring["min_steps"])) {
    min_steps = as<IntegerVector>(scoring["min_steps"]);
  }
  double concavity = as<double>(scoring["concavity"]);
  bool xpiwe = as<bool>(scoring["xpiwe"]);
  double xpiwe_r = as<double>(scoring["xpiwe_r"]);
  double xpiwe_max_f = as<double>(scoring["xpiwe_max_f"]);
  IntegerVector obs_count;
  if (scoring.containsElementNamed("obs_count") &&
      !Rf_isNull(scoring["obs_count"])) {
    obs_count = as<IntegerVector>(scoring["obs_count"]);
  }

  // Use SEXP to avoid Nullable assignment ambiguity
  SEXP ia_sexp = R_NilValue;
  if (scoring.containsElementNamed("infoAmounts") &&
      !Rf_isNull(scoring["infoAmounts"])) {
    ia_sexp = scoring["infoAmounts"];
  }
  Nullable<NumericMatrix> infoAmounts(ia_sexp);

  return make_dataset(contrast, tip_data, weight, levels,
                      min_steps, concavity, infoAmounts,
                      xpiwe, xpiwe_r, xpiwe_max_f, obs_count);
}

// Unpack constraintConfig list into ConstraintData.
// Uses SEXP intermediaries to avoid Nullable<T> assignment ambiguity.
static ts::ConstraintData unpack_constraint(int n_tips,
                                            Nullable<List> constraintConfig) {
  if (constraintConfig.isNotNull()) {
    List cc(constraintConfig.get());

    SEXP csm = cc.containsElementNamed("consSplitMatrix")
                   ? SEXP(cc["consSplitMatrix"]) : R_NilValue;
    SEXP cco = cc.containsElementNamed("consContrast")
                   ? SEXP(cc["consContrast"]) : R_NilValue;
    SEXP ctd = cc.containsElementNamed("consTipData")
                   ? SEXP(cc["consTipData"]) : R_NilValue;
    SEXP cw  = cc.containsElementNamed("consWeight")
                   ? SEXP(cc["consWeight"]) : R_NilValue;
    SEXP cl  = cc.containsElementNamed("consLevels")
                   ? SEXP(cc["consLevels"]) : R_NilValue;
    int consExpectedScore = 0;
    if (cc.containsElementNamed("consExpectedScore"))
      consExpectedScore = as<int>(cc["consExpectedScore"]);

    return build_constraint_from_r(
        n_tips,
        Nullable<IntegerMatrix>(csm),
        Nullable<NumericMatrix>(cco),
        Nullable<IntegerMatrix>(ctd),
        Nullable<IntegerVector>(cw),
        Nullable<CharacterVector>(cl),
        consExpectedScore);
  }
  return ts::ConstraintData{};
}

// Unpack hsjConfig list into DataSet fields.
static void unpack_hsj(Nullable<List> hsjConfig, ts::DataSet& ds) {
  if (hsjConfig.isNotNull()) {
    List hc(hsjConfig.get());
    int hsjAbsentState = 0;
    if (hc.containsElementNamed("hsjAbsentState"))
      hsjAbsentState = as<int>(hc["hsjAbsentState"]);

    List hb_list = as<List>(hc["hierarchyBlocks"]);
    for (int b = 0; b < hb_list.size(); ++b) {
      List rb = hb_list[b];
      ts::HierarchyBlock block;
      block.primary_char = as<int>(rb["primary"]);
      block.secondary_chars = as<std::vector<int>>(rb["secondaries"]);
      block.n_secondaries = static_cast<int>(block.secondary_chars.size());
      block.absent_state = hsjAbsentState;
      ds.hierarchy_blocks.push_back(block);
    }
    ds.hsj_alpha = as<double>(hc["hsjAlpha"]);
    ds.scoring_mode = ts::ScoringMode::HSJ;

    if (hc.containsElementNamed("hsjTipLabels") &&
        !Rf_isNull(hc["hsjTipLabels"])) {
      IntegerMatrix tl = as<IntegerMatrix>(hc["hsjTipLabels"]);
      int n_t = tl.nrow();
      int n_c = tl.ncol();
      ds.n_orig_chars = n_c;
      ds.tip_labels.resize(n_t * n_c);
      for (int t = 0; t < n_t; ++t) {
        for (int c = 0; c < n_c; ++c) {
          ds.tip_labels[t * n_c + c] = tl(t, c);
        }
      }
    }
  }
}

// Unpack xformConfig list into DataSet fields.
static void unpack_xform(Nullable<List> xformConfig,
                          IntegerMatrix tip_data,
                          ts::DataSet& ds) {
  if (xformConfig.isNotNull()) {
    List xc(xformConfig.get());
    List xf_list = as<List>(xc["xformChars"]);
    int n_xf = xf_list.size();
    int max_ns = 0;
    std::vector<int> ns_vec(n_xf);
    std::vector<int> fr_vec(n_xf);

    for (int ch = 0; ch < n_xf; ++ch) {
      List rc = xf_list[ch];
      ns_vec[ch] = as<int>(rc["n_states"]);
      fr_vec[ch] = as<int>(rc["forced_root_state"]);
      if (ns_vec[ch] > max_ns) max_ns = ns_vec[ch];
    }

    ds.sankoff_n_chars = n_xf;
    ds.sankoff_max_states = max_ns;
    ds.sankoff_n_states = ns_vec;
    ds.sankoff_forced_root = fr_vec;

    // Cost matrices: flat [n_chars * max_ns * max_ns]
    ds.sankoff_cost_matrices.assign(
        static_cast<size_t>(n_xf) * max_ns * max_ns, 0.0);
    for (int ch = 0; ch < n_xf; ++ch) {
      List rc = xf_list[ch];
      NumericMatrix cm = as<NumericMatrix>(rc["cost_matrix"]);
      int ns = ns_vec[ch];
      double* dst = ds.sankoff_cost_matrices.data() +
          static_cast<size_t>(ch) * max_ns * max_ns;
      for (int r = 0; r < ns; ++r)
        for (int c = 0; c < ns; ++c)
          dst[r * max_ns + c] = cm(r, c);
    }

    // Tip costs: flat [n_tips * stride], stride = n_chars * max_ns
    int n_t = tip_data.nrow();
    int stride = n_xf * max_ns;
    const double INF = std::numeric_limits<double>::infinity();
    ds.sankoff_tip_costs.assign(static_cast<size_t>(n_t) * stride, INF);
    for (int ch = 0; ch < n_xf; ++ch) {
      List rc = xf_list[ch];
      IntegerVector ts_r = as<IntegerVector>(rc["tip_states"]);
      int ns = ns_vec[ch];
      for (int t = 0; t < n_t; ++t) {
        int state = ts_r[t];
        double* tip_ptr = ds.sankoff_tip_costs.data() +
            t * stride + ch * max_ns;
        if (state == -1) {
          for (int s = 0; s < ns; ++s) tip_ptr[s] = 0.0;
        } else if (state == -2) {
          for (int s = 1; s < ns; ++s) tip_ptr[s] = 0.0;
        } else if (state >= 0 && state < ns) {
          tip_ptr[state] = 0.0;
        }
      }
    }

    ds.scoring_mode = ts::ScoringMode::XFORM;
  }
}

// [[Rcpp::export]]
List ts_driven_search(
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    List searchControl,
    List runtimeConfig,
    List scoringConfig,
    Nullable<List> constraintConfig = R_NilValue,
    Nullable<List> hsjConfig = R_NilValue,
    Nullable<List> xformConfig = R_NilValue)
{
  // Build dataset from flat data matrices + scoring config
  ts::DataSet ds = unpack_scoring(contrast, tip_data, weight, levels,
                                  scoringConfig);

  // HSJ hierarchy scoring setup
  unpack_hsj(hsjConfig, ds);

  // Xform (step-matrix) scoring setup
  unpack_xform(xformConfig, tip_data, ds);

  // Build constraint if provided
  int n_tips = tip_data.nrow();
  ts::ConstraintData cd = unpack_constraint(n_tips, constraintConfig);
  ts::ConstraintData* cd_ptr = cd.active ? &cd : nullptr;

  // Populate DrivenParams from searchControl + runtimeConfig
  ts::DrivenParams params;
  unpack_search_control(searchControl, params);
  int nThreads = unpack_runtime(runtimeConfig, params);

  ts::TreePool pool(params.pool_max_size, params.pool_suboptimal);
  ts::DrivenResult result;
  if (nThreads > 1) {
    result = ts::parallel_driven_search(pool, ds, params, cd_ptr, nThreads);
  } else {
    result = ts::driven_search(pool, ds, params, cd_ptr);
  }

  // Build timings as a NumericVector (lighter than List)
  NumericVector timings = NumericVector::create(
    Named("wagner_ms")    = result.timings.wagner_ms,
    Named("nni_ms")       = result.timings.nni_ms,
    Named("tbr_ms")       = result.timings.tbr_ms,
    Named("xss_ms")       = result.timings.xss_ms,
    Named("rss_ms")       = result.timings.rss_ms,
    Named("css_ms")       = result.timings.css_ms,
    Named("ratchet_ms")   = result.timings.ratchet_ms,
    Named("nni_perturb_ms") = result.timings.nni_perturb_ms,
    Named("drift_ms")     = result.timings.drift_ms,
    Named("anneal_ms")    = result.timings.anneal_ms,
    Named("prune_reinsert_ms") = result.timings.prune_reinsert_ms,
    Named("final_tbr_ms") = result.timings.final_tbr_ms,
    Named("fuse_ms")      = result.timings.fuse_ms
  );

  // Per-strategy diagnostics (T-190)
  List strategy_diag = R_NilValue;
  if (params.adaptive_start || nThreads > 1) {
    CharacterVector sn(ts::N_STRAT);
    IntegerVector sa(ts::N_STRAT), ss(ts::N_STRAT);
    for (int i = 0; i < ts::N_STRAT; ++i) {
      sn[i] = ts::strategy_name(static_cast<ts::StartStrategy>(i));
      sa[i] = result.strategy_attempts[i];
      ss[i] = result.strategy_successes[i];
    }
    sa.names() = sn;
    ss.names() = sn;
    strategy_diag = List::create(Named("attempts") = sa, Named("successes") = ss);
  }

  // Per-replicate scores for ScoreSpectrum() coverage estimation
  NumericVector rep_scores(result.replicate_scores.begin(),
                           result.replicate_scores.end());

  if (result.pool_size == 0) {
    return List::create(
      Named("trees") = List::create(),
      Named("scores") = NumericVector::create(),
      Named("best_score") = result.best_score,
      Named("replicates") = result.replicates_completed,
      Named("hits_to_best") = result.hits_to_best,
      Named("pool_size") = 0,
      Named("n_topologies") = 0,
      Named("last_improved_rep") = result.last_improved_rep,
      Named("timed_out") = result.timed_out,
      Named("consensus_stable") = result.consensus_stable,
      Named("perturb_stop") = result.perturb_stop,
      Named("timings") = timings,
      Named("strategy_diagnostics") = strategy_diag,
      Named("replicate_scores") = rep_scores,
      Named("candidates_evaluated") = (double) result.candidates_evaluated
    );
  }

  // Return all pool trees as a list of edge matrices
  const auto& entries = pool.all();
  List tree_list(entries.size());
  NumericVector score_vec(entries.size());
  for (size_t i = 0; i < entries.size(); ++i) {
    tree_list[i] = tree_to_edge(entries[i].tree);
    score_vec[i] = entries[i].score;
  }

  return List::create(
    Named("trees") = tree_list,
    Named("scores") = score_vec,
    Named("best_score") = result.best_score,
    Named("replicates") = result.replicates_completed,
    Named("hits_to_best") = result.hits_to_best,
    Named("pool_size") = result.pool_size,
    Named("n_topologies") = result.n_topologies_at_best,
    Named("last_improved_rep") = result.last_improved_rep,
    Named("timed_out") = result.timed_out,
    Named("consensus_stable") = result.consensus_stable,
    Named("perturb_stop") = result.perturb_stop,
    Named("timings") = timings,
    Named("strategy_diagnostics") = strategy_diag,
    Named("replicate_scores") = rep_scores,
    Named("candidates_evaluated") = (double) result.candidates_evaluated
  );
}

// Collapse a pool of binary trees to TNT-style "zero-length-branch" polytopies
// and deduplicate on the resulting collapsed topologies, entirely in C++.  Used
// by MaximizeParsimony(collapse = TRUE) (the default) so the output path avoids
// per-tree R surgery (RootTree/RenumberTips/CollapseNode/write.tree).
//
// Input `edges` are binary edge matrices straight from tree_to_edge(), so tip i
// already corresponds to tip_data row i — no R rerooting (which would permute
// tips and mis-score against tip_data).  Pass only the best-score trees: a
// collapsed topology has a unique min-resolution length, so the returned set is
// the distinct MPTs at the best score.
//
// Per tree: rebuild it, re-root on tip 0 IN C++ (preserves tip indices), so the
// only root-adjacent edges are trivial splits the collapse criterion skips —
// this is what makes the contraction rooting-invariant.  Then full_rescore,
// flag aggressive (TNT `collapse 3`) min-length-0 internal edges in the search's
// real scoring mode (falls back to the conservative criterion for inapplicable
// data), and dedup on the collapsed split set (canonical bipartitions, so
// rooting-invariant given the shared tip-0 rooting).
//
// Returns: list `trees` of collapsed (multifurcating) edge matrices for the
// distinct topologies, and `n_topologies` (their count).
// [[Rcpp::export]]
List ts_collapse_pool(
    List edges,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    List scoringConfig,
    Nullable<List> hsjConfig = R_NilValue,
    Nullable<List> xformConfig = R_NilValue,
    Nullable<IntegerMatrix> consSplitMatrix = R_NilValue)
{
  ts::DataSet ds = unpack_scoring(contrast, tip_data, weight, levels,
                                  scoringConfig);
  unpack_hsj(hsjConfig, ds);
  unpack_xform(xformConfig, tip_data, ds);

  // "Show the enforced clade": a topological constraint encodes external
  // evidence for a grouping the matrix does not capture, so an enforced clade
  // must stay visible even when its branch is zero-length (it would otherwise be
  // contracted, leaving the result looking unconstrained).  We protect every
  // constraint split from collapse — the unsupported NON-constraint branches
  // still collapse.  Store each constraint split as a canonical (tip-0-excluded)
  // bitset: trees are re-rooted on tip 0 below, so every internal node's
  // descendant set excludes tip 0 and is directly comparable to these.
  const int n_tip = tip_data.nrow();
  const int wps = (n_tip + 63) / 64;
  std::vector<std::vector<uint64_t>> cons_canon;
  if (consSplitMatrix.isNotNull()) {
    IntegerMatrix cs(consSplitMatrix.get());
    for (int r = 0; r < cs.nrow(); ++r) {
      std::vector<uint64_t> b(wps, 0);
      for (int c = 0; c < n_tip && c < cs.ncol(); ++c) {
        if (cs(r, c)) b[c >> 6] |= (1ULL << (c & 63));
      }
      if (b[0] & 1ULL) {                       // canonicalize: exclude tip 0
        for (int w = 0; w < wps; ++w) b[w] = ~b[w];
        int rem = n_tip & 63;
        if (rem) b[wps - 1] &= ((1ULL << rem) - 1);  // clear padding bits
      }
      cons_canon.push_back(std::move(b));
    }
  }

  std::vector<IntegerMatrix> reps;          // representative collapsed edges
  std::vector<uint64_t> rep_hash;           // collapsed-split hash per rep
  std::vector<ts::SplitSet> rep_splits;     // for exact dedup on hash collision
  std::vector<uint8_t> flags;

  for (R_xlen_t t = 0; t < edges.size(); ++t) {
    IntegerMatrix edge = edges[t];
    ts::TreeState tree;
    tree.init_from_edge(&edge(0, 0), &edge(0, 1), edge.nrow(), ds);

    // Root on tip 0 so root-adjacent edges are trivial (rooting-invariant
    // collapse), then refresh state arrays for the flag computation.
    ts::reroot_at_tip(tree, 0);
    tree.reset_states(ds);
    ts::score_tree(tree, ds);

    ts::compute_collapsed_flags_aggressive(tree, ds, flags);

    // Protect constraint splits: clear the collapse flag of any internal edge
    // whose bipartition realises a constraint (keeps the enforced clade
    // visible).  Per-node descendant tip sets via a postorder OR; rooted on
    // tip 0, so every internal set excludes tip 0 == the canonical form above.
    if (!cons_canon.empty()) {
      std::vector<uint64_t> tb(static_cast<size_t>(tree.n_node) * wps, 0);
      for (int tp = 0; tp < n_tip; ++tp) {
        tb[static_cast<size_t>(tp) * wps + (tp >> 6)] = 1ULL << (tp & 63);
      }
      for (size_t pi = 0; pi < tree.postorder.size(); ++pi) {
        int node = tree.postorder[pi];
        if (node < n_tip) continue;
        int ni = node - n_tip;
        uint64_t* dst = &tb[static_cast<size_t>(node) * wps];
        const uint64_t* L = &tb[static_cast<size_t>(tree.left[ni]) * wps];
        const uint64_t* R = &tb[static_cast<size_t>(tree.right[ni]) * wps];
        for (int w = 0; w < wps; ++w) dst[w] = L[w] | R[w];
      }
      for (int v = n_tip + 1; v < tree.n_node; ++v) {
        if (v >= static_cast<int>(flags.size()) || !flags[v]) continue;
        const uint64_t* nb = &tb[static_cast<size_t>(v) * wps];
        for (const auto& cb : cons_canon) {
          bool eq = true;
          for (int w = 0; w < wps; ++w) {
            if (nb[w] != cb[w]) { eq = false; break; }
          }
          if (eq) { flags[v] = 0; break; }
        }
      }
    }

    // Dedup on the collapsed split set (skips the flagged zero-length edges).
    ts::SplitSet css = ts::compute_collapsed_splits(tree, flags);
    uint64_t h = ts::hash_splits(css);
    bool dup = false;
    for (size_t i = 0; i < rep_hash.size(); ++i) {
      if (rep_hash[i] == h && ts::splits_equal(rep_splits[i], css)) {
        dup = true;
        break;
      }
    }
    if (dup) continue;

    rep_hash.push_back(h);
    rep_splits.push_back(std::move(css));
    reps.push_back(tree_to_collapsed_edge(tree, flags));
  }

  List tree_list(reps.size());
  for (size_t i = 0; i < reps.size(); ++i) tree_list[i] = reps[i];
  return List::create(
    Named("trees") = tree_list,
    Named("n_topologies") = static_cast<int>(reps.size())
  );
}

// [[Rcpp::export]]
List ts_resample_search(
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    bool bootstrap = false,
    double jackProportion = 2.0 / 3.0,
    int maxReplicates = 5,
    int targetHits = 2,
    int tbrMaxHits = 1,
    int ratchetCycles = 3,
    double ratchetPerturbProb = 0.04,
    int driftCycles = 0,
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0,
    Nullable<IntegerMatrix> consSplitMatrix = R_NilValue,
    Nullable<NumericMatrix> consContrast = R_NilValue,
    Nullable<IntegerMatrix> consTipData = R_NilValue,
    Nullable<IntegerVector> consWeight = R_NilValue,
    Nullable<CharacterVector> consLevels = R_NilValue,
    int consExpectedScore = 0,
    Nullable<NumericMatrix> infoAmounts = R_NilValue,
    bool xpiwe = false,
    double xpiwe_r = 0.5,
    double xpiwe_max_f = 5.0,
    IntegerVector obs_count = IntegerVector())
{
  if (concavity < 0) concavity = HUGE_VAL;
  int n_tips = tip_data.nrow();
  int n_patterns = tip_data.ncol();
  int n_tokens = contrast.nrow();
  int n_states = contrast.ncol();

  // Validate parallel-vector lengths at the boundary: build_dataset() indexes
  // these by raw pointer over [0, n_patterns)/[0, n_states), so a short vector
  // is an out-of-bounds read. Public wrappers always size them correctly; this
  // turns a direct internal call with mismatched lengths into a clean error
  // rather than undefined behaviour.
  if (weight.size() != n_patterns) {
    Rcpp::stop("`weight` length (%d) must equal the number of characters (%d)",
               static_cast<int>(weight.size()), n_patterns);
  }
  if (levels.size() != n_states) {
    Rcpp::stop("`levels` length (%d) must equal ncol(contrast) (%d)",
               static_cast<int>(levels.size()), n_states);
  }
  if (min_steps.size() > 0 && min_steps.size() != n_patterns) {
    Rcpp::stop("`min_steps` length (%d) must equal the number of characters "
               "(%d)", static_cast<int>(min_steps.size()), n_patterns);
  }
  if (xpiwe && obs_count.size() > 0 && obs_count.size() != n_patterns) {
    Rcpp::stop("`obs_count` length (%d) must equal the number of characters "
               "(%d)", static_cast<int>(obs_count.size()), n_patterns);
  }

  std::vector<std::string> level_strs(n_states);
  std::vector<const char*> level_ptrs(n_states);
  for (int i = 0; i < n_states; ++i) {
    level_strs[i] = as<std::string>(levels[i]);
    level_ptrs[i] = level_strs[i].c_str();
  }

  const int* min_steps_ptr = (min_steps.size() > 0) ? INTEGER(min_steps)
                                                     : nullptr;
  const int* obs_count_ptr = (xpiwe && obs_count.size() > 0)
                                 ? INTEGER(obs_count)
                                 : nullptr;

  // Profile parsimony
  const double* info_amounts_ptr;
  int info_max_steps;
  extract_info_amounts(infoAmounts, info_amounts_ptr, info_max_steps);

  // Constraints
  ts::ConstraintData cd = build_constraint_from_r(
      n_tips, consSplitMatrix, consContrast, consTipData,
      consWeight, consLevels, consExpectedScore);
  ts::ConstraintData* cd_ptr = cd.active ? &cd : nullptr;

  ts::ResampleParams params;
  params.bootstrap = bootstrap;
  params.jack_proportion = jackProportion;
  params.search.max_replicates = maxReplicates;
  params.search.target_hits = targetHits;
  params.search.tbr_max_hits = tbrMaxHits;
  params.search.ratchet_cycles = ratchetCycles;
  params.search.ratchet_perturb_prob = ratchetPerturbProb;
  params.search.drift_cycles = driftCycles;

  ts::ResampleResult result = ts::resample_search(
      REAL(contrast), n_tokens, n_states,
      INTEGER(tip_data), n_tips, n_patterns,
      INTEGER(weight),
      level_ptrs.data(),
      min_steps_ptr,
      concavity,
      params,
      info_amounts_ptr,
      info_max_steps,
      cd_ptr,
      xpiwe,
      xpiwe_r,
      xpiwe_max_f,
      obs_count_ptr);

  if (result.edge_parent.empty()) {
    return List::create(
      Named("edge") = IntegerMatrix(0, 2),
      Named("score") = result.score
    );
  }

  int n_edge = static_cast<int>(result.edge_parent.size());
  IntegerMatrix edge(n_edge, 2);
  for (int i = 0; i < n_edge; ++i) {
    edge(i, 0) = result.edge_parent[i];
    edge(i, 1) = result.edge_child[i];
  }

  return List::create(
    Named("edge") = edge,
    Named("score") = result.score
  );
}

// [[Rcpp::export]]
List ts_parallel_resample(
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int nReplicates = 1,
    int nThreads = 1,
    bool bootstrap = false,
    double jackProportion = 2.0 / 3.0,
    int maxReplicates = 5,
    int targetHits = 2,
    int tbrMaxHits = 1,
    int ratchetCycles = 3,
    double ratchetPerturbProb = 0.04,
    int driftCycles = 0,
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0,
    Nullable<IntegerMatrix> consSplitMatrix = R_NilValue,
    Nullable<NumericMatrix> consContrast = R_NilValue,
    Nullable<IntegerMatrix> consTipData = R_NilValue,
    Nullable<IntegerVector> consWeight = R_NilValue,
    Nullable<CharacterVector> consLevels = R_NilValue,
    int consExpectedScore = 0,
    Nullable<NumericMatrix> infoAmounts = R_NilValue,
    bool xpiwe = false,
    double xpiwe_r = 0.5,
    double xpiwe_max_f = 5.0,
    IntegerVector obs_count = IntegerVector())
{
  if (concavity < 0) concavity = HUGE_VAL;
  int n_tips = tip_data.nrow();
  int n_patterns = tip_data.ncol();
  int n_tokens = contrast.nrow();
  int n_states = contrast.ncol();

  // Validate parallel-vector lengths at the boundary: build_dataset() indexes
  // these by raw pointer over [0, n_patterns)/[0, n_states), so a short vector
  // is an out-of-bounds read. Public wrappers always size them correctly; this
  // turns a direct internal call with mismatched lengths into a clean error
  // rather than undefined behaviour.
  if (weight.size() != n_patterns) {
    Rcpp::stop("`weight` length (%d) must equal the number of characters (%d)",
               static_cast<int>(weight.size()), n_patterns);
  }
  if (levels.size() != n_states) {
    Rcpp::stop("`levels` length (%d) must equal ncol(contrast) (%d)",
               static_cast<int>(levels.size()), n_states);
  }
  if (min_steps.size() > 0 && min_steps.size() != n_patterns) {
    Rcpp::stop("`min_steps` length (%d) must equal the number of characters "
               "(%d)", static_cast<int>(min_steps.size()), n_patterns);
  }
  if (xpiwe && obs_count.size() > 0 && obs_count.size() != n_patterns) {
    Rcpp::stop("`obs_count` length (%d) must equal the number of characters "
               "(%d)", static_cast<int>(obs_count.size()), n_patterns);
  }

  std::vector<std::string> level_strs(n_states);
  std::vector<const char*> level_ptrs(n_states);
  for (int i = 0; i < n_states; ++i) {
    level_strs[i] = as<std::string>(levels[i]);
    level_ptrs[i] = level_strs[i].c_str();
  }

  const int* min_steps_ptr = (min_steps.size() > 0) ? INTEGER(min_steps)
                                                     : nullptr;
  const int* obs_count_ptr = (xpiwe && obs_count.size() > 0)
                                 ? INTEGER(obs_count)
                                 : nullptr;

  const double* info_amounts_ptr;
  int info_max_steps;
  extract_info_amounts(infoAmounts, info_amounts_ptr, info_max_steps);

  ts::ConstraintData cd = build_constraint_from_r(
      n_tips, consSplitMatrix, consContrast, consTipData,
      consWeight, consLevels, consExpectedScore);
  ts::ConstraintData* cd_ptr = cd.active ? &cd : nullptr;

  ts::ResampleParams params;
  params.bootstrap = bootstrap;
  params.jack_proportion = jackProportion;
  params.search.max_replicates = maxReplicates;
  params.search.target_hits = targetHits;
  params.search.tbr_max_hits = tbrMaxHits;
  params.search.ratchet_cycles = ratchetCycles;
  params.search.ratchet_perturb_prob = ratchetPerturbProb;
  params.search.drift_cycles = driftCycles;

  if (nReplicates < 1) nReplicates = 1;

  std::vector<ts::ResampleResult> results;
  if (nThreads > 1 && nReplicates > 1) {
    results = ts::parallel_resample(
        REAL(contrast), n_tokens, n_states,
        INTEGER(tip_data), n_tips, n_patterns,
        INTEGER(weight), level_ptrs.data(), min_steps_ptr,
        concavity, params, nReplicates, nThreads,
        info_amounts_ptr, info_max_steps, cd_ptr,
        xpiwe, xpiwe_r, xpiwe_max_f, obs_count_ptr);
  } else {
    // Serial path: run each replicate sequentially
    results.resize(nReplicates);
    for (int r = 0; r < nReplicates; ++r) {
      results[r] = ts::resample_search(
          REAL(contrast), n_tokens, n_states,
          INTEGER(tip_data), n_tips, n_patterns,
          INTEGER(weight), level_ptrs.data(), min_steps_ptr,
          concavity, params,
          info_amounts_ptr, info_max_steps, cd_ptr,
          xpiwe, xpiwe_r, xpiwe_max_f, obs_count_ptr);
    }
  }

  // Package results as list of edge matrices + score vector
  List edges(nReplicates);
  NumericVector scores(nReplicates);
  for (int r = 0; r < nReplicates; ++r) {
    const auto& res = results[r];
    scores[r] = res.score;
    if (res.edge_parent.empty()) {
      edges[r] = IntegerMatrix(0, 2);
    } else {
      int n_edge = static_cast<int>(res.edge_parent.size());
      IntegerMatrix em(n_edge, 2);
      for (int i = 0; i < n_edge; ++i) {
        em(i, 0) = res.edge_parent[i];
        em(i, 1) = res.edge_child[i];
      }
      edges[r] = em;
    }
  }

  return List::create(
    Named("edges") = edges,
    Named("scores") = scores,
    Named("n_replicates") = nReplicates
  );
}

// [[Rcpp::export]]
List ts_successive_approx(
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    double saK = 3.0,
    int maxSAIter = 20,
    int maxReplicates = 10,
    int targetHits = 3,
    int tbrMaxHits = 1,
    int ratchetCycles = 5,
    double ratchetPerturbProb = 0.04,
    int driftCycles = 0,
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0,
    Nullable<IntegerMatrix> consSplitMatrix = R_NilValue,
    Nullable<NumericMatrix> consContrast = R_NilValue,
    Nullable<IntegerMatrix> consTipData = R_NilValue,
    Nullable<IntegerVector> consWeight = R_NilValue,
    Nullable<CharacterVector> consLevels = R_NilValue,
    int consExpectedScore = 0,
    Nullable<NumericMatrix> infoAmounts = R_NilValue,
    bool xpiwe = false,
    double xpiwe_r = 0.5,
    double xpiwe_max_f = 5.0,
    IntegerVector obs_count = IntegerVector())
{
  if (concavity < 0) concavity = HUGE_VAL;
  int n_tips = tip_data.nrow();
  int n_patterns = tip_data.ncol();
  int n_tokens = contrast.nrow();
  int n_states = contrast.ncol();

  // Validate parallel-vector lengths at the boundary: build_dataset() indexes
  // these by raw pointer over [0, n_patterns)/[0, n_states), so a short vector
  // is an out-of-bounds read. Public wrappers always size them correctly; this
  // turns a direct internal call with mismatched lengths into a clean error
  // rather than undefined behaviour.
  if (weight.size() != n_patterns) {
    Rcpp::stop("`weight` length (%d) must equal the number of characters (%d)",
               static_cast<int>(weight.size()), n_patterns);
  }
  if (levels.size() != n_states) {
    Rcpp::stop("`levels` length (%d) must equal ncol(contrast) (%d)",
               static_cast<int>(levels.size()), n_states);
  }
  if (min_steps.size() > 0 && min_steps.size() != n_patterns) {
    Rcpp::stop("`min_steps` length (%d) must equal the number of characters "
               "(%d)", static_cast<int>(min_steps.size()), n_patterns);
  }
  if (xpiwe && obs_count.size() > 0 && obs_count.size() != n_patterns) {
    Rcpp::stop("`obs_count` length (%d) must equal the number of characters "
               "(%d)", static_cast<int>(obs_count.size()), n_patterns);
  }

  std::vector<std::string> level_strs(n_states);
  std::vector<const char*> level_ptrs(n_states);
  for (int i = 0; i < n_states; ++i) {
    level_strs[i] = as<std::string>(levels[i]);
    level_ptrs[i] = level_strs[i].c_str();
  }

  const int* min_steps_ptr = (min_steps.size() > 0) ? INTEGER(min_steps)
                                                     : nullptr;
  const int* obs_count_ptr = (xpiwe && obs_count.size() > 0)
                                 ? INTEGER(obs_count)
                                 : nullptr;

  // Profile parsimony
  const double* info_amounts_ptr;
  int info_max_steps;
  extract_info_amounts(infoAmounts, info_amounts_ptr, info_max_steps);

  // Constraints
  ts::ConstraintData cd = build_constraint_from_r(
      n_tips, consSplitMatrix, consContrast, consTipData,
      consWeight, consLevels, consExpectedScore);
  ts::ConstraintData* cd_ptr = cd.active ? &cd : nullptr;

  ts::SAParams params;
  params.k = saK;
  params.max_sa_iter = maxSAIter;
  params.search.max_replicates = maxReplicates;
  params.search.target_hits = targetHits;
  params.search.tbr_max_hits = tbrMaxHits;
  params.search.ratchet_cycles = ratchetCycles;
  params.search.ratchet_perturb_prob = ratchetPerturbProb;
  params.search.drift_cycles = driftCycles;

  ts::SAResult result = ts::successive_approximations(
      REAL(contrast), n_tokens, n_states,
      INTEGER(tip_data), n_tips, n_patterns,
      INTEGER(weight),
      level_ptrs.data(),
      min_steps_ptr,
      concavity,
      params,
      info_amounts_ptr,
      info_max_steps,
      cd_ptr,
      xpiwe,
      xpiwe_r,
      xpiwe_max_f,
      obs_count_ptr);

  if (result.edge_parent.empty()) {
    return List::create(
      Named("edge") = IntegerMatrix(0, 2),
      Named("score") = result.score,
      Named("sa_iterations") = result.sa_iterations,
      Named("converged") = result.converged
    );
  }

  int n_edge = static_cast<int>(result.edge_parent.size());
  IntegerMatrix edge(n_edge, 2);
  for (int i = 0; i < n_edge; ++i) {
    edge(i, 0) = result.edge_parent[i];
    edge(i, 1) = result.edge_child[i];
  }

  return List::create(
    Named("edge") = edge,
    Named("score") = result.score,
    Named("sa_iterations") = result.sa_iterations,
    Named("converged") = result.converged
  );
}

// --- Phase 3D: Benchmarking diagnostics ---

// [[Rcpp::export]]
List ts_bench_tbr_phases(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0)
{
  using Clock = std::chrono::high_resolution_clock;
  using Us = std::chrono::microseconds;

  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(&edge(0, 0), &edge(0, 1), edge.nrow(), ds);

  // --- Phase A: Full rescore (baseline) ---
  auto t0 = Clock::now();
  tree.reset_states(ds);
  double score = ts::score_tree(tree, ds);
  auto t1 = Clock::now();
  double time_full_rescore_us =
      std::chrono::duration_cast<Us>(t1 - t0).count();

  // Check for NA blocks
  bool has_na = false;
  for (int b = 0; b < ds.n_blocks; ++b) {
    if (ds.blocks[b].has_inapplicable) { has_na = true; break; }
  }
  bool use_iw = std::isfinite(ds.concavity);

  // Seed RNG
  std::mt19937 rng = ts::make_rng();

  // Build clip candidates (non-root)
  std::vector<int> clip_candidates;
  for (int node = 0; node < tree.n_node; ++node) {
    if (node == tree.n_tip) continue;
    clip_candidates.push_back(node);
  }
  std::shuffle(clip_candidates.begin(), clip_candidates.end(), rng);

  // Subtree sizes for filtering
  std::vector<int> subtree_sizes(tree.n_node, 0);
  {
    // Compute subtree sizes (same as TBR search)
    for (int i = 0; i < tree.n_tip; ++i) subtree_sizes[i] = 1;
    for (int node : tree.postorder) {
      int ni = node - tree.n_tip;
      subtree_sizes[node] = subtree_sizes[tree.left[ni]]
                           + subtree_sizes[tree.right[ni]];
    }
  }

  // Timing accumulators (microseconds)
  double time_clip_incr_us = 0;
  double time_indirect_us = 0;
  double time_unclip_us = 0;
  int n_clips = 0;
  int n_candidates_total = 0;

  // Snapshot timing
  double time_snapshot_save_us = 0;
  double time_snapshot_restore_us = 0;
  int n_snapshot_ops = 0;

  // IW buffers
  std::vector<int> divided_steps;
  std::vector<double> iw_delta;
  if (use_iw) {
    divided_steps.resize(ds.n_patterns, 0);
    iw_delta.resize(ds.n_patterns, 0.0);
  }

  // Process all clip candidates (one pass, no moves applied)
  std::vector<std::pair<int,int>> main_edges;
  std::vector<std::pair<int,int>> sub_edges;
  std::vector<uint64_t> from_above(
      static_cast<size_t>(tree.n_node) * tree.total_words, 0);
  std::vector<uint64_t> virtual_prelim(tree.total_words);
  std::vector<uint64_t> vroot_cache;

  for (int clip_node : clip_candidates) {
    if (tree.parent[clip_node] == tree.n_tip) continue;
    int clip_size = subtree_sizes[clip_node];
    if (clip_size > tree.n_tip / 2) continue;

    // Save clip actives for NA
    std::vector<uint64_t> clip_actives_buf;
    const uint64_t* clip_actives = nullptr;
    if (has_na) {
      size_t clip_sa_base =
          static_cast<size_t>(clip_node) * tree.total_words;
      clip_actives_buf.assign(
          tree.subtree_actives.begin() + clip_sa_base,
          tree.subtree_actives.begin() + clip_sa_base + tree.total_words);
      clip_actives = clip_actives_buf.data();
    }

    // --- Time clip + incremental scoring ---
    auto tc0 = Clock::now();

    tree.spr_clip(clip_node);
    tree.build_postorder();

    int nz = tree.clip_state.clip_grandpar;

    if (has_na) {
      ts::fitch_na_incremental_downpass(tree, ds, nz);
      ts::fitch_na_incremental_uppass(tree, ds, nz);
      (void)ts::fitch_na_pass3_score(tree, ds);
    } else {
      ts::fitch_incremental_downpass(tree, ds, nz);
      ts::fitch_incremental_uppass(tree, ds, nz);
    }

    // IW base
    double base_iw = 0.0;
    if (use_iw) {
      ts::extract_char_steps(tree, ds, divided_steps);
      base_iw = ts::compute_weighted_score(ds, divided_steps);
      ts::precompute_weighted_delta(ds, divided_steps, iw_delta);
    }

    auto tc1 = Clock::now();
    time_clip_incr_us +=
        std::chrono::duration_cast<Us>(tc1 - tc0).count();

    // --- Time indirect evaluation ---
    auto ti0 = Clock::now();

    // Collect main edges
    main_edges.clear();
    for (int node : tree.postorder) {
      int ni = node - tree.n_tip;
      main_edges.push_back({node, tree.left[ni]});
      main_edges.push_back({node, tree.right[ni]});
    }

    int ns = tree.clip_state.clip_sibling;
    size_t clip_base = static_cast<size_t>(clip_node) * tree.total_words;
    const uint64_t* clip_prelim = &tree.prelim[clip_base];
    int n_spr_candidates = 0;

    // SPR candidates (unbounded for timing)
    for (auto& [above, below] : main_edges) {
      if (above == nz && below == ns) continue;
      if (has_na) {
        if (use_iw) {
          ts::indirect_na_iw_length_bounded(clip_prelim, clip_actives,
              tree, ds, above, below, base_iw, iw_delta, HUGE_VAL);
        } else {
          ts::fitch_na_indirect_length_bounded(clip_prelim, clip_actives,
              tree, ds, above, below, INT_MAX);
        }
      } else if (use_iw) {
        ts::indirect_iw_length_bounded(clip_prelim, tree, ds, above, below,
            base_iw, iw_delta, HUGE_VAL);
      } else {
        ts::fitch_indirect_length_bounded(clip_prelim, tree, ds,
            above, below, INT_MAX);
      }
      ++n_spr_candidates;
    }

    // TBR candidates (with vroot cache)
    int n_tbr_candidates = 0;
    if (clip_node >= tree.n_tip) {
      // Precompute vroot cache
      int n_main = static_cast<int>(main_edges.size());
      vroot_cache.resize(static_cast<size_t>(n_main) * tree.total_words);
      for (int ei = 0; ei < n_main; ++ei) {
        int a = main_edges[ei].first;
        int d = main_edges[ei].second;
        size_t a_base = static_cast<size_t>(a) * tree.total_words;
        size_t d_base = static_cast<size_t>(d) * tree.total_words;
        size_t out_base = static_cast<size_t>(ei) * tree.total_words;
        for (int s = 0; s < tree.total_words; ++s) {
          vroot_cache[out_base + s] = tree.final_[a_base + s]
                                    | tree.final_[d_base + s];
        }
      }

      // Collect subtree edges
      sub_edges.clear();
      std::vector<int> sub_stack;
      sub_stack.push_back(clip_node);
      while (!sub_stack.empty()) {
        int nd = sub_stack.back();
        sub_stack.pop_back();
        if (nd >= tree.n_tip) {
          int ni = nd - tree.n_tip;
          int lc = tree.left[ni], rc = tree.right[ni];
          if (lc >= 0) { sub_edges.push_back({nd, lc}); sub_stack.push_back(lc); }
          if (rc >= 0) { sub_edges.push_back({nd, rc}); sub_stack.push_back(rc); }
        }
      }

      // Compute from_above for subtree
      // (simplified: use final_ of clip_node's parent as pseudo-above)
      // For correct timing we don't need exact from_above, just measure
      // the iteration cost. Use clip_prelim as virtual_prelim placeholder.
      for (auto& [sp, sc] : sub_edges) {
        if (sp == clip_node) continue;
        // Quick virtual prelim (just use prelim as-is for timing)
        size_t sc_base = static_cast<size_t>(sc) * tree.total_words;
        const uint64_t* vp = &tree.prelim[sc_base];

        for (int ei = 0; ei < n_main; ++ei) {
          auto& [above, below] = main_edges[ei];
          if (above == nz && below == ns) continue;
          if (use_iw) {
            ts::indirect_iw_length_cached(
                vp, &vroot_cache[static_cast<size_t>(ei) * tree.total_words],
                ds, base_iw, iw_delta, HUGE_VAL);
          } else {
            ts::fitch_indirect_length_cached(
                vp, &vroot_cache[static_cast<size_t>(ei) * tree.total_words],
                ds, INT_MAX);
          }
          ++n_tbr_candidates;
        }
      }
    }

    auto ti1 = Clock::now();
    time_indirect_us +=
        std::chrono::duration_cast<Us>(ti1 - ti0).count();
    n_candidates_total += n_spr_candidates + n_tbr_candidates;

    // --- Time unclip ---
    auto tu0 = Clock::now();
    tree.spr_unclip();
    tree.build_postorder();
    auto tu1 = Clock::now();
    time_unclip_us +=
        std::chrono::duration_cast<Us>(tu1 - tu0).count();

    ++n_clips;

    ts::check_interrupt();
  }

  // --- Time snapshot save/restore (separate measurement) ---
  {
    // Allocate snapshot
    struct SnapBench {
      std::vector<uint64_t> prelim, final_, local_cost, down2, sub_act;
      void alloc(const ts::TreeState& t, bool na) {
        size_t ssz = static_cast<size_t>(t.n_node) * t.total_words;
        size_t csz = static_cast<size_t>(t.n_node) * t.n_blocks;
        prelim.resize(ssz); final_.resize(ssz); local_cost.resize(csz);
        if (na) { down2.resize(ssz); sub_act.resize(ssz); }
      }
    } snap_bench;
    snap_bench.alloc(tree, has_na);

    size_t state_bytes = snap_bench.prelim.size() * sizeof(uint64_t);
    size_t cost_bytes = snap_bench.local_cost.size() * sizeof(uint64_t);

    // Warm up
    std::memcpy(snap_bench.prelim.data(), tree.prelim.data(), state_bytes);
    std::memcpy(snap_bench.final_.data(), tree.final_.data(), state_bytes);

    int n_snap_iters = std::max(100, 10000 / std::max(1, tree.n_node));
    auto ts0 = Clock::now();
    for (int i = 0; i < n_snap_iters; ++i) {
      std::memcpy(snap_bench.prelim.data(), tree.prelim.data(), state_bytes);
      std::memcpy(snap_bench.final_.data(), tree.final_.data(), state_bytes);
      std::memcpy(snap_bench.local_cost.data(), tree.local_cost.data(), cost_bytes);
      if (has_na) {
        std::memcpy(snap_bench.down2.data(), tree.down2.data(), state_bytes);
        std::memcpy(snap_bench.sub_act.data(), tree.subtree_actives.data(), state_bytes);
      }
    }
    auto ts1 = Clock::now();
    time_snapshot_save_us =
        static_cast<double>(std::chrono::duration_cast<Us>(ts1 - ts0).count())
        / n_snap_iters;

    auto tr0 = Clock::now();
    for (int i = 0; i < n_snap_iters; ++i) {
      std::memcpy(tree.prelim.data(), snap_bench.prelim.data(), state_bytes);
      std::memcpy(tree.final_.data(), snap_bench.final_.data(), state_bytes);
      std::memcpy(tree.local_cost.data(), snap_bench.local_cost.data(), cost_bytes);
      if (has_na) {
        std::memcpy(tree.down2.data(), snap_bench.down2.data(), state_bytes);
        std::memcpy(tree.subtree_actives.data(), snap_bench.sub_act.data(), state_bytes);
      }
    }
    auto tr1 = Clock::now();
    time_snapshot_restore_us =
        static_cast<double>(std::chrono::duration_cast<Us>(tr1 - tr0).count())
        / n_snap_iters;
    n_snapshot_ops = n_snap_iters;
  }

  // Dataset info
  int total_chars = 0;
  std::vector<int> block_n_states(ds.n_blocks);
  for (int b = 0; b < ds.n_blocks; ++b) {
    total_chars += ds.blocks[b].n_chars;
    block_n_states[b] = ds.blocks[b].n_states;
  }

  // Snapshot size in bytes
  size_t snap_bytes = static_cast<size_t>(tree.n_node) * tree.total_words
                    * sizeof(uint64_t);
  size_t snap_total = snap_bytes * 2 + // prelim + final_
      static_cast<size_t>(tree.n_node) * tree.n_blocks * sizeof(uint64_t);
  if (has_na) snap_total += snap_bytes * 2; // down2 + subtree_actives

  return List::create(
    Named("n_tips") = tree.n_tip,
    Named("n_node") = tree.n_node,
    Named("n_blocks") = ds.n_blocks,
    Named("total_words") = tree.total_words,
    Named("total_chars") = total_chars,
    Named("block_n_states") = wrap(block_n_states),
    Named("has_na") = has_na,
    Named("use_iw") = use_iw,
    Named("score") = score,
    // Timing (microseconds)
    Named("time_full_rescore_us") = time_full_rescore_us,
    Named("time_clip_incr_us") = time_clip_incr_us,
    Named("time_indirect_us") = time_indirect_us,
    Named("time_unclip_us") = time_unclip_us,
    Named("time_snapshot_save_us") = time_snapshot_save_us,
    Named("time_snapshot_restore_us") = time_snapshot_restore_us,
    Named("snapshot_bytes") = static_cast<double>(snap_total),
    // Counts
    Named("n_clips") = n_clips,
    Named("n_candidates") = n_candidates_total,
    Named("n_snapshot_iters") = n_snapshot_ops
  );
}

// [[Rcpp::export]]
List ts_simplify_diag(
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels)
{
  int n_tokens = contrast.nrow();
  int n_states = contrast.ncol();
  int n_tips = tip_data.nrow();
  int n_patterns = tip_data.ncol();

  // Identify inapp_state
  int inapp_state = -1;
  for (int s = 0; s < n_states; ++s) {
    if (std::string(levels[s]) == "-") {
      inapp_state = s;
      break;
    }
  }

  // Build token_states
  std::vector<uint32_t> token_states(n_tokens, 0);
  for (int t = 0; t < n_tokens; ++t) {
    for (int s = 0; s < n_states; ++s) {
      if (contrast(t, s) > 0.5) {
        token_states[t] |= (1u << s);
      }
    }
  }

  ts::SimplificationResult simpl = ts::simplify_patterns(
      token_states, &tip_data(0, 0), n_tips, n_patterns,
      &weight[0], n_states, inapp_state);

  // Build return: per-pattern info
  IntegerVector precomputed(n_patterns);
  LogicalVector informative(n_patterns);
  IntegerVector n_states_remaining(n_patterns);
  for (int p = 0; p < n_patterns; ++p) {
    precomputed[p] = simpl.patterns[p].precomputed_steps;
    informative[p] = simpl.patterns[p].informative;
    n_states_remaining[p] = simpl.patterns[p].n_states_remaining;
  }

  // Compute ew_offset
  int ew_offset = 0;
  for (int p = 0; p < n_patterns; ++p) {
    ew_offset += simpl.patterns[p].precomputed_steps * weight[p];
  }

  return List::create(
    Named("n_patterns_removed") = simpl.n_patterns_removed,
    Named("n_states_reduced") = simpl.n_states_reduced,
    Named("total_offset_steps") = simpl.total_offset_steps,
    Named("ew_offset") = ew_offset,
    Named("precomputed_steps") = precomputed,
    Named("informative") = informative,
    Named("n_states_remaining") = n_states_remaining
  );
}

// [[Rcpp::export]]
double ts_hsj_score(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    List hierarchy_blocks_r,
    double alpha,
    IntegerMatrix tip_labels_r,
    int absent_state)
{
  // Build DataSet for non-hierarchy characters (weight already adjusted)
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);

  // Build tree
  ts::TreeState tree;
  tree.init_from_edge(&edge(0, 0), &edge(0, 1), edge.nrow(), ds);

  // Convert hierarchy blocks from R list to C++ vector
  std::vector<ts::HierarchyBlock> blocks;
  for (int b = 0; b < hierarchy_blocks_r.size(); ++b) {
    List rb = hierarchy_blocks_r[b];
    ts::HierarchyBlock block;
    block.primary_char = as<int>(rb["primary"]);     // 0-based
    block.secondary_chars = as<std::vector<int>>(rb["secondaries"]); // 0-based
    block.n_secondaries = static_cast<int>(block.secondary_chars.size());
    block.absent_state = absent_state;
    blocks.push_back(block);
  }

  // Flatten tip_labels matrix (n_tips x n_orig_chars, row-major)
  int n_tips = tip_labels_r.nrow();
  int n_orig_chars = tip_labels_r.ncol();
  std::vector<int> tip_labels(n_tips * n_orig_chars);
  for (int t = 0; t < n_tips; ++t) {
    for (int c = 0; c < n_orig_chars; ++c) {
      tip_labels[t * n_orig_chars + c] = tip_labels_r(t, c);
    }
  }

  return ts::hsj_score(tree, ds, blocks, alpha, tip_labels, n_orig_chars);
}

// =========================================================================
// Sankoff parsimony scoring — test bridge
// =========================================================================

#include "ts_sankoff.h"

// Build tree topology from R edge matrix (no DataSet needed).
static void build_topo_from_edge(
    const int* edge_parent, const int* edge_child, int n_edge,
    int n_tip,
    std::vector<int>& left_out, std::vector<int>& right_out,
    std::vector<int>& postorder_out)
{
  int n_internal = n_tip - 1;
  left_out.assign(n_internal, -1);
  right_out.assign(n_internal, -1);

  for (int i = 0; i < n_edge; ++i) {
    int p = edge_parent[i] - 1;
    int c = edge_child[i] - 1;
    int pi = p - n_tip;
    if (left_out[pi] == -1) left_out[pi] = c;
    else                     right_out[pi] = c;
  }

  // Two-stack postorder (internal nodes only, leaves-to-root)
  postorder_out.clear();
  postorder_out.reserve(n_internal);
  std::vector<int> stk;
  stk.push_back(n_tip);
  while (!stk.empty()) {
    int nd = stk.back(); stk.pop_back();
    if (nd >= n_tip) {
      postorder_out.push_back(nd);
      int ni = nd - n_tip;
      if (left_out[ni]  >= 0) stk.push_back(left_out[ni]);
      if (right_out[ni] >= 0) stk.push_back(right_out[ni]);
    }
  }
  std::reverse(postorder_out.begin(), postorder_out.end());
}

// [[Rcpp::export]]
List ts_sankoff_test(
    IntegerMatrix edge,
    IntegerVector n_states_r,
    List cost_matrices_r,
    IntegerMatrix tip_states_r,
    IntegerVector forced_root_r)
{
  int n_edge = edge.nrow();
  int n_tip  = (n_edge / 2) + 1;
  int n_chars = n_states_r.size();
  int n_node = 2 * n_tip - 1;
  const double INF = std::numeric_limits<double>::infinity();

  // Build topology
  std::vector<int> left_v, right_v, postorder;
  build_topo_from_edge(&edge(0, 0), &edge(0, 1), n_edge, n_tip,
                       left_v, right_v, postorder);
  int n_internal = static_cast<int>(postorder.size());

  // Build SankoffData
  ts::SankoffData sd;
  sd.n_tips = n_tip;
  sd.n_chars = n_chars;
  sd.max_states = 0;
  sd.chars.resize(n_chars);

  for (int ch = 0; ch < n_chars; ++ch) {
    int ns = n_states_r[ch];
    sd.chars[ch].n_states = ns;
    sd.chars[ch].forced_root_state = forced_root_r[ch];
    if (ns > sd.max_states) sd.max_states = ns;

    NumericMatrix cm = as<NumericMatrix>(cost_matrices_r[ch]);
    sd.chars[ch].cost_matrix.resize(ns * ns);
    for (int r = 0; r < ns; ++r)
      for (int c = 0; c < ns; ++c)
        sd.chars[ch].cost_matrix[r * ns + c] = cm(r, c);
  }

  // Build tip costs
  int stride = sd.stride();
  sd.tip_costs.assign(static_cast<size_t>(n_tip) * stride, INF);
  for (int t = 0; t < n_tip; ++t) {
    for (int ch = 0; ch < n_chars; ++ch) {
      int state = tip_states_r(t, ch);
      const int ns_ch = sd.chars[ch].n_states;
      double* tip_ptr = sd.tip_costs.data() +
          static_cast<size_t>(t) * stride + static_cast<size_t>(ch) * sd.max_states;
      // Mirror the live xform path (ts_driven_search): -1 = fully ambiguous
      // ("?" in a controlling character), -2 = present but in an unknown
      // secondary combination (any present state). Previously these sentinels
      // were skipped, leaving every state at INF, so any "?" inflated the
      // hierarchy score to Inf.
      if (state == -1) {
        for (int s = 0; s < ns_ch; ++s) tip_ptr[s] = 0.0;
      } else if (state == -2) {
        for (int s = 1; s < ns_ch; ++s) tip_ptr[s] = 0.0;
      } else if (state >= 0 && state < ns_ch) {
        tip_ptr[state] = 0.0;
      }
    }
  }

  // Score all characters
  double total = ts::sankoff_score(
      left_v.data(), right_v.data(),
      postorder.data(), n_internal, n_tip, sd);

  // Per-character scores + node costs + uppass
  NumericVector per_char(n_chars);
  IntegerMatrix opt_states(n_node, n_chars);

  for (int ch = 0; ch < n_chars; ++ch) {
    const ts::SankoffChar& sc = sd.chars[ch];
    const double* ch_tip = sd.tip_costs.data() + ch * sd.max_states;

    std::vector<double> node_costs(static_cast<size_t>(n_node) * sc.n_states);
    per_char[ch] = ts::sankoff_score_char(
        left_v.data(), right_v.data(),
        postorder.data(), n_internal, n_tip,
        sc, ch_tip, stride, node_costs.data());

    std::vector<int> opt(n_node);
    ts::sankoff_uppass(
        left_v.data(), right_v.data(),
        postorder.data(), n_internal, n_tip,
        sc, node_costs.data(), opt.data());

    for (int nd = 0; nd < n_node; ++nd)
      opt_states(nd, ch) = opt[nd];
  }

  return List::create(
    Named("score") = total,
    Named("per_char") = per_char,
    Named("optimal_states") = opt_states);
}


// --- Wagner bias benchmark ---
//
// For each of n_reps random seeds, builds a Wagner tree under the specified
// biasing criterion and optionally runs TBR to the local optimum.  Returns
// per-replicate Wagner scores (and TBR scores if run_tbr = TRUE) so that
// callers can compare average starting-tree quality across criteria.
//
// bias:        0 = RANDOM, 1 = GOLOBOFF, 2 = ENTROPY
// temperature: softmax temperature (0 = greedy; applied to [0,1]-normalised
//              scores so the parameter is dataset-independent)
// n_reps:      number of trees to build
// run_tbr:     if TRUE, run TBR convergence and record its score too

// [[Rcpp::export]]
List ts_wagner_bias_bench(
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    IntegerVector min_steps,
    double concavity,
    int    bias,
    double temperature,
    int    n_reps,
    bool   run_tbr)
{
  if (concavity < 0) concavity = HUGE_VAL;
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::BiasedWagnerParams wp;
  wp.bias        = static_cast<ts::WagnerBias>(bias);
  wp.temperature = temperature;

  NumericVector wagner_scores(n_reps, NA_REAL);
  NumericVector tbr_scores(n_reps, NA_REAL);
  // Per-tip Goloboff and entropy scores (computed once)
  NumericVector goloboff_scores_r(ds.n_tips, NA_REAL);
  NumericVector entropy_scores_r(ds.n_tips, NA_REAL);
  {
    auto gs = ts::wagner_goloboff_scores(ds);
    auto es = ts::wagner_entropy_scores(ds);
    for (int t = 0; t < ds.n_tips; ++t) {
      goloboff_scores_r[t] = gs[t];
      entropy_scores_r[t]  = es[t];
    }
  }

  for (int rep = 0; rep < n_reps; ++rep) {
    ts::TreeState tree;
    ts::biased_wagner_tree(tree, ds, wp, nullptr);
    wagner_scores[rep] = ts::score_tree(tree, ds);

    if (run_tbr) {
      ts::TBRParams tp;
      ts::tbr_search(tree, ds, tp, nullptr, nullptr, nullptr, nullptr);
      tbr_scores[rep] = ts::score_tree(tree, ds);
    }
  }

  return List::create(
    Named("wagner_score")    = wagner_scores,
    Named("tbr_score")       = tbr_scores,
    Named("goloboff_scores") = goloboff_scores_r,
    Named("entropy_scores")  = entropy_scores_r
  );
}


// Parallel tempering functions (ts_stochastic_tbr, ts_parallel_temper)
// removed — live on feature/parallel-temper branch.

// [[Rcpp::export]]
List ts_test_strategy_tracker(int seed, int n_draws) {
  using ts::StrategyTracker;
  using ts::StartStrategy;
  using ts::N_STRAT;

  StrategyTracker tracker;
  std::mt19937 rng(seed);

  // 1. Draw `n_draws` strategies and count selections
  IntegerVector counts(N_STRAT, 0);
  for (int i = 0; i < n_draws; ++i) {
    auto s = tracker.select(rng);
    counts[static_cast<int>(s)]++;
  }

  // 2. Record initial alpha/beta
  NumericVector alpha_init(N_STRAT), beta_init(N_STRAT);
  for (int i = 0; i < N_STRAT; ++i) {
    alpha_init[i] = tracker.alpha(static_cast<StartStrategy>(i));
    beta_init[i] = tracker.beta_param(static_cast<StartStrategy>(i));
  }

  // 3. Update: arm 0 gets 5 successes, arm 1 gets 5 failures
  for (int i = 0; i < 5; ++i) {
    tracker.update(StartStrategy::WAGNER_RANDOM, true);
    tracker.update(StartStrategy::WAGNER_GOLOBOFF, false);
  }

  NumericVector alpha_after_update(N_STRAT), beta_after_update(N_STRAT);
  for (int i = 0; i < N_STRAT; ++i) {
    alpha_after_update[i] = tracker.alpha(static_cast<StartStrategy>(i));
    beta_after_update[i] = tracker.beta_param(static_cast<StartStrategy>(i));
  }

  // 4. Decay
  tracker.decay(0.5);
  NumericVector alpha_after_decay(N_STRAT), beta_after_decay(N_STRAT);
  for (int i = 0; i < N_STRAT; ++i) {
    alpha_after_decay[i] = tracker.alpha(static_cast<StartStrategy>(i));
    beta_after_decay[i] = tracker.beta_param(static_cast<StartStrategy>(i));
  }

  // 5. Post-update selection distribution (arm 0 should dominate)
  IntegerVector counts_biased(N_STRAT, 0);
  for (int i = 0; i < n_draws; ++i) {
    auto s = tracker.select(rng);
    counts_biased[static_cast<int>(s)]++;
  }

  // 6. Round-robin
  auto rr = StrategyTracker::round_robin(12);
  IntegerVector round_robin_seq(12);
  for (int i = 0; i < 12; ++i) {
    round_robin_seq[i] = static_cast<int>(rr[i]);
  }

  // 7. Strategy names
  CharacterVector names(N_STRAT);
  for (int i = 0; i < N_STRAT; ++i) {
    names[i] = ts::strategy_name(static_cast<StartStrategy>(i));
  }

  return List::create(
    Named("n_strategies") = N_STRAT,
    Named("strategy_names") = names,
    Named("initial_counts") = counts,
    Named("alpha_init") = alpha_init,
    Named("beta_init") = beta_init,
    Named("alpha_after_update") = alpha_after_update,
    Named("beta_after_update") = beta_after_update,
    Named("alpha_after_decay") = alpha_after_decay,
    Named("beta_after_decay") = beta_after_decay,
    Named("biased_counts") = counts_biased,
    Named("round_robin") = round_robin_seq
  );
}

// Run TBR to convergence with per-pass diagnostic instrumentation.
// Returns the final tree, scalar summary statistics, and a data frame
// of TBRPassRecords (one row per pass of the outer while loop).
// clipOrder: 0=RANDOM (baseline), 1=INV_WEIGHT, 2=TIPS_FIRST, 3=BUCKET
// [[Rcpp::export]]
List ts_tbr_diagnostics(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int maxHits = 1,
    bool acceptEqual = false,
    int maxChanges = 0,
    IntegerVector min_steps = IntegerVector(),
    double concavity = -1.0,
    int clipOrder = 0,
    bool unrooted = true)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  ts::TBRParams params;
  params.max_hits = maxHits;
  params.accept_equal = acceptEqual;
  params.max_accepted_changes = maxChanges;
  params.diagnostics = true;
  params.clip_order = static_cast<ts::ClipOrder>(clipOrder);
  params.unrooted = unrooted;

  ts::TBRResult result = ts::tbr_search(tree, ds, params);

  // Unpack pass records into parallel vectors for a data frame
  int n_passes = static_cast<int>(result.pass_records.size());
  IntegerVector pass_index(n_passes);
  LogicalVector productive(n_passes);
  IntegerVector accepted_clip_size(n_passes);
  IntegerVector n_clips_tried(n_passes);
  IntegerVector n_candidates_evaluated(n_passes);

  for (int i = 0; i < n_passes; ++i) {
    const ts::TBRPassRecord& rec = result.pass_records[i];
    pass_index[i]            = rec.pass_index;
    productive[i]            = rec.productive;
    accepted_clip_size[i]    = rec.accepted_clip_size;
    n_clips_tried[i]         = rec.n_clips_tried;
    n_candidates_evaluated[i]= rec.n_candidates_evaluated;
  }

  DataFrame passes = DataFrame::create(
    Named("pass_index")             = pass_index,
    Named("productive")             = productive,
    Named("accepted_clip_size")     = accepted_clip_size,
    Named("n_clips_tried")          = n_clips_tried,
    Named("n_candidates_evaluated") = n_candidates_evaluated
  );

  return List::create(
    Named("edge")          = tree_to_edge(tree),
    Named("score")         = result.best_score,
    Named("n_accepted")    = result.n_accepted,
    Named("n_evaluated")   = result.n_evaluated,
    Named("n_zero_skipped")= result.n_zero_skipped,
    Named("converged")     = result.converged,
    Named("passes")        = passes
  );
}

// Regression probe for the exact_verify_sweep optimum cache (test-ts-na-evcache.R).
// Returns the EXACT 64-bit cache key the NA convergence certifier would use for
// this (topology, dataset, weighting-regime) triple, as a 16-hex-digit string,
// by calling ts::exact_verify_cache_key — the same helper the cache itself uses.
//
// The three perturbation flags reproduce the ratchet mutating the live DataSet
// in place (ts_ratchet.cpp): `zero_active` clears one active_mask bit (ZERO_ONLY,
// the DEFAULT NA ratchet strategy), `set_upweight` sets one upweight_mask bit
// (UPWEIGHT/MIXED modes), `bump_pattern_freq` increments pattern_freq (IW).  Each
// must change the returned key; if it does not, a base-regime "optimal" verdict
// would leak into a perturbed pass and silently skip the improving moves the
// ratchet exists to find.  None of these flags touch ds_fingerprint's inputs
// (n_tips/n_blocks/tip_states), so any observed key change is attributable to the
// weighting-regime term alone — pinning that weight_fingerprint is XORed into the
// composite key, not merely that it exists.
// [[Rcpp::export]]
std::string ts_ev_cache_key_probe(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    double concavity = -1.0,
    bool zero_active = false,
    bool set_upweight = false,
    bool bump_pattern_freq = false)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                IntegerVector(), concavity);
  ts::TreeState tree;
  tree.init_from_edge(&edge(0, 0), &edge(0, 1), edge.nrow(), ds);

  if (zero_active) {
    // Clear the lowest set active bit in the first block that has one, keeping
    // the FlatBlock cache in sync (exactly perturb_zero's bookkeeping).
    for (int b = 0; b < ds.n_blocks; ++b) {
      uint64_t m = ds.blocks[b].active_mask;
      if (m) {
        ds.blocks[b].active_mask = m & (m - 1);          // clear lowest set bit
        ds.flat_blocks[b].active_mask = ds.blocks[b].active_mask;
        break;
      }
    }
  }
  if (set_upweight) {
    // Set the lowest active bit in upweight_mask (must be a subset of active).
    for (int b = 0; b < ds.n_blocks; ++b) {
      uint64_t a = ds.blocks[b].active_mask;
      if (a) {
        ds.blocks[b].upweight_mask |= (a & (~a + 1));    // lowest set active bit
        break;
      }
    }
  }
  if (bump_pattern_freq && !ds.pattern_freq.empty()) {
    ds.pattern_freq[0] += 1;
  }

  uint64_t key = ts::exact_verify_cache_key(tree, ds);
  char buf[17];
  std::snprintf(buf, sizeof(buf), "%016llx", (unsigned long long)key);
  return std::string(buf);
}
