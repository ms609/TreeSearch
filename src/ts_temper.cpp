#include "ts_temper.h"
#include "ts_collapsed.h"
#include "ts_constraint.h"
#include "ts_fitch.h"
#include "ts_tbr.h"
#include "ts_rng.h"
#include <algorithm>
#include <random>
#include <vector>
#include <cmath>
#include <climits>

#include <Rcpp.h>
#include <R.h>
#include <Rinternals.h>

namespace ts {

// --- Helpers (file-local) ---

static double temper_full_rescore(TreeState& tree, const DataSet& ds) {
  tree.reset_states(ds);
  return score_tree(tree, ds);
}

// Topology snapshot for undo after failed move application.
struct TemperTopoSnapshot {
  std::vector<int> parent;
  std::vector<int> left;
  std::vector<int> right;
};

static void temper_save_topology(const TreeState& tree,
                                 TemperTopoSnapshot& snap) {
  snap.parent = tree.parent;
  snap.left = tree.left;
  snap.right = tree.right;
}

static void temper_restore_topology(TreeState& tree,
                                    const TemperTopoSnapshot& snap) {
  tree.parent = snap.parent;
  tree.left = snap.left;
  tree.right = snap.right;
}

static bool temper_validate_topology(const TreeState& tree) {
  int root = tree.n_tip;
  std::vector<bool> visited(tree.n_node, false);
  std::vector<int> stack;
  stack.push_back(root);
  int n_visited = 0;

  while (!stack.empty()) {
    int node = stack.back();
    stack.pop_back();
    if (node < 0 || node >= tree.n_node) return false;
    if (visited[node]) return false;
    visited[node] = true;
    ++n_visited;
    if (node >= tree.n_tip) {
      int ni = node - tree.n_tip;
      stack.push_back(tree.left[ni]);
      stack.push_back(tree.right[ni]);
    }
  }
  return (n_visited == tree.n_node);
}

// Collect all parent->child edges in the (divided) tree.
static void temper_collect_edges(
    const TreeState& tree,
    std::vector<std::pair<int, int>>& edges) {
  edges.clear();
  std::vector<int> stack;
  stack.push_back(tree.n_tip);  // root
  while (!stack.empty()) {
    int node = stack.back();
    stack.pop_back();
    if (node < tree.n_tip) continue;
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    edges.push_back({node, lc});
    edges.push_back({node, rc});
    stack.push_back(lc);
    stack.push_back(rc);
  }
}

// Apply SPR move to full (unclipped) tree.
// Detaches clip_node's subtree, regrafts onto edge (above, below).
static bool temper_apply_spr_move(
    TreeState& tree, int clip_node,
    int above, int below) {
  int nx = tree.parent[clip_node];
  int nz = tree.parent[nx];
  int nxi = nx - tree.n_tip;
  int ns = (tree.left[nxi] == clip_node)
               ? tree.right[nxi] : tree.left[nxi];

  // Step 1: Detach — bypass nx, connect sibling to grandparent
  tree.parent[ns] = nz;
  if (nz >= tree.n_tip) {
    int nzi = nz - tree.n_tip;
    if (tree.left[nzi] == nx)
      tree.left[nzi] = ns;
    else
      tree.right[nzi] = ns;
  }

  // Step 2: Regraft — insert nx between above and below
  if (above >= tree.n_tip) {
    int ai = above - tree.n_tip;
    if (tree.left[ai] == below)
      tree.left[ai] = nx;
    else
      tree.right[ai] = nx;
  }
  tree.parent[nx] = above;
  tree.left[nxi] = clip_node;
  tree.right[nxi] = below;
  tree.parent[clip_node] = nx;
  tree.parent[below] = nx;

  return true;
}

// --- Main stochastic TBR phase ---

TemperResult stochastic_tbr_phase(
    TreeState& tree, const DataSet& ds,
    const TemperParams& params,
    ConstraintData* cd,
    std::function<bool()> check_timeout) {

  bool constrained = cd && cd->active;
  if (constrained) update_constraint(tree, *cd);

  double score = temper_full_rescore(tree, ds);
  double best_score = score;
  const bool use_iw = std::isfinite(ds.concavity);
  const double eps = use_iw ? 1e-10 : 0.0;
  const double temperature = params.temperature;

  bool has_na = false;
  for (int b = 0; b < ds.n_blocks; ++b) {
    if (ds.blocks[b].has_inapplicable) { has_na = true; break; }
  }

  std::mt19937 rng = ts::make_rng();

  // Build initial clip candidate list (all nodes except root and its children)
  std::vector<int> clip_candidates;
  clip_candidates.reserve(tree.n_node);
  for (int node = 0; node < tree.n_node; ++node) {
    if (node == tree.n_tip) continue;
    clip_candidates.push_back(node);
  }
  int n_cand = static_cast<int>(clip_candidates.size());

  // Collapsed flags
  std::vector<uint8_t> collapsed;
  compute_collapsed_flags(tree, ds, collapsed);

  // Edge collection buffer
  std::vector<std::pair<int, int>> main_edges;
  main_edges.reserve(2 * tree.n_tip);

  // IW buffers
  std::vector<int> divided_steps;
  std::vector<double> iw_delta;
  if (use_iw) {
    divided_steps.resize(ds.n_patterns, 0);
    iw_delta.resize(ds.n_patterns, 0.0);
  }

  // NA clip-actives buffer
  std::vector<uint64_t> clip_actives_buf(has_na ? tree.total_words : 0);

  // Pre-allocated undo stack
  TreeState::PreallocUndo fast_undo;
  fast_undo.init(3 * tree.n_node, tree.total_words, tree.n_blocks, has_na);
  tree.prealloc_undo = &fast_undo;

  // Work buffer for postorder rebuild
  std::vector<int> work_stack;
  work_stack.reserve(tree.n_node * 2);

  // Save postorder for restore after each clip
  std::vector<int> saved_postorder = tree.postorder;

  TemperTopoSnapshot snap;
  TemperResult result = {best_score, score, 0, 0, 0};

  for (int step = 0; step < params.n_moves; ++step) {
    ++result.n_attempted;

    // --- Pick random clip node ---
    int clip_node = clip_candidates[
        std::uniform_int_distribution<int>(0, n_cand - 1)(rng)];

    // Guard: root's children can't be clipped (would disconnect tree)
    if (tree.parent[clip_node] == tree.n_tip) continue;

    // Skip collapsed edges (can't change score)
    if (!collapsed.empty() && collapsed[clip_node]) continue;

    // --- Clip ---
    const uint64_t* clip_actives = nullptr;
    if (has_na) {
      size_t clip_sa_base =
          static_cast<size_t>(clip_node) * tree.total_words;
      std::copy(
          tree.subtree_actives.begin() + clip_sa_base,
          tree.subtree_actives.begin() + clip_sa_base + tree.total_words,
          clip_actives_buf.begin());
      clip_actives = clip_actives_buf.data();
    }

    fast_undo.clear();
    tree.spr_clip(clip_node);
    tree.build_postorder_prealloc(work_stack);

    int ns = tree.clip_state.clip_sibling;
    int nz = tree.clip_state.clip_grandpar;
    int nx = tree.clip_state.clip_parent;

    // --- Incremental scoring on divided tree ---
    double divided_length;
    if (has_na) {
      fitch_na_incremental_downpass(tree, ds, nz);
      fitch_na_incremental_uppass(tree, ds, nz);
      divided_length = static_cast<double>(fitch_na_pass3_score(tree, ds));
    } else {
      int delta = fitch_incremental_downpass(tree, ds, nz);
      fitch_incremental_uppass(tree, ds, nz);
      // Subtract the removed node's cost
      int nx_cost = 0;
      for (int b = 0; b < ds.n_blocks; ++b) {
        uint64_t lc =
            tree.local_cost[static_cast<size_t>(nx) * tree.n_blocks + b];
        int nu = popcount64(lc);
        if (ds.blocks[b].upweight_mask)
          nu += popcount64(lc & ds.blocks[b].upweight_mask);
        nx_cost += ds.blocks[b].weight * nu;
      }
      divided_length = score + delta - nx_cost;
    }

    // IW precomputation for this clip
    double base_iw = 0.0;
    if (use_iw) {
      std::fill(divided_steps.begin(), divided_steps.end(), 0);
      extract_char_steps(tree, ds, divided_steps);
      base_iw = compute_weighted_score(ds, divided_steps);
      precompute_weighted_delta(ds, divided_steps, iw_delta);
    }

    // Constraint classification for this clip
    if (constrained) classify_clip_constraints(tree, clip_node, *cd);

    // --- Collect edges and pick ONE random regraft ---
    temper_collect_edges(tree, main_edges);
    int n_edges = static_cast<int>(main_edges.size());

    bool found_candidate = false;
    double candidate_score = HUGE_VAL;
    int target_above = -1, target_below = -1;

    // Retry a few times if we hit the original position or a constraint
    for (int attempt = 0; attempt < 5 && !found_candidate; ++attempt) {
      int edge_idx =
          std::uniform_int_distribution<int>(0, n_edges - 1)(rng);
      auto [above, below] = main_edges[edge_idx];

      // Skip original position
      if (above == nz && below == ns) continue;
      // Skip collapsed regraft targets
      if (!collapsed.empty() && collapsed[below]) continue;
      // Skip constraint violations
      if (constrained && regraft_violates_constraint(below, *cd)) continue;

      // --- Evaluate candidate via indirect scoring ---
      size_t clip_base =
          static_cast<size_t>(clip_node) * tree.total_words;
      const uint64_t* clip_prelim = &tree.prelim[clip_base];

      if (has_na) {
        if (use_iw) {
          candidate_score = indirect_na_iw_length_bounded(
              clip_prelim, clip_actives, tree, ds, above, below,
              base_iw, iw_delta, HUGE_VAL);
        } else {
          int indirect = fitch_na_indirect_length(
              clip_prelim, clip_actives, tree, ds, above, below);
          candidate_score = divided_length + indirect;
        }
      } else if (use_iw) {
        candidate_score = indirect_iw_length(
            clip_prelim, tree, ds, above, below, base_iw, iw_delta);
      } else {
        int indirect = fitch_indirect_length(
            clip_prelim, tree, ds, above, below);
        candidate_score = divided_length + indirect;
      }

      target_above = above;
      target_below = below;
      found_candidate = true;
    }

    // --- Restore divided tree ---
    tree.restore_prealloc_undo();
    tree.spr_unclip();
    tree.postorder.assign(saved_postorder.begin(), saved_postorder.end());

    if (!found_candidate) continue;

    // --- Boltzmann acceptance ---
    double delta_score = candidate_score - score;
    bool accept = false;

    if (delta_score < -eps) {
      // Improvement: always accept
      accept = true;
    } else if (std::fabs(delta_score) <= eps) {
      // Equal score: always accept
      accept = true;
    } else if (temperature > 0.0) {
      // Suboptimal: accept with probability exp(-delta/T)
      double p = std::exp(-delta_score / temperature);
      std::uniform_real_distribution<double> unif(0.0, 1.0);
      accept = (unif(rng) < p);
    }
    // temperature == 0 and delta > 0: reject (strict hill-climbing)

    if (!accept) continue;

    // --- Apply move ---
    temper_save_topology(tree, snap);
    bool ok = temper_apply_spr_move(tree, clip_node,
                                     target_above, target_below);

    if (!ok || !temper_validate_topology(tree)) {
      temper_restore_topology(tree, snap);
      tree.build_postorder();
      temper_full_rescore(tree, ds);
      continue;
    }

    tree.build_postorder();
    score = temper_full_rescore(tree, ds);
    saved_postorder = tree.postorder;
    ++result.n_accepted;

    if (delta_score < -eps) ++result.n_improved;
    if (score < best_score) best_score = score;

    // Recompute collapsed flags
    compute_collapsed_flags(tree, ds, collapsed);

    if (constrained) update_constraint(tree, *cd);

    if (ts::check_interrupt()) break;
    if (check_timeout && check_timeout()) break;
  }

  tree.prealloc_undo = nullptr;

  result.best_score = best_score;
  result.final_score = score;
  return result;
}

// --- Layer 2: Multi-chain parallel tempering ---

// Deep-copy a TreeState (topology + all state arrays).
static TreeState copy_tree_state(const TreeState& src) {
  TreeState dst;
  dst.n_tip = src.n_tip;
  dst.n_internal = src.n_internal;
  dst.n_node = src.n_node;
  dst.total_words = src.total_words;
  dst.n_blocks = src.n_blocks;
  dst.parent = src.parent;
  dst.left = src.left;
  dst.right = src.right;
  dst.prelim = src.prelim;
  dst.final_ = src.final_;
  dst.down2 = src.down2;
  dst.subtree_actives = src.subtree_actives;
  dst.local_cost = src.local_cost;
  dst.postorder = src.postorder;
  return dst;
}

// Swap two chains' trees and scores (topology only — state arrays are
// rebuilt by full_rescore after swap).
static void swap_chains(TreeState& a, double& score_a,
                        TreeState& b, double& score_b) {
  std::swap(a.parent, b.parent);
  std::swap(a.left, b.left);
  std::swap(a.right, b.right);
  std::swap(a.postorder, b.postorder);
  std::swap(score_a, score_b);
}

PTResult parallel_temper_search(
    TreeState& tree, const DataSet& ds,
    const PTParams& params,
    TreePool* pool,
    ConstraintData* cd,
    std::function<bool()> check_timeout) {

  int nc = params.n_chains;
  if (nc < 2) nc = 2;

  // Build temperature ladder (default: geometric spacing)
  std::vector<double> temps = params.temperatures;
  if (static_cast<int>(temps.size()) < nc) {
    temps.resize(nc);
    temps[0] = 0.0;
    for (int i = 1; i < nc; ++i) {
      temps[i] = std::pow(3.0, i);  // 0, 3, 9, 27, ...
    }
  }

  int moves_per_round = params.moves_per_round;
  if (moves_per_round <= 0) moves_per_round = tree.n_tip;

  // Initialize chains: chain 0 is cold, others are hot copies
  std::vector<TreeState> chains(nc);
  std::vector<double> scores(nc);

  tree.reset_states(ds);
  chains[0] = copy_tree_state(tree);
  scores[0] = score_tree(chains[0], ds);

  for (int i = 1; i < nc; ++i) {
    chains[i] = copy_tree_state(tree);
    scores[i] = scores[0];
  }

  double pool_best = (pool && pool->best_score() < 1e17)
                         ? pool->best_score()
                         : scores[0];
  double overall_best = scores[0];

  std::mt19937 rng = ts::make_rng();

  PTResult result = {};
  result.best_score = overall_best;
  result.cold_final_score = scores[0];

  for (int round = 0; round < params.rounds; ++round) {

    // --- Cold chain: standard TBR to convergence ---
    {
      TBRParams tp;
      tp.accept_equal = false;
      tp.max_accepted_changes = 0;
      tp.max_hits = 1;
      TBRResult tr = tbr_search(chains[0], ds, tp, cd,
                                 nullptr, nullptr, check_timeout);
      scores[0] = tr.best_score;
    }

    // --- Hot chains: stochastic TBR ---
    for (int i = 1; i < nc; ++i) {
      TemperParams tp;
      tp.temperature = temps[i];
      tp.n_moves = moves_per_round;
      TemperResult tr = stochastic_tbr_phase(
          chains[i], ds, tp, cd, check_timeout);
      scores[i] = tr.final_score;

      // If a hot chain found a new best score, record it
      if (tr.best_score < overall_best) {
        overall_best = tr.best_score;
        ++result.hot_discoveries;
      }
    }

    // --- Attempt adjacent-temperature swaps (from hottest to coldest) ---
    for (int i = nc - 1; i > 0; --i) {
      ++result.total_swaps_attempted;

      // Metropolis swap criterion:
      //   accept with prob = min(1, exp((1/T_lo - 1/T_hi) * (s_lo - s_hi)))
      // where T_lo = temps[i-1], T_hi = temps[i]
      // For cold chain (T=0): accept only if hot chain score <= cold score.
      double s_lo = scores[i - 1];
      double s_hi = scores[i];

      bool accept_swap = false;
      if (temps[i - 1] == 0.0) {
        // Cold chain: only accept if hot chain is at least as good
        accept_swap = (s_hi <= s_lo + 1e-10);
      } else {
        // Both chains have T > 0: standard Metropolis
        double inv_t_lo = 1.0 / temps[i - 1];
        double inv_t_hi = 1.0 / temps[i];
        double log_alpha = (inv_t_lo - inv_t_hi) * (s_lo - s_hi);
        if (log_alpha >= 0.0) {
          accept_swap = true;
        } else {
          std::uniform_real_distribution<double> unif(0.0, 1.0);
          accept_swap = (std::log(unif(rng)) < log_alpha);
        }
      }

      if (accept_swap) {
        swap_chains(chains[i - 1], scores[i - 1],
                    chains[i], scores[i]);
        ++result.total_swaps_accepted;

        // If cold chain improved, full rescore to refresh state arrays
        if (i == 1) {
          chains[0].reset_states(ds);
          scores[0] = score_tree(chains[0], ds);
        }
      }
    }

    // --- Pool insertion for cold chain improvements ---
    if (pool && scores[0] < pool_best) {
      pool_best = scores[0];
      // Pool insertion will happen via caller (driven_search); we just
      // track the best score here.
    }

    // Track overall best
    for (int i = 0; i < nc; ++i) {
      if (scores[i] < overall_best) overall_best = scores[i];
    }

    if (ts::check_interrupt()) break;
    if (check_timeout && check_timeout()) break;
  }

  // Return cold chain's tree to caller
  tree = std::move(chains[0]);
  tree.reset_states(ds);
  score_tree(tree, ds);

  result.best_score = overall_best;
  result.cold_final_score = scores[0];
  return result;
}

} // namespace ts
