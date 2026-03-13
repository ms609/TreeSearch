#include "ts_drift.h"
#include "ts_fitch.h"
#include "ts_tbr.h"
#include <algorithm>
#include <random>
#include <vector>
#include <climits>
#include <cmath>

#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

namespace ts {

// --- Helpers (file-local, mirrored from ts_tbr.cpp) ---

static double drift_full_rescore(TreeState& tree, const DataSet& ds) {
  tree.reset_states(ds);
  return score_tree(tree, ds);
}

// Extract per-character step counts from local_cost only (standard Fitch).
// Same as TBR's extract_divided_steps helper.
static void drift_extract_divided_steps(
    const TreeState& tree, const DataSet& ds,
    std::vector<int>& char_steps) {
  std::fill(char_steps.begin(), char_steps.end(), 0);
  for (int node : tree.postorder) {
    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      uint64_t mask =
          tree.local_cost[static_cast<size_t>(node) * tree.n_blocks + b];
      while (mask) {
        int c = ts::ctz64(mask);
        char_steps[blk.pattern_index[c]] += 1;
        mask &= mask - 1;
      }
    }
  }
}

static void drift_collect_main_edges(
    const TreeState& tree,
    std::vector<std::pair<int,int>>& edges)
{
  edges.clear();
  std::vector<int> stack;
  stack.push_back(tree.n_tip);

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

static void drift_collect_subtree_edges(
    const TreeState& tree, int subtree_root,
    std::vector<std::pair<int,int>>& edges)
{
  edges.clear();
  if (subtree_root < tree.n_tip) return;

  std::vector<int> stack;
  stack.push_back(subtree_root);

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

static void drift_fitch_join_states(
    const uint64_t* state_a,
    const uint64_t* state_b,
    uint64_t* out,
    const DataSet& ds)
{
  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    int offset = ds.block_word_offset[b];

    uint64_t any_isect = 0;
    for (int s = 0; s < blk.n_states; ++s) {
      any_isect |= (state_a[offset + s] & state_b[offset + s]);
    }
    uint64_t no_isect = ~any_isect & blk.active_mask;

    for (int s = 0; s < blk.n_states; ++s) {
      uint64_t isect = state_a[offset + s] & state_b[offset + s];
      uint64_t uni   = state_a[offset + s] | state_b[offset + s];
      out[offset + s] = (isect & any_isect) | (uni & no_isect);
    }
  }
}

static void drift_compute_from_above(
    const TreeState& tree, const DataSet& ds,
    int subtree_root,
    std::vector<uint64_t>& from_above)
{
  int tw = tree.total_words;

  std::vector<int> preorder;
  {
    std::vector<int> stack;
    stack.push_back(subtree_root);
    while (!stack.empty()) {
      int node = stack.back();
      stack.pop_back();
      if (node < tree.n_tip) continue;
      preorder.push_back(node);
      int ni = node - tree.n_tip;
      stack.push_back(tree.right[ni]);
      stack.push_back(tree.left[ni]);
    }
  }

  if (preorder.empty()) return;

  int root = preorder[0];
  int ri = root - tree.n_tip;
  int lc = tree.left[ri];
  int rc = tree.right[ri];

  for (int w = 0; w < tw; ++w) {
    from_above[static_cast<size_t>(lc) * tw + w] =
        tree.prelim[static_cast<size_t>(rc) * tw + w];
    from_above[static_cast<size_t>(rc) * tw + w] =
        tree.prelim[static_cast<size_t>(lc) * tw + w];
  }

  for (size_t i = 1; i < preorder.size(); ++i) {
    int node = preorder[i];
    int ni = node - tree.n_tip;
    int lc_i = tree.left[ni];
    int rc_i = tree.right[ni];

    drift_fitch_join_states(
        &from_above[static_cast<size_t>(node) * tw],
        &tree.prelim[static_cast<size_t>(rc_i) * tw],
        &from_above[static_cast<size_t>(lc_i) * tw], ds);
    drift_fitch_join_states(
        &from_above[static_cast<size_t>(node) * tw],
        &tree.prelim[static_cast<size_t>(lc_i) * tw],
        &from_above[static_cast<size_t>(rc_i) * tw], ds);
  }
}

// --- Topology snapshot for safe undo ---

struct DriftTopoSnapshot {
  std::vector<int> parent;
  std::vector<int> left;
  std::vector<int> right;
};

static void drift_save_topology(const TreeState& tree,
                                DriftTopoSnapshot& snap) {
  snap.parent = tree.parent;
  snap.left = tree.left;
  snap.right = tree.right;
}

static void drift_restore_topology(TreeState& tree,
                                   const DriftTopoSnapshot& snap) {
  tree.parent = snap.parent;
  tree.left = snap.left;
  tree.right = snap.right;
}

// --- Topology validation ---

static bool drift_validate_topology(const TreeState& tree) {
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

// Apply a TBR move directly to the tree topology.
// Mirrors apply_tbr_move from ts_tbr.cpp.
static bool drift_apply_tbr_move(
    TreeState& tree,
    int clip_node,
    int reroot_parent, int reroot_child,
    int above, int below)
{
  int nx = tree.parent[clip_node];
  int nz = tree.parent[nx];
  int nxi = nx - tree.n_tip;
  int ns;
  ns = (tree.left[nxi] == clip_node) ? tree.right[nxi] : tree.left[nxi];

  // Step 1: Detach
  tree.parent[ns] = nz;
  if (nz >= tree.n_tip) {
    int nzi = nz - tree.n_tip;
    if (tree.left[nzi] == nx) {
      tree.left[nzi] = ns;
    } else {
      tree.right[nzi] = ns;
    }
  }

  // Step 2: Reroot the clipped subtree if needed
  int new_subtree_root = clip_node;
  if (reroot_parent >= 0 && reroot_parent != clip_node) {
    std::vector<int> path;
    {
      std::vector<int> dfs_stack;
      std::vector<int> sub_parent(tree.n_node, -1);
      dfs_stack.push_back(clip_node);

      while (!dfs_stack.empty()) {
        int node = dfs_stack.back();
        dfs_stack.pop_back();
        if (node == reroot_parent) break;
        if (node < tree.n_tip) continue;
        int ni = node - tree.n_tip;
        int lc = tree.left[ni];
        int rc = tree.right[ni];
        sub_parent[lc] = node;
        sub_parent[rc] = node;
        dfs_stack.push_back(lc);
        dfs_stack.push_back(rc);
      }

      int cur = reroot_parent;
      while (cur != clip_node && cur >= 0) {
        path.push_back(cur);
        cur = sub_parent[cur];
      }
      if (cur < 0) return false;
      path.push_back(clip_node);
      std::reverse(path.begin(), path.end());
    }

    if (path.size() < 2) return false;

    for (size_t j = 0; j + 1 < path.size(); ++j) {
      int A = path[j];
      int B = path[j + 1];

      int ai = A - tree.n_tip;
      int bi = B - tree.n_tip;

      int B_off_path;
      if (j + 2 < path.size()) {
        int next_on_path = path[j + 2];
        B_off_path = (tree.left[bi] == next_on_path)
                     ? tree.right[bi] : tree.left[bi];
      } else {
        B_off_path = (tree.left[bi] == reroot_child)
                     ? tree.right[bi] : tree.left[bi];
      }

      if (tree.left[ai] == B) {
        tree.left[ai] = B_off_path;
      } else {
        tree.right[ai] = B_off_path;
      }
      tree.parent[B_off_path] = A;

      if (tree.left[bi] == B_off_path) {
        tree.left[bi] = A;
      } else {
        tree.right[bi] = A;
      }
      tree.parent[A] = B;
    }

    new_subtree_root = reroot_parent;
  }

  // Step 3: Regraft
  if (above >= tree.n_tip) {
    int ai = above - tree.n_tip;
    if (tree.left[ai] == below) {
      tree.left[ai] = nx;
    } else {
      tree.right[ai] = nx;
    }
  }
  tree.parent[nx] = above;

  tree.left[nxi] = new_subtree_root;
  tree.right[nxi] = below;
  tree.parent[new_subtree_root] = nx;
  tree.parent[below] = nx;

  return true;
}

// --- Drift phase ---
//
// Modified TBR loop that accepts suboptimal moves based on AFD/RFD criteria.
// Returns the number of accepted moves.
static int drift_phase(TreeState& tree, const DataSet& ds,
                       int afd_limit, double rfd_limit,
                       int max_changes, std::mt19937& rng) {
  double score = drift_full_rescore(tree, ds);
  int n_accepted = 0;
  const bool use_iw = std::isfinite(ds.concavity);
  const double eps = use_iw ? 1e-10 : 0.0;

  std::vector<int> clip_candidates;
  for (int node = 0; node < tree.n_node; ++node) {
    if (node == tree.n_tip) continue;
    clip_candidates.push_back(node);
  }

  std::vector<std::pair<int,int>> main_edges;
  std::vector<std::pair<int,int>> sub_edges;

  std::vector<uint64_t> from_above(
      static_cast<size_t>(tree.n_node) * tree.total_words, 0);
  std::vector<uint64_t> virtual_prelim(tree.total_words);

  // IW buffers
  std::vector<int> divided_steps;
  std::vector<double> iw_delta;
  if (use_iw) {
    divided_steps.resize(ds.n_patterns, 0);
    iw_delta.resize(ds.n_patterns, 0.0);
  }

  DriftTopoSnapshot snap;
  std::vector<uint64_t> old_local_cost;

  std::shuffle(clip_candidates.begin(), clip_candidates.end(), rng);

  for (int clip_node : clip_candidates) {
    if (tree.parent[clip_node] == tree.n_tip) continue;

    // --- Phase 1: Clip + indirect evaluation ---
    tree.spr_clip(clip_node);
    tree.build_postorder();

    int ns = tree.clip_state.clip_sibling;
    int nz = tree.clip_state.clip_grandpar;
    int nx = tree.clip_state.clip_parent;
    int delta = fitch_incremental_downpass(tree, ds, nz);
    fitch_incremental_uppass(tree, ds, nz);

    int nx_cost = 0;
    for (int b = 0; b < ds.n_blocks; ++b) {
      nx_cost += ds.blocks[b].weight * popcount64(
          tree.local_cost[static_cast<size_t>(nx) * tree.n_blocks + b]);
    }
    double divided_length = score + delta - nx_cost;

    // IW: precompute base IW and marginal deltas
    double base_iw = 0.0;
    if (use_iw) {
      drift_extract_divided_steps(tree, ds, divided_steps);
      base_iw = compute_iw(ds, divided_steps);
      precompute_iw_delta(ds, divided_steps, iw_delta);
    }

    drift_collect_main_edges(tree, main_edges);

    // Find best candidate via indirect evaluation
    double best_candidate = HUGE_VAL;
    int best_above = -1, best_below = -1;
    int best_reroot_parent = -1, best_reroot_child = -1;

    size_t clip_base = static_cast<size_t>(clip_node) * tree.total_words;
    const uint64_t* clip_prelim = &tree.prelim[clip_base];

    // SPR candidates
    for (auto& [above, below] : main_edges) {
      if (above == nz && below == ns) continue;
      double candidate;
      if (use_iw) {
        candidate = indirect_iw_length(clip_prelim, tree, ds,
                                       above, below, base_iw, iw_delta);
      } else {
        candidate = divided_length +
            fitch_indirect_length(clip_prelim, tree, ds, above, below);
      }
      if (candidate < best_candidate) {
        best_candidate = candidate;
        best_above = above;
        best_below = below;
        best_reroot_parent = -1;
        best_reroot_child = -1;
      }
    }

    // TBR candidates (rerooting)
    if (clip_node >= tree.n_tip) {
      drift_compute_from_above(tree, ds, clip_node, from_above);
      drift_collect_subtree_edges(tree, clip_node, sub_edges);

      for (auto& [sp, sc] : sub_edges) {
        if (sp == clip_node) continue;

        drift_fitch_join_states(
            &from_above[static_cast<size_t>(sc) * tree.total_words],
            &tree.prelim[static_cast<size_t>(sc) * tree.total_words],
            virtual_prelim.data(), ds);

        for (auto& [above, below] : main_edges) {
          if (above == nz && below == ns) continue;
          double candidate;
          if (use_iw) {
            candidate = indirect_iw_length(virtual_prelim.data(), tree, ds,
                                           above, below, base_iw, iw_delta);
          } else {
            candidate = divided_length +
                fitch_indirect_length(virtual_prelim.data(), tree, ds,
                                      above, below);
          }
          if (candidate < best_candidate) {
            best_candidate = candidate;
            best_above = above;
            best_below = below;
            best_reroot_parent = sp;
            best_reroot_child = sc;
          }
        }
      }
    }

    // --- Phase 2: Restore and decide ---
    tree.spr_unclip();
    tree.build_postorder();

    if (best_candidate >= HUGE_VAL || best_above < 0) continue;

    double delta_score = best_candidate - score;

    if (delta_score > afd_limit + eps) {
      continue;
    }

    // Save topology for potential undo / RFD
    drift_save_topology(tree, snap);

    bool ok = drift_apply_tbr_move(tree, clip_node,
                                    best_reroot_parent, best_reroot_child,
                                    best_above, best_below);

    if (!ok || !drift_validate_topology(tree)) {
      drift_restore_topology(tree, snap);
      tree.build_postorder();
      drift_full_rescore(tree, ds);
      continue;
    }

    if (delta_score < -eps) {
      // Improvement: always accept
      tree.build_postorder();
      score = drift_full_rescore(tree, ds);
      ++n_accepted;
    } else if (std::fabs(delta_score) <= eps) {
      // Equal: always accept
      tree.build_postorder();
      score = drift_full_rescore(tree, ds);
      ++n_accepted;
    } else {
      // Suboptimal but within AFD limit: check RFD
      tree.build_postorder();
      double new_score = drift_full_rescore(tree, ds);

      if (use_iw) {
        // Under IW, use score-based RFD: (worsening - improving) / worsening
        // Simplify to score delta ratio
        double rfd = (new_score > score && score > 0.0)
            ? (new_score - score) / new_score : 0.0;

        if (rfd <= rfd_limit) {
          score = new_score;
          ++n_accepted;
        } else {
          drift_restore_topology(tree, snap);
          tree.build_postorder();
          score = drift_full_rescore(tree, ds);
        }
      } else {
        // EW: original local_cost-based RFD
        std::vector<uint64_t> new_local_cost = tree.local_cost;

        drift_restore_topology(tree, snap);
        tree.build_postorder();
        drift_full_rescore(tree, ds);
        old_local_cost = tree.local_cost;

        double F = 0, C = 0;
        for (int node = tree.n_tip; node < tree.n_node; ++node) {
          for (int b = 0; b < ds.n_blocks; ++b) {
            size_t idx = static_cast<size_t>(node) * tree.n_blocks + b;
            int old_steps = popcount64(old_local_cost[idx]);
            int new_steps = popcount64(new_local_cost[idx]);
            int d = (new_steps - old_steps) * ds.blocks[b].weight;
            if (d > 0) F += d;
            if (d < 0) C += (-d);
          }
        }
        double rfd = (F == 0.0) ? 0.0 : (F - C) / F;

        if (rfd <= rfd_limit) {
          ok = drift_apply_tbr_move(tree, clip_node,
                                     best_reroot_parent, best_reroot_child,
                                     best_above, best_below);
          if (!ok || !drift_validate_topology(tree)) {
            tree.build_postorder();
            drift_full_rescore(tree, ds);
            continue;
          }
          tree.build_postorder();
          score = drift_full_rescore(tree, ds);
          ++n_accepted;
        } else {
          score = drift_full_rescore(tree, ds);
        }
      }
    }

    if (n_accepted >= max_changes) break;

    R_CheckUserInterrupt();
  }

  return n_accepted;
}

// --- Main drift search ---

DriftResult drift_search(TreeState& tree, const DataSet& ds,
                         const DriftParams& params) {
  double best_score = drift_full_rescore(tree, ds);
  int total_drift_moves = 0;
  int total_tbr_moves = 0;

  // Seed from R's RNG for reproducibility with set.seed()
  GetRNGstate();
  std::mt19937 rng(static_cast<unsigned>(unif_rand() * 4294967295.0));
  PutRNGstate();

  // Save the best tree topology
  DriftTopoSnapshot best_snap;
  drift_save_topology(tree, best_snap);

  int max_drift_changes = std::max(20,
      std::min(200, tree.n_tip / 8));

  for (int cycle = 1; cycle <= params.n_cycles; ++cycle) {

    // --- Perturbation phase ---
    if (cycle % 2 == 1) {
      // Suboptimal drift phase
      int drift_moves = drift_phase(tree, ds,
                                     params.afd_limit, params.rfd_limit,
                                     max_drift_changes, rng);
      total_drift_moves += drift_moves;
    } else {
      // Equal-score drift phase
      TBRParams eq_params;
      eq_params.accept_equal = true;
      eq_params.max_accepted_changes = tree.n_tip / 8;
      eq_params.max_hits = 100;  // generous for equal-score exploration

      TBRResult eq_result = tbr_search(tree, ds, eq_params);
      total_drift_moves += eq_result.n_accepted;
    }

    // --- Search phase: standard TBR to converge ---
    TBRParams search_params;
    search_params.accept_equal = false;
    search_params.max_accepted_changes = 0;  // run to convergence
    search_params.max_hits = params.max_hits;

    TBRResult search_result = tbr_search(tree, ds, search_params);
    total_tbr_moves += search_result.n_accepted;

    // Update best if improved
    if (search_result.best_score < best_score) {
      best_score = search_result.best_score;
      drift_save_topology(tree, best_snap);
    } else {
      // Restore best tree for next cycle
      drift_restore_topology(tree, best_snap);
      tree.build_postorder();
    }

    R_CheckUserInterrupt();
  }

  // Ensure tree is the best found
  drift_restore_topology(tree, best_snap);
  tree.build_postorder();
  drift_full_rescore(tree, ds);

  return DriftResult{best_score, params.n_cycles,
                     total_drift_moves, total_tbr_moves};
}

} // namespace ts
