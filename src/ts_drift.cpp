#include "ts_drift.h"
#include "ts_collapsed.h"
#include "ts_constraint.h"
#include "ts_fitch.h"
#include "ts_tbr.h"
#include "ts_rng.h"
#include <algorithm>
#include <random>
#include <vector>
#include <climits>
#include <cmath>
#include <cstdlib>

#include <Rcpp.h>
#include <R.h>
#include <Rinternals.h>

namespace ts {

// Scan-accuracy census for the pure-EW drift candidate scorer (env
// TS_DRIFT_SCANCHK).  Quantifies how far the deployed union-of-finals
// approximation (fitch_indirect_length_bounded, reads the incrementally
// maintained tree.final_) diverges from the EXACT directional insertion cost
// (fitch_indirect_length_cached over compute_insertion_edge_sets), the same
// exact quantity tbr_search migrated to in 2b299e4b.
//
// The union error is TWO-SIDED on the deployed path (empirically): it both
// over- and under-counts by up to tens of steps -- NOT the strict over-count
// the exactness gate measured with FRESH finals (dev/profiling/exactness-gate.md
// P2).  The decision-relevant harm comes from the UNDER-counts via an
// optimizer's-curse effect: selecting argmin over (divided_length + u_cost)
// preferentially picks candidates whose u_cost < e_cost (the moves the estimator
// is most optimistic about), so the SELECTED move's error is biased negative --
// a truly-expensive move can look cheap, pass the AFD gate, and be applied
// outside drift's intended drift envelope.  The exact scan (error == 0) has no
// such bias.  Recorded scores stay exact regardless (every accept is
// drift_full_rescore'd) -- this measures the scan/decision error, not a leak.
struct DriftEwCensus {
  long long clips = 0;             // EW drift clips the census examined
  long long candidates = 0;        // EW candidates scored (SPR + reroot)
  long long overcount_cands = 0;   // candidates where union > exact
  long long overcount_sum = 0;     // sum of (union - exact) over over-counts
  int overcount_max = 0;           // largest single (union - exact)
  long long undercount_cands = 0;  // candidates where union < exact (dangerous)
  long long undercount_sum = 0;    // sum of (exact - union) over under-counts
  int undercount_max = 0;          // largest single (exact - union)
  long long select_flip = 0;       // clips where union-best move != exact-best move
  double select_regret_sum = 0.0;  // sum over select_flip clips of exact-cost regret
  long long env_violation = 0;     // clips where union admits a move whose TRUE
                                   // delta exceeds afd_limit (drift walks outside
                                   // its intended AFD envelope)
  // TS_IW_SCANCHK analog: predicted best_candidate vs post-apply full_rescore of
  // the applied move.  On the exact path this MUST be 0 (validates the wiring).
  long long applied_checked = 0;
  long long applied_mismatch = 0;
  long long applied_under = 0;     // applied moves truly WORSE than the scan said
                                   // (optimizer's-curse signature)
  double applied_maxdiff = 0.0;
};

// --- Helpers (file-local, mirrored from ts_tbr.cpp) ---

static double drift_full_rescore(TreeState& tree, const DataSet& ds) {
  tree.reset_states(ds);
  return score_tree(tree, ds);
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
                       int max_changes, std::mt19937& rng,
                       ConstraintData* cd = nullptr,
                       const std::vector<bool>* sector_mask = nullptr,
                       bool ew_exact = false, bool scanchk = false,
                       DriftEwCensus* census = nullptr) {
  bool constrained = cd && cd->active;
  if (constrained) update_constraint(tree, *cd);
  double score = drift_full_rescore(tree, ds);
  int n_accepted = 0;
  const bool use_iw = std::isfinite(ds.concavity);
  const double eps = use_iw ? 1e-10 : 0.0;

  bool has_na = false;
  for (int b = 0; b < ds.n_blocks; ++b) {
    if (ds.blocks[b].has_inapplicable) { has_na = true; break; }
  }

  // Compute subtree sizes for smaller-subtree filter
  std::vector<int> subtree_sizes(tree.n_node, 0);
  for (int i = 0; i < tree.n_tip; ++i) subtree_sizes[i] = 1;
  for (int node : tree.postorder) {
    int ni = node - tree.n_tip;
    subtree_sizes[node] = subtree_sizes[tree.left[ni]]
                        + subtree_sizes[tree.right[ni]];
  }
  int half_n = tree.n_tip / 2;

  std::vector<int> clip_candidates;
  for (int node = 0; node < tree.n_node; ++node) {
    if (node == tree.n_tip) continue;
    // Skip clips where subtree > n/2 (same optimization as tbr_search)
    if (subtree_sizes[node] > half_n) continue;
    // Sector-mask: only clip inside the masked clade (pins out-of-mask nodes,
    // e.g. a sector's HTU pseudo-tip and the synthetic root edge).
    if (sector_mask && !(*sector_mask)[node]) continue;
    clip_candidates.push_back(node);
  }

  // Collapsed flags: edges that provably cannot yield an improvement
  // (clip skipping + regraft merging).
  std::vector<uint8_t> collapsed;
  compute_collapsed_flags(tree, ds, collapsed);

  std::vector<std::pair<int,int>> main_edges;
  std::vector<std::pair<int,int>> sub_edges;

  std::vector<uint64_t> from_above(
      static_cast<size_t>(tree.n_node) * tree.total_words, 0);
  std::vector<uint64_t> virtual_prelim(tree.total_words);

  // Exact directional insertion edge sets for the pure-EW path, mirroring
  // tbr_search (ts_tbr.cpp:1749).  Computed once per clip; edge_set[below] is
  // the exact edge set above node `below`, used by both the SPR and reroot
  // scans.  Caller-owned scratch, size-ensured (non-zeroing) and reused across
  // clips.  Only touched on the EW path when the exact scorer or the census is
  // active; NA/IW drift never allocate these.
  const int tw = tree.total_words;
  const bool ew_path = !has_na && !use_iw;
  const bool need_edge_set = ew_path && (ew_exact || scanchk);
  std::vector<uint64_t> edge_set_buf;
  std::vector<uint64_t> edge_set_up;
  std::vector<int> edge_set_pre;

  // IW buffers
  std::vector<int> divided_steps;
  std::vector<double> iw_delta;
  if (use_iw) {
    divided_steps.resize(ds.n_patterns, 0);
    iw_delta.resize(ds.n_patterns, 0.0);
  }

  // Buffer for clip subtree's actives (NA indirect length)
  std::vector<uint64_t> clip_actives_buf(has_na ? tree.total_words : 0);

  DriftTopoSnapshot snap;
  std::vector<uint64_t> old_local_cost;

  // Pre-allocated undo stack (eliminates heap allocs in save_node_state)
  TreeState::PreallocUndo fast_undo;
  // Capacity must cover downpass + uppass + tips: up to 3 * n_node saves per clip
  fast_undo.init(3 * tree.n_node, tree.total_words, tree.n_blocks, has_na);
  tree.prealloc_undo = &fast_undo;

  // Pre-allocated work buffer for build_postorder_prealloc
  std::vector<int> work_stack;
  work_stack.reserve(tree.n_node * 2);

  // Save postorder for restore after unclip (avoids O(n) rebuild)
  std::vector<int> saved_postorder = tree.postorder;

  std::shuffle(clip_candidates.begin(), clip_candidates.end(), rng);

  for (int clip_node : clip_candidates) {
    if (tree.parent[clip_node] == tree.n_tip) continue;

    // Skip collapsed edges (zero-length, provably unimprovable).
    if (!collapsed.empty() && collapsed[clip_node])
      continue;

    // --- Phase 1: Clip + indirect evaluation ---

    // Save clip subtree actives before clipping
    const uint64_t* clip_actives = nullptr;
    if (has_na) {
      size_t clip_sa_base =
          static_cast<size_t>(clip_node) * tree.total_words;
      std::copy(tree.subtree_actives.begin() + clip_sa_base,
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

    double divided_length;
    if (has_na) {
      fitch_na_incremental_downpass(tree, ds, nz);
      fitch_na_incremental_uppass(tree, ds, nz);
      divided_length = static_cast<double>(fitch_na_pass3_score(tree, ds));
    } else {
      int delta = fitch_incremental_downpass(tree, ds, nz);
      fitch_incremental_uppass(tree, ds, nz);

      int nx_cost = 0;
      for (int b = 0; b < ds.n_blocks; ++b) {
        uint64_t lc = tree.local_cost[static_cast<size_t>(nx) * tree.n_blocks + b];
        int nu = popcount64(lc);
        if (ds.blocks[b].upweight_mask) nu += popcount64(lc & ds.blocks[b].upweight_mask);
        nx_cost += ds.blocks[b].weight * nu;
      }
      divided_length = score + delta - nx_cost;
    }

    // Exact directional insertion edge sets for the pure-EW path (mirrors
    // ts_tbr.cpp:1749): computed once per clip from the current clipped-tree
    // downpass, then indexed by `below` in both the SPR and reroot scans.
    if (need_edge_set) {
      compute_insertion_edge_sets(tree, ds, edge_set_buf,
                                  edge_set_up, edge_set_pre);
    }

    // Weighted scoring (IW or profile): precompute base score and deltas
    double base_iw = 0.0;
    if (use_iw) {
      std::fill(divided_steps.begin(), divided_steps.end(), 0);
      extract_char_steps(tree, ds, divided_steps);
      base_iw = compute_weighted_score(ds, divided_steps);
      precompute_weighted_delta(ds, divided_steps, iw_delta);
    }

    drift_collect_main_edges(tree, main_edges);
    // Partial shuffle: seed bound with diverse sample
    {
      int ne = static_cast<int>(main_edges.size());
      int k = std::min(20, ne);
      for (int i = 0; i < k; ++i) {
        std::uniform_int_distribution<int> dist(i, ne - 1);
        std::swap(main_edges[i], main_edges[dist(rng)]);
      }
    }

    // Constraint: classify this clip
    if (constrained) classify_clip_constraints(tree, clip_node, *cd);

    // Find best candidate via indirect evaluation
    double best_candidate = HUGE_VAL;
    int best_above = -1, best_below = -1;
    int best_reroot_parent = -1, best_reroot_child = -1;

    size_t clip_base = static_cast<size_t>(clip_node) * tree.total_words;
    const uint64_t* clip_prelim = &tree.prelim[clip_base];

    // Per-clip census bests (only used under scanchk): the union-selected best
    // move (cen_bc_union) with the EXACT cost of that same move (cen_e_at_union),
    // and the exact-selected best move (cen_bc_exact).  All are full candidate
    // scores (divided_length + insertion cost).
    double cen_bc_union = HUGE_VAL, cen_bc_exact = HUGE_VAL,
           cen_e_at_union = HUGE_VAL;

    // Pure-EW candidate scorer: returns divided_length + insertion cost.  When
    // ew_exact, selects on the EXACT directional edge set (edge_set_buf[below]),
    // matching tbr_search; otherwise the deployed union-of-finals approximation.
    // Under scanchk it computes BOTH (unbounded -- byte-identical selection to
    // the bounded form, since the bail only triggers once the value already
    // exceeds the cutoff) and feeds the decision-flip census.  `above` is only
    // needed by the union scorer; the exact scorer indexes by `below`.
    auto score_ew = [&](const uint64_t* prelim, int above, int below,
                        int ew_cutoff) -> double {
      if (!scanchk) {
        int extra = ew_exact
            ? fitch_indirect_length_cached(
                  prelim, &edge_set_buf[static_cast<size_t>(below) * tw],
                  ds, ew_cutoff)
            : fitch_indirect_length_bounded(prelim, tree, ds,
                                            above, below, ew_cutoff);
        return divided_length + extra;
      }
      int u_cost = fitch_indirect_length_bounded(prelim, tree, ds,
                                                 above, below, INT_MAX);
      int e_cost = fitch_indirect_length_cached(
          prelim, &edge_set_buf[static_cast<size_t>(below) * tw], ds, INT_MAX);
      ++census->candidates;
      int over = u_cost - e_cost;  // union - exact; TWO-SIDED on the deployed path
      if (over > 0) {
        ++census->overcount_cands;
        census->overcount_sum += over;
        if (over > census->overcount_max) census->overcount_max = over;
      } else if (over < 0) {
        ++census->undercount_cands;
        census->undercount_sum += (-over);
        if (-over > census->undercount_max) census->undercount_max = -over;
      }
      double cu = divided_length + u_cost, ce = divided_length + e_cost;
      if (cu < cen_bc_union) { cen_bc_union = cu; cen_e_at_union = ce; }
      if (ce < cen_bc_exact) { cen_bc_exact = ce; }
      return divided_length + (ew_exact ? e_cost : u_cost);
    };

    // SPR candidates (bounded to skip losing positions early)
    for (auto& [above, below] : main_edges) {
      if (above == nz && below == ns) continue;
      if (constrained && regraft_violates_constraint(below, *cd)) continue;
      // Sector-mask: only regraft onto edges inside the masked clade.
      if (sector_mask && !(*sector_mask)[below]) continue;
      // Collapsed-region regraft merging: skip interior collapsed edges.
      if (!collapsed.empty() && collapsed[below])
        continue;
      double candidate;
      if (has_na) {
        if (use_iw) {
          candidate = indirect_na_iw_length_bounded(clip_prelim, clip_actives,
              tree, ds, above, below, base_iw, iw_delta, best_candidate);
        } else {
          int ew_cutoff = (best_candidate < HUGE_VAL)
              ? static_cast<int>(best_candidate - divided_length) : INT_MAX;
          candidate = divided_length +
              fitch_na_indirect_length_bounded(clip_prelim, clip_actives,
                  tree, ds, above, below, ew_cutoff);
        }
      } else if (use_iw) {
        candidate = indirect_iw_length_bounded(clip_prelim, tree, ds,
                                       above, below, base_iw, iw_delta,
                                       best_candidate);
      } else {
        int ew_cutoff = (best_candidate < HUGE_VAL)
            ? static_cast<int>(best_candidate - divided_length) : INT_MAX;
        candidate = score_ew(clip_prelim, above, below, ew_cutoff);
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
          if (constrained && regraft_violates_constraint(below, *cd))
            continue;
          // Sector-mask: only regraft onto edges inside the masked clade.
          if (sector_mask && !(*sector_mask)[below]) continue;
          // Collapsed-region regraft merging (same as SPR loop).
          if (!collapsed.empty() && collapsed[below])
            continue;
          double candidate;
          if (has_na) {
            if (use_iw) {
              candidate = indirect_na_iw_length_bounded(
                  virtual_prelim.data(),
                  clip_actives, tree, ds, above, below,
                  base_iw, iw_delta, best_candidate);
            } else {
              int ew_cutoff = (best_candidate < HUGE_VAL)
                  ? static_cast<int>(best_candidate - divided_length)
                  : INT_MAX;
              candidate = divided_length +
                  fitch_na_indirect_length_bounded(virtual_prelim.data(),
                      clip_actives, tree, ds, above, below, ew_cutoff);
            }
          } else if (use_iw) {
            candidate = indirect_iw_length_bounded(
                virtual_prelim.data(), tree, ds,
                above, below, base_iw, iw_delta, best_candidate);
          } else {
            int ew_cutoff = (best_candidate < HUGE_VAL)
                ? static_cast<int>(best_candidate - divided_length)
                : INT_MAX;
            candidate = score_ew(virtual_prelim.data(), above, below, ew_cutoff);
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

    // Census: per-clip decision accounting.  select_flip counts clips where the
    // union-chosen move is EXACT-suboptimal (the exact scorer would pick a
    // strictly better move).  env_violation counts clips where the union scan
    // GATES IN its chosen move (union delta within the AFD limit) but that
    // move's TRUE cost exceeds the AFD limit -- drift then walks outside its
    // intended drift envelope, the direct consequence of the under-count.
    if (scanchk && cen_bc_union < HUGE_VAL) {
      ++census->clips;
      if (cen_e_at_union > cen_bc_exact + eps) {
        ++census->select_flip;
        census->select_regret_sum += (cen_e_at_union - cen_bc_exact);
      }
      bool union_gated_in = (cen_bc_union - score) <= afd_limit + eps;
      bool true_out_of_env = (cen_e_at_union - score) > afd_limit + eps;
      if (union_gated_in && true_out_of_env) ++census->env_violation;
    }

    // --- Phase 2: Restore and decide ---
    tree.restore_prealloc_undo();
    tree.spr_unclip();
    tree.postorder.assign(saved_postorder.begin(), saved_postorder.end());

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

    // Post-hoc constraint validation: TBR rerooting can break
    // splits classified as UNCONSTRAINED during clip phase.
    if (constrained) {
      tree.build_postorder();
      map_constraint_nodes(tree, *cd);
      bool violation = false;
      for (int _s = 0; _s < cd->n_splits; ++_s) {
        if (cd->constraint_node[_s] < 0) {
          violation = true;
          break;
        }
      }
      if (violation) {
        drift_restore_topology(tree, snap);
        tree.build_postorder();
        drift_full_rescore(tree, ds);
        update_constraint(tree, *cd);
        continue;
      }
    }

    // TS_IW_SCANCHK analog: the move has been applied and validated; compare the
    // scan's predicted best_candidate against the authoritative full_rescore.
    // One extra rescore per applied clip (diagnostic only); the branch logic
    // below recomputes score idempotently.
    if (scanchk) {
      tree.build_postorder();
      double applied_exact = drift_full_rescore(tree, ds);
      ++census->applied_checked;
      double signed_err = applied_exact - best_candidate;  // >0 => truth worse
      double d = std::fabs(signed_err);
      if (d > 1e-6) {
        ++census->applied_mismatch;
        if (signed_err > 0) ++census->applied_under;  // optimizer's-curse signature
        if (d > census->applied_maxdiff) census->applied_maxdiff = d;
      }
    }

    int n_before = n_accepted;

    if (delta_score < -eps) {
      // Improvement: always accept
      tree.build_postorder();
      score = drift_full_rescore(tree, ds);
      ++n_accepted;
      if (constrained) update_constraint(tree, *cd);
    } else if (std::fabs(delta_score) <= eps) {
      // Equal: always accept
      tree.build_postorder();
      score = drift_full_rescore(tree, ds);
      ++n_accepted;
      if (constrained) update_constraint(tree, *cd);
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
          if (constrained) update_constraint(tree, *cd);
        } else {
          drift_restore_topology(tree, snap);
          tree.build_postorder();
          score = drift_full_rescore(tree, ds);
          if (constrained) update_constraint(tree, *cd);
        }
      } else {
        // EW: original local_cost-based RFD
        std::vector<uint64_t> new_local_cost = tree.local_cost;

        drift_restore_topology(tree, snap);
        tree.build_postorder();
        score = drift_full_rescore(tree, ds);
        old_local_cost = tree.local_cost;

        double F = 0, C = 0;
        for (int node = tree.n_tip; node < tree.n_node; ++node) {
          for (int b = 0; b < ds.n_blocks; ++b) {
            size_t idx = static_cast<size_t>(node) * tree.n_blocks + b;
            int old_nu = popcount64(old_local_cost[idx]);
            if (ds.blocks[b].upweight_mask) old_nu += popcount64(old_local_cost[idx] & ds.blocks[b].upweight_mask);
            int new_nu = popcount64(new_local_cost[idx]);
            if (ds.blocks[b].upweight_mask) new_nu += popcount64(new_local_cost[idx] & ds.blocks[b].upweight_mask);
            int d = (new_nu - old_nu) * ds.blocks[b].weight;
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
            drift_restore_topology(tree, snap);
            tree.build_postorder();
            drift_full_rescore(tree, ds);
            continue;
          }
          tree.build_postorder();
          score = drift_full_rescore(tree, ds);
          ++n_accepted;
          if (constrained) update_constraint(tree, *cd);
        } else {
          // score already set when topology was restored above
          if (constrained) update_constraint(tree, *cd);
        }
      }
    }

    // Recompute collapsed regions after any accepted move.
    if (n_accepted > n_before) {
      compute_collapsed_flags(tree, ds, collapsed);
    }

    if (n_accepted >= max_changes) break;

    if (ts::check_interrupt()) break;
  }

  tree.prealloc_undo = nullptr;

  // Ensure postorder matches current topology. saved_postorder is only set
  // once (before any moves); after accepted moves + rejected unclip restore,
  // it can be stale.  The caller (tbr_search via drift_search) relies on
  // a correct postorder for full_rescore.
  if (n_accepted > 0) {
    tree.build_postorder();
  }
  return n_accepted;
}

// --- Main drift search ---

DriftResult drift_search(TreeState& tree, const DataSet& ds,
                         const DriftParams& params,
                         ConstraintData* cd,
                         std::function<bool()> check_timeout,
                         const std::vector<bool>* sector_mask) {
  double best_score = drift_full_rescore(tree, ds);

  // No informative characters: all trees have the same score.
  if (ds.total_words == 0) return {best_score, 0, 0};

  int total_drift_moves = 0;
  int total_tbr_moves = 0;

  // Env flags read ONCE per drift_search (never per candidate: std::getenv is
  // ~2.4us on ucrt and would confound matched-wall timing -- see memory node
  // getenv-ucrt-cost).  TS_DRIFT_EXACT routes the pure-EW drift scorer to the
  // exact directional edge set (opt-in; union-of-finals remains the default so
  // committed behaviour is unchanged until the matched-wall gate clears it).
  // TS_DRIFT_SCANCHK enables the divergence census below.
  const bool ew_exact = std::getenv("TS_DRIFT_EXACT") != nullptr;
  const bool scanchk  = std::getenv("TS_DRIFT_SCANCHK") != nullptr;
  DriftEwCensus census;

  // Seed RNG (from R in serial mode, from thread-local in parallel mode)
  std::mt19937 rng = ts::make_rng();

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
                                     max_drift_changes, rng, cd, sector_mask,
                                     ew_exact, scanchk,
                                     scanchk ? &census : nullptr);
      total_drift_moves += drift_moves;
    } else {
      // Equal-score drift phase
      TBRParams eq_params;
      eq_params.accept_equal = true;
      eq_params.max_accepted_changes = tree.n_tip / 8;
      eq_params.max_hits = 100;  // generous for equal-score exploration
      eq_params.tabu_size = params.tabu_size;

      TBRResult eq_result = tbr_search(tree, ds, eq_params, cd,
                                        sector_mask, nullptr, check_timeout);
      total_drift_moves += eq_result.n_accepted;
    }

    // --- Search phase: standard TBR to converge ---
    TBRParams search_params;
    search_params.accept_equal = false;
    search_params.max_accepted_changes = 0;  // run to convergence
    search_params.max_hits = params.max_hits;
    search_params.tabu_size = params.tabu_size;

    TBRResult search_result = tbr_search(tree, ds, search_params, cd,
                                          sector_mask, nullptr, check_timeout);
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

    if (ts::check_interrupt()) break;
    if (check_timeout && check_timeout()) break;
  }

  // Ensure tree is the best found
  drift_restore_topology(tree, best_snap);
  tree.build_postorder();
  drift_full_rescore(tree, ds);

  if (scanchk) {
    // over/under: two-sided per-candidate scan error vs the exact directional
    // cost.  select_flip/env_violation: per-clip decision harm on the union
    // path (property of the candidate set; on the exact path they measure what
    // union WOULD have done).  applied_mismatch/under: scan-vs-full_rescore of
    // the APPLIED move -- MUST be 0 on the exact path (ew_exact=1), validating
    // the wiring; on the union path applied_under exposes the optimizer's curse.
    REprintf("DRIFT-SCANCHK ew_exact=%d clips=%lld cands=%lld "
             "over=%lld(%.1f%%,max%d) under=%lld(%.1f%%,max%d) "
             "select_flip=%lld(%.1f%%) regret=%.0f env_violation=%lld(%.1f%%) "
             "applied=%lld mism=%lld under=%lld maxdiff=%.1f\n",
             ew_exact ? 1 : 0,
             census.clips, census.candidates,
             census.overcount_cands,
             census.candidates
                 ? 100.0 * census.overcount_cands / census.candidates : 0.0,
             census.overcount_max,
             census.undercount_cands,
             census.candidates
                 ? 100.0 * census.undercount_cands / census.candidates : 0.0,
             census.undercount_max,
             census.select_flip,
             census.clips ? 100.0 * census.select_flip / census.clips : 0.0,
             census.select_regret_sum,
             census.env_violation,
             census.clips ? 100.0 * census.env_violation / census.clips : 0.0,
             census.applied_checked, census.applied_mismatch,
             census.applied_under, census.applied_maxdiff);
  }

  return DriftResult{best_score, params.n_cycles,
                     total_drift_moves, total_tbr_moves};
}

} // namespace ts
