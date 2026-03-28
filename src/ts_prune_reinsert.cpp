#include "ts_prune_reinsert.h"
#include "ts_fitch.h"
#include "ts_tbr.h"
#include "ts_wagner.h"
#include "ts_pool.h"
#include "ts_splits.h"
#include "ts_rng.h"

#include <R.h>
#include <algorithm>
#include <climits>
#include <cstring>
#include <numeric>
#include <vector>

namespace ts {

// -----------------------------------------------------------------------
// Internal helpers
// -----------------------------------------------------------------------

namespace {

// Select which tips to drop.
// Returns a sorted vector of 0-based tip indices.
std::vector<int> select_tips_to_drop(
    const TreeState& tree,
    const DataSet& ds,
    const PruneReinsertParams& params,
    const SplitFrequencyTable* split_freq)
{
  int n_tip = tree.n_tip;
  int k = static_cast<int>(params.drop_fraction * n_tip);
  k = std::max(k, params.min_drop);
  if (params.max_drop > 0) k = std::min(k, params.max_drop);
  // Never drop so many that the reduced tree has < 4 tips
  k = std::min(k, n_tip - 4);
  if (k <= 0) return {};

  std::vector<int> candidates(n_tip);
  std::iota(candidates.begin(), candidates.end(), 0);

  std::vector<int> dropped;
  dropped.reserve(k);

  if (params.selection == PruneSelection::INSTABILITY && split_freq &&
      split_freq->n_trees >= 2)
  {
    // Instability-weighted selection: tips whose parent-edge split is
    // rare in the pool are more likely to be dropped.
    //
    // For each tip, compute the hash of the split at its parent edge
    // (the bipartition separating parent(tip)'s subtree from the rest),
    // then look up its frequency in the pool.  Instability = 1 - freq.

    // Build tip descendants for each internal node (bitset per node)
    int nw = (n_tip + 63) / 64;
    std::vector<uint64_t> node_tips(
        static_cast<size_t>(tree.n_node) * nw, 0ULL);

    // Tips: set their own bit
    for (int t = 0; t < n_tip; ++t) {
      node_tips[static_cast<size_t>(t) * nw + t / 64] |=
          (1ULL << (t % 64));
    }

    // Internal nodes in postorder: union of children
    for (int po_idx = 0;
         po_idx < static_cast<int>(tree.postorder.size()); ++po_idx) {
      int node = tree.postorder[po_idx];
      if (node < n_tip) continue;
      int ni = node - n_tip;
      int lc = tree.left[ni];
      int rc = tree.right[ni];
      uint64_t* dst = &node_tips[static_cast<size_t>(node) * nw];
      const uint64_t* lp = &node_tips[static_cast<size_t>(lc) * nw];
      const uint64_t* rp = &node_tips[static_cast<size_t>(rc) * nw];
      for (int w = 0; w < nw; ++w) dst[w] = lp[w] | rp[w];
    }

    // Compute instability score for each tip.
    // Use the split at the edge above parent(tip).  Canonicalize so
    // that bit 0 (tip 0) is always 0, matching the pool's hashing.
    std::vector<double> instability(n_tip, 1.0);
    std::vector<uint64_t> canon(nw);  // scratch buffer for canonicalization

    for (int t = 0; t < n_tip; ++t) {
      int par = tree.parent[t];
      if (par == tree.n_tip) {
        // Tip is child of root — root splits are trivial
        instability[t] = 1.0;
        continue;
      }

      // The split defined by the edge (grandparent → parent): separates
      // parent(t)'s subtree from the rest of the tree.
      const uint64_t* raw =
          &node_tips[static_cast<size_t>(par) * nw];

      // Canonicalize: ensure bit 0 is 0 (flip if needed)
      bool flip = (raw[0] & 1ULL) != 0;
      for (int w = 0; w < nw; ++w) {
        canon[w] = flip ? ~raw[w] : raw[w];
      }
      // Mask trailing bits in the last word
      int rem = n_tip % 64;
      if (rem > 0) canon[nw - 1] &= (1ULL << rem) - 1;

      uint64_t h = hash_single_split(canon.data(), nw);

      auto it = split_freq->freq.find(h);
      if (it != split_freq->freq.end()) {
        double freq =
            static_cast<double>(it->second) / split_freq->n_trees;
        instability[t] = 1.0 - freq;
      }
      // else: split not in pool → instability = 1.0 (maximally unstable)
    }

    // Weighted sampling without replacement
    std::vector<double> weights = instability;
    for (int j = 0; j < k; ++j) {
      // Compute cumulative weights
      double total = 0.0;
      for (int i = 0; i < n_tip; ++i) total += weights[i];
      if (total <= 0.0) break;

      double r = ts::thread_safe_unif() * total;
      double cum = 0.0;
      int pick = n_tip - 1;
      for (int i = 0; i < n_tip; ++i) {
        cum += weights[i];
        if (cum >= r) { pick = i; break; }
      }
      dropped.push_back(pick);
      weights[pick] = 0.0;  // remove from future picks
    }
  } else {
    // Random selection: Fisher-Yates partial shuffle
    for (int j = 0; j < k; ++j) {
      int idx = j + static_cast<int>(ts::thread_safe_unif() * (n_tip - j));
      if (idx >= n_tip) idx = n_tip - 1;
      std::swap(candidates[j], candidates[idx]);
      dropped.push_back(candidates[j]);
    }
  }

  std::sort(dropped.begin(), dropped.end());
  return dropped;
}

// Extract a pruned tree topology by removing specified tips.
// Returns the edge matrix (1-based, R format) for the reduced tree,
// plus the mapping from reduced tip indices to original tip indices.
//
// The reduced tree has m = n_tip - |dropped| tips, numbered 0..m-1,
// and m-1 internal nodes numbered m..2m-2.
struct PrunedTopology {
  int n_red_tip;
  std::vector<int> edge_parent;  // 1-based
  std::vector<int> edge_child;   // 1-based
  std::vector<int> tip_map;      // red_tip_i → original_tip_index
};

// Recursive helper: prune a subtree, returning the reduced node index
// or -1 if the subtree contains no kept tips.
int prune_subtree(int node, const TreeState& tree,
                  const std::vector<bool>& keep,
                  std::vector<int>& tip_map,
                  std::vector<std::pair<int,int>>& edges,  // 0-based
                  int n_red_tip,
                  int& next_internal)
{
  if (node < tree.n_tip) {
    if (keep[node]) {
      int red_tip = static_cast<int>(tip_map.size());
      tip_map.push_back(node);
      return red_tip;
    }
    return -1;
  }

  int ni = node - tree.n_tip;
  int left_red = prune_subtree(tree.left[ni], tree, keep, tip_map,
                               edges, n_red_tip, next_internal);
  int right_red = prune_subtree(tree.right[ni], tree, keep, tip_map,
                                edges, n_red_tip, next_internal);

  if (left_red < 0 && right_red < 0) return -1;
  if (left_red < 0) return right_red;   // collapse: bypass this node
  if (right_red < 0) return left_red;

  // Both children have kept descendants → create reduced internal node
  int red_int = next_internal++;
  edges.push_back({red_int, left_red});
  edges.push_back({red_int, right_red});
  return red_int;
}

PrunedTopology extract_pruned_topology(
    const TreeState& tree,
    const std::vector<int>& dropped)
{
  int n_tip = tree.n_tip;
  std::vector<bool> keep(n_tip, true);
  for (int t : dropped) keep[t] = false;

  int m = n_tip - static_cast<int>(dropped.size());

  PrunedTopology result;
  result.n_red_tip = m;
  result.tip_map.reserve(m);

  std::vector<std::pair<int,int>> edges_0;  // 0-based edges
  edges_0.reserve(2 * (m - 1));

  int next_internal = m;  // first reduced internal node
  int root_node = prune_subtree(
      tree.n_tip, tree, keep, result.tip_map, edges_0, m, next_internal);

  // The root of the reduced tree should be 'm' (the TreeState convention).
  // If prune_subtree returned a different root, renumber.
  // Actually, the root is the last internal node created by the postorder
  // recursion. It might not be m. We need to relabel so that root = m.
  if (root_node != m) {
    // Swap root_node and m in all edges
    for (auto& e : edges_0) {
      if (e.first == root_node) e.first = m;
      else if (e.first == m) e.first = root_node;
      if (e.second == root_node) e.second = m;
      else if (e.second == m) e.second = root_node;
    }
    // Also swap in tip_map if one of them is a tip (shouldn't be, both >= m)
    // root_node and m are both >= m so they are internal nodes — no tip_map fix needed
  }

  // Convert to 1-based edge matrix
  result.edge_parent.resize(edges_0.size());
  result.edge_child.resize(edges_0.size());
  for (size_t i = 0; i < edges_0.size(); ++i) {
    result.edge_parent[i] = edges_0[i].first + 1;
    result.edge_child[i] = edges_0[i].second + 1;
  }

  return result;
}

// Build a reduced DataSet containing only the kept tips.
// Copies scoring metadata verbatim; only tip_states is subsetted.
DataSet build_reduced_dataset(const DataSet& ds,
                              const std::vector<int>& tip_map)
{
  int m = static_cast<int>(tip_map.size());
  DataSet red;

  // Copy all scoring metadata
  red.n_tips = m;
  red.n_blocks = ds.n_blocks;
  red.total_words = ds.total_words;
  red.blocks = ds.blocks;
  red.block_word_offset = ds.block_word_offset;
  red.flat_blocks = ds.flat_blocks;
  red.all_weight_one = ds.all_weight_one;
  red.n_patterns = ds.n_patterns;
  red.min_steps = ds.min_steps;
  red.pattern_freq = ds.pattern_freq;
  red.concavity = ds.concavity;
  red.eff_k = ds.eff_k;
  red.phi = ds.phi;
  red.scoring_mode = ds.scoring_mode;
  red.ew_offset = ds.ew_offset;
  red.precomputed_steps = ds.precomputed_steps;
  red.inapp_state = ds.inapp_state;

  // Subset tip_states
  int tw = ds.total_words;
  red.tip_states.resize(static_cast<size_t>(m) * tw);
  for (int i = 0; i < m; ++i) {
    int orig = tip_map[i];
    std::memcpy(&red.tip_states[static_cast<size_t>(i) * tw],
                &ds.tip_states[static_cast<size_t>(orig) * tw],
                tw * sizeof(uint64_t));
  }

  return red;
}

// Expand the optimized reduced tree back into a full-size TreeState,
// then Wagner-insert each dropped tip at its best position.
//
// The reduced tree's tips are mapped to original tip indices via tip_map.
// After expansion, the full tree has all n_tip tips, with the backbone
// from the reduced tree and the dropped tips greedily reinserted.
void expand_and_reinsert(
    TreeState& tree,
    const DataSet& ds,
    const TreeState& red_tree,
    const std::vector<int>& tip_map,
    const std::vector<int>& dropped,
    ConstraintData* cd)
{
  int n_tip = ds.n_tips;
  int m = red_tree.n_tip;
  int tw = ds.total_words;

  // 1. Allocate full-size tree
  init_wagner_state(tree, ds);

  // 2. Map reduced tree topology into full tree's index space.
  //    Reduced tips 0..m-1 → original tips tip_map[0..m-1]
  //    Reduced internal m..2m-2 → full internal n_tip..n_tip+m-2
  //    Root: reduced m → full n_tip

  // Build node mapping: reduced_node → full_node
  std::vector<int> node_map(red_tree.n_node);
  for (int i = 0; i < m; ++i) {
    node_map[i] = tip_map[i];
  }
  for (int i = m; i < red_tree.n_node; ++i) {
    node_map[i] = n_tip + (i - m);
  }

  // Set topology
  for (int ri = 0; ri < red_tree.n_internal; ++ri) {
    int red_node = m + ri;
    int full_node = node_map[red_node];
    int full_ni = full_node - n_tip;

    int red_lc = red_tree.left[ri];
    int red_rc = red_tree.right[ri];
    int full_lc = node_map[red_lc];
    int full_rc = node_map[red_rc];

    tree.left[full_ni] = full_lc;
    tree.right[full_ni] = full_rc;
    tree.parent[full_lc] = full_node;
    tree.parent[full_rc] = full_node;
  }
  tree.parent[n_tip] = n_tip;  // root

  // 3. Build postorder and score the backbone
  tree.build_postorder();
  score_tree(tree, ds);  // sets prelim/final_ for all backbone nodes

  // 4. Wagner-insert each dropped tip.
  // Internal nodes for new insertions: n_tip + (m-1), n_tip + m, ...
  int next_internal = n_tip + (m - 1);

  // Randomize reinsertion order for stochasticity
  std::vector<int> reinsert_order = dropped;
  for (int i = static_cast<int>(reinsert_order.size()) - 1; i > 0; --i) {
    int j = static_cast<int>(ts::thread_safe_unif() * (i + 1));
    if (j > i) j = i;
    std::swap(reinsert_order[i], reinsert_order[j]);
  }

  // Pre-allocate DFS stack once; cleared per tip.
  // wagner_incremental_rescore walks up via tree.parent (no postorder needed),
  // and the DFS below uses tree.left/right only — so build_postorder() is NOT
  // required between insertions.  One call after the loop suffices.
  std::vector<int> stack;
  stack.reserve(tree.n_node);

  for (int tip : reinsert_order) {
    int new_internal = next_internal++;

    const uint64_t* tip_prelim =
        &ds.tip_states[static_cast<size_t>(tip) * tw];

    // Find best insertion edge via DFS from root
    int best_above = -1, best_below = -1;
    int best_extra = INT_MAX;

    stack.clear();
    stack.push_back(n_tip);  // root

    while (!stack.empty()) {
      int node = stack.back();
      stack.pop_back();

      if (node < n_tip) continue;

      int ni = node - n_tip;
      int lc = tree.left[ni];
      int rc = tree.right[ni];

      // Skip unconnected nodes (left == -1 means unused internal)
      if (lc < 0 || rc < 0) continue;

      // Evaluate edge (node, lc)
      int extra = fitch_indirect_length_bounded(
          tip_prelim, tree, ds, node, lc, best_extra);
      if (extra < best_extra) {
        best_extra = extra;
        best_above = node;
        best_below = lc;
      }

      // Evaluate edge (node, rc)
      extra = fitch_indirect_length_bounded(
          tip_prelim, tree, ds, node, rc, best_extra);
      if (extra < best_extra) {
        best_extra = extra;
        best_above = node;
        best_below = rc;
      }

      stack.push_back(lc);
      stack.push_back(rc);
    }

    // Fallback: if no valid edge found, insert at root edge
    if (best_above < 0 || best_below < 0) {
      best_above = n_tip;
      best_below = tree.left[0];
    }

    insert_tip_at_edge(tree, tip, new_internal, best_above, best_below);
    wagner_incremental_rescore(tree, ds, new_internal);
  }

  // Rebuild postorder once after all insertions — required by the subsequent
  // TBR polish in prune_reinsert_search (and by score_tree).
  tree.build_postorder();
}

} // anonymous namespace

// -----------------------------------------------------------------------
// Public API
// -----------------------------------------------------------------------

PruneReinsertResult prune_reinsert_search(
    TreeState& tree,
    DataSet& ds,
    const PruneReinsertParams& params,
    ConstraintData* cd,
    const SplitFrequencyTable* split_freq,
    std::function<bool()> check_timeout)
{
  PruneReinsertResult result;
  result.n_improvements = 0;

  // build_reduced_dataset() copies EW/IW/XPIWE metadata completely.
  // For PROFILE (info_amounts), HSJ (hierarchy_blocks/tip_labels), and
  // XFORM (sankoff_* fields), the reduced dataset would be missing
  // scoring-specific fields, causing incorrect reduced-tree scores (T-275).
  // Guard: skip prune-reinsert for unsupported scoring modes.
  if (ds.scoring_mode == ScoringMode::PROFILE ||
      ds.scoring_mode == ScoringMode::HSJ ||
      ds.scoring_mode == ScoringMode::XFORM) {
    return result;
  }

  double current_score = score_tree(tree, ds);
  result.best_score = current_score;

  for (int cyc = 0; cyc < params.n_cycles; ++cyc) {
    if (check_timeout && check_timeout()) break;

    // 1. Select tips to drop
    std::vector<int> dropped = select_tips_to_drop(
        tree, ds, params, split_freq);
    if (dropped.empty()) break;

    // Save current tree for revert
    TreeState backup = tree;

    // 2. Extract pruned topology
    PrunedTopology pt = extract_pruned_topology(tree, dropped);

    // 3. Build reduced dataset
    DataSet red_ds = build_reduced_dataset(ds, pt.tip_map);

    // 4. Init reduced tree and run TBR
    TreeState red_tree;
    red_tree.init_from_edge(pt.edge_parent.data(), pt.edge_child.data(),
                            static_cast<int>(pt.edge_parent.size()),
                            red_ds);
    {
      TBRParams tp;
      tp.accept_equal = false;
      tp.max_accepted_changes = params.tbr_max_moves;
      tp.max_hits = params.tbr_max_hits;
      tp.tabu_size = params.tabu_size;
      tbr_search(red_tree, red_ds, tp, nullptr, nullptr, nullptr,
                 check_timeout);
    }

    if (check_timeout && check_timeout()) break;

    // 5. Expand reduced tree and reinsert dropped tips
    expand_and_reinsert(tree, ds, red_tree, pt.tip_map, dropped, cd);

    // Re-sync constraint metadata: expand_and_reinsert completely rebuilds
    // the topology (init_wagner_state + node mapping + tip reinsertion),
    // so cd->constraint_node and DFS timestamps are stale.
    // Same bug class as T-278 (TBR), T-279 (drift), F-015 (ratchet),
    // F-016 (NNI-perturb).
    if (cd) update_constraint(tree, *cd);

    // 6. TBR polish on full tree
    {
      TBRParams tp;
      tp.accept_equal = false;
      tp.max_accepted_changes = 0;  // converge
      tp.max_hits = params.tbr_max_hits;
      tp.tabu_size = params.tabu_size;
      tbr_search(tree, ds, tp, cd, nullptr, nullptr, check_timeout);
    }

    // 7. Accept or revert
    double new_score = score_tree(tree, ds);
    if (new_score < current_score - 1e-10) {
      current_score = new_score;
      result.best_score = new_score;
      ++result.n_improvements;
    } else {
      tree = backup;  // revert
      // Re-sync constraint metadata after topology revert.
      // Same bug class as F-015 (ratchet), F-016 (NNI-perturb).
      if (cd) update_constraint(tree, *cd);
    }
  }

  return result;
}

} // namespace ts
