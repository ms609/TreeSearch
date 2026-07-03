#include "ts_wagner.h"
#include "ts_constraint.h"
#include "ts_fitch.h"
#include "ts_rng.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <climits>
#include <vector>
#include <R.h>

namespace ts {

// Guard: Wagner tree construction requires at least 3 taxa
static void check_wagner_precondition(int n_tip) {
  if (n_tip < 3) {
    Rf_error("Wagner tree requires at least 3 taxa (got %d)", n_tip);
  }
}

// Initialize a full-sized TreeState for Wagner construction.
// Allocates all arrays but does not set topology (caller does that).
void init_wagner_state(TreeState& tree, const DataSet& ds) {
  tree.n_tip = ds.n_tips;
  tree.n_internal = ds.n_tips - 1;
  tree.n_node = 2 * ds.n_tips - 1;
  tree.total_words = ds.total_words;
  tree.n_blocks = ds.n_blocks;

  tree.parent.assign(tree.n_node, -1);
  tree.left.assign(tree.n_internal, -1);
  tree.right.assign(tree.n_internal, -1);

  size_t state_size = static_cast<size_t>(tree.n_node) * tree.total_words;
  tree.prelim.assign(state_size, 0ULL);
  tree.final_.assign(state_size, 0ULL);
  tree.down2.assign(state_size, 0ULL);
  tree.subtree_actives.assign(state_size, 0ULL);
  tree.local_cost.assign(
      static_cast<size_t>(tree.n_node) * tree.n_blocks, 0ULL);

  tree.load_tip_states(ds);
}

/* Build the initial 3-taxon tree:
 *
 *        root (= n_tip)
 *       /    \
 *    int1     t2
 *  (n_tip+1)
 *   /    \
 * t0      t1
 */
void build_three_taxon_tree(TreeState& tree, int t0, int t1, int t2) {
  int root = tree.n_tip;
  int int1 = tree.n_tip + 1;

  tree.parent[root] = root;
  tree.parent[int1] = root;
  tree.parent[t2] = root;
  tree.parent[t0] = int1;
  tree.parent[t1] = int1;

  tree.left[0] = int1;  // root's left
  tree.right[0] = t2;   // root's right
  tree.left[1] = t0;    // int1's left
  tree.right[1] = t1;   // int1's right
}

/* Insert tip at edge (above, below), creating new_internal between them.
 *
 * Before:     above          After:     above
 *               |                        |
 *             below                   new_internal
 *                                     /        \
 *                                   tip       below
 */
void insert_tip_at_edge(TreeState& tree, int tip, int new_internal,
                        int above, int below) {
  int ni = new_internal - tree.n_tip;

  tree.parent[new_internal] = above;
  tree.left[ni] = tip;
  tree.right[ni] = below;
  tree.parent[tip] = new_internal;
  tree.parent[below] = new_internal;

  // Update above's child pointer
  int ai = above - tree.n_tip;
  if (tree.left[ai] == below) {
    tree.left[ai] = new_internal;
  } else {
    tree.right[ai] = new_internal;
  }
}

// Incremental two-pass Fitch scoring after Wagner insertion.
//
// After inserting a tip at edge (above, below) with new_internal between them,
// recompute downpass states from new_internal to root, then uppass from root
// back down through changed nodes. Returns the score delta (positive = score
// increased, which it always will during Wagner construction).
//
// Downpass: O(depth×C) — walk from new_internal to root, stop when prelim
// stabilizes.
// Uppass: DFS from root with early termination — O(affected_region × C),
// typically much less than a full uppass.
int wagner_incremental_rescore(TreeState& tree, const DataSet& ds,
                               int new_internal) {
  int n_tip = tree.n_tip;
  int tw = tree.total_words;
  int nb = tree.n_blocks;
  int root = n_tip;
  int score_delta = 0;

  // --- Phase 1: Incremental downpass (new_internal → root) ---
  // new_internal's children are the inserted tip and `below` (which already
  // has valid prelim states). Compute prelim at new_internal, then walk
  // upward recomputing ancestors whose children changed.

  int node = new_internal;
  while (true) {
    int ni = node - n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    size_t node_base = static_cast<size_t>(node) * tw;
    size_t lb = static_cast<size_t>(lc) * tw;
    size_t rb = static_cast<size_t>(rc) * tw;
    bool changed = false;

    for (int b = 0; b < nb; ++b) {
      const CharBlock& blk = ds.blocks[b];
      int offset = ds.block_word_offset[b];

      const uint64_t* left_state = &tree.prelim[lb + offset];
      const uint64_t* right_state = &tree.prelim[rb + offset];
      uint64_t* node_state = &tree.prelim[node_base + offset];

      // Subtract old local cost for this node (0 for newly created nodes)
      size_t cost_idx = static_cast<size_t>(node) * nb + b;
      uint64_t old_cost = tree.local_cost[cost_idx];
      int old_nu = popcount64(old_cost);
      if (blk.upweight_mask) old_nu += popcount64(old_cost & blk.upweight_mask);
      score_delta -= blk.weight * old_nu;

      // Fitch downpass: intersection/union
      uint64_t any_intersect = 0;
      for (int s = 0; s < blk.n_states; ++s) {
        any_intersect |= (left_state[s] & right_state[s]);
      }
      uint64_t needs_union = ~any_intersect & blk.active_mask;

      int new_nu = popcount64(needs_union);
      if (blk.upweight_mask) new_nu += popcount64(needs_union & blk.upweight_mask);
      score_delta += blk.weight * new_nu;

      tree.local_cost[cost_idx] = needs_union;

      for (int s = 0; s < blk.n_states; ++s) {
        uint64_t isect = left_state[s] & right_state[s];
        uint64_t uni = left_state[s] | right_state[s];
        uint64_t new_val = (isect & any_intersect) | (uni & needs_union);
        if (new_val != node_state[s]) changed = true;
        node_state[s] = new_val;
      }
    }

    if (node == root) break;
    // Early termination: if prelim didn't change, ancestors are unaffected
    if (!changed) break;
    node = tree.parent[node];
  }

  // --- Phase 2: Uppass from root with early termination ---
  // Set root final = prelim, then DFS through internal nodes, computing
  // final from parent's final + own prelim. Skip subtrees where final
  // didn't change.

  size_t root_base = static_cast<size_t>(root) * tw;
  for (int w = 0; w < tw; ++w) {
    tree.final_[root_base + w] = tree.prelim[root_base + w];
  }

  // DFS from root's children. Each internal node whose final_ changes
  // has its children pushed for processing.
  std::vector<int> up_stack;
  {
    int ri = root - n_tip;
    int lc = tree.left[ri];
    int rc = tree.right[ri];
    if (lc >= 0) up_stack.push_back(lc);
    if (rc >= 0) up_stack.push_back(rc);
  }

  while (!up_stack.empty()) {
    int n = up_stack.back();
    up_stack.pop_back();
    if (n < n_tip) continue;  // tip: final = prelim, already correct

    int anc = tree.parent[n];
    size_t n_base = static_cast<size_t>(n) * tw;
    size_t a_base = static_cast<size_t>(anc) * tw;
    bool changed = false;

    for (int b = 0; b < nb; ++b) {
      const CharBlock& blk = ds.blocks[b];
      int offset = ds.block_word_offset[b];

      uint64_t any_isect = 0;
      for (int s = 0; s < blk.n_states; ++s) {
        any_isect |= (tree.final_[a_base + offset + s]
                    & tree.prelim[n_base + offset + s]);
      }
      uint64_t no_isect = ~any_isect & blk.active_mask;

      for (int s = 0; s < blk.n_states; ++s) {
        uint64_t isect = tree.final_[a_base + offset + s]
                       & tree.prelim[n_base + offset + s];
        uint64_t new_val = (isect & any_isect)
                         | (tree.prelim[n_base + offset + s] & no_isect);
        if (new_val != tree.final_[n_base + offset + s]) changed = true;
        tree.final_[n_base + offset + s] = new_val;
      }
    }

    // Propagate to children only if final changed
    if (changed) {
      int ni_idx = n - n_tip;
      int lc = tree.left[ni_idx];
      int rc = tree.right[ni_idx];
      if (lc >= 0) up_stack.push_back(lc);
      if (rc >= 0) up_stack.push_back(rc);
    }
  }

  return score_delta;
}

// Expose the is_ancestor_or_equal helper for constraint checking.
static inline bool is_ancestor_or_equal(
    int u, int v,
    const std::vector<int>& entry, const std::vector<int>& exit)
{
  return entry[u] <= entry[v] && exit[u] >= exit[v];
}

// Wagner-specific constraint node mapping.
// During incremental construction, the full split may not be present yet.
// Instead of requiring an exact match, find the LCA of added inside tips:
// the smallest internal node whose subtree contains all added inside tips.
// This correctly constrains placement even when only some inside tips
// have been added so far.
static void wagner_map_constraint_nodes(
    const TreeState& tree, ConstraintData& cd,
    const std::vector<uint64_t>& added_tips)
{
  if (!cd.active) return;

  int n_tip = tree.n_tip;
  int nw = cd.n_words;

  // Build per-node subtree tip bitmasks via postorder traversal.
  std::vector<uint64_t> node_tips(
      static_cast<size_t>(tree.n_node) * nw, 0ULL);

  for (int t = 0; t < n_tip; ++t) {
    int w = t / 64;
    int b = t % 64;
    node_tips[static_cast<size_t>(t) * nw + w] = (1ULL << b);
  }

  for (int node : tree.postorder) {
    int ni = node - n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    uint64_t* nd = &node_tips[static_cast<size_t>(node) * nw];
    const uint64_t* lt = &node_tips[static_cast<size_t>(lc) * nw];
    const uint64_t* rt = &node_tips[static_cast<size_t>(rc) * nw];
    for (int w = 0; w < nw; ++w) {
      nd[w] = lt[w] | rt[w];
    }
  }

  // For each split, find the smallest node (tip or internal) whose subtree
  // contains all added inside tips.  Including tips handles the case where
  // only 1 inside tip has been added: the constraint node is that tip itself,
  // forcing the next inside tip to be placed adjacent to it.
  for (int s = 0; s < cd.n_splits; ++s) {
    const uint64_t* split =
        &cd.split_tips[static_cast<size_t>(s) * nw];

    // Compute which inside tips have been added
    bool any_added = false;
    std::vector<uint64_t> needed(nw, 0ULL);
    for (int w = 0; w < nw; ++w) {
      needed[w] = split[w] & added_tips[w];
      if (needed[w]) any_added = true;
    }
    if (!any_added) {
      cd.constraint_node[s] = -1;
      continue;
    }

    // Check tips first: if exactly 1 inside tip is added, use it directly.
    // (A tip's "subtree" is just itself — size 1.)
    int best_node = -1;
    int best_size = tree.n_node + 1;

    for (int t = 0; t < n_tip; ++t) {
      int w = t / 64;
      int b = t % 64;
      if (!((added_tips[w] >> b) & 1)) continue;  // tip not added yet
      // Does this tip's bitmask contain all needed?
      // A tip's bitmask is just itself, so this is true only if needed = {t}.
      const uint64_t* nd = &node_tips[static_cast<size_t>(t) * nw];
      bool superset = true;
      for (int w2 = 0; w2 < nw; ++w2) {
        if ((nd[w2] & needed[w2]) != needed[w2]) {
          superset = false;
          break;
        }
      }
      if (superset && 1 < best_size) {
        best_size = 1;
        best_node = t;
      }
    }

    // Then check internal nodes for a smaller clade.
    for (int node : tree.postorder) {
      const uint64_t* nd = &node_tips[static_cast<size_t>(node) * nw];
      bool superset = true;
      for (int w = 0; w < nw; ++w) {
        if ((nd[w] & needed[w]) != needed[w]) {
          superset = false;
          break;
        }
      }
      if (superset) {
        int sz = 0;
        for (int w = 0; w < nw; ++w) {
          sz += popcount64(nd[w]);
        }
        if (sz < best_size) {
          best_size = sz;
          best_node = node;
        }
      }
    }
    cd.constraint_node[s] = best_node;
  }
}

// Wagner-specific: remap constraint nodes using LCA, then recompute DFS.
static void wagner_update_constraint(
    const TreeState& tree, ConstraintData& cd,
    const std::vector<uint64_t>& added_tips)
{
  if (!cd.active) return;
  wagner_map_constraint_nodes(tree, cd, added_tips);
  compute_dfs_timestamps(tree, cd);
}

// Check if an edge (above, below) is legal under the constraint.
// `added_tips`: bitmask of tips already in the tree.
// For each constraint split where both sides have previously-added tips:
//   - If the new tip is "inside" the split, the insertion must be inside
//     the LCA clade of already-added inside tips.
//   - If the new tip is "outside", the insertion must be outside.
// Uses DFS timestamps for O(1) descendant test per constraint.
static bool wagner_edge_violates_constraint(
    const TreeState& tree, int below, int tip,
    const ConstraintData& cd,
    const std::vector<uint64_t>& added_tips)
{
  for (int s = 0; s < cd.n_splits; ++s) {
    const uint64_t* split =
        &cd.split_tips[static_cast<size_t>(s) * cd.n_words];

    // Is the new tip inside or outside this split?
    int tw = tip / 64;
    int tb = tip % 64;
    bool tip_inside = (split[tw] >> tb) & 1;

    // Split constrains placement when the opposite side of the new tip
    // has at least one previously-added tip.  An inside tip is only
    // constrained when outside tips already exist (and vice versa).
    bool has_prev_inside = false, has_prev_outside = false;
    for (int w = 0; w < cd.n_words; ++w) {
      uint64_t prev = added_tips[w];
      if (prev & split[w]) has_prev_inside = true;
      uint64_t outside_mask = ~split[w];
      if (w == cd.n_words - 1) {
        int rem = tree.n_tip % 64;
        if (rem > 0) outside_mask &= (1ULL << rem) - 1;
      }
      if (prev & outside_mask) has_prev_outside = true;
    }
    if (tip_inside && !has_prev_outside) continue;
    if (!tip_inside && !has_prev_inside) continue;

    // The constraint is active. constraint_node[s] is the LCA of
    // added inside tips (set by wagner_map_constraint_nodes).
    int cn = cd.constraint_node[s];
    if (cn < 0) continue;

    // If the LCA is the root, inside tips span both sides of the root.
    // The inside/outside distinction is meaningless (every node is a
    // "descendant" of root), so skip this constraint for this insertion.
    if (cn == tree.n_tip) continue;

    bool below_inside =
        is_ancestor_or_equal(cn, below, cd.dfs_entry, cd.dfs_exit);

    if (tip_inside && !below_inside) return true;
    // Exclude the boundary edge (above_cn, cn): inserting an outside tip there
    // makes it sibling of the entire constraint clade, which preserves
    // monophyly.  Only reject if the tip would go *strictly inside* the clade.
    if (!tip_inside && below_inside && below != cn) return true;
  }
  return false;
}

WagnerResult wagner_tree(TreeState& tree, const DataSet& ds,
                         const std::vector<int>& addition_order,
                         ConstraintData* cd) {
  int n_tip = ds.n_tips;
  check_wagner_precondition(n_tip);

  bool constrained = cd && cd->active;

  // Validate or create addition order
  std::vector<int> order;
  if (addition_order.empty()) {
    order.resize(n_tip);
    std::iota(order.begin(), order.end(), 0);
  } else {
    order = addition_order;
  }

  // Initialize full-sized TreeState
  init_wagner_state(tree, ds);

  // Build initial 3-taxon tree and do full two-pass scoring.
  // This sets up prelim/final_ for indirect length evaluation of the first
  // insertion candidate. Subsequent insertions use incremental scoring.
  build_three_taxon_tree(tree, order[0], order[1], order[2]);
  tree.build_postorder();
  int ew_score = fitch_score(tree, ds);

  // Track which tips have been added (bitmask)
  int n_words = constrained ? cd->n_words : 0;
  std::vector<uint64_t> added_tips(n_words, 0ULL);
  if (constrained) {
    for (int j = 0; j < 3; ++j) {
      int t = order[j];
      added_tips[t / 64] |= (1ULL << (t % 64));
    }
    wagner_update_constraint(tree, *cd, added_tips);
  }

  // Set if the constraint filter ever exhausted every legal edge and we fell
  // back to the unconstrained root edge (see the guard below); warned once
  // after construction so a constraint-violating tree is never returned mutely.
  bool constraint_fallback = false;

  // Scratch buffer for the exact directional insertion edge sets, recomputed
  // each step from the current downpass (prelim, kept current by the
  // incremental rescore below).  edge_set_up / edge_set_pre are caller-owned
  // scratch reused across the insertion loop (size-ensured, non-zeroing) so
  // compute_insertion_edge_sets avoids per-step allocation and zero-fill.
  std::vector<uint64_t> edge_set;
  std::vector<uint64_t> edge_set_up;
  std::vector<int> edge_set_pre;

  // Add remaining taxa one at a time
  for (int i = 3; i < n_tip; ++i) {
    int tip = order[i];
    int new_internal = n_tip + i - 1;

    const uint64_t* tip_prelim =
        &ds.tip_states[static_cast<size_t>(tip) * ds.total_words];

    // Exact insertion cost: edge_set[D] = combine(prelim[D], up[D]).  Replaces
    // the union-of-finals (final_[node] | final_[child]) approximation that
    // undercut insertion cost and produced ~+30% Wagner trees.
    compute_insertion_edge_sets(tree, ds, edge_set, edge_set_up, edge_set_pre);
    const int tw = tree.total_words;

    // Find best insertion edge via DFS from root
    int best_above = -1, best_below = -1;
    int best_extra = INT_MAX;

    std::vector<int> stack;
    stack.push_back(n_tip);  // root

    while (!stack.empty()) {
      int node = stack.back();
      stack.pop_back();

      if (node < n_tip) continue;  // tip — no children to enumerate

      int ni = node - n_tip;
      int lc = tree.left[ni];
      int rc = tree.right[ni];

      // Evaluate edge (node, lc)
      if (!constrained ||
          !wagner_edge_violates_constraint(tree, lc, tip, *cd,
                                            added_tips)) {
        int extra = fitch_indirect_length_cached(
            tip_prelim, &edge_set[static_cast<size_t>(lc) * tw],
            ds, best_extra);
        if (extra < best_extra) {
          best_extra = extra;
          best_above = node;
          best_below = lc;
        }
      }

      // Evaluate edge (node, rc)
      if (!constrained ||
          !wagner_edge_violates_constraint(tree, rc, tip, *cd,
                                            added_tips)) {
        int extra = fitch_indirect_length_cached(
            tip_prelim, &edge_set[static_cast<size_t>(rc) * tw],
            ds, best_extra);
        if (extra < best_extra) {
          best_extra = extra;
          best_above = node;
          best_below = rc;
        }
      }

      stack.push_back(lc);
      stack.push_back(rc);
    }

    // Guard: if constraint filtered every edge, fall back to root edge.
    // This should not happen after the cn==root fix, but protects against
    // any remaining edge case in constraint logic. The fallback edge is not
    // constraint-checked, so flag it and warn after construction.
    if (best_above < 0 || best_below < 0) {
      best_above = n_tip;
      best_below = tree.left[0];
      if (constrained) constraint_fallback = true;
    }

    // Insert tip at the best edge
    insert_tip_at_edge(tree, tip, new_internal, best_above, best_below);

    // Incremental rescore: update only the insertion-to-root path
    int delta = wagner_incremental_rescore(tree, ds, new_internal);
    ew_score += delta;

    if (constrained) {
      added_tips[tip / 64] |= (1ULL << (tip % 64));
      // Rebuild postorder before updating constraint mapping — the previous
      // postorder is stale (doesn't include newly created internal nodes).
      tree.build_postorder();
      wagner_update_constraint(tree, *cd, added_tips);
    }
  }

  // Build postorder (needed by subsequent TBR search) and compute final score.
  // The incremental EW score tracked during construction is exact for standard
  // Fitch; for NA datasets or IW, score_tree gives the authoritative result.
  tree.build_postorder();
  double score = score_tree(tree, ds);

  if (constraint_fallback) {
    Rf_warning(
      "AdditionTree(): constraint could not be honoured for at least one "
      "taxon insertion; the returned tree may violate the constraint. "
      "Consider supplying a `sequence` that adds constrained taxa earlier.");
  }

  WagnerResult result;
  result.score = score;
  return result;
}

// ---------------------------------------------------------------------------
// Biased addition-order scoring
// ---------------------------------------------------------------------------

// Goloboff "informative" score (Goloboff 2014, "Hide and vanish").
// score[t] = number of characters for which tip t has a specific
// (non-ambiguous) state.  Ambiguous = all n_states bits set for that
// character.  Characters are counted once per block position regardless
// of per-pattern frequency weighting, so the score reflects the number of
// independently coded characters.
std::vector<double> wagner_goloboff_scores(const DataSet& ds) {
  int n_tip = ds.n_tips;
  int tw    = ds.total_words;
  std::vector<double> scores(n_tip, 0.0);

  for (int t = 0; t < n_tip; ++t) {
    double score = 0.0;
    const uint64_t* tip_base = &ds.tip_states[static_cast<size_t>(t) * tw];

    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      int offset = ds.block_word_offset[b];

      // tip_ambiguous: bit c is set when tip has ALL n_states for char c.
      // Compute as the AND of all per-state words masked to active chars.
      uint64_t tip_ambiguous = blk.active_mask;
      for (int s = 0; s < blk.n_states; ++s) {
        tip_ambiguous &= tip_base[offset + s];
      }
      // Non-ambiguous characters = active chars minus those that are ambiguous.
      score += blk.n_chars - popcount64(tip_ambiguous);
    }
    scores[t] = score;
  }
  return scores;
}

// Information-theoretic (entropy) score.
// score[t] = Σ_c (n_states_c - |state_set of t at c|)
//          = Σ_b (n_states_b * n_chars_b - Σ_s popcount(tip_word[s] & active))
//
// A taxon with entirely specific single-state codings scores highest.
// A fully ambiguous taxon scores 0.  The maximum is (n_states-1)*n_chars.
std::vector<double> wagner_entropy_scores(const DataSet& ds) {
  int n_tip = ds.n_tips;
  int tw    = ds.total_words;
  std::vector<double> scores(n_tip, 0.0);

  for (int t = 0; t < n_tip; ++t) {
    double score = 0.0;
    const uint64_t* tip_base = &ds.tip_states[static_cast<size_t>(t) * tw];

    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      int offset = ds.block_word_offset[b];

      int total_state_bits = 0;
      for (int s = 0; s < blk.n_states; ++s) {
        total_state_bits += popcount64(tip_base[offset + s] & blk.active_mask);
      }
      // Each character contributes (n_states - actual_state_count).
      // Ambiguous chars contribute 0; single-state chars contribute n_states-1.
      score += blk.n_states * blk.n_chars - total_state_bits;
    }
    scores[t] = score;
  }
  return scores;
}

// Softmax-weighted sampling without replacement.
// Scores are normalised to [0, 1] before applying temperature so that
// the parameter is dataset-independent.  temperature == 0 → greedy argmax.
// Returns a permutation of 0..n_tip-1 in the sampled addition order.
static std::vector<int> softmax_sample_order(
    const std::vector<double>& scores,
    double temperature)
{
  int n = static_cast<int>(scores.size());
  std::vector<int> order;
  order.reserve(n);

  if (temperature <= 0.0) {
    // Greedy: sort descending by score, break ties arbitrarily
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(),
      [&scores](int a, int b){ return scores[a] > scores[b]; });
    return idx;
  }

  // Normalise scores to [0, 1] so temperature is dataset-independent
  double mn = *std::min_element(scores.begin(), scores.end());
  double mx = *std::max_element(scores.begin(), scores.end());
  double rng = mx - mn;

  std::vector<double> norm(n);
  if (rng < 1e-12) {
    // All scores identical → uniform random; temperature irrelevant
    std::fill(norm.begin(), norm.end(), 1.0);
  } else {
    for (int i = 0; i < n; ++i) {
      norm[i] = (scores[i] - mn) / rng;  // in [0, 1]
    }
  }

  // Softmax weights: w[i] = exp(norm[i] / temperature).
  // Subtract 1.0 (max of norm) for numerical stability.
  std::vector<double> w(n);
  for (int i = 0; i < n; ++i) {
    w[i] = std::exp((norm[i] - 1.0) / temperature);
  }

  // Sample without replacement: at each step draw from remaining
  // taxa proportional to their weights.
  std::vector<bool> used(n, false);

  ts::rng_state_begin();
  for (int step = 0; step < n; ++step) {
    double total = 0.0;
    for (int i = 0; i < n; ++i) {
      if (!used[i]) total += w[i];
    }

    double r = ts::thread_safe_unif() * total;
    double cum = 0.0;
    int chosen = -1;
    for (int i = 0; i < n; ++i) {
      if (used[i]) continue;
      cum += w[i];
      if (r <= cum) { chosen = i; break; }
    }
    // Numerical safety: pick last unused taxon
    if (chosen < 0) {
      for (int i = n - 1; i >= 0; --i) {
        if (!used[i]) { chosen = i; break; }
      }
    }
    used[chosen] = true;
    order.push_back(chosen);
  }
  ts::rng_state_end();

  return order;
}

WagnerResult biased_wagner_tree(TreeState& tree, const DataSet& ds,
                                const BiasedWagnerParams& params,
                                ConstraintData* cd) {
  if (params.bias == WagnerBias::RANDOM) {
    return random_wagner_tree(tree, ds, cd);
  }

  const std::vector<double> scores =
    (params.bias == WagnerBias::GOLOBOFF)
      ? wagner_goloboff_scores(ds)
      : wagner_entropy_scores(ds);

  const std::vector<int> order =
    softmax_sample_order(scores, params.temperature);

  WagnerResult result = wagner_tree(tree, ds, order, cd);

  // Constraint post-hoc check + retry (mirrors random_wagner_tree)
  bool constrained = cd && cd->active && cd->has_posthoc;
  if (constrained) {
    for (int attempt = 1; attempt < 100; ++attempt) {
      if (!violates_constraint_posthoc(tree, *cd)) break;
      const std::vector<int> retry_order =
        softmax_sample_order(scores, params.temperature);
      result = wagner_tree(tree, ds, retry_order, cd);
    }
  }

  return result;
}

// ---------------------------------------------------------------------------

WagnerResult random_wagner_tree(TreeState& tree, const DataSet& ds,
                                ConstraintData* cd) {
  int n_tips = ds.n_tips;
  std::vector<int> order(n_tips);
  std::iota(order.begin(), order.end(), 0);

  bool constrained = cd && cd->active && cd->has_posthoc;

  // Fisher-Yates shuffle
  ts::rng_state_begin();
  for (int i = n_tips - 1; i > 0; --i) {
    int j = static_cast<int>(ts::thread_safe_unif() * (i + 1));
    if (j > i) j = i;
    std::swap(order[i], order[j]);
  }
  ts::rng_state_end();

  WagnerResult result = wagner_tree(tree, ds, order, cd);

  // For constrained search: verify the Wagner tree satisfies the constraint.
  // If not, retry with different random orders. The constraint-aware edge
  // filter in wagner_tree makes violation increasingly unlikely, but can't
  // guarantee satisfaction for all split configurations.
  if (constrained) {
    for (int attempt = 1; attempt < 100; ++attempt) {
      if (!violates_constraint_posthoc(tree, *cd)) break;
      // Reshuffle and rebuild
      ts::rng_state_begin();
      for (int i = n_tips - 1; i > 0; --i) {
        int j = static_cast<int>(ts::thread_safe_unif() * (i + 1));
        if (j > i) j = i;
        std::swap(order[i], order[j]);
      }
      ts::rng_state_end();
      result = wagner_tree(tree, ds, order, cd);
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// Purely random tree topology (no scoring, no character-based placement).
// Builds a tree by inserting tips at uniformly random edges.
void random_topology_tree(TreeState& tree, const DataSet& ds) {
  int n_tip = ds.n_tips;
  check_wagner_precondition(n_tip);
  init_wagner_state(tree, ds);

  // Random tip insertion order (Fisher-Yates)
  std::vector<int> order(n_tip);
  std::iota(order.begin(), order.end(), 0);
  ts::rng_state_begin();
  for (int i = n_tip - 1; i > 0; --i) {
    int j = static_cast<int>(ts::thread_safe_unif() * (i + 1));
    if (j > i) j = i;
    std::swap(order[i], order[j]);
  }

  // Build initial 3-taxon tree
  build_three_taxon_tree(tree, order[0], order[1], order[2]);

  // Track all edge children (every node except root)
  std::vector<int> edge_children;
  edge_children.reserve(2 * n_tip - 2);
  edge_children.push_back(order[0]);
  edge_children.push_back(order[1]);
  edge_children.push_back(order[2]);
  edge_children.push_back(n_tip + 1);  // first internal node below root

  // Insert remaining tips at random edges
  for (int k = 3; k < n_tip; ++k) {
    int tip = order[k];
    int new_internal = n_tip + k - 1;

    int n_edges = static_cast<int>(edge_children.size());
    int edge_idx = static_cast<int>(ts::thread_safe_unif() * n_edges);
    if (edge_idx >= n_edges) edge_idx = n_edges - 1;
    int below = edge_children[edge_idx];
    int above = tree.parent[below];

    insert_tip_at_edge(tree, tip, new_internal, above, below);

    // New edges: new_internal (child of above) and tip (child of new_internal)
    edge_children.push_back(new_internal);
    edge_children.push_back(tip);
  }
  ts::rng_state_end();

  tree.build_postorder();
}

// =========================================================================
// Random constrained tree
// =========================================================================
//
// Algorithm:
// 1. Identify constraint splits ordered from largest to smallest (by
//    popcount of the "inside" set). Larger splits enclose smaller ones.
// 2. Assign each tip to its tightest (smallest) enclosing constraint
//    split, or "root level" if unconstrained.
// 3. Build the tree bottom-up: for each constraint split (smallest first),
//    randomly wire all its direct children (tips + smaller split roots)
//    into a binary subtree via random edge insertion.
// 4. Finally, wire all root-level items (unconstrained tips + top-level
//    split roots) into the tree.
//
// The result is a uniformly random binary tree among those that satisfy
// all constraint splits. (Uniform conditional on the split nesting
// structure, which determines the partition of items across polytomy
// resolution steps.)

namespace {

// Randomly resolve a set of items into a binary subtree.
// `items` are node indices (tips or internal subtree roots).
// Returns the root node of the resolved subtree.
// `next_internal` is incremented as internal nodes are consumed.
// Requires items.size() >= 1.
int resolve_randomly(TreeState& tree, std::vector<int>& items,
                     int& next_internal) {
  if (items.size() == 1) return items[0];

  // Shuffle items
  for (int i = static_cast<int>(items.size()) - 1; i > 0; --i) {
    int j = static_cast<int>(ts::thread_safe_unif() * (i + 1));
    if (j > i) j = i;
    std::swap(items[i], items[j]);
  }

  if (items.size() == 2) {
    int nd = next_internal++;
    int ni = nd - tree.n_tip;
    tree.left[ni] = items[0];
    tree.right[ni] = items[1];
    tree.parent[items[0]] = nd;
    tree.parent[items[1]] = nd;
    return nd;
  }

  // Build initial pair from first two items
  int nd = next_internal++;
  int ni = nd - tree.n_tip;
  tree.left[ni] = items[0];
  tree.right[ni] = items[1];
  tree.parent[items[0]] = nd;
  tree.parent[items[1]] = nd;

  // Track edges available for insertion (below-node of each edge)
  std::vector<int> edge_children;
  edge_children.push_back(items[0]);
  edge_children.push_back(items[1]);

  int subtree_root = nd;

  // Insert remaining items at random edges within this subtree
  for (size_t k = 2; k < items.size(); ++k) {
    int item = items[k];
    int new_nd = next_internal++;

    int n_edges = static_cast<int>(edge_children.size());
    int edge_idx = static_cast<int>(ts::thread_safe_unif() * n_edges);
    if (edge_idx >= n_edges) edge_idx = n_edges - 1;
    int below = edge_children[edge_idx];
    int above = tree.parent[below];

    // Insert: above -> new_nd -> {item, below}
    int new_ni = new_nd - tree.n_tip;
    tree.parent[new_nd] = above;
    tree.left[new_ni] = item;
    tree.right[new_ni] = below;
    tree.parent[item] = new_nd;
    tree.parent[below] = new_nd;

    // Update above's child pointer
    if (above >= tree.n_tip) {
      int ai = above - tree.n_tip;
      if (tree.left[ai] == below) {
        tree.left[ai] = new_nd;
      } else {
        tree.right[ai] = new_nd;
      }
    }

    // Update subtree root if we inserted above it
    if (below == subtree_root) {
      subtree_root = new_nd;
    }

    edge_children.push_back(new_nd);
    edge_children.push_back(item);
  }

  return subtree_root;
}

// Count set bits in a bitmask span
int popcount_span(const uint64_t* mask, int n_words) {
  int count = 0;
  for (int w = 0; w < n_words; ++w) {
    count += popcount64(mask[w]);
  }
  return count;
}

// Check if split a is a strict subset of split b (a ⊂ b)
bool is_strict_subset(const uint64_t* a, const uint64_t* b, int n_words) {
  bool proper = false;
  for (int w = 0; w < n_words; ++w) {
    if (a[w] & ~b[w]) return false;  // a has bits not in b
    if (b[w] & ~a[w]) proper = true; // b has bits not in a
  }
  return proper;
}

// Check if tip t is in split mask
bool tip_in_split(int t, const uint64_t* mask) {
  return (mask[t / 64] >> (t % 64)) & 1ULL;
}

} // anonymous namespace


void random_constrained_tree(TreeState& tree, const DataSet& ds,
                             ConstraintData& cd) {
  if (!cd.active || cd.n_splits == 0) {
    random_topology_tree(tree, ds);
    return;
  }

  int n_tip = ds.n_tips;
  check_wagner_precondition(n_tip);
  init_wagner_state(tree, ds);

  int n_words = cd.n_words;
  int n_splits = cd.n_splits;

  // --- Step 1: Order splits by popcount descending (largest first) ---
  // Then we'll process smallest first for bottom-up construction.
  std::vector<int> split_order(n_splits);
  std::iota(split_order.begin(), split_order.end(), 0);
  std::vector<int> split_size(n_splits);
  for (int s = 0; s < n_splits; ++s) {
    split_size[s] = popcount_span(
        &cd.split_tips[static_cast<size_t>(s) * n_words], n_words);
  }
  // Sort ascending by size (process smallest first)
  std::sort(split_order.begin(), split_order.end(),
    [&](int a, int b) { return split_size[a] < split_size[b]; });

  // --- Step 2: Find each split's parent (tightest enclosing split) ---
  // parent_split[s] = index into split_order of the parent, or -1 for root.
  std::vector<int> parent_split(n_splits, -1);
  for (int i = 0; i < n_splits; ++i) {
    int si = split_order[i];
    const uint64_t* si_mask =
        &cd.split_tips[static_cast<size_t>(si) * n_words];
    // Find smallest enclosing split (next larger split that contains si)
    for (int j = i + 1; j < n_splits; ++j) {
      int sj = split_order[j];
      const uint64_t* sj_mask =
          &cd.split_tips[static_cast<size_t>(sj) * n_words];
      if (is_strict_subset(si_mask, sj_mask, n_words)) {
        parent_split[i] = j;
        break;
      }
    }
  }

  // --- Step 3: Assign each tip to its tightest split ---
  // tip_owner[t] = index into split_order, or -1 for root level.
  std::vector<int> tip_owner(n_tip, -1);
  for (int t = 0; t < n_tip; ++t) {
    for (int i = 0; i < n_splits; ++i) {
      int si = split_order[i];
      const uint64_t* mask =
          &cd.split_tips[static_cast<size_t>(si) * n_words];
      if (tip_in_split(t, mask)) {
        tip_owner[t] = i;
        break;  // split_order is ascending, so first hit is tightest
      }
    }
  }

  // --- Step 4: Build bottom-up ---
  // For each split, collect its direct children (tips + child split roots)
  // and resolve them randomly.
  // Internal node allocation: root is n_tip (index 0 in left/right).
  // Sub-splits consume nodes n_tip+1, n_tip+2, ...
  // Root-level wiring uses the root node directly (no resolve_randomly).
  int root = n_tip;
  int next_internal = n_tip + 1;
  tree.parent[root] = root;

  // split_root[i] = node index of the subtree root for split_order[i]
  std::vector<int> split_root(n_splits, -1);

  ts::rng_state_begin();

  for (int i = 0; i < n_splits; ++i) {
    // Collect items that belong directly to this split
    std::vector<int> items;

    // Tips owned by this split
    for (int t = 0; t < n_tip; ++t) {
      if (tip_owner[t] == i) items.push_back(t);
    }

    // Child splits whose parent is this split
    for (int j = 0; j < i; ++j) {
      if (parent_split[j] == i && split_root[j] >= 0) {
        items.push_back(split_root[j]);
      }
    }

    if (items.empty()) {
      split_root[i] = -1;
      continue;
    }

    split_root[i] = resolve_randomly(tree, items, next_internal);
  }

  // --- Step 5: Wire root level ---
  // Collect unconstrained tips + top-level split roots, then build
  // directly onto the root node (avoiding extra node allocation).
  std::vector<int> root_items;

  for (int t = 0; t < n_tip; ++t) {
    if (tip_owner[t] == -1) root_items.push_back(t);
  }
  for (int i = 0; i < n_splits; ++i) {
    if (parent_split[i] == -1 && split_root[i] >= 0) {
      root_items.push_back(split_root[i]);
    }
  }

  // Shuffle root items
  for (int i = static_cast<int>(root_items.size()) - 1; i > 0; --i) {
    int j = static_cast<int>(ts::thread_safe_unif() * (i + 1));
    if (j > i) j = i;
    std::swap(root_items[i], root_items[j]);
  }

  if (root_items.size() >= 2) {
    // Wire first two items as root's children
    tree.left[0] = root_items[0];
    tree.right[0] = root_items[1];
    tree.parent[root_items[0]] = root;
    tree.parent[root_items[1]] = root;

    // Track edges for subsequent insertions
    std::vector<int> edge_children;
    edge_children.push_back(root_items[0]);
    edge_children.push_back(root_items[1]);

    // Insert remaining root items at random edges
    for (size_t k = 2; k < root_items.size(); ++k) {
      int item = root_items[k];
      int new_nd = next_internal++;
      int new_ni = new_nd - n_tip;

      int n_edges = static_cast<int>(edge_children.size());
      int edge_idx = static_cast<int>(ts::thread_safe_unif() * n_edges);
      if (edge_idx >= n_edges) edge_idx = n_edges - 1;
      int below = edge_children[edge_idx];
      int above = tree.parent[below];

      tree.parent[new_nd] = above;
      tree.left[new_ni] = item;
      tree.right[new_ni] = below;
      tree.parent[item] = new_nd;
      tree.parent[below] = new_nd;

      int ai = above - n_tip;
      if (tree.left[ai] == below) {
        tree.left[ai] = new_nd;
      } else {
        tree.right[ai] = new_nd;
      }

      edge_children.push_back(new_nd);
      edge_children.push_back(item);
    }
  } else if (root_items.size() == 1) {
    // Single top-level item must be an internal node; adopt its children
    int sub = root_items[0];
    if (sub >= n_tip) {
      int si = sub - n_tip;
      tree.left[0] = tree.left[si];
      tree.right[0] = tree.right[si];
      tree.parent[tree.left[0]] = root;
      tree.parent[tree.right[0]] = root;
    }
  }

  ts::rng_state_end();

  tree.build_postorder();
  update_constraint(tree, cd);
}

} // namespace ts
