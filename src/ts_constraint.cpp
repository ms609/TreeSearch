#include "ts_constraint.h"
#include "ts_fitch.h"
#include <algorithm>
#include <cstring>
#include <vector>

namespace ts {

// =========================================================================
// Build constraint from R-side split matrix
// =========================================================================

ConstraintData build_constraint(
    const int* split_matrix, int n_splits, int n_tips)
{
  ConstraintData cd;
  if (n_splits == 0) return cd;

  cd.active = true;
  cd.n_splits = n_splits;
  cd.n_words = (n_tips + 63) / 64;

  cd.split_tips.resize(
      static_cast<size_t>(n_splits) * cd.n_words, 0ULL);
  cd.constraint_node.assign(n_splits, -1);

  // Pack split_matrix rows into bitmasks
  for (int s = 0; s < n_splits; ++s) {
    uint64_t* mask = &cd.split_tips[static_cast<size_t>(s) * cd.n_words];
    for (int t = 0; t < n_tips; ++t) {
      if (split_matrix[s * n_tips + t]) {
        int w = t / 64;
        int b = t % 64;
        mask[w] |= (1ULL << b);
      }
    }
    // Canonicalize: tip 0 must be on the "outside" (bit 0 = 0).
    // If bit 0 is set, flip the entire mask.
    if (mask[0] & 1ULL) {
      for (int w = 0; w < cd.n_words; ++w) {
        mask[w] = ~mask[w];
      }
      // Clear bits beyond n_tips
      int remainder = n_tips % 64;
      if (remainder > 0) {
        mask[cd.n_words - 1] &= (1ULL << remainder) - 1;
      }
    }
  }

  // Allocate DFS timestamp arrays
  int n_node = 2 * n_tips - 1;
  cd.dfs_entry.assign(n_node, 0);
  cd.dfs_exit.assign(n_node, 0);

  // Allocate per-clip workspace
  cd.clip_zones.resize(n_splits, ClipZone::UNCONSTRAINED);
  cd.clip_tip_mask.resize(cd.n_words, 0ULL);

  return cd;
}

// =========================================================================
// Post-hoc fallback: build a DataSet from constraint phyDat
// =========================================================================

void build_constraint_posthoc(
    ConstraintData& cd,
    const double* contrast_r, int n_tokens, int n_states,
    const int* tip_data_r, int n_tips, int n_patterns,
    const int* weight_r,
    const char** levels_r,
    int expected_score)
{
  cd.posthoc_data = build_dataset(
      contrast_r, n_tokens, n_states,
      tip_data_r, n_tips, n_patterns,
      weight_r, levels_r);
  cd.expected_score = expected_score;
  cd.has_posthoc = true;
}

// =========================================================================
// Map constraint nodes: find which internal node holds each split
// =========================================================================

void map_constraint_nodes(const TreeState& tree, ConstraintData& cd)
{
  if (!cd.active) return;

  int n_tip = tree.n_tip;

  // Build per-node subtree tip bitmasks via postorder traversal.
  // For tips: bit[t] = 1.
  // For internal nodes: OR of children's masks.
  std::vector<uint64_t> node_tips(
      static_cast<size_t>(tree.n_node) * cd.n_words, 0ULL);

  // Initialize tips
  for (int t = 0; t < n_tip; ++t) {
    int w = t / 64;
    int b = t % 64;
    node_tips[static_cast<size_t>(t) * cd.n_words + w] = (1ULL << b);
  }

  // Postorder: compute internal node tip masks bottom-up
  for (int node : tree.postorder) {
    int ni = node - n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    uint64_t* nd = &node_tips[static_cast<size_t>(node) * cd.n_words];
    const uint64_t* lt = &node_tips[static_cast<size_t>(lc) * cd.n_words];
    const uint64_t* rt = &node_tips[static_cast<size_t>(rc) * cd.n_words];
    for (int w = 0; w < cd.n_words; ++w) {
      nd[w] = lt[w] | rt[w];
    }
  }

  // For each constraint split, find the internal node whose subtree
  // tip mask matches (after canonicalization with tip 0 outside).
  // A node matches split S if its subtree tips == S.
  for (int s = 0; s < cd.n_splits; ++s) {
    const uint64_t* split = &cd.split_tips[static_cast<size_t>(s) * cd.n_words];
    cd.constraint_node[s] = -1;

    for (int node : tree.postorder) {
      const uint64_t* nd = &node_tips[static_cast<size_t>(node) * cd.n_words];
      bool match = true;
      for (int w = 0; w < cd.n_words; ++w) {
        if (nd[w] != split[w]) { match = false; break; }
      }
      if (match) {
        cd.constraint_node[s] = node;
        break;
      }
    }
    // If constraint_node[s] == -1, the tree doesn't display this split.
    // This is fine during search — moves that break a constraint will be
    // rejected, and the node will be remapped after the next valid move.
  }
}

// =========================================================================
// DFS timestamps for O(1) descendant queries
// =========================================================================

void compute_dfs_timestamps(const TreeState& tree, ConstraintData& cd)
{
  if (!cd.active) return;

  int n_tip = tree.n_tip;
  int root = n_tip;
  int counter = 0;

  // Iterative DFS using an explicit stack.
  // Stack entries: (node, is_exit).
  struct DFSEntry { int node; bool is_exit; };
  std::vector<DFSEntry> stack;
  stack.push_back({root, false});

  while (!stack.empty()) {
    auto [node, is_exit] = stack.back();
    stack.pop_back();

    if (is_exit) {
      cd.dfs_exit[node] = counter++;
      continue;
    }

    cd.dfs_entry[node] = counter++;

    if (node < n_tip) {
      // Tip: exit immediately
      cd.dfs_exit[node] = counter++;
      continue;
    }

    int ni = node - n_tip;
    // Push exit marker first (will be processed last)
    stack.push_back({node, true});
    // Push children (right first so left is processed first)
    stack.push_back({tree.right[ni], false});
    stack.push_back({tree.left[ni], false});
  }
}

// =========================================================================
// Combined update
// =========================================================================

void update_constraint(const TreeState& tree, ConstraintData& cd)
{
  if (!cd.active) return;
  map_constraint_nodes(tree, cd);
  compute_dfs_timestamps(tree, cd);
}

// =========================================================================
// Per-clip: compute clipped subtree tip mask
// =========================================================================

void compute_clip_tip_mask(const TreeState& tree, int clip_node,
                           std::vector<uint64_t>& mask)
{
  std::fill(mask.begin(), mask.end(), 0ULL);

  // DFS through clipped subtree
  std::vector<int> stack;
  stack.push_back(clip_node);
  while (!stack.empty()) {
    int node = stack.back();
    stack.pop_back();
    if (node < tree.n_tip) {
      int w = node / 64;
      int b = node % 64;
      mask[w] |= (1ULL << b);
    } else {
      int ni = node - tree.n_tip;
      stack.push_back(tree.left[ni]);
      stack.push_back(tree.right[ni]);
    }
  }
}

// =========================================================================
// Per-clip: classify each constraint split
// =========================================================================

void classify_clip_constraints(const TreeState& tree, int clip_node,
                               ConstraintData& cd)
{
  if (!cd.active) return;

  compute_clip_tip_mask(tree, clip_node, cd.clip_tip_mask);

  for (int s = 0; s < cd.n_splits; ++s) {
    const uint64_t* split =
        &cd.split_tips[static_cast<size_t>(s) * cd.n_words];

    bool any_inside = false;
    bool any_outside = false;
    for (int w = 0; w < cd.n_words; ++w) {
      if (cd.clip_tip_mask[w] & split[w]) any_inside = true;
      if (cd.clip_tip_mask[w] & ~split[w]) any_outside = true;
      if (any_inside && any_outside) break;
    }

    if (any_inside && any_outside) {
      // Clip subtree straddles the split — the split boundary is
      // internal to the clipped subtree, so regraft position doesn't
      // affect it.
      cd.clip_zones[s] = ClipZone::UNCONSTRAINED;
    } else if (any_inside) {
      cd.clip_zones[s] = ClipZone::MUST_INSIDE;
    } else if (any_outside) {
      cd.clip_zones[s] = ClipZone::MUST_OUTSIDE;
    } else {
      cd.clip_zones[s] = ClipZone::UNCONSTRAINED;
    }
  }
}

// =========================================================================
// Per-candidate: check regraft legality
// =========================================================================

// Helper: is node `u` an ancestor of node `v` (or equal to v)?
static inline bool is_ancestor_or_equal(
    int u, int v,
    const std::vector<int>& entry, const std::vector<int>& exit)
{
  return entry[u] <= entry[v] && exit[u] >= exit[v];
}

bool regraft_violates_constraint(int below,
                                 const ConstraintData& cd)
{
  if (!cd.active) return false;

  for (int s = 0; s < cd.n_splits; ++s) {
    if (cd.clip_zones[s] == ClipZone::UNCONSTRAINED) continue;

    int cn = cd.constraint_node[s];
    if (cn < 0) {
      // Constraint not currently displayed — tree already violates.
      // Reject all moves to avoid entrenching a bad state.
      // (This shouldn't happen if the starting tree is valid and
      // we only accept valid moves.)
      return true;
    }

    // Is `below` a descendant of cn (= inside the constraint clade)?
    bool inside = is_ancestor_or_equal(cn, below,
                                        cd.dfs_entry, cd.dfs_exit);

    if (cd.clip_zones[s] == ClipZone::MUST_INSIDE && !inside) {
      return true;
    }
    if (cd.clip_zones[s] == ClipZone::MUST_OUTSIDE && inside) {
      return true;
    }
  }

  return false;
}

// =========================================================================
// Post-hoc: full Fitch check
// =========================================================================

bool violates_constraint_posthoc(const TreeState& tree,
                                 const ConstraintData& cd)
{
  if (!cd.active || !cd.has_posthoc) return false;

  // Build a temporary TreeState for the constraint dataset using
  // the current tree's topology.
  TreeState ctree;
  ctree.n_tip = tree.n_tip;
  ctree.n_internal = tree.n_internal;
  ctree.n_node = tree.n_node;
  ctree.total_words = cd.posthoc_data.total_words;
  ctree.n_blocks = cd.posthoc_data.n_blocks;

  // Copy topology
  ctree.parent = tree.parent;
  ctree.left = tree.left;
  ctree.right = tree.right;

  // Allocate state arrays
  size_t state_size =
      static_cast<size_t>(ctree.n_node) * ctree.total_words;
  ctree.prelim.assign(state_size, 0ULL);
  ctree.final_.assign(state_size, 0ULL);
  ctree.down2.assign(state_size, 0ULL);
  ctree.subtree_actives.assign(state_size, 0ULL);
  ctree.local_cost.assign(
      static_cast<size_t>(ctree.n_node) * ctree.n_blocks, 0ULL);

  ctree.load_tip_states(cd.posthoc_data);
  ctree.build_postorder();

  int score = fitch_score(ctree, cd.posthoc_data);
  return score != cd.expected_score;
}

} // namespace ts
