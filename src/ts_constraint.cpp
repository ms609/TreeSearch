#include "ts_constraint.h"
#include "ts_fitch.h"
#include "ts_rng.h"
#include <algorithm>
#include <cstring>
#include <numeric>
#include <random>
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

  // Pack split_matrix rows into bitmasks.
  // split_matrix is column-major (from R): element [s, t] is at
  // index s + n_splits * t.
  for (int s = 0; s < n_splits; ++s) {
    uint64_t* mask = &cd.split_tips[static_cast<size_t>(s) * cd.n_words];
    for (int t = 0; t < n_tips; ++t) {
      if (split_matrix[s + n_splits * t]) {
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
// Build constraint from pre-canonicalized split bitsets
// =========================================================================

ConstraintData build_constraint_from_bitsets(
    const uint64_t* split_bits, int n_splits,
    int words_per_split, int n_tips) {
  ConstraintData cd;
  if (n_splits == 0) return cd;

  cd.active = true;
  cd.n_splits = n_splits;
  cd.n_words = words_per_split;

  // Copy split data
  size_t total = static_cast<size_t>(n_splits) * words_per_split;
  cd.split_tips.assign(split_bits, split_bits + total);
  cd.constraint_node.assign(n_splits, -1);

  int n_node = 2 * n_tips - 1;
  cd.dfs_entry.assign(n_node, 0);
  cd.dfs_exit.assign(n_node, 0);

  cd.clip_zones.resize(n_splits, ClipZone::UNCONSTRAINED);
  cd.clip_tip_mask.resize(words_per_split, 0ULL);

  // No posthoc fallback — sector/fuse won't enforce these constraints
  cd.has_posthoc = false;

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
// Compute per-node subtree tip bitmasks
// =========================================================================

std::vector<uint64_t> compute_node_tips(const TreeState& tree, int n_words)
{
  std::vector<uint64_t> node_tips(
      static_cast<size_t>(tree.n_node) * n_words, 0ULL);

  // Initialize tips
  for (int t = 0; t < tree.n_tip; ++t) {
    int w = t / 64;
    int b = t % 64;
    node_tips[static_cast<size_t>(t) * n_words + w] = (1ULL << b);
  }

  // Postorder: compute internal node tip masks bottom-up
  for (int node : tree.postorder) {
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    uint64_t* nd = &node_tips[static_cast<size_t>(node) * n_words];
    const uint64_t* lt = &node_tips[static_cast<size_t>(lc) * n_words];
    const uint64_t* rt = &node_tips[static_cast<size_t>(rc) * n_words];
    for (int w = 0; w < n_words; ++w) {
      nd[w] = lt[w] | rt[w];
    }
  }

  return node_tips;
}

// =========================================================================
// Map constraint nodes: find which internal node holds each split
// =========================================================================

void map_constraint_nodes(const TreeState& tree, ConstraintData& cd)
{
  if (!cd.active) return;

  auto node_tips = compute_node_tips(tree, cd.n_words);

  // For each constraint split, find the internal node whose subtree
  // tip mask matches (after canonicalization with tip 0 outside).
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
      // Clip subtree straddles the split.  Check whether the rest of
      // the tree also has tips on both sides.  If so, no single edge
      // can separate IN from OUT after any regraft → FORBIDDEN.
      // If the clip contains ALL of one side, the split boundary could
      // be internal to the clip subtree → UNCONSTRAINED.
      bool rest_has_in = false;
      bool rest_has_out = false;
      for (int w = 0; w < cd.n_words; ++w) {
        uint64_t rest = ~cd.clip_tip_mask[w];
        // Mask out bits beyond n_tips in the last word
        if (w == cd.n_words - 1) {
          int remainder = tree.n_tip % 64;
          if (remainder > 0)
            rest &= (1ULL << remainder) - 1;
        }
        if (rest & split[w]) rest_has_in = true;
        if (rest & ~split[w]) {
          uint64_t out_bits = ~split[w];
          if (w == cd.n_words - 1) {
            int remainder = tree.n_tip % 64;
            if (remainder > 0)
              out_bits &= (1ULL << remainder) - 1;
          }
          if (rest & out_bits) rest_has_out = true;
        }
      }
      if (rest_has_in && rest_has_out) {
        cd.clip_zones[s] = ClipZone::FORBIDDEN;
      } else {
        cd.clip_zones[s] = ClipZone::UNCONSTRAINED;
      }
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

    // Clip straddles the split AND rest also straddles: no regraft
    // can preserve this split — reject unconditionally.
    if (cd.clip_zones[s] == ClipZone::FORBIDDEN) return true;

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
    // Exclude the boundary edge (above_cn, cn): regrafting an outside-only
    // clade just above the constraint clade makes it a sibling of that clade,
    // preserving monophyly.  Only reject if the clade would land *strictly
    // inside* the constraint clade.
    if (cd.clip_zones[s] == ClipZone::MUST_OUTSIDE && inside && below != cn) {
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

// =========================================================================
// Post-hoc repair: impose constraints via minimal SPR moves
// =========================================================================

namespace {

// Topology-only SPR: move `clip` to the edge (above, below).
// Unlike spr_clip/spr_regraft, this handles root-child clips correctly
// and doesn't save/restore state (caller must rebuild postorder and rescore).
void topology_spr(TreeState& tree, int clip, int above, int below) {
  const int root = tree.n_tip;
  int nx = tree.parent[clip];
  int ns = (tree.left[nx - root] == clip)
               ? tree.right[nx - root]
               : tree.left[nx - root];

  if (nx != root) {
    // --- Normal case: detach nx, connect ns to grandparent ---
    int nz = tree.parent[nx];
    tree.parent[ns] = nz;
    if (nz >= tree.n_tip) {
      int nzi = nz - tree.n_tip;
      if (tree.left[nzi] == nx)
        tree.left[nzi] = ns;
      else
        tree.right[nzi] = ns;
    }

    // Insert nx between above and below
    if (above >= tree.n_tip) {
      int ai = above - tree.n_tip;
      if (tree.left[ai] == below)
        tree.left[ai] = nx;
      else
        tree.right[ai] = nx;
    }
    tree.parent[nx] = above;
    int nxi = nx - tree.n_tip;
    tree.left[nxi] = clip;
    tree.right[nxi] = below;
    tree.parent[clip] = nx;
    tree.parent[below] = nx;
  } else {
    // --- Root-child case: clip is a direct child of root ---
    // Can't float root (identity is fixed at n_tip).
    // Absorb ns into root and repurpose ns as the insertion node.
    if (ns < tree.n_tip) return;  // ns is a tip — degenerate, bail out

    int nsi = ns - tree.n_tip;
    int ns_left = tree.left[nsi];
    int ns_right = tree.right[nsi];

    // Root absorbs ns's children
    tree.left[0] = ns_left;
    tree.right[0] = ns_right;
    tree.parent[ns_left] = root;
    tree.parent[ns_right] = root;

    // Insert ns between above and below, with clip as its other child
    if (above >= tree.n_tip) {
      int ai = above - tree.n_tip;
      if (tree.left[ai] == below)
        tree.left[ai] = ns;
      else
        tree.right[ai] = ns;
    }
    tree.parent[ns] = above;
    tree.left[nsi] = clip;
    tree.right[nsi] = below;
    tree.parent[clip] = ns;
    tree.parent[below] = ns;
  }
}

// Collect (above, below) edge pairs within the subtree rooted at node.
// Iterative DFS; does NOT include the edge above `node` itself.
void collect_edges_in_subtree(const TreeState& tree, int sub_root,
                              std::vector<std::pair<int,int>>& edges) {
  std::vector<int> stack;
  stack.push_back(sub_root);
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

// Collect (above, below) edge pairs NOT in the subtree rooted at
// `exclude_root`. DFS from tree root, skipping the excluded subtree.
// Includes the edge leading to exclude_root (a valid outside target:
// regrafting there makes the moved subtree a sibling of exclude_root).
void collect_edges_outside_subtree(
    const TreeState& tree, int exclude_root,
    std::vector<std::pair<int,int>>& edges) {
  std::vector<int> stack;
  stack.push_back(tree.n_tip); // tree root
  while (!stack.empty()) {
    int node = stack.back();
    stack.pop_back();
    if (node < tree.n_tip) continue;
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    edges.push_back({node, lc});
    edges.push_back({node, rc});
    // Descend into children, but skip the excluded subtree
    if (lc != exclude_root) stack.push_back(lc);
    if (rc != exclude_root) stack.push_back(rc);
  }
}

// Iterative DFS to find maximal subtrees whose tips are entirely
// within `mask`. Searches from `search_root`, optionally skipping
// `exclude` (-1 to skip nothing).
void find_maximal_subtrees(const TreeState& tree, int search_root,
                           int exclude,
                           const std::vector<uint64_t>& node_tips,
                           const std::vector<uint64_t>& mask,
                           int n_words,
                           std::vector<int>& out) {
  std::vector<int> stack;
  stack.push_back(search_root);
  while (!stack.empty()) {
    int node = stack.back();
    stack.pop_back();
    if (node == exclude) continue;
    const uint64_t* nt =
        &node_tips[static_cast<size_t>(node) * n_words];
    bool all_in = true;
    bool any_in = false;
    for (int w = 0; w < n_words; ++w) {
      if (nt[w] & ~mask[w]) all_in = false;
      if (nt[w] & mask[w]) any_in = true;
    }
    if (all_in && any_in) {
      out.push_back(node); // Maximal: don't recurse
      continue;
    }
    if (any_in && node >= tree.n_tip) {
      int ni = node - tree.n_tip;
      stack.push_back(tree.left[ni]);
      stack.push_back(tree.right[ni]);
    }
  }
}

// Full structural validity of a rooted binary TreeState (tips 0..n_tip-1,
// internal n_tip..2*n_tip-2, root == n_tip).  Returns true iff the tree is a
// well-formed arborescence AND parent[] is consistent with left/right:
//   (1) left/right slots in range, left != right;
//   (2) in-degree via left/right is exactly 1 for every non-root node and 0
//       for the root (no double-reference, no orphan);
//   (3) a DFS from the root over left/right visits every node exactly once —
//       acyclic AND fully connected (catches a disjoint in-degree-1 cycle,
//       which the in-degree pass alone cannot);
//   (4) parent[root] == root and parent[] is the exact inverse of left/right,
//       so a parent-ASCENDING consumer (reroot_at_tip's
//       `while (cur != root) cur = parent[cur]`, ts_tbr.cpp) can neither loop
//       nor wander, and the next topology_spr's `parent[clip]` read is sound.
// O(n).  This is the complete check the T-327 backstop only approximated:
// build_postorder counts only reachable INTERNAL nodes, so a net-zero
// corruption (one node double-referenced +1, one orphaned -1 — topology_spr's
// root-child case) lands on exactly n_internal and slips a postorder.size()
// test.  Exhaustively characterised in dev/red-team/heavy-tests/impose_validity/.
bool structurally_valid(const TreeState& tree) {
  const int nt = tree.n_tip, ni = tree.n_internal, nn = tree.n_node,
            root = nt;
  if (static_cast<int>(tree.left.size())   != ni) return false;
  if (static_cast<int>(tree.right.size())  != ni) return false;
  if (static_cast<int>(tree.parent.size()) != nn) return false;

  // (1) child slots in range and distinct.
  for (int i = 0; i < ni; ++i) {
    const int l = tree.left[i], r = tree.right[i];
    if (l < 0 || l >= nn || r < 0 || r >= nn) return false;
    if (l == r) return false;
  }

  // (2) in-degree: exactly one parent per non-root node, none for the root.
  std::vector<int> indeg(nn, 0);
  for (int i = 0; i < ni; ++i) { ++indeg[tree.left[i]]; ++indeg[tree.right[i]]; }
  if (indeg[root] != 0) return false;
  for (int v = 0; v < nn; ++v) {
    if (v != root && indeg[v] != 1) return false;
  }

  // (3) root DFS reaches every node exactly once (acyclic + connected).
  std::vector<char> seen(nn, 0);
  std::vector<int> stack;
  stack.reserve(nn);
  stack.push_back(root);
  int visited = 0;
  while (!stack.empty()) {
    const int node = stack.back();
    stack.pop_back();
    if (seen[node]) return false;           // revisit => double-reference/cycle
    seen[node] = 1;
    ++visited;
    if (node >= nt) {
      const int i = node - nt;
      stack.push_back(tree.left[i]);
      stack.push_back(tree.right[i]);
    }
  }
  if (visited != nn) return false;          // orphaned node(s)

  // (4) parent[] is the exact inverse of left/right.
  if (tree.parent[root] != root) return false;
  for (int i = 0; i < ni; ++i) {
    const int node = nt + i;
    if (tree.parent[tree.left[i]]  != node) return false;
    if (tree.parent[tree.right[i]] != node) return false;
  }
  return true;
}

} // anonymous namespace


// Single pass: fix all currently-violated splits (smallest first).
// Returns number of SPR moves performed.
static int impose_one_pass(TreeState& tree, ConstraintData& cd,
                           std::mt19937& rng) {
  const int n_words = cd.n_words;
  const int root = tree.n_tip;

  tree.build_postorder();
  auto node_tips = compute_node_tips(tree, n_words);

  // --- Identify violated splits ---
  std::vector<int> violated;
  for (int s = 0; s < cd.n_splits; ++s) {
    const uint64_t* split =
        &cd.split_tips[static_cast<size_t>(s) * n_words];
    bool found = false;
    for (int node : tree.postorder) {
      const uint64_t* nd =
          &node_tips[static_cast<size_t>(node) * n_words];
      bool match = true;
      for (int w = 0; w < n_words; ++w) {
        if (nd[w] != split[w]) { match = false; break; }
      }
      if (match) { found = true; break; }
    }
    if (!found) violated.push_back(s);
  }

  if (violated.empty()) return 0;

  // Sort by popcount ascending (smallest first)
  std::sort(violated.begin(), violated.end(),
    [&](int a, int b) {
      int pa = 0, pb = 0;
      const uint64_t* sa =
          &cd.split_tips[static_cast<size_t>(a) * n_words];
      const uint64_t* sb =
          &cd.split_tips[static_cast<size_t>(b) * n_words];
      for (int w = 0; w < n_words; ++w) {
        pa += popcount64(sa[w]);
        pb += popcount64(sb[w]);
      }
      return pa < pb;
    });

  int total_moves = 0;

  for (size_t vi = 0; vi < violated.size(); ++vi) {
    int s = violated[vi];
    const uint64_t* split =
        &cd.split_tips[static_cast<size_t>(s) * n_words];

    // Rebuild bitmasks after previous split's moves
    if (vi > 0) {
      tree.build_postorder();
      node_tips = compute_node_tips(tree, n_words);
    }

    // --- Find best candidate node (min symmetric difference) ---
    int best_node = -1;
    int best_cost = tree.n_tip + 1;
    for (int node : tree.postorder) {
      const uint64_t* nd =
          &node_tips[static_cast<size_t>(node) * n_words];
      int cost = 0;
      for (int w = 0; w < n_words; ++w) {
        cost += popcount64(nd[w] ^ split[w]);
      }
      if (cost < best_cost) {
        best_cost = cost;
        best_node = node;
      }
    }

    if (best_cost == 0) continue; // Already satisfied

    // --- Compute misplaced tip masks ---
    std::vector<uint64_t> move_out_mask(n_words);
    std::vector<uint64_t> move_in_mask(n_words);
    const uint64_t* best_nt =
        &node_tips[static_cast<size_t>(best_node) * n_words];
    for (int w = 0; w < n_words; ++w) {
      move_out_mask[w] = best_nt[w] & ~split[w];
      move_in_mask[w]  = split[w] & ~best_nt[w];
    }

    // --- Find maximal subtrees to move ---
    std::vector<int> move_out_roots;
    find_maximal_subtrees(tree, best_node, -1, node_tips,
                          move_out_mask, n_words, move_out_roots);

    std::vector<int> move_in_roots;
    find_maximal_subtrees(tree, root, best_node, node_tips,
                          move_in_mask, n_words, move_in_roots);

    // Safety cap: abandon this pass if the repair is unexpectedly large.
    int n_moves = static_cast<int>(
        move_out_roots.size() + move_in_roots.size());
    if (total_moves + n_moves > tree.n_tip / 4 + 2) {
      return -1;  // Distinguish "bailed out" from "no violations" (0)
    }

    // --- Execute topology moves ---
    // Uses topology_spr() which handles root-child moves correctly
    // (unlike spr_clip which can't detach root children).
    //
    // T-327: topology_spr can relocate the node id `best_node` — when a
    // move root's parent IS best_node it detaches that internal node and
    // reuses it at the far graft point, so the captured `best_node` then
    // points at the wrong region. The stale node makes collect_edges_*()
    // enumerate the wrong edges, and a graft target adjacent to (or inside)
    // the moved subtree splices a node under its own descendant → a
    // cyclic / double-parented tree. Any later DFS on that tree
    // (collect_edges_*, compute_node_tips, build_postorder) then walks the
    // cycle and allocates without bound → std::bad_alloc.
    //
    // Two defences make each move independently safe:
    //  (a) re-anchor best_node against the split before every move (cheap;
    //      not the hot path) so edge enumeration is relative to the current
    //      topology; and
    //  (b) snapshot the topology, apply the move, and validate that the
    //      result is still an acyclic tree (postorder visits every internal
    //      node exactly once). If not, revert the move and skip it. This is
    //      complete without having to enumerate every topology_spr corner
    //      case, and guarantees no corrupt tree ever reaches a DFS helper.
    auto reanchor_best_node = [&]() -> int {
      tree.build_postorder();
      auto nt = compute_node_tips(tree, n_words);
      int bn = -1;
      int bc = tree.n_tip + 1;
      for (int node : tree.postorder) {
        const uint64_t* nd = &nt[static_cast<size_t>(node) * n_words];
        int cost = 0;
        for (int w = 0; w < n_words; ++w) {
          cost += popcount64(nd[w] ^ split[w]);
        }
        if (cost < bc) { bc = cost; bn = node; }
      }
      return bn;
    };

    // Apply one SPR, reverting if it corrupts the tree. Returns true if the
    // move was applied (tree still valid), false if reverted.
    auto try_move = [&](int M, bool outside) -> bool {
      int bn = reanchor_best_node();  // rebuilds postorder on a valid tree
      std::vector<std::pair<int,int>> targets;
      if (outside) collect_edges_outside_subtree(tree, bn, targets);
      else         collect_edges_in_subtree(tree, bn, targets);
      if (targets.empty()) return false;
      auto [above, below] =
          targets[std::uniform_int_distribution<int>(
              0, static_cast<int>(targets.size()) - 1)(rng)];
      // Snapshot topology (topology_spr mutates only parent/left/right).
      std::vector<int> save_parent = tree.parent;
      std::vector<int> save_left = tree.left;
      std::vector<int> save_right = tree.right;
      topology_spr(tree, M, above, below);
      // T-333: validate full structure, not a node count.  topology_spr's
      // root-child case can emit a net-zero corruption (one internal node
      // double-referenced, one orphaned) that lands on exactly n_internal and
      // would slip a postorder.size() != n_internal test — yet still feed a
      // cyclic/dangling tree to the next reanchor DFS and to reroot_at_tip's
      // parent walk.  structurally_valid() is O(n), the same complexity class
      // as the build_postorder it replaces here, and executes ONLY on this
      // constraint-repair path (never in unconstrained search); validating
      // before the rebuild keeps build_postorder off corrupt trees entirely.
      // (TreeState::build_postorder keeps its own > n_internal cap as the
      // bad_alloc backstop for its ~96 other, hot-path call sites.)
      if (!structurally_valid(tree)) {
        tree.parent = std::move(save_parent);
        tree.left = std::move(save_left);
        tree.right = std::move(save_right);
        tree.build_postorder();  // restore postorder for the (valid) saved tree
        return false;
      }
      tree.build_postorder();    // valid: refresh postorder for downstream use
      return true;
    };

    for (int M : move_out_roots) {
      if (try_move(M, /*outside=*/true)) ++total_moves;
    }
    for (int M : move_in_roots) {
      if (try_move(M, /*outside=*/false)) ++total_moves;
    }
  }

  return total_moves;
}


int impose_constraint(TreeState& tree, ConstraintData& cd)
{
  if (!cd.active) return 0;

  std::mt19937 rng = ts::make_rng();
  int total_moves = 0;

  // Iterate: fixing one split can break another (e.g. moving a tip
  // outside a small split may land it outside a larger enclosing split).
  // Each pass fixes at least the smallest violated split, so convergence
  // is bounded by n_splits.  Cap at n_splits + 1 for safety.
  for (int pass = 0; pass <= cd.n_splits; ++pass) {
    int moves = impose_one_pass(tree, cd, rng);
    if (moves < 0) break;  // Bailed out — too many moves needed
    if (moves == 0) break;  // No violations found — done
    total_moves += moves;
  }

  tree.build_postorder();
  update_constraint(tree, cd);
  return total_moves;
}

} // namespace ts
