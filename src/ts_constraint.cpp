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
  cd.split_zeros.resize(
      static_cast<size_t>(n_splits) * cd.n_words, 0ULL);
  cd.constraint_node.assign(n_splits, -1);
  cd.constraint_node_hi.assign(n_splits, -1);
  cd.constraint_complement.assign(n_splits, 0);

  // Pack split_matrix rows into a pair of bitmasks.
  // split_matrix is column-major (from R): element [s, t] is at
  // index s + n_splits * t.  1 -> "together" group, 0 -> "apart" group,
  // anything else (NA_INTEGER) -> free, in neither mask (#54).
  for (int s = 0; s < n_splits; ++s) {
    uint64_t* ones = &cd.split_tips[static_cast<size_t>(s) * cd.n_words];
    uint64_t* zeros = &cd.split_zeros[static_cast<size_t>(s) * cd.n_words];
    for (int t = 0; t < n_tips; ++t) {
      const int v = split_matrix[s + n_splits * t];
      if (v != 1 && v != 0) continue;               // free tip
      int w = t / 64;
      int b = t % 64;
      (v == 1 ? ones : zeros)[w] |= (1ULL << b);
    }
    // Canonicalize: tip 0 must be outside split_tips (bit 0 = 0).  The two
    // groups of a bipartition are interchangeable, so SWAP them rather than
    // complementing either — complementing would swallow the free tips into
    // the "apart" group and reinstate the exact-clade reading.
    if (ones[0] & 1ULL) {
      for (int w = 0; w < cd.n_words; ++w) {
        std::swap(ones[w], zeros[w]);
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
  // These splits come from pool bipartitions, which partition every tip: there
  // are no free tips, so the "apart" group is exactly the complement and the
  // free-taxa machinery collapses back to the exact-clade test (#54).
  cd.split_zeros.assign(total, 0ULL);
  {
    const int rem = n_tips % 64;
    const uint64_t top = rem ? ((1ULL << rem) - 1ULL) : ~0ULL;
    for (int s = 0; s < n_splits; ++s) {
      const size_t off = static_cast<size_t>(s) * words_per_split;
      for (int w = 0; w < words_per_split; ++w) {
        cd.split_zeros[off + w] = ~cd.split_tips[off + w];
        if (w == words_per_split - 1) cd.split_zeros[off + w] &= top;
      }
    }
  }
  cd.constraint_node.assign(n_splits, -1);
  cd.constraint_node_hi.assign(n_splits, -1);
  cd.constraint_complement.assign(n_splits, 0);

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
// Negative (converse) constraints
// =========================================================================

void add_negative_constraint(
    ConstraintData& cd, const int* neg_matrix, int n_neg, int n_tips)
{
  if (n_neg <= 0) return;

  int n_words = (n_tips + 63) / 64;
  if (!cd.active) {
    // Initialize an otherwise-empty ConstraintData so the search's
    // `constrained` path runs (it will be a no-op over zero positive splits).
    cd.active = true;
    cd.n_words = n_words;
    int n_node = 2 * n_tips - 1;
    cd.dfs_entry.assign(n_node, 0);
    cd.dfs_exit.assign(n_node, 0);
    cd.clip_tip_mask.resize(n_words, 0ULL);
  }

  cd.neg_active = true;
  cd.n_neg_splits = n_neg;
  cd.neg_split_tips.assign(
      static_cast<size_t>(n_neg) * cd.n_words, 0ULL);

  // Pack rows into bitmasks (column-major from R: [s, t] at s + n_neg * t),
  // canonicalizing so tip 0 is always outside (bit 0 = 0), exactly as
  // build_constraint() does for positive splits.
  for (int s = 0; s < n_neg; ++s) {
    uint64_t* mask = &cd.neg_split_tips[static_cast<size_t>(s) * cd.n_words];
    for (int t = 0; t < n_tips; ++t) {
      if (neg_matrix[s + n_neg * t]) {
        mask[t / 64] |= (1ULL << (t % 64));
      }
    }
    if (mask[0] & 1ULL) {
      for (int w = 0; w < cd.n_words; ++w) {
        mask[w] = ~mask[w];
      }
      int remainder = n_tips % 64;
      if (remainder > 0) {
        mask[cd.n_words - 1] &= (1ULL << remainder) - 1;
      }
    }
  }
}

bool displays_forbidden_clade(const TreeState& tree, const ConstraintData& cd)
{
  if (!cd.neg_active || cd.n_neg_splits == 0) return false;

  auto node_tips = compute_node_tips(tree, cd.n_words);

  // The internal tree is not guaranteed to be rooted on tip 0, so a bipartition
  // can be displayed as EITHER side.  Canonicalize each node's subtree mask to
  // the tip-0-outside form (flipping if tip 0 is inside) before comparing to the
  // canonical forbidden split; this matches a clade regardless of which side is
  // the rooted subtree.
  int n_tips = tree.n_tip;
  int rem = n_tips % 64;
  uint64_t last_mask = rem ? ((1ULL << rem) - 1) : ~0ULL;

  std::vector<uint64_t> canon(cd.n_words);
  for (int node : tree.postorder) {
    const uint64_t* nd = &node_tips[static_cast<size_t>(node) * cd.n_words];
    bool flip = (nd[0] & 1ULL) != 0;
    for (int w = 0; w < cd.n_words; ++w) {
      canon[w] = flip ? ~nd[w] : nd[w];
    }
    if (flip) canon[cd.n_words - 1] &= last_mask;

    for (int s = 0; s < cd.n_neg_splits; ++s) {
      const uint64_t* split =
          &cd.neg_split_tips[static_cast<size_t>(s) * cd.n_words];
      bool match = true;
      for (int w = 0; w < cd.n_words; ++w) {
        if (canon[w] != split[w]) { match = false; break; }
      }
      if (match) return true;
    }
  }
  return false;
}

// =========================================================================
// Map constraint nodes: find which internal node holds each split
// =========================================================================

// node_displays_split() — the shared "does this node display the split"
// predicate — lives in ts_constraint.h, so the Wagner build and the collapse
// pass answer the question with the same code rather than a lookalike.

// Tightest and highest node displaying `together` | `apart`, or {-1, -1}.
//
// Every node that displays the split covers `together`, so all of them are
// ancestors of LCA(together) and they form one unbroken upward chain: each
// step up adds tips, and the moment a step adds a tip of `apart` the chain
// ends (higher nodes keep it).  So the tight end is the first match in an
// order that visits descendants before ancestors, and the high end is found by
// walking parents from there.  Tips are candidates too, for the single-taxon
// group whose "clade" is the tip itself — tree.postorder holds only internal
// nodes, so scanning it alone left those splits unmapped, which
// regraft_violates_constraint() reads as "already violating".
static void find_displaying_chain(
    const TreeState& tree, const std::vector<uint64_t>& node_tips,
    const uint64_t* together, const uint64_t* apart, int n_words,
    int& lo, int& hi)
{
  lo = -1;
  hi = -1;

  // Tip candidates without scanning the tips: a tip's set is the singleton
  // {t}, so it can only cover `together` when `together` is {t} itself (or,
  // degenerately, empty — then the lowest tip outside `apart` wins).
  int n_together = 0, lone_together = -1;
  for (int w = 0; w < n_words; ++w) {
    if (together[w]) {
      n_together += popcount64(together[w]);
      lone_together = w * 64 + ctz64(together[w]);
    }
  }
  if (n_together == 1) {
    const uint64_t* nd = &node_tips[static_cast<size_t>(lone_together) * n_words];
    if (node_displays_split(nd, together, apart, n_words)) lo = lone_together;
  } else if (n_together == 0) {
    for (int w = 0; w < n_words && lo < 0; ++w) {
      uint64_t free_here = ~apart[w];
      const int lim = tree.n_tip - w * 64;
      if (lim < 64) free_here &= (1ULL << lim) - 1ULL;
      if (free_here) lo = w * 64 + ctz64(free_here);
    }
  }
  if (lo < 0) {
    for (int node : tree.postorder) {
      const uint64_t* nd = &node_tips[static_cast<size_t>(node) * n_words];
      if (node_displays_split(nd, together, apart, n_words)) { lo = node; break; }
    }
  }
  if (lo < 0) return;

  // Walk to the top of the chain.  Bounded by n_node rather than trusting the
  // root to be reachable: impose_one_pass() calls this on trees it is midway
  // through repairing, and a parent-ascending loop over a corrupt parent[] is
  // exactly the hang T-327/T-333 had to be defended against elsewhere.
  hi = lo;
  const int root = tree.n_tip;
  for (int guard = 0; guard < tree.n_node && hi != root; ++guard) {
    const int up = tree.parent[hi];
    if (up < 0 || up >= tree.n_node || up == hi) break;
    const uint64_t* nd = &node_tips[static_cast<size_t>(up) * n_words];
    if (!node_displays_split(nd, together, apart, n_words)) break;
    hi = up;
  }
}

void map_constraint_nodes(const TreeState& tree, ConstraintData& cd)
{
  if (!cd.active) return;

  auto node_tips = compute_node_tips(tree, cd.n_words);

  // For each constraint split, find the node that displays it.
  //
  // T-384: a constraint split is an *unrooted* bipartition A|B, but a clade is
  // a rooted subtree, so the split is displayed whenever EITHER side is a
  // clade.  Exactly one of the two is, except when the split is the root's own
  // bipartition (then both are): for an edge (parent(v), v) with v != root the
  // two sides are desc(v) and its complement.  build_constraint() canonicalises
  // A so that tip 0 is outside it, which makes A the clade side only when tip 0
  // sits on the root's own edge -- true of a tip-0-rooted tree and of nothing
  // else.  Testing the complement as well is what makes this mapping
  // rooting-agnostic, and it costs one extra scan only for splits that used to
  // map to -1 (which regraft_violates_constraint reads as "tree already
  // violates", rejecting every move).  Phase 1 is run to completion first so
  // that every tree which mapped successfully before maps to exactly the same
  // node now.
  //
  // #54: "is a clade" is the free-taxa reading, not set equality -- see
  // node_displays_split().  The chain of displaying nodes is recorded at both
  // ends, because a regraft that must land INSIDE the constrained group may use
  // the whole chain while one that must land outside may not; see
  // regraft_violates_constraint().
  for (int s = 0; s < cd.n_splits; ++s) {
    const uint64_t* ones = &cd.split_tips[static_cast<size_t>(s) * cd.n_words];
    const uint64_t* zeros = &cd.split_zeros[static_cast<size_t>(s) * cd.n_words];
    cd.constraint_complement[s] = 0;

    int lo = -1, hi = -1;
    find_displaying_chain(tree, node_tips, ones, zeros, cd.n_words, lo, hi);
    if (lo < 0) {
      // Phase 2: the tip-0 side is the clade in this rooting.
      find_displaying_chain(tree, node_tips, zeros, ones, cd.n_words, lo, hi);
      if (lo >= 0) cd.constraint_complement[s] = 1;
    }
    cd.constraint_node[s] = lo;
    cd.constraint_node_hi[s] = hi;
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

  // "Inside" means the clip holds a tip of the group that must stay together;
  // "outside", a tip of the group that must stay apart from it.  A clip made
  // only of FREE tips is in neither, and lands here as UNCONSTRAINED — the
  // whole point of the free-taxa reading, and what lets such a clip be
  // regrafted anywhere without breaking the separating edge (#54).  Reading
  // "outside" as ~split, which is what this did before free tips existed,
  // pinned every free tip to the far side of the constraint.
  for (int s = 0; s < cd.n_splits; ++s) {
    const uint64_t* ones =
        &cd.split_tips[static_cast<size_t>(s) * cd.n_words];
    const uint64_t* zeros =
        &cd.split_zeros[static_cast<size_t>(s) * cd.n_words];

    // The clip carries the constraint with it: its own tip set covers one
    // group and holds none of the other.  Wherever it is regrafted, the node
    // at the attachment point has exactly the clip's tip set, so the split
    // stays displayed — and TBR's rerooting of the clip cannot change that,
    // since the set is the same however the subtree hangs.  Testing this
    // first is what unpins a clip that contains the whole displaying chain:
    // the anchor is then inside the clipped subtree, no surviving `below` can
    // be its descendant, and the MUST_INSIDE test below would reject every
    // regraft of a subtree that is in fact free to go anywhere.
    if (node_displays_split(cd.clip_tip_mask.data(), ones, zeros,
                            cd.n_words) ||
        node_displays_split(cd.clip_tip_mask.data(), zeros, ones,
                            cd.n_words)) {
      cd.clip_zones[s] = ClipZone::UNCONSTRAINED;
      continue;
    }

    bool any_inside = false;
    bool any_outside = false;
    for (int w = 0; w < cd.n_words; ++w) {
      if (cd.clip_tip_mask[w] & ones[w]) any_inside = true;
      if (cd.clip_tip_mask[w] & zeros[w]) any_outside = true;
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
        const uint64_t rest = ~cd.clip_tip_mask[w];
        // No width mask needed: ones/zeros carry zeros above tip n_tip - 1.
        if (rest & ones[w]) rest_has_in = true;
        if (rest & zeros[w]) rest_has_out = true;
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

    const int cn_lo = cd.constraint_node[s];
    if (cn_lo < 0) {
      // Constraint genuinely not displayed by the current tree (both sides
      // tested — see map_constraint_nodes).  Reject all moves to avoid
      // entrenching a bad state.
      return true;
    }
    const int cn_hi = cd.constraint_node_hi[s];

    // Which side of the split does cn's subtree hold?  Under the canonical
    // orientation it is the split itself; in a rooting where only the tip-0
    // side is a clade, map_constraint_nodes() maps that side instead (T-384)
    // and the two zones swap: a clip whose tips are all OUTSIDE the split then
    // has to land INSIDE cn's subtree, and vice versa.
    const ClipZone zone_in  = cd.constraint_complement[s]
                            ? ClipZone::MUST_OUTSIDE : ClipZone::MUST_INSIDE;
    const ClipZone zone_out = cd.constraint_complement[s]
                            ? ClipZone::MUST_INSIDE : ClipZone::MUST_OUTSIDE;

    // The two ends of the displaying chain answer two different questions, and
    // each wants the end that permits most (#54; with no free tips the chain
    // is one node long and both reduce to the pre-#54 test):
    //
    //  * a clip that must land INSIDE carries tips of the together-group but
    //    none of the apart-group, so anywhere within the HIGHEST displaying
    //    node keeps that node covering the group and free of the other.  Only
    //    above cn_hi does the enclosing node pick up an apart-group tip.
    //  * a clip that must land OUTSIDE carries apart-group tips, so it may go
    //    anywhere that leaves some displaying node intact — and the TIGHTEST
    //    is the one hardest to contaminate, so it forbids least.
    if (cd.clip_zones[s] == zone_in &&
        !is_ancestor_or_equal(cn_hi, below, cd.dfs_entry, cd.dfs_exit)) {
      return true;
    }
    // Exclude the boundary edge (above_cn, cn): regrafting an outside-only
    // clade just above the constraint clade makes it a sibling of that clade,
    // preserving monophyly.  Only reject if the clade would land *strictly
    // inside* the constraint clade.
    if (cd.clip_zones[s] == zone_out && below != cn_lo &&
        is_ancestor_or_equal(cn_lo, below, cd.dfs_entry, cd.dfs_exit)) {
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
  // A split is violated only when NEITHER side of the bipartition is a clade
  // (T-384): testing the canonical side alone made this "repair" a tree that
  // already displayed every constraint, spending up to n_tip / 4 + 2 arbitrary
  // SPR moves on it.  That mattered most at ts_nni_perturb.cpp's unconditional
  // impose_constraint() call, which runs after every perturbation cycle.
  // "Is a clade" is the free-taxa reading (#54), so a tree that satisfies
  // what `constraint` documents is likewise left alone.  The repair below aims
  // at making the canonical side a clade, which displays the split either way.
  std::vector<int> violated;
  for (int s = 0; s < cd.n_splits; ++s) {
    const uint64_t* ones =
        &cd.split_tips[static_cast<size_t>(s) * n_words];
    const uint64_t* zeros =
        &cd.split_zeros[static_cast<size_t>(s) * n_words];
    int lo = -1, hi = -1;
    find_displaying_chain(tree, node_tips, ones, zeros, n_words, lo, hi);
    if (lo < 0) {
      find_displaying_chain(tree, node_tips, zeros, ones, n_words, lo, hi);
    }
    if (lo < 0) violated.push_back(s);
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

  // Tips that keep node `nd` from displaying the split: those of the
  // together-group it is missing, plus those of the apart-group it holds.
  // Free tips appear in neither, so the repair never moves one (#54) — they
  // may sit on whichever side they already do.
  auto repair_cost = [&](const uint64_t* nd, const uint64_t* ones,
                         const uint64_t* zeros) {
    int cost = 0;
    for (int w = 0; w < n_words; ++w) {
      cost += popcount64(ones[w] & ~nd[w]) + popcount64(zeros[w] & nd[w]);
    }
    return cost;
  };

  for (size_t vi = 0; vi < violated.size(); ++vi) {
    int s = violated[vi];
    const uint64_t* split =
        &cd.split_tips[static_cast<size_t>(s) * n_words];
    const uint64_t* split_out =
        &cd.split_zeros[static_cast<size_t>(s) * n_words];

    // Rebuild bitmasks after previous split's moves
    if (vi > 0) {
      tree.build_postorder();
      node_tips = compute_node_tips(tree, n_words);
    }

    // --- Find best candidate node (fewest misplaced tips) ---
    int best_node = -1;
    int best_cost = tree.n_tip + 1;
    for (int node : tree.postorder) {
      const uint64_t* nd =
          &node_tips[static_cast<size_t>(node) * n_words];
      int cost = repair_cost(nd, split, split_out);
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
      move_out_mask[w] = best_nt[w] & split_out[w];
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
        int cost = repair_cost(nd, split, split_out);
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
