#ifndef TS_CONSTRAINT_H
#define TS_CONSTRAINT_H

// Topological constraint enforcement (TNT-style locked nodes).
//
// A constraint is a set of splits.  Each split names two disjoint groups of
// tips — the "1" group and the "0" group of one constraint character — and any
// remaining tips are FREE: coded `?`, or absent from the constraint phyDat
// altogether.  A tree satisfies the split iff some edge separates the 1 group
// from the 0 group, free tips falling on either side.  That is the contract
// `?MaximizeParsimony`'s `constraint` argument documents, and since
// agent-issues/TreeSearch#54 it is the one enforced here: a node DISPLAYS the
// split when its descendant tip set is a superset of one group and disjoint
// from the other.
//
// Requiring the tip set to EQUAL a group (the pre-#54 test) is strictly
// stronger.  It never returned a wrong answer, but a start tree that satisfied
// the documented contract without making either group an exact clade mapped to
// no node at all, which regraft_violates_constraint() reads as "the tree
// already violates" — freezing the replicate on its start.
//
// Implementation:
//   1. At init: store constraint splits as tip bitmasks.
//   2. Map each split to the internal node that currently represents it.
//   3. Compute DFS timestamps for O(1) descendant queries.
//   4. Per TBR/SPR clip: classify each constraint as MUST_INSIDE,
//      MUST_OUTSIDE, or UNCONSTRAINED (based on clip tip set).
//   5. Per candidate regraft edge: O(1) check per active constraint.
//   6. After each accepted move: remap constraint nodes (cheap).

#include "ts_data.h"
#include "ts_tree.h"
#include <vector>
#include <cstdint>

namespace ts {

enum class ClipZone : int {
  UNCONSTRAINED = 0,  // clip tips don't interact with this split
  MUST_INSIDE   = 1,  // all clip tips inside split → regraft inside
  MUST_OUTSIDE  = 2,  // all clip tips outside split → regraft outside
  FORBIDDEN     = 3   // clip straddles split AND rest straddles → no valid regraft
};

struct ConstraintData {
  bool active = false;
  int n_splits = 0;
  int n_words = 0;           // ceil(n_tips / 64)

  // Tip bitmasks: split_tips[i * n_words .. (i+1) * n_words - 1].
  // The tips that must end up TOGETHER, on one side of some edge.
  // Canonical: bit 0 (tip 0) is always on the "outside" (= 0).
  std::vector<uint64_t> split_tips;

  // The tips that must end up on the OTHER side of that edge, same layout.
  // Disjoint from split_tips; the two need NOT be complements — a tip in
  // neither mask is free to fall on either side (#54).  When the caller
  // supplies no free tips this is exactly ~split_tips, and every check below
  // reduces to the pre-#54 exact-clade test.
  std::vector<uint64_t> split_zeros;

  // Current mapping: constraint_node[i] = the TIGHTEST node (tip or internal)
  // that displays split i in the current tree — its descendant tip set covers
  // one of the two groups and avoids the other.
  // -1 if not yet mapped (or the tree does not display split i).
  std::vector<int> constraint_node;

  // The HIGHEST node that displays split i, in the same polarity as
  // constraint_node[i]; equal to it when no free tip sits directly above.
  // The displaying nodes form an unbroken chain from constraint_node[i] up to
  // this one (each step adds only free tips), so the two ends are all a
  // regraft test needs — see regraft_violates_constraint(), which uses this
  // end for "must land inside" and the tight end for "must land outside".
  // -1 exactly when constraint_node[i] is.
  std::vector<int> constraint_node_hi;

  // Polarity of constraint_node[i] (T-384).  0: the node's descendant tip set
  // covers split_tips[i] and avoids split_zeros[i].  1: the other way round —
  // it covers split_zeros[i], the tip-0 side of the bipartition, and avoids
  // split_tips[i].  A constraint split is an UNROOTED
  // bipartition, so a tree displays it whenever EITHER side is a rooted clade,
  // and which side that is depends on the rooting alone -- see
  // map_constraint_nodes().  Consumers that treat constraint_node[i] as "the
  // inside clade" must swap MUST_INSIDE/MUST_OUTSIDE when this is 1.
  // Written by every map_constraint_nodes() call, alongside constraint_node.
  std::vector<char> constraint_complement;

  // DFS timestamps for O(1) descendant checks.
  // Node u is ancestor of v iff dfs_entry[u] <= dfs_entry[v]
  //                           && dfs_exit[u]  >= dfs_exit[v].
  std::vector<int> dfs_entry;
  std::vector<int> dfs_exit;

  // Post-hoc fallback: constraint as a DataSet + expected Fitch score.
  // Used for sector/fuse where full topology reconstruction makes
  // the locked-node approach impractical.
  DataSet posthoc_data;
  int expected_score = 0;
  bool has_posthoc = false;

  // Per-clip workspace (reused across clips, sized at init)
  std::vector<ClipZone> clip_zones;          // [n_splits]
  std::vector<uint64_t> clip_tip_mask;       // [n_words]
};

// Build ConstraintData from R-side split membership matrix.
// split_matrix: n_splits x n_tips, column-major.  Element [s, t] is
//   1  tip t is in split s's "together" group;
//   0  tip t is in split s's "apart" group;
//   anything else (NA_INTEGER, as .PrepareConstraint() writes for a `?`-coded
//      or unconstrained taxon) — tip t is FREE, and may fall on either side.
// A pure 0/1 matrix therefore means "no free tips", i.e. the exact-clade
// reading that predates #54; callers that build one by hand keep it.
// The two groups are swapped where needed so that tip 0 is never in
// split_tips ("outside" the canonical side) — the same invariant the Wagner
// and pool paths have always relied on, and harmless because a split is an
// unrooted bipartition whose two groups are interchangeable.
ConstraintData build_constraint(
    const int* split_matrix, int n_splits, int n_tips);

// Also set up the post-hoc fallback DataSet from R-side phyDat components.
void build_constraint_posthoc(
    ConstraintData& cd,
    const double* contrast_r, int n_tokens, int n_states,
    const int* tip_data_r, int n_tips, int n_patterns,
    const int* weight_r,
    const char** levels_r,
    int expected_score);

// --- Node mapping and DFS timestamps ---

// Does the edge above a node whose descendant tip set is `nd` separate
// `together` from `apart`?  It does when the set covers every tip of
// `together` and holds none of `apart`; the tips in neither group are free and
// are not looked at.  With `apart` the exact complement of `together` the two
// conditions force set equality, which is the exact-clade test this replaced.
//
// THE definition of "displays a constraint split", shared by every entry point
// that has to decide it: the search/TBR mapping (map_constraint_nodes), the
// Wagner build's own check (wagner_tree_displays_constraint, ts_wagner.cpp),
// and the collapse pass's branch protection (ts_collapse_pool, ts_rcpp.cpp).
// They must not drift apart: the stricter of any two would reject trees
// another searches happily, or accept ones it will not move from.
inline bool node_displays_split(
    const uint64_t* nd, const uint64_t* together, const uint64_t* apart,
    int n_words)
{
  for (int w = 0; w < n_words; ++w) {
    if ((together[w] & ~nd[w]) != 0ULL) return false;  // a required tip missing
    if ((apart[w] & nd[w]) != 0ULL) return false;      // an excluded tip present
  }
  return true;
}

// Find which internal node holds each constraint split in the current tree.
// Must be called after each accepted move and at search init.
void map_constraint_nodes(const TreeState& tree, ConstraintData& cd);

// Compute DFS entry/exit timestamps for the current tree.
// Must be called after map_constraint_nodes (or any topology change).
void compute_dfs_timestamps(const TreeState& tree, ConstraintData& cd);

// Combined: remap + recompute DFS. Convenience function.
void update_constraint(const TreeState& tree, ConstraintData& cd);

// --- Per-clip classification ---

// Compute the tip bitmask of the subtree rooted at clip_node.
void compute_clip_tip_mask(const TreeState& tree, int clip_node,
                           std::vector<uint64_t>& mask);

// Classify each constraint split for this clip.
// Populates cd.clip_zones[].
void classify_clip_constraints(const TreeState& tree, int clip_node,
                               ConstraintData& cd);

// --- Per-candidate check ---

// Returns true if regrafting onto the edge whose child endpoint is `below`
// would violate any active constraint, given the current clip_zones
// classification.  (Only `below` is needed: the parent endpoint of the target
// edge never changes which side of a constraint clade the clip lands on.)
// Uses DFS timestamps for O(1) descendant test per constraint.
//
// Screening only.  It classifies against the PRE-clip mapping, so a TBR
// rerooting can still break a split it passed; every caller must re-verify the
// applied move (map_constraint_nodes + constraint_node[s] >= 0) before
// accepting it.
bool regraft_violates_constraint(int below,
                                 const ConstraintData& cd);

// Build ConstraintData directly from pre-canonicalized split bitsets.
// `split_bits` is contiguous: n_splits * words_per_split uint64_t values.
// Splits must already be canonicalized (bit 0 clear).
// No posthoc DataSet is built (has_posthoc = false).
ConstraintData build_constraint_from_bitsets(
    const uint64_t* split_bits, int n_splits,
    int words_per_split, int n_tips);

// --- Post-hoc check (for sector/fuse) ---

// Full Fitch check: score the tree against the constraint DataSet.
// Returns true if constraint is violated.
bool violates_constraint_posthoc(const TreeState& tree,
                                 const ConstraintData& cd);

// --- Post-hoc repair ---

// Compute per-node subtree tip bitmasks via postorder traversal.
// Returns array of size n_node * n_words.
// For tips: bit[t] = 1. For internal nodes: OR of children.
std::vector<uint64_t> compute_node_tips(const TreeState& tree, int n_words);

// Repair constraint violations by minimal SPR moves.
//
// HEURISTIC, and it can fail: a pass bails out when the repair needs more than
// n_tip / 4 + 2 moves, individual moves are skipped when they would corrupt the
// tree (see try_move), and nothing guarantees the fixed-point is reached within
// the n_splits + 1 pass cap.  The return value does NOT distinguish "repaired"
// from "gave up" -- so EVERY caller must re-verify (map_constraint_nodes, then
// constraint_node[s] >= 0 for all s) and discard the tree if it still violates.
//
// update_constraint() has been called on return. Caller must rescore.
// Returns the number of SPR moves performed (0 if tree was already valid).
int impose_constraint(TreeState& tree, ConstraintData& cd);

} // namespace ts

#endif // TS_CONSTRAINT_H
