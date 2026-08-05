#ifndef TS_CONSTRAINT_H
#define TS_CONSTRAINT_H

// Topological constraint enforcement (TNT-style locked nodes).
//
// A constraint is a set of splits (bipartitions). A tree satisfies the
// constraint iff every constraint split is displayed — i.e., for each
// split there is an internal node whose subtree tip set matches (after
// accounting for unconstrained taxa that may sit on either side).
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

  // Tip bitmasks: split_tips[i * n_words .. (i+1) * n_words - 1]
  // Canonical: bit 0 (tip 0) is always on the "outside" (= 0).
  std::vector<uint64_t> split_tips;

  // Current mapping: constraint_node[i] = the internal node whose
  // subtree tips match split i in the current tree.
  // -1 if not yet mapped (or the tree does not display split i).
  std::vector<int> constraint_node;

  // Polarity of constraint_node[i] (T-384).  0: the node's descendant tip set
  // is split_tips[i] itself.  1: it is the *complement* of split_tips[i], i.e.
  // the tip-0 side of the bipartition.  A constraint split is an UNROOTED
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

// Build ConstraintData from R-side split bitmask matrix.
// split_matrix: n_splits x n_tips, each row is 0/1 indicating split membership.
// The matrix is canonicalized so tip 0 is always "outside" (= 0).
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
