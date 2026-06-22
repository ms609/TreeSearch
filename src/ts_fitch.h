#ifndef TS_FITCH_H
#define TS_FITCH_H

// Bit-packed Fitch parsimony scoring.
//
// Standard (non-inapplicable) Fitch downpass + uppass for equal weights:
// processes 64 characters per block using bitwise AND/OR + popcount.

#include "ts_data.h"
#include "ts_simd.h"
#include "ts_tree.h"
#include <cmath>
#include <vector>

namespace ts {

// Full Fitch downpass: compute prelim state sets and local_cost at all nodes.
// Returns the total EW parsimony score.
int fitch_downpass(TreeState& tree, const DataSet& ds);

// Full Fitch uppass: compute final_ state sets at all internal nodes.
// Must be called after fitch_downpass(). Tips retain their observed states.
void fitch_uppass(TreeState& tree, const DataSet& ds);

// Full two-pass score: downpass + uppass. Returns the EW score.
// After this, tree.prelim, tree.final_, and tree.local_cost are all current.
int fitch_score(TreeState& tree, const DataSet& ds);

// Score a single block at a single internal node.
// Returns the number of characters (out of 64) that needed a union (= step).
// Writes the result state into node_state[0..n_states-1].
int fitch_downpass_node(
    const uint64_t* left_state,
    const uint64_t* right_state,
    uint64_t* node_state,
    int n_states,
    uint64_t active_mask);

// --- Incremental scoring for SPR ---

// Incremental downpass after clipping: walk from start_node rootward,
// recomputing prelim and local_cost. Stops when prelim stabilizes.
// Returns the length delta (new_score - old_score for the main tree).
// Saves old states to tree.clip_undo_stack for restoration.
// cs_delta (optional, IW dirty-region): if non-null, accumulates the per-pattern
// char_steps change over the visited path (new local_cost bits +1, old -1); the
// function zeroes it first. Lets the IW path derive divided_steps incrementally
// (full_char_steps + cs_delta - nx) instead of an O(n_node) extract_char_steps.
int fitch_incremental_downpass(TreeState& tree, const DataSet& ds,
                               int start_node,
                               std::vector<int>* cs_delta = nullptr);

// Incremental uppass after clip: recompute final_ for nodes whose
// ancestor's final states changed. Propagates from root downward,
// visiting only nodes in the changed region.
void fitch_incremental_uppass(TreeState& tree, const DataSet& ds,
                              int start_node);

// Dirty-set rescore after an SPR move (T-300).
//
// Recomputes prelim and local_cost for every node on the union of paths
// start_a -> root and start_b -> root, visiting each node exactly once in
// postorder.  start_a and start_b are the two clip endpoints whose children
// changed after apply_tbr_move (typically nz = clip grandparent and
// nx = regraft point).
//
// Caller must call tree.build_postorder_prealloc() first so that
// tree.postorder reflects the post-move topology.
//
// Returns the EW length delta: actual = prior_score + delta.
// For IW/profile, ignore the return value and use extract_char_steps +
// compute_weighted_score after this call (local_cost is correct).
int fitch_dirty_downpass(TreeState& tree, const DataSet& ds,
                         int start_a, int start_b);

// Companion uppass for fitch_dirty_downpass.  Recomputes final_ for nodes
// whose ancestor's final_ may have changed, seeded from the same start
// points.  Propagates downward.
void fitch_dirty_uppass(TreeState& tree, const DataSet& ds,
                        int start_a, int start_b);

// --- NA-aware dirty-set incremental rescore (T-300 NA variant) ---
//
// Same dirty-set approach as fitch_dirty_downpass / fitch_dirty_uppass but
// handles inapplicable-bearing blocks via the NA-aware Pass 1 / Pass 2
// logic.  Used for the SPR accept path under has_inapplicable to avoid
// full_rescore.  The return value is the EW length delta for standard
// blocks only — NA block step counts require Pass 3, so call
// fitch_na_pass3_score(tree, ds) on the updated state to obtain the
// authoritative score.
int fitch_na_dirty_downpass(TreeState& tree, const DataSet& ds,
                             int start_a, int start_b);

void fitch_na_dirty_uppass(TreeState& tree, const DataSet& ds,
                            int start_a, int start_b);

// Indirect tree length calculation: given the clipped subtree's basal
// state set (prelim of clip_node) and a candidate destination edge (A, D),
// compute the length increase from joining them.
// Returns the step count increase (0 = free join, positive = extra steps).
// `clip_prelim` points to total_words entries for the clip node.
int fitch_indirect_length(const uint64_t* clip_prelim,
                          const TreeState& tree,
                          const DataSet& ds,
                          int node_a, int node_d);

// Early-termination variant: returns as soon as extra_steps >= cutoff.
// Returns exact result if below cutoff, or a value >= cutoff otherwise.
int fitch_indirect_length_bounded(const uint64_t* clip_prelim,
                                  const TreeState& tree,
                                  const DataSet& ds,
                                  int node_a, int node_d,
                                  int cutoff);

// Precomputed-vroot variant: uses pre-computed virtual root states instead
// of reading final_ from tree. vroot points to total_words entries.
// Supports early termination via cutoff (INT_MAX = no cutoff).
int fitch_indirect_length_cached(const uint64_t* clip_prelim,
                                 const uint64_t* vroot,
                                 const DataSet& ds,
                                 int cutoff);

// Compute the EXACT per-node insertion edge set for every non-root node D:
//   edge_set[D] = combine(prelim[D], up[D])   (per-character intersect-else-union)
// where up[D] is the directional Fitch up-message
//   up[D] = combine(up[parent(D)], prelim[sibling(D)])   (root degree-2 vertex:
//   up[child] = prelim[other child]).
// Inserting a clip with downpass set T on the edge above D then costs exactly
// #chars where T & edge_set[D] == 0 (use fitch_indirect_length_cached with
// edge_set[D] as the vroot).  This is the CORRECT replacement for the
// union-of-finals (final_[A] | final_[D]) approximation, which undercounts.
// Requires a current downpass (prelim).  `edge_set` is sized n_node*total_words;
// the root entry is left unspecified.  `up` (n_node*total_words scratch up-
// message buffer) and `pre` (preorder node list) are caller-owned and reused
// across calls; pass the same vectors on every call to avoid per-call
// allocation and zero-fill.  They are size-ensured (non-zeroing) and fully
// overwritten internally, so the caller need not initialize them.
void compute_insertion_edge_sets(const TreeState& tree, const DataSet& ds,
                                 std::vector<uint64_t>& edge_set,
                                 std::vector<uint64_t>& up,
                                 std::vector<int>& pre);

// --- Flat EW specializations (skip weight/upweight overhead) ---
//
// These use FlatBlock metadata (1 cache line for all blocks) instead of
// the full CharBlock array. Valid only when all blocks have weight==1
// and no upweight_mask is set (normal search, not ratchet).

// 3-operand bounded: for SPR candidates (reads final_ directly).
int fitch_indirect_bounded_flat(const uint64_t* clip_prelim,
                                const TreeState& tree,
                                const DataSet& ds,
                                int node_a, int node_d,
                                int cutoff);

// 2-operand cached: for TBR rerooting candidates (pre-computed vroot).
int fitch_indirect_cached_flat(const uint64_t* clip_prelim,
                               const uint64_t* vroot,
                               const DataSet& ds,
                               int cutoff);

// NA-aware bounded: for SPR candidates (reads final_ + subtree_actives).
int fitch_na_indirect_bounded_flat(const uint64_t* clip_prelim,
                                   const uint64_t* clip_actives,
                                   const TreeState& tree,
                                   const DataSet& ds,
                                   int node_a, int node_d,
                                   int cutoff);

// NA-aware cached: for TBR rerooting candidates (pre-computed vroot).
int fitch_na_indirect_cached_flat(const uint64_t* clip_prelim,
                                  const uint64_t* clip_actives,
                                  const uint64_t* vroot,
                                  const uint64_t* below_actives,
                                  const DataSet& ds,
                                  int cutoff);

// --- 4-wide TBR rerooting batch (T-245) ---
//
// Process 4 regraft candidates simultaneously in the TBR rerooting inner loop.
// Each block iteration issues 4 independent loads from vroot_cache, letting
// the out-of-order CPU serve them concurrently and hide L2 latency.
//
// out[i]: extra steps for candidate i; may equal or exceed cutoff.
// Exits early when ALL 4 exceed cutoff after any block.
//
// Requires: ds.all_weight_one, no block has upweight_mask (normal EW search,
//           not ratchet). Use the use_flat guard in the caller.

// EW flat 4-wide batch (no inapplicable characters).
void fitch_indirect_cached_flat_x4(
    const uint64_t* clip_prelim,
    const uint64_t* vroot0, const uint64_t* vroot1,
    const uint64_t* vroot2, const uint64_t* vroot3,
    const DataSet& ds, int cutoff, int out[4]);

// NA-aware flat 4-wide batch (mixed standard + inapplicable blocks).
// ba0..ba3: below_actives_cache rows (1 uint64 per block) for each candidate.
void fitch_na_indirect_cached_flat_x4(
    const uint64_t* clip_prelim,
    const uint64_t* clip_actives,
    const uint64_t* vroot0, const uint64_t* vroot1,
    const uint64_t* vroot2, const uint64_t* vroot3,
    const uint64_t* ba0, const uint64_t* ba1,
    const uint64_t* ba2, const uint64_t* ba3,
    const DataSet& ds, int cutoff, int out[4]);

// --- Inapplicable (NA) three-pass scoring ---

// Full three-pass score for datasets with inapplicable characters.
// Handles both standard blocks (one-pass Fitch) and inapplicable blocks
// (Brazeau et al. three-pass algorithm). Returns the total EW score.
int fitch_na_score(TreeState& tree, const DataSet& ds);

// --- Incremental NA-aware scoring for SPR/TBR ---

// NA-aware incremental first downpass. Walks rootward from start_node,
// computing prelim with NA-aware logic for inapplicable blocks and
// standard Fitch for standard blocks. Also maintains subtree_actives.
// Returns the length delta for standard blocks. NA blocks require
// fitch_na_pass3_score() for exact step counts.
int fitch_na_incremental_downpass(TreeState& tree, const DataSet& ds,
                                   int start_node);

// NA-aware incremental first uppass. Recomputes final_ for nodes in
// the dirty region, using NA-aware uppass logic for inapplicable blocks.
// Also updates tip down2 and subtree_actives for affected tips.
void fitch_na_incremental_uppass(TreeState& tree, const DataSet& ds,
                                  int clip_ancestor);

// Full Pass 3 (second downpass) on a divided tree. Computes down2 for
// all internal nodes, counts steps for both standard and NA blocks.
// Requires Passes 1+2 to be current (from full or incremental scoring).
// Returns the total EW score.
int fitch_na_pass3_score(TreeState& tree, const DataSet& ds);

// NA-aware indirect length calculation. For standard blocks, identical to
// fitch_indirect_length. For NA blocks, suppresses steps where either the
// clip subtree or the edge-below subtree has no applicable tips.
// clip_actives: subtree_actives for the clip subtree (total_words entries).
int fitch_na_indirect_length(
    const uint64_t* clip_prelim,
    const uint64_t* clip_actives,
    const TreeState& tree,
    const DataSet& ds,
    int node_a, int node_d);

// NA-aware bounded indirect length (early termination at cutoff).
int fitch_na_indirect_length_bounded(
    const uint64_t* clip_prelim,
    const uint64_t* clip_actives,
    const TreeState& tree,
    const DataSet& ds,
    int node_a, int node_d,
    int cutoff);

// NA-aware cached indirect length (pre-computed vroot + below_actives).
// below_actives: per-edge OR of applicable subtree_actives[D] (1 uint64 per block).
int fitch_na_indirect_length_cached(
    const uint64_t* clip_prelim,
    const uint64_t* clip_actives,
    const uint64_t* vroot,
    const uint64_t* below_actives,
    const DataSet& ds,
    int cutoff);

// NA-aware indirect IW length. Same NA-suppression as above, for IW scoring.
double indirect_na_iw_length(
    const uint64_t* clip_prelim,
    const uint64_t* clip_actives,
    const TreeState& tree, const DataSet& ds,
    int node_a, int node_d,
    double base_iw,
    const std::vector<double>& iw_delta);

// NA-aware bounded IW indirect length.
double indirect_na_iw_length_bounded(
    const uint64_t* clip_prelim,
    const uint64_t* clip_actives,
    const TreeState& tree, const DataSet& ds,
    int node_a, int node_d,
    double base_iw,
    const std::vector<double>& iw_delta,
    double cutoff);

// NA-aware cached IW indirect length.
double indirect_na_iw_length_cached(
    const uint64_t* clip_prelim,
    const uint64_t* clip_actives,
    const uint64_t* vroot,
    const uint64_t* below_actives,
    const DataSet& ds,
    double base_iw,
    const std::vector<double>& iw_delta,
    double cutoff);

// --- Per-character step extraction ---

// Extract per-pattern step counts from local_cost masks (standard blocks)
// and down2 states (NA blocks) after a full scoring pass.
// char_steps must be pre-sized to ds.n_patterns and zero-initialized.
void extract_char_steps(const TreeState& tree, const DataSet& ds,
                        std::vector<int>& char_steps);

// --- Unified scoring ---

// Score non-hierarchy characters via Fitch (EW, IW, or profile).
// For HSJ mode, treats the DataSet as EW (non-hierarchy chars only).
// Does NOT include HSJ hierarchy block contributions.
double fitch_score_ew(TreeState& tree, const DataSet& ds);

// Score the tree using the appropriate scoring mode.
// For HSJ: calls fitch_score_ew() + HSJ DP on hierarchy blocks.
// For EW/IW/PROFILE: delegates to fitch_score_ew().
double score_tree(TreeState& tree, const DataSet& ds);

// Compute IW score from per-character step counts.
double compute_iw(const DataSet& ds, const std::vector<int>& char_steps);

// --- IW indirect calculation ---

// Evaluate a regraft candidate under IW.
// base_iw: precomputed IW of the divided tree (without reconnection cost).
// iw_delta: precomputed marginal IW cost if pattern p gains one more step.
// Returns the total IW score of the candidate tree.
double indirect_iw_length(
    const uint64_t* clip_prelim,
    const TreeState& tree, const DataSet& ds,
    int node_a, int node_d,
    double base_iw,
    const std::vector<double>& iw_delta);

// Early-termination IW variant: returns as soon as candidate >= cutoff.
double indirect_iw_length_bounded(
    const uint64_t* clip_prelim,
    const TreeState& tree, const DataSet& ds,
    int node_a, int node_d,
    double base_iw,
    const std::vector<double>& iw_delta,
    double cutoff);

// Precomputed-vroot IW variant with early termination.
double indirect_iw_length_cached(
    const uint64_t* clip_prelim,
    const uint64_t* vroot,
    const DataSet& ds,
    double base_iw,
    const std::vector<double>& iw_delta,
    double cutoff);

// 4-wide IW batch for TBR rerooting (pure-IW, no inapplicables): the IW analog
// of fitch_indirect_cached_flat_x4. 4 independent any_hit_reduce load streams +
// 4 independent double accumulators (each candidate's weighted sum keeps its own
// scalar add order => bit-identical to indirect_iw_length_cached per candidate;
// the batch only shares the batch-start cutoff, which affects early-exit not the
// final < best comparison => byte-identical search outcome).
void indirect_iw_cached_flat_x4(
    const uint64_t* clip_prelim,
    const uint64_t* vroot0, const uint64_t* vroot1,
    const uint64_t* vroot2, const uint64_t* vroot3,
    const DataSet& ds, double base_iw,
    const std::vector<double>& iw_delta,
    double cutoff, double out[4]);

// 4-wide NA-IW batch for TBR rerooting (inapplicable characters + implied
// weights): the NA analog of indirect_iw_cached_flat_x4. Fuses the NA
// active-mask candidate logic of fitch_na_indirect_cached_flat_x4 with the
// per-candidate iw_delta ctz-gather of indirect_iw_cached_flat_x4, so each
// es{k} reproduces indirect_na_iw_length_cached's scalar add order exactly
// (bit-identical per candidate; the shared all-4-exceed-cutoff bail only
// affects early-exit, not the final < best comparison). ba0..ba3 are
// below_actives_cache rows (1 uint64 per block) for each candidate.
void indirect_na_iw_cached_flat_x4(
    const uint64_t* clip_prelim,
    const uint64_t* clip_actives,
    const uint64_t* vroot0, const uint64_t* vroot1,
    const uint64_t* vroot2, const uint64_t* vroot3,
    const uint64_t* ba0, const uint64_t* ba1,
    const uint64_t* ba2, const uint64_t* ba3,
    const DataSet& ds, double base_iw,
    const std::vector<double>& iw_delta,
    double cutoff, double out[4]);

// Precompute iw_delta[p] = marginal cost of one additional step for pattern p.
// divided_steps: per-pattern step counts of the divided tree.
void precompute_iw_delta(const DataSet& ds,
                         const std::vector<int>& divided_steps,
                         std::vector<double>& iw_delta);

// --- Profile parsimony scoring ---

// Compute profile parsimony score from per-character step counts.
// Looks up info_amounts[total_steps, pattern] for each pattern.
double compute_profile(const DataSet& ds, const std::vector<int>& char_steps);

// Precompute marginal profile cost of one additional step per pattern.
void precompute_profile_delta(const DataSet& ds,
                               const std::vector<int>& divided_steps,
                               std::vector<double>& delta);

// --- Weighted scoring dispatch (IW or profile) ---

// Dispatch to compute_iw or compute_profile based on ds.scoring_mode.
double compute_weighted_score(const DataSet& ds,
                               const std::vector<int>& char_steps);

// Dispatch to precompute_iw_delta or precompute_profile_delta.
void precompute_weighted_delta(const DataSet& ds,
                                const std::vector<int>& divided_steps,
                                std::vector<double>& delta);

} // namespace ts

#endif // TS_FITCH_H
