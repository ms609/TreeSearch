#ifndef TS_FITCH_H
#define TS_FITCH_H

// Bit-packed Fitch parsimony scoring.
//
// Standard (non-inapplicable) Fitch downpass + uppass for equal weights:
// processes 64 characters per block using bitwise AND/OR + popcount.

#include "ts_data.h"
#include "ts_tree.h"

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
int fitch_incremental_downpass(TreeState& tree, const DataSet& ds,
                               int start_node);

// Incremental uppass after clip: recompute final_ for nodes whose
// ancestor's final states changed. Propagates from root downward,
// visiting only nodes in the changed region.
void fitch_incremental_uppass(TreeState& tree, const DataSet& ds,
                              int start_node);

// Indirect tree length calculation: given the clipped subtree's basal
// state set (prelim of clip_node) and a candidate destination edge (A, D),
// compute the length increase from joining them.
// Returns the step count increase (0 = free join, positive = extra steps).
// `clip_prelim` points to total_words entries for the clip node.
int fitch_indirect_length(const uint64_t* clip_prelim,
                          const TreeState& tree,
                          const DataSet& ds,
                          int node_a, int node_d);

} // namespace ts

#endif // TS_FITCH_H
