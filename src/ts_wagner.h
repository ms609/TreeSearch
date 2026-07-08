#ifndef TS_WAGNER_H
#define TS_WAGNER_H

// Wagner tree construction: greedy taxon addition for parsimony.
//
// Builds a tree by adding taxa one at a time at the position that
// minimizes the parsimony score. Used to produce starting trees for
// heuristic search (TBR, ratchet, etc.).

#include "ts_data.h"
#include "ts_tree.h"
#include "ts_constraint.h"
#include <vector>

namespace ts {

struct WagnerResult {
  double score;
};

// Build a Wagner tree by greedy addition.
// `tree` is populated in place (overwritten).
// `addition_order`: tip indices in insertion order (size n_tips).
//   If empty, uses sequential order 0..n_tips-1.
WagnerResult wagner_tree(TreeState& tree, const DataSet& ds,
                         const std::vector<int>& addition_order,
                         ConstraintData* cd = nullptr);

// Build a random-addition-sequence Wagner tree.
// Uses R's RNG (respects set.seed).
WagnerResult random_wagner_tree(TreeState& tree, const DataSet& ds,
                                ConstraintData* cd = nullptr);

// Criterion for biasing taxon addition order.
//
// RANDOM      — uniform random (existing behaviour)
// GOLOBOFF    — prioritise taxa with more non-ambiguous parsimony-informative
//               characters (Goloboff 2014 "informative" addition sequence)
// ENTROPY     — prioritise taxa with more specific (lower-entropy) state
//               codings; score(t) = Σ_c (n_states_c - |state_set of t at c|)
// CLOSEST     — TNT-style "cas": greedy nearest-neighbour addition. At each
//               step, taxa close (in pairwise Fitch character distance) to
//               the already-placed set are preferred.
// FURTHEST    — TNT-style "fas": greedy furthest-neighbour addition (the
//               mirror image of CLOSEST).
// INFORMATIVE — TNT-style "ias": at each step, prefer the not-yet-placed
//               taxon that renders the largest number of characters
//               parsimony-informative given the taxa already placed.
//               Targets missing-data-heavy / non-overlapping-block matrices
//               (Goloboff's stated use case for `ias`).
//
// NOTE: TNT also has `sas` ("largest best/worst insertion-cost difference"),
// which requires evaluating exact insertion cost against the *actual*
// growing tree at each step (not just a prefix-set property) — a materially
// different construction than the other four. It is deliberately NOT
// implemented here; see dev notes for T-addseq before adding it, so it is
// not silently approximated by a distance-based proxy and mislabelled.
//
// GOLOBOFF/ENTROPY use softmax-weighted sampling WITHOUT replacement over a
// single static per-tip score vector (computed once, independent of
// placement order). CLOSEST/FURTHEST/INFORMATIVE are inherently sequential
// (their criterion depends on which taxa are already placed), so they are
// built step-by-step instead, softmax-sampling one taxon at a time from a
// freshly recomputed score over the remaining candidates (see
// greedy_addseq_order()). All modes share the same `temperature` knob:
// 0 → greedy argmax; large → near-uniform. Temperature is applied to scores
// normalised to [0, 1] so that behaviour is consistent across datasets
// regardless of character count.
enum class WagnerBias {
  RANDOM      = 0,
  GOLOBOFF    = 1,
  ENTROPY     = 2,
  CLOSEST     = 3,
  FURTHEST    = 4,
  INFORMATIVE = 5
};

struct BiasedWagnerParams {
  WagnerBias bias        = WagnerBias::RANDOM;
  double     temperature = 1.0;   // in [0, ∞); 0 → greedy
};

// Build a biased-addition-sequence Wagner tree.
// Falls through to random_wagner_tree() when bias == RANDOM.
WagnerResult biased_wagner_tree(TreeState& tree, const DataSet& ds,
                                const BiasedWagnerParams& params,
                                ConstraintData* cd = nullptr);

// Compute per-tip Goloboff informative scores (exported for diagnostics).
// score[t] = number of non-ambiguous characters for tip t.
std::vector<double> wagner_goloboff_scores(const DataSet& ds);

// Compute per-tip entropy scores (exported for diagnostics).
// score[t] = Σ_c (n_states_c - |state_set of t at c|)
std::vector<double> wagner_entropy_scores(const DataSet& ds);

// Pairwise inter-taxon character distance matrix (n_tips × n_tips, symmetric,
// zero diagonal), flattened row-major. distance(i, j) = the two-taxon Fitch
// parsimony score if i and j were forced to be sisters in isolation (i.e. the
// number of characters where i and j share no state, weighted like normal
// scoring). Used by CLOSEST/FURTHEST addition-sequence ordering.
// Returns all-zero when ds.total_words == 0 (no informative characters).
std::vector<double> wagner_pairwise_distances(const DataSet& ds);

// Build a full addition order (permutation of 0..n_tips-1) for the
// CLOSEST / FURTHEST / INFORMATIVE criteria via sequential greedy
// construction: the first taxon is chosen uniformly at random (no prefix
// yet to compare against); each subsequent taxon is softmax-sampled from a
// criterion recomputed against the current prefix (nearest/furthest
// neighbour distance, or parsimony-informativeness gain). `temperature`
// has the same meaning as in BiasedWagnerParams (0 = greedy argmax).
// Falls back to a uniform random permutation when ds.total_words == 0.
// Uses R's RNG (respects set.seed) via ts::thread_safe_unif, exactly like
// softmax_sample_order().
std::vector<int> greedy_addseq_order(const DataSet& ds, WagnerBias bias,
                                      double temperature);

// Build a purely random tree topology (no character data used).
// Inserts tips in random order at random edges. The resulting tree has
// valid topology + tip states loaded but is NOT scored (prelim/final
// for internal nodes are zeroed). Caller should call score_tree().
//
// Goloboff (2014) found that random starting trees sometimes reach basins
// inaccessible to Wagner trees, justifying inclusion in a strategy mix.
void random_topology_tree(TreeState& tree, const DataSet& ds);

// Build a random tree topology that satisfies topological constraints.
// Constructs the constraint backbone (one node per constraint split),
// then randomly resolves all multifurcations by uniform random binary
// insertion.  Like random_topology_tree(), the result is NOT scored.
//
// Falls back to random_topology_tree() if no constraints are active.
void random_constrained_tree(TreeState& tree, const DataSet& ds,
                             ConstraintData& cd);

// --- Low-level helpers (used by prune-reinsert and Wagner) ---

// Allocate a full-sized TreeState for n_tips taxa and load tip states.
void init_wagner_state(TreeState& tree, const DataSet& ds);

// Insert a new tip at edge (above, below), creating new_internal between them.
void insert_tip_at_edge(TreeState& tree, int tip, int new_internal,
                        int above, int below);

// Incremental two-pass Fitch rescore after insert_tip_at_edge().
// Returns the score delta (always positive during construction).
int wagner_incremental_rescore(TreeState& tree, const DataSet& ds,
                               int new_internal);

} // namespace ts

#endif // TS_WAGNER_H
