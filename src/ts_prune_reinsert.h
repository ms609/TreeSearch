#ifndef TS_PRUNE_REINSERT_H
#define TS_PRUNE_REINSERT_H

// Taxon pruning-reinsertion perturbation for parsimony search.
//
// Perturbation strategy complementary to the ratchet (weight-space) and
// NNI-perturbation (topology-space).  Operates in taxon-composition space:
// drop a subset of leaves → TBR on the reduced tree → greedily re-add
// the dropped taxa via Wagner insertion → TBR polish on the full tree.
//
// Some taxa can be trapped in suboptimal positions that no single TBR move
// fixes (coordinated relocation required).  Temporarily removing them
// lets the backbone restructure, and greedy re-addition may land them in
// a different basin of attraction.

#include "ts_data.h"
#include "ts_tree.h"
#include "ts_constraint.h"
#include <functional>
#include <vector>

namespace ts {

// Forward declarations (defined in ts_pool.h)
struct SplitFrequencyTable;

// Tip selection strategy.
// RANDOM      — uniform random (baseline)
// INSTABILITY — sample weighted by positional instability in the pool
//               (tips whose parent-edge split is rare across pool trees
//               are more likely to be dropped)
// MISSING     — sample weighted by uninformative character count: characters
//               where the tip is fully ambiguous or inapplicable.
//               High-missingness taxa are hardest to score and most likely
//               to be trapped in suboptimal positions.
// COMBINED    — product of INSTABILITY and MISSING scores (normalised):
//               w(t) = instability(t) * (1 + miss_fraction(t)).
//               Targets taxa that are both unstably placed and data-poor.
//               Falls back to INSTABILITY when pool has < 2 trees, and to
//               MISSING when there is no missingness variation.
enum class PruneSelection { RANDOM = 0, INSTABILITY = 1, MISSING = 2, COMBINED = 3 };

struct PruneReinsertParams {
  int n_cycles = 1;
  double drop_fraction = 0.10;      // fraction of tips to drop
  int min_drop = 3;                 // floor: always drop at least this many
  int max_drop = 0;                 // 0 = no cap
  PruneSelection selection = PruneSelection::RANDOM;
  int tbr_max_moves = 0;           // TBR on reduced tree: 0 = converge
  int tbr_full_max_moves = 0;      // TBR on full tree after reinsert: 0 = converge
  bool nni_full = false;           // use NNI (not TBR) for full-tree polish
  int tbr_max_hits = 1;            // TBR max equal-score hits
  int tabu_size = 100;             // TBR tabu list size
};

struct PruneReinsertResult {
  double best_score;
  int n_improvements;              // cycles that improved the score
};

// Run one or more prune-reinsert cycles on the current tree.
// Modifies `tree` in-place.  Returns the best score achieved.
//
// The pool's SplitFrequencyTable is optional; only needed for
// instability-guided selection (ignored in RANDOM mode).
PruneReinsertResult prune_reinsert_search(
    TreeState& tree,
    DataSet& ds,
    const PruneReinsertParams& params,
    ConstraintData* cd = nullptr,
    const SplitFrequencyTable* split_freq = nullptr,
    std::function<bool()> check_timeout = nullptr);

} // namespace ts

#endif // TS_PRUNE_REINSERT_H
