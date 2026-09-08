#ifndef TS_NNI_PERTURB_H
#define TS_NNI_PERTURB_H

// Stochastic NNI-perturbation: escape local optima by randomly applying
// NNI swaps to a fraction of internal branches, then re-optimizing.
//
// Complementary to the weight-perturbation ratchet (ts_ratchet.h):
// - Ratchet perturbs the objective function (character weights), then
//   re-optimizes the perturbed landscape before restoring weights.
// - NNI-perturbation perturbs the topology directly (random NNI swaps),
//   then re-optimizes on the original landscape.
//
// Inspired by IQ-TREE's doRandomNNIs() (Nguyen et al. 2015).
//
// Algorithm:
// 1. Collect all internal NNI edges, shuffle randomly.
// 2. For each edge (with probability perturb_fraction), apply a random
//    NNI swap — but skip edges that conflict with already-applied swaps
//    (two NNI operations conflict if their edges are adjacent in the tree).
// 3. Rebuild postorder and full rescore.
// 4. TBR to a new local optimum.
// 5. If improved, keep; otherwise revert to best topology.

#include "ts_data.h"
#include "ts_tree.h"
#include "ts_constraint.h"
#include <functional>

namespace ts {

struct NNIPerturbParams {
  int n_cycles = 5;              // perturbation + TBR cycles
  double perturb_fraction = 0.5; // fraction of internal branches to perturb
  int max_hits = 1;              // for the post-perturbation TBR
  int tabu_size = 0;             // tabu list for TBR (0 = disabled)
};

struct NNIPerturbResult {
  double best_score;
  int n_cycles_completed;
  int total_tbr_moves;
  int n_escapes;               // cycles that improved the best score
};

// Apply random NNI-perturbation cycles on `tree` with dataset `ds`.
// Modifies `tree` in place to the best tree found across all cycles.
NNIPerturbResult nni_perturb_search(
    TreeState& tree, const DataSet& ds,
    const NNIPerturbParams& params,
    ConstraintData* cd = nullptr,
    std::function<bool()> check_timeout = nullptr);

// Low-level: apply a single batch of random compatible NNI swaps.
// Returns the number of swaps applied.
int random_nni_perturb(TreeState& tree, double fraction);

} // namespace ts

#endif // TS_NNI_PERTURB_H
