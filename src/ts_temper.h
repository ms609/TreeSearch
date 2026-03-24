#ifndef TS_TEMPER_H
#define TS_TEMPER_H

// Parallel tempering for parsimony search.
//
// Two layers:
// 1. Stochastic TBR: sample random SPR moves, accept/reject with
//    Boltzmann probability exp(-delta/T).  T=0 = strict hill-climbing.
// 2. Multi-chain framework: N chains at different temperatures, with
//    periodic Metropolis swaps.  The cold chain (T=0) uses standard TBR;
//    hot chains use stochastic TBR.

#include "ts_data.h"
#include "ts_tree.h"
#include "ts_constraint.h"
#include "ts_pool.h"
#include <functional>
#include <vector>

namespace ts {

// --- Layer 1: Stochastic TBR ---

struct TemperParams {
  double temperature = 1.0;  // Boltzmann temperature; 0 = strict
  int n_moves = 100;         // stochastic moves to attempt
};

struct TemperResult {
  double best_score;    // best score seen during this phase
  double final_score;   // score at end (chain may have wandered)
  int n_accepted;       // accepted moves (improvements + suboptimal)
  int n_improved;       // moves that strictly improved score
  int n_attempted;      // total attempts (may be < n_moves if interrupted)
};

// Run stochastic SPR with Boltzmann acceptance on `tree`.
// Modifies `tree` in place: at exit, the tree is in whatever state the
// chain wandered to (NOT necessarily the best-scoring tree seen).
TemperResult stochastic_tbr_phase(
    TreeState& tree, const DataSet& ds,
    const TemperParams& params,
    ConstraintData* cd = nullptr,
    std::function<bool()> check_timeout = nullptr);

// --- Layer 2: Multi-chain parallel tempering ---

struct PTParams {
  int n_chains = 4;                          // total chains (including cold)
  std::vector<double> temperatures;          // per-chain; [0]=0 (cold)
  int rounds = 5;                            // swap rounds
  int moves_per_round = 0;                   // stochastic moves per hot chain
                                             // per round; 0 = n_tip
};

struct PTResult {
  double best_score;          // best score found by any chain
  double cold_final_score;    // cold chain's final score
  int total_swaps_accepted;   // accepted chain swaps
  int total_swaps_attempted;  // total swap attempts
  int hot_discoveries;        // times a hot chain beat pool best
};

// Run parallel tempering search.
// `tree` is the starting tree (used for cold chain; hot chains get copies).
// `pool` receives any trees that beat the current pool best score.
// On return, `tree` holds the cold chain's final state.
PTResult parallel_temper_search(
    TreeState& tree, const DataSet& ds,
    const PTParams& params,
    TreePool* pool = nullptr,
    ConstraintData* cd = nullptr,
    std::function<bool()> check_timeout = nullptr);

} // namespace ts

#endif // TS_TEMPER_H
