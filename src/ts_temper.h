#ifndef TS_TEMPER_H
#define TS_TEMPER_H

// Stochastic TBR with Boltzmann acceptance (Layer 1) and simulated
// annealing schedule (Layer 3).
//
// Layer 1 (stochastic_tbr_phase) ported from Agent C's T-198 on
// feature/parallel-temper.  Layer 2 (multi-chain parallel tempering)
// lives on that branch and is NOT included here.
//
// Layer 1: sample random SPR moves on a clipped tree, accept/reject
//   with Boltzmann probability exp(-delta/T).  T=0 = strict hill-climbing.
// Layer 3: linear cooling schedule calling stochastic_tbr_phase() at
//   decreasing temperatures.

#include "ts_data.h"
#include "ts_tree.h"
#include "ts_constraint.h"
#include <functional>

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

// --- Layer 3: Simulated annealing (single-chain scheduled cooling) ---

struct AnnealParams {
  double t_start = 20.0;   // initial Boltzmann temperature
  double t_end = 0.0;      // final temperature (0 = strict at end)
  int n_phases = 5;        // number of temperature steps (linear schedule)
  int moves_per_phase = 0; // stochastic moves per phase; 0 = n_tip
};

struct AnnealResult {
  double best_score;     // best score seen across all phases
  double final_score;    // score at end of last phase
  int total_accepted;    // total accepted moves across all phases
  int total_improved;    // total strictly-improving moves
  int total_attempted;   // total move attempts
};

// Run simulated annealing: stochastic SPR at linearly decreasing
// temperature from t_start to t_end over n_phases steps.  Each step
// calls stochastic_tbr_phase().  The tree is NOT restored to its best
// state — it follows the annealing trajectory.  Caller should run TBR
// polish afterwards to converge to the nearest local optimum.
AnnealResult anneal_search(
    TreeState& tree, const DataSet& ds,
    const AnnealParams& params,
    ConstraintData* cd = nullptr,
    std::function<bool()> check_timeout = nullptr);

} // namespace ts

#endif // TS_TEMPER_H
