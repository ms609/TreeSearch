#ifndef TS_TEMPER_H
#define TS_TEMPER_H

// Parallel tempering building blocks for parsimony search.
//
// Stochastic TBR: instead of exhaustive TBR (O(n^2) per pass, picking the
// best move), sample random clip+regraft moves and accept/reject with
// Boltzmann probability exp(-delta/T).  Temperature T=0 means strict
// hill-climbing; higher T accepts worse moves more readily.
//
// Designed as the "hot chain" operator for a parallel tempering framework.
// The cold chain (T=0) uses standard TBR search (ts_tbr.h).

#include "ts_data.h"
#include "ts_tree.h"
#include "ts_constraint.h"
#include <functional>

namespace ts {

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
// Caller should snapshot the tree before calling if the best tree is needed.
TemperResult stochastic_tbr_phase(
    TreeState& tree, const DataSet& ds,
    const TemperParams& params,
    ConstraintData* cd = nullptr,
    std::function<bool()> check_timeout = nullptr);

} // namespace ts

#endif // TS_TEMPER_H
