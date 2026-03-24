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

TemperResult stochastic_tbr_phase(
    TreeState& tree, const DataSet& ds,
    const TemperParams& params,
    ConstraintData* cd = nullptr,
    std::function<bool()> check_timeout = nullptr);

// --- Layer 2: Multi-chain parallel tempering ---

struct PTParams {
  int n_chains = 4;
  std::vector<double> temperatures;   // per-chain; [0]=0 (cold)
  int rounds = 5;
  int moves_per_round = 0;            // 0 = n_tip
};

struct PTResult {
  double best_score;
  double cold_final_score;
  int total_swaps_accepted;
  int total_swaps_attempted;
  int hot_discoveries;
};

// --- Diagnostics ---

struct SwapDiag {
  int round;
  int pair_lo;
  int pair_hi;
  double score_lo;
  double score_hi;
  double metropolis_prob;
  bool accepted;
};

struct ChainRoundDiag {
  int round;
  int chain_idx;
  double temperature;
  double score_before;
  double score_after;
  int n_accepted;
  int n_improved;
  int n_attempted;
  double elapsed_ms;
};

struct PTDiagnostics {
  std::vector<SwapDiag> swap_log;
  std::vector<ChainRoundDiag> chain_log;
  std::vector<double> pair_acceptance_rates;  // length = n_chains - 1
  std::vector<double> cold_scores;            // per-round after swaps
  double cold_tbr_total_ms = 0.0;
  double hot_stochastic_total_ms = 0.0;
  double total_pt_ms = 0.0;
  int cold_improvements_from_swaps = 0;
};

// If diag is non-null, detailed diagnostics are collected.
PTResult parallel_temper_search(
    TreeState& tree, const DataSet& ds,
    const PTParams& params,
    TreePool* pool = nullptr,
    ConstraintData* cd = nullptr,
    std::function<bool()> check_timeout = nullptr,
    PTDiagnostics* diag = nullptr);

} // namespace ts

#endif // TS_TEMPER_H
