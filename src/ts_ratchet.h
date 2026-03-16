#ifndef TS_RATCHET_H
#define TS_RATCHET_H

// Parsimony ratchet: escape local optima by temporarily perturbing character
// weights (via active_mask zeroing and/or upweight_mask doubling), then
// searching the perturbed landscape with TBR before reverting to original
// weights.
//
// Supports three perturbation modes: zero-only (original Nixon/Goloboff),
// upweight-only, and mixed (zero some + double others).
//
// Optional adaptive tuning adjusts perturbation intensity based on escape rate.

#include "ts_data.h"
#include "ts_tree.h"
#include "ts_constraint.h"

namespace ts {

enum class PerturbMode : int {
  ZERO_ONLY    = 0,  // Zero each char with probability p
  UPWEIGHT_ONLY = 1, // Double each char with probability p (keep all active)
  MIXED        = 2   // Zero with prob p, upweight with prob p (disjoint sets)
};

struct RatchetParams {
  int n_cycles = 10;           // number of perturbation+search cycles
  double perturb_prob = 0.04;  // probability per character
  int max_hits = 1;            // for the search-phase TBR
  PerturbMode perturb_mode = PerturbMode::ZERO_ONLY;

  // Inner search intensity (configurable; was hardcoded)
  int perturb_max_moves = 0;         // 0 = auto: max(20, min(200, n_tip/8))
  bool perturb_accept_equal = true;  // Accept equal-score moves during perturb?

  // Adaptive perturbation
  bool adaptive = false;             // Auto-tune perturb_prob?
  double target_escape_rate = 0.3;   // Fraction of cycles that find improvement
  double adapt_min_prob = 0.02;      // Lower bound for adaptive tuning
  double adapt_max_prob = 0.50;      // Upper bound for adaptive tuning

  int tabu_size = 0;                 // Tabu list size for TBR calls (0 = disabled)
};

struct RatchetResult {
  double best_score;
  int n_cycles_completed;
  int total_tbr_moves;         // across all cycles
  int n_escapes;               // cycles that found a new local optimum
  double final_perturb_prob;   // final perturbation probability (adaptive)
};

RatchetResult ratchet_search(TreeState& tree, DataSet& ds,
                             const RatchetParams& params,
                             ConstraintData* cd = nullptr);

} // namespace ts

#endif // TS_RATCHET_H
