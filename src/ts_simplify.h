#ifndef TS_SIMPLIFY_H
#define TS_SIMPLIFY_H

// Character simplification: precompute topology-independent steps.
//
// Three transforms applied to each pattern (non-NA characters only):
//   1. Remove parsimony-uninformative characters (constant steps on all trees)
//   2. Remove singleton states (unique unambiguous states → wildcard + offset)
//   3. Remove redundant ambiguity states (never unambiguous, always co-occur
//      with another state in every token)
//
// After simplification, each pattern has:
//   - A possibly-reduced set of tip tokens (fewer states)
//   - A precomputed step offset (topology-independent steps removed)
//   - An informativeness flag (false → exclude from scoring entirely)
//
// The offsets are stored in DataSet and applied:
//   - EW: ds.ew_offset added to score_tree() return value
//   - IW: min_steps[p] reduced by offset[p], so extra_p is unchanged
//   - Profile: offset[p] added back to char_steps[p] before lookup

#include <cstdint>
#include <vector>

namespace ts {

struct SimplifiedPattern {
  int original_index;                   // pattern index (into original arrays)
  int precomputed_steps;                // topology-independent step offset
  bool informative;                     // false → skip entirely
  bool has_genuine_inapp;               // true → some tip has genuine "-" (not just "?")
  std::vector<uint32_t> tip_tokens;     // simplified token bitmask per tip
  int n_states_remaining;               // number of applicable states after simplification
};

struct SimplificationResult {
  std::vector<SimplifiedPattern> patterns;
  int n_patterns_removed;               // uninformative patterns removed
  int n_states_reduced;                 // total state-count reduction across patterns
  int total_offset_steps;               // sum of precomputed_steps (unweighted)
};

// Simplify all patterns. Skips patterns with inapplicable tokens (Phase 1).
//
// token_states: per-token bitmask of possible states (from contrast matrix)
// tip_data_r:   n_tips × n_patterns, 1-based token indices (column-major)
// weight_r:     per-pattern weights
// inapp_state:  column index of "-" state, or -1 if none
SimplificationResult simplify_patterns(
    const std::vector<uint32_t>& token_states,
    const int* tip_data_r, int n_tips, int n_patterns,
    const int* weight_r,
    int n_states, int inapp_state);

} // namespace ts

#endif // TS_SIMPLIFY_H
