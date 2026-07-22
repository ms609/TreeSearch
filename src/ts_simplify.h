#ifndef TS_SIMPLIFY_H
#define TS_SIMPLIFY_H

// Character simplification: precompute topology-independent steps.
//
// Every reduction MUST be topology-independent (score-preserving on ALL trees).
// The optimal resolution of an ambiguous token is topology-DEPENDENT, so the
// transforms may only touch a state that appears in <= 1 tip (a lone leaf
// option cannot form a cost-free clade). Transforms applied to each pattern
// (non-NA characters only):
//   1. Remove parsimony-uninformative characters. Exact for fully-unambiguous
//      characters (classical criterion); for characters with genuine ambiguity
//      only the provable case is used — a state shared by every tip (score 0).
//      (A prior caterpillar-sampling heuristic mis-classified ambiguous
//      characters and over-/under-counted; it was removed.)
//   2. Remove singleton states: a state unambiguous in exactly one tip AND
//      absent from every other tip (a genuine autapomorphy) → wildcard + 1.
//   3. Remove a redundant ambiguity state that occurs in exactly one tip.
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
