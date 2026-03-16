#include "ts_simplify.h"
#include <algorithm>

namespace ts {

// Count how many tips have state s as their ONLY state (unambiguous).
// Also count how many tips have state s in ANY token (including ambiguous).
static void count_state_occurrences(
    const std::vector<uint32_t>& tips, int n_tips, int n_states,
    int inapp_state,
    std::vector<int>& unambig_count,    // [n_states]
    std::vector<int>& any_count,        // [n_states]
    std::vector<int>& unambig_tip_idx)  // [n_states] last tip with unambig
{
  std::fill(unambig_count.begin(), unambig_count.end(), 0);
  std::fill(any_count.begin(), any_count.end(), 0);
  std::fill(unambig_tip_idx.begin(), unambig_tip_idx.end(), -1);

  for (int tip = 0; tip < n_tips; ++tip) {
    uint32_t tok = tips[tip];
    // Count applicable states in this token
    int n_set = 0;
    int single_state = -1;
    for (int s = 0; s < n_states; ++s) {
      if (s == inapp_state) continue;
      if (tok & (1u << s)) {
        ++n_set;
        single_state = s;
        any_count[s]++;
      }
    }
    if (n_set == 1 && single_state >= 0) {
      unambig_count[single_state]++;
      unambig_tip_idx[single_state] = tip;
    }
  }
}

// Check if a pattern is parsimony-uninformative: at most one state appears
// in 2+ tips (unambiguously).
static bool is_uninformative(const std::vector<int>& unambig_count,
                              int n_states, int inapp_state) {
  int informative_states = 0;
  for (int s = 0; s < n_states; ++s) {
    if (s == inapp_state) continue;
    if (unambig_count[s] >= 2) ++informative_states;
  }
  return informative_states <= 1;
}

// Compute fixed step count for an uninformative pattern.
// = number of states with unambig_count >= 1, minus 1 (if any states present).
// This is the minimum number of state changes on any tree.
static int compute_fixed_steps(const std::vector<int>& unambig_count,
                                int n_states, int inapp_state) {
  int distinct = 0;
  for (int s = 0; s < n_states; ++s) {
    if (s == inapp_state) continue;
    if (unambig_count[s] >= 1) ++distinct;
  }
  return (distinct > 0) ? distinct - 1 : 0;
}

// Get the "all applicable states" mask (excluding inapp_state).
static uint32_t all_applicable_mask(const std::vector<uint32_t>& tips,
                                      int n_tips, int n_states,
                                      int inapp_state) {
  uint32_t all = 0;
  for (int tip = 0; tip < n_tips; ++tip) {
    all |= tips[tip];
  }
  if (inapp_state >= 0) all &= ~(1u << inapp_state);
  return all;
}

SimplificationResult simplify_patterns(
    const std::vector<uint32_t>& token_states,
    const int* tip_data_r, int n_tips, int n_patterns,
    const int* weight_r,
    int n_states, int inapp_state)
{
  SimplificationResult result;
  result.n_patterns_removed = 0;
  result.n_states_reduced = 0;
  result.total_offset_steps = 0;
  result.patterns.resize(n_patterns);

  std::vector<int> unambig_count(n_states);
  std::vector<int> any_count(n_states);
  std::vector<int> unambig_tip_idx(n_states);

  for (int p = 0; p < n_patterns; ++p) {
    auto& sp = result.patterns[p];
    sp.original_index = p;
    sp.precomputed_steps = 0;
    sp.informative = true;

    // Skip zero-weight patterns (they'll be removed by build_dataset anyway)
    if (weight_r[p] == 0) {
      sp.informative = false;
      sp.n_states_remaining = 0;
      sp.tip_tokens.resize(n_tips, 0);
      result.n_patterns_removed++;
      continue;
    }

    // Build per-tip token bitmasks from the original data
    sp.tip_tokens.resize(n_tips);
    bool has_inapp = false;
    for (int tip = 0; tip < n_tips; ++tip) {
      int token = tip_data_r[tip + n_tips * p] - 1;  // 1-based to 0-based
      sp.tip_tokens[tip] = token_states[token];
      if (inapp_state >= 0 &&
          (token_states[token] & (1u << inapp_state))) {
        has_inapp = true;
      }
    }

    // Phase 1: skip inapplicable characters entirely
    if (has_inapp) {
      // Count states for metadata only
      uint32_t all = all_applicable_mask(sp.tip_tokens, n_tips, n_states,
                                          inapp_state);
      int nc = 0;
      for (int s = 0; s < n_states; ++s) {
        if (s == inapp_state) continue;
        if (all & (1u << s)) ++nc;
      }
      sp.n_states_remaining = nc;
      continue;
    }

    // --- Transforms 2 & 3: iterate until stable ---
    // Transform 2 (singleton removal) can create situations where
    // Transform 3 (ambiguity removal) applies, and vice versa.
    bool changed_outer = true;
    while (changed_outer) {
      changed_outer = false;

      // Transform 3: Remove redundant ambiguity states
      // A state is redundant if it never appears unambiguously AND every
      // token containing it also contains at least one other state that
      // appears unambiguously somewhere.
      bool changed_t3 = true;
      while (changed_t3) {
        changed_t3 = false;
        count_state_occurrences(sp.tip_tokens, n_tips, n_states,
                                inapp_state, unambig_count, any_count,
                                unambig_tip_idx);
        for (int s = 0; s < n_states; ++s) {
          if (s == inapp_state) continue;
          if (any_count[s] == 0) continue;
          if (unambig_count[s] > 0) continue;

          bool safe_to_remove = true;
          for (int tip = 0; tip < n_tips; ++tip) {
            if (!(sp.tip_tokens[tip] & (1u << s))) continue;
            uint32_t other = sp.tip_tokens[tip] & ~(1u << s);
            if (inapp_state >= 0) other &= ~(1u << inapp_state);
            bool has_other_unambig = false;
            for (int s2 = 0; s2 < n_states; ++s2) {
              if (s2 == inapp_state || s2 == s) continue;
              if ((other & (1u << s2)) && unambig_count[s2] > 0) {
                has_other_unambig = true;
                break;
              }
            }
            if (!has_other_unambig) {
              safe_to_remove = false;
              break;
            }
          }

          if (safe_to_remove) {
            uint32_t mask = ~(1u << s);
            for (int tip = 0; tip < n_tips; ++tip) {
              sp.tip_tokens[tip] &= mask;
            }
            changed_t3 = true;
            changed_outer = true;
            result.n_states_reduced++;
          }
        }
      }

      // Transform 2: Remove singleton states
      // A state is a removable singleton if it appears unambiguously in
      // exactly 1 tip, AND every ambiguity token containing it also has
      // at least one other state with unambig_count >= 2 (so the
      // optimizer never needs to resolve the ambiguity to this state).
      bool changed_t2 = true;
      while (changed_t2) {
        changed_t2 = false;
        count_state_occurrences(sp.tip_tokens, n_tips, n_states,
                                inapp_state, unambig_count, any_count,
                                unambig_tip_idx);

        for (int s = 0; s < n_states; ++s) {
          if (s == inapp_state) continue;
          if (unambig_count[s] != 1) continue;

          // Check all ambiguity tokens containing s
          bool safe = true;
          for (int tip = 0; tip < n_tips; ++tip) {
            uint32_t tok = sp.tip_tokens[tip];
            if (!(tok & (1u << s))) continue;
            // Count applicable states in this token
            int n_set = 0;
            for (int s2 = 0; s2 < n_states; ++s2) {
              if (s2 != inapp_state && (tok & (1u << s2))) ++n_set;
            }
            if (n_set == 1) continue;  // unambiguous token — the singleton tip itself

            // Ambiguity token: must contain another state with unambig >= 2
            bool has_dominant = false;
            for (int s2 = 0; s2 < n_states; ++s2) {
              if (s2 == inapp_state || s2 == s) continue;
              if ((tok & (1u << s2)) && unambig_count[s2] >= 2) {
                has_dominant = true;
                break;
              }
            }
            if (!has_dominant) { safe = false; break; }
          }
          if (!safe) continue;

          int tip = unambig_tip_idx[s];

          uint32_t all = all_applicable_mask(sp.tip_tokens, n_tips,
                                              n_states, inapp_state);
          uint32_t wildcard = all & ~(1u << s);
          if (wildcard == 0) continue;

          sp.tip_tokens[tip] = wildcard;
          sp.precomputed_steps += 1;
          changed_t2 = true;
          changed_outer = true;
          result.n_states_reduced++;
          break;
        }
      }
    }

    // --- Transform 1: Check informativeness ---
    count_state_occurrences(sp.tip_tokens, n_tips, n_states,
                            inapp_state, unambig_count, any_count,
                            unambig_tip_idx);

    if (is_uninformative(unambig_count, n_states, inapp_state)) {
      // Add the fixed steps of the reduced character to the accumulated
      // singleton steps from Transform 2
      sp.precomputed_steps += compute_fixed_steps(unambig_count, n_states,
                                                   inapp_state);
      sp.informative = false;
      sp.n_states_remaining = 0;
      result.n_patterns_removed++;
      result.total_offset_steps += sp.precomputed_steps;
    } else {
      // Count remaining states
      uint32_t all = all_applicable_mask(sp.tip_tokens, n_tips, n_states,
                                          inapp_state);
      int nc = 0;
      for (int s = 0; s < n_states; ++s) {
        if (s == inapp_state) continue;
        if (all & (1u << s)) ++nc;
      }
      sp.n_states_remaining = nc;
      result.total_offset_steps += sp.precomputed_steps;
    }
  }

  return result;
}

} // namespace ts
